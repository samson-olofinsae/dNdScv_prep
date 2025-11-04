#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, glob, shutil, argparse, shlex
import subprocess as sp
from pathlib import Path
from typing import List
import pandas as pd

# --------------------------- helpers (shell + fs) ---------------------------

def qq(p) -> str:
    return shlex.quote(str(p))

def sh(cmd, check=True):
    print(f"[cmd] {cmd}")
    return sp.run(cmd, shell=True, check=check)

def exists_or_die(p: Path, what="file"):
    if not p.exists():
        sys.exit(f"ERROR: {what} not found: {p}")
    return p

def which_or_die(tool):
    if shutil.which(tool) is None:
        sys.exit(f"ERROR: '{tool}' not in PATH. Please install or load your env.")
    return tool

def tabix_index(vcf_gz: Path):
    sh(f"tabix -p vcf {qq(vcf_gz)}")

def bcftools_query_to_tsv(vcf_gz: Path, out_tsv: Path):
    qfmt = r'[%CHROM\t%POS\t%REF\t%ALT\n]'
    sh(f"bcftools query -f '{qfmt}' {qq(vcf_gz)} > {qq(out_tsv)}")

# --------------------------- helpers (uniform prompts) ---------------------------

def prompt_value(label, default=None, current=None, cast=None, validator=None):
    shown_default = current if (current not in [None, ""]) else default
    prompt = f"{label}"
    if shown_default is not None:
        prompt += f" [{shown_default}]"
    prompt += ": "
    entered = input(prompt).strip()
    value = entered if entered != "" else shown_default
    if cast:
        try:
            value = cast(value)
        except Exception:
            sys.exit(f"ERROR: invalid value for '{label}': {value}")
    if validator:
        value = validator(value)
    return value

def prompt_yes_no(label, default="Y", assume_yes=False):
    default = default.upper()
    suffix = "Y/n" if default == "Y" else "y/N"
    if assume_yes:
        print(f"{label} [{suffix}]: {default.lower()}")
        return default == "Y"
    resp = input(f"{label} [{suffix}]: ").strip().lower()
    if not resp:
        return default == "Y"
    return resp.startswith("y")

# --------------------------- discovery ---------------------------

def discover_r1s_from_source(source: str) -> List[Path]:
    p = Path(source).expanduser()
    pattern = str(p / "*_R1.fastq.gz") if p.is_dir() else source
    r1s = [Path(x) for x in glob.glob(pattern)]
    return sorted(r1s)

# --------------------------- sanitizer helper ---------------------------

def _sanitize_for_dndscv(df: pd.DataFrame, sample_name: str):
    '''
    Make sure ref/mut contain only A/C/G/T/N or '-' and drop spanning deletions (ALT='*').
    Returns (clean_df, dropped_star, dropped_invalid).
    '''
    if df.empty:
        return df, 0, 0

    df["ref"] = df["ref"].astype(str).str.strip().str.upper()
    df["mut"] = df["mut"].astype(str).str.strip().str.upper()

    mask_star = (df["mut"] == "*")
    dropped_star = int(mask_star.sum())
    df = df[~mask_star].copy()

    valid_nt = set(list("ACGTN") + ["-"])

    def _ok(s: str) -> bool:
        return set(s) <= valid_nt

    mask_valid = df["ref"].apply(_ok) & df["mut"].apply(_ok)
    dropped_invalid = int((~mask_valid).sum())

    if dropped_star or dropped_invalid:
        print(f"[warn] {sample_name}: dropped {dropped_star} spanning-deletion rows (ALT='*') "
              f"and {dropped_invalid} invalid rows (non-ACGTN/-) before merge.")

    return df[mask_valid].copy(), dropped_star, dropped_invalid

# --------------------------- core runner ---------------------------

def run_pipeline(ref: Path, r1_glob_or_dir: str, outdir: Path, threads: int, ploidy: int,
                 min_qual: int, min_dp: int, auto_index: bool, assume_yes: bool):

    for tool in ["bwa", "samtools", "bcftools", "bgzip", "tabix"]:
        which_or_die(tool)

    wd = Path.cwd()
    print("Accepted format for FASTQs is *_R1.fastq.gz and *_R2.fastq.gz")
    print(f"Your working directory: {wd}")

    (outdir/"bam").mkdir(parents=True, exist_ok=True)
    (outdir/"bcf").mkdir(parents=True, exist_ok=True)
    (outdir/"vcf"/"snvs").mkdir(parents=True, exist_ok=True)
    (outdir/"vcf"/"indels").mkdir(parents=True, exist_ok=True)
    (outdir/"dndscv").mkdir(parents=True, exist_ok=True)

    ref = exists_or_die(ref, "reference FASTA").resolve()
    print(f"Your ref genome file path: {ref}")

    if not Path(str(ref)+".bwt").exists():
        if auto_index or assume_yes:
            sh(f"bwa index {qq(ref)}")
        else:
            if prompt_yes_no(f"bwa index not found for {ref}. Create it now?", "Y", False):
                sh(f"bwa index {qq(ref)}")
            else:
                sys.exit("Cannot continue without BWA index.")
    if not Path(str(ref)+".fai").exists():
        if auto_index or assume_yes:
            sh(f"samtools faidx {qq(ref)}")
        else:
            if prompt_yes_no(f"samtools faidx not found for {ref}. Create it now?", "Y", False):
                sh(f"samtools faidx {qq(ref)}")
            else:
                sys.exit("Cannot continue without FASTA index (.fai).")

    r1s = discover_r1s_from_source(r1_glob_or_dir)
    if not r1s:
        sys.exit(f"ERROR: No R1 FASTQs found with source: {r1_glob_or_dir}")

    print(f"Found {len(r1s)} R1 files (source: {r1_glob_or_dir})")
    for p in r1s[:10]:
        print(f"  - {p}")
    if len(r1s) > 10:
        print(f"  ... and {len(r1s)-10} more")
    print(f"FASTQ discovery pattern/source: {r1_glob_or_dir}")

    print("\nSummary:")
    print(f"  Working dir : {wd}")
    print(f"  Reference   : {ref}")
    print(f"  Outputs to  : {outdir.resolve()}")
    print(f"  Threads     : {threads}")
    print(f"  Ploidy      : {ploidy}")
    print(f"  Filters     : QUAL>={min_qual}, DP>={min_dp}")
    print(f"  Samples     : {len(r1s)}")
    if not prompt_yes_no("Proceed with processing?", "Y", assume_yes):
        sys.exit("Aborted by user.")
    if assume_yes:
        print("Proceeding (auto-accepted with --yes).")

    snv_rows, indel_rows = [], []

    for r1 in r1s:
        base = r1.name.split("_R1.fastq.gz")[0]
        r2 = r1.with_name(f"{base}_R2.fastq.gz")
        if not r2.exists():
            print(f"[WARN] Skipping {base} (missing {r2.name})")
            continue

        print(f"\n=== Processing sample: {base} ===")
        bam_sorted = outdir/"bam"/f"{base}.aligned.sorted.bam"
        raw_bcf    = outdir/"bcf"/f"{base}_raw.bcf"
        vcf_all    = outdir/"vcf"/f"{base}_final_variants.vcf.gz"
        vcf_snvs   = outdir/"vcf"/"snvs"/f"{base}_snvs.vcf.gz"
        vcf_indels = outdir/"vcf"/"indels"/f"{base}_indels.vcf.gz"

        sh(f"bwa mem -t {threads} {qq(ref)} {qq(r1)} {qq(r2)} | samtools sort -@ {threads} -o {qq(bam_sorted)}")
        sh(f"samtools index -@ {threads} {qq(bam_sorted)}")
        sh(f"bcftools mpileup -O b -o {qq(raw_bcf)} -f {qq(ref)} {qq(bam_sorted)} && "
           f"bcftools call --ploidy {ploidy} -m -v -o - {qq(raw_bcf)} | "
           f"bcftools filter -e 'QUAL<{min_qual} || DP<{min_dp}' - | bgzip -c > {qq(vcf_all)}")
        tabix_index(vcf_all)

        sh(f"bcftools view -v snps   {qq(vcf_all)} -Oz -o {qq(vcf_snvs)}")
        tabix_index(vcf_snvs)
        sh(f"bcftools view -v indels {qq(vcf_all)} -Oz -o {qq(vcf_indels)}")
        tabix_index(vcf_indels)

        snv_tsv, indel_tsv = outdir/"vcf"/"snvs"/f"{base}_snvs.tsv", outdir/"vcf"/"indels"/f"{base}_indels.tsv"
        bcftools_query_to_tsv(vcf_snvs, snv_tsv)
        bcftools_query_to_tsv(vcf_indels, indel_tsv)

        def load_tbl(p: Path):
            if not p.exists() or p.stat().st_size == 0:
                return pd.DataFrame(columns=["sampleID","chr","pos","ref","mut"])
            df = pd.read_csv(p, sep="\t", header=None, names=["chr","pos","ref","mut"])
            df.insert(0, "sampleID", base)
            return df

        snv_df = load_tbl(snv_tsv)
        snv_df, star_snv, bad_snv = _sanitize_for_dndscv(snv_df, base + " (SNVs)")

        indel_df = load_tbl(indel_tsv)
        indel_df, star_indel, bad_indel = _sanitize_for_dndscv(indel_df, base + " (INDELs)")

        snv_rows.append(snv_df)
        indel_rows.append(indel_df)

    ddir = outdir/"dndscv"
    snv_all = pd.concat(snv_rows, ignore_index=True) if snv_rows else pd.DataFrame(columns=["sampleID","chr","pos","ref","mut"])
    indel_all = pd.concat(indel_rows, ignore_index=True) if indel_rows else pd.DataFrame(columns=["sampleID","chr","pos","ref","mut"])
    dndscv_all = pd.concat([snv_all, indel_all], ignore_index=True)

    required_cols = ["sampleID","chr","pos","ref","mut"]
    if list(dndscv_all.columns) != required_cols:
        sys.exit(f"ERROR: dNdScv schema mismatch. Got columns {list(dndscv_all.columns)} expected {required_cols}")
    if not dndscv_all["pos"].apply(lambda x: str(x).isdigit() and int(x) > 0).all():
        sys.exit("ERROR: Invalid POS values detected.")

    ddir.mkdir(parents=True, exist_ok=True)
    (ddir/"combined_snv_variants.csv").write_text(snv_all.to_csv(index=False))
    (ddir/"combined_indels_variants.csv").write_text(indel_all.to_csv(index=False))
    (ddir/"dndscv_input.csv").write_text(dndscv_all.to_csv(index=False))

    print("\nDone.")
    print(f"  dNdScv input : {ddir/'dndscv_input.csv'}")
    print(f"  SNVs   : {ddir/'combined_snv_variants.csv'}")
    print(f"  INDELs : {ddir/'combined_indels_variants.csv'}")

# --------------------------- entrypoint ---------------------------

def build_parser():
    p = argparse.ArgumentParser(description="FASTQ→BAM→VCF→dNdScv prep workflow (interactive or non-interactive)")
    p.add_argument("--non-interactive", action="store_true", help="Run fully non-interactive (for CI/CD).")
    p.add_argument("--yes", "--assume-yes", dest="assume_yes", action="store_true",
                   help="Auto-accept prompts but still allow typing.")
    p.add_argument("--ref", help="Path to reference FASTA")
    p.add_argument("--r1-source", default="*_R1.fastq.gz", help="Directory or glob for R1 FASTQs.")
    p.add_argument("--outdir", default="results", help="Output directory.")
    p.add_argument("--threads", type=int, default=4, help="Number of threads.")
    p.add_argument("--ploidy", type=int, default=2, help="Ploidy for bcftools call.")
    p.add_argument("--min-qual", type=int, default=20, help="Minimum QUAL.")
    p.add_argument("--min-dp", type=int, default=10, help="Minimum depth.")
    p.add_argument("--auto-index", action="store_true", help="Auto-create FASTA indexes if missing.")
    return p

def main():
    args = build_parser().parse_args()
    if args.non_interactive:
        if not args.ref:
            sys.exit("ERROR: --non-interactive requires --ref <fasta>")
        run_pipeline(Path(args.ref).expanduser().resolve(), args.r1_source,
                     Path(args.outdir), args.threads, args.ploidy,
                     args.min_qual, args.min_dp, args.auto_index, True)
        return

    print("Accepted format for FASTQs is *_R1.fastq.gz and *_R2.fastq.gz")
    print("TIP: Reference FASTA path example: ./demo_ref/demo.fa")
    print("TIP: FASTQ glob (R1) example:      ./demo_fastq/*_R1.fastq.gz\n")

    wd = Path.cwd()
    default_ref = "./demo_ref/demo.fa" if (wd/"demo_ref"/"demo.fa").exists() else None
    ref_input = prompt_value("Enter path to reference FASTA", default=default_ref, current=args.ref)
    ref_path = Path(ref_input).expanduser()
    if not ref_path.is_absolute():
        ref_path = (wd / ref_path).resolve()

    r1_source = prompt_value("FASTQ location (directory OR glob)", "*_R1.fastq.gz", args.r1_source)
    outdir = Path(prompt_value("Choose output directory", args.outdir, args.outdir))
    threads = prompt_value("Threads to use", str(args.threads), str(args.threads), int)
    ploidy  = prompt_value("Ploidy for bcftools call", str(args.ploidy), str(args.ploidy), int)

    run_pipeline(ref_path, r1_source, outdir, threads, ploidy,
                 args.min_qual, args.min_dp, args.auto_index, args.assume_yes)

if __name__ == "__main__":
    main()
