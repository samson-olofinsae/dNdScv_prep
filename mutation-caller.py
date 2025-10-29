#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, glob, shutil, argparse, shlex
import subprocess as sp
from pathlib import Path
from typing import List
import pandas as pd

# --------------------------- helpers ---------------------------

def qq(p) -> str:
    """Shell-safe quote a path/arg."""
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

def prompt_yn(msg, default="y"):
    d = default.lower()
    while True:
        ans = input(f"{msg} [{'Y/n' if d=='y' else 'y/N'}]: ").strip().lower()
        if ans == "" and d in ("y","n"):
            return d == "y"
        if ans in ("y","yes"): return True
        if ans in ("n","no"):  return False
        print("Please answer y or n.")

def prompt_str(msg, default=None, allow_blank=False):
    while True:
        ans = input(f"{msg}{f' [{default}]' if default else ''}: ").strip()
        if ans == "" and default is not None:
            return default
        if ans == "" and allow_blank:
            return ""
        if ans != "":
            return ans
        print("Please enter a value.")

def prompt_int(msg, default=None, minval=None, maxval=None):
    while True:
        raw = input(f"{msg}{f' [{default}]' if default is not None else ''}: ").strip()
        if raw == "" and default is not None:
            val = default
        else:
            try:
                val = int(raw)
            except ValueError:
                print("Please enter a whole number.")
                continue
        if minval is not None and val < minval:
            print(f"Must be ≥ {minval}."); continue
        if maxval is not None and val > maxval:
            print(f"Must be ≤ {maxval}."); continue
        return val

def tabix_index(vcf_gz: Path):
    sh(f"tabix -p vcf {qq(vcf_gz)}")

def bcftools_query_to_tsv(vcf_gz: Path, out_tsv: Path):
    qfmt = r'[%CHROM\t%POS\t%REF\t%ALT\n]'
    sh(f"bcftools query -f '{qfmt}' {qq(vcf_gz)} > {qq(out_tsv)}")

def discover_r1s_from_source(source: str) -> List[Path]:
    """
    If 'source' is a directory, use '<dir>/*_R1.fastq.gz'.
    Otherwise treat 'source' as a glob (absolute or relative).
    """
    p = Path(source).expanduser()
    if p.is_dir():
        pattern = str(p / "*_R1.fastq.gz")
    else:
        pattern = source
    r1s = [Path(x) for x in glob.glob(pattern)]
    return sorted(r1s)

# --------------------------- core runner ---------------------------

def run_pipeline(ref: Path, r1_glob_or_dir: str, outdir: Path, threads: int, ploidy: int,
                 min_qual: int, min_dp: int, auto_index: bool, assume_yes: bool,
                 run_dndscv: bool=False, assembly: str="auto", refdb_rds: str=None):

    # tool checks
    for tool in ["bwa", "samtools", "bcftools", "bgzip", "tabix"]:
        which_or_die(tool)

    wd = Path.cwd()
    print("Accepted format for FASTQs is *_R1.fastq.gz and *_R2.fastq.gz")
    print(f"Your working directory: {wd}")

    # Prepare dirs
    (outdir/"bam").mkdir(parents=True, exist_ok=True)
    (outdir/"bcf").mkdir(parents=True, exist_ok=True)
    (outdir/"vcf"/"snvs").mkdir(parents=True, exist_ok=True)
    (outdir/"vcf"/"indels").mkdir(parents=True, exist_ok=True)
    (outdir/"dndscv").mkdir(parents=True, exist_ok=True)

    # Reference
    ref = exists_or_die(ref, "reference FASTA").resolve()
    print(f"Your ref genome file path: {ref}")

    # Reference indexes (bwa + faidx)
    if not Path(str(ref)+".bwt").exists():
        if auto_index or assume_yes:
            sh(f"bwa index {qq(ref)}")
        else:
            if prompt_yn(f"bwa index not found for {ref}. Create it now?", default="y"):
                sh(f"bwa index {qq(ref)}")
            else:
                sys.exit("Cannot continue without BWA index.")
    if not Path(str(ref)+".fai").exists():
        if auto_index or assume_yes:
            sh(f"samtools faidx {qq(ref)}")
        else:
            if prompt_yn(f"samtools faidx not found for {ref}. Create it now?", default="y"):
                sh(f"samtools faidx {qq(ref)}")
            else:
                sys.exit("Cannot continue without FASTA index (.fai).")

    # Detect FASTQs anywhere
    r1s = discover_r1s_from_source(r1_glob_or_dir)
    if not r1s:
        sys.exit(f"ERROR: No R1 FASTQs found with source: {r1_glob_or_dir}")
    print(f"Found {len(r1s)} R1 files (source: {r1_glob_or_dir})")
    for p in r1s[:10]:
        print(f"  - {p}")
    if len(r1s) > 10:
        print(f"  ... and {len(r1s)-10} more")

    print(f"FASTQ discovery pattern/source: {r1_glob_or_dir}")

    # Confirm summary (only in interactive mode)
    if not assume_yes:
        print("\nSummary:")
        print(f"  Working dir : {wd}")
        print(f"  Reference   : {ref}")
        print(f"  Outputs to  : {outdir.resolve()}")
        print(f"  Threads     : {threads}")
        print(f"  Ploidy      : {ploidy}")
        print(f"  Filters     : QUAL>={min_qual}, DP>={min_dp}")
        print(f"  Samples     : {len(r1s)}")
        if run_dndscv:
            print(f"  dNdScv      : enabled (assembly={assembly}{', custom refdb' if (assembly=='custom' and refdb_rds) else ''})")
        if not prompt_yn("Proceed with processing?", default="y"):
            sys.exit("Aborted by user.")

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

        # Align → sort → index
        sh(f"bwa mem -t {threads} {qq(ref)} {qq(r1)} {qq(r2)} | "
           f"samtools sort -@ {threads} -o {qq(bam_sorted)}")
        sh(f"samtools index -@ {threads} {qq(bam_sorted)}")

        # Call → filter → compress/index
        sh(
            f"bcftools mpileup -O b -o {qq(raw_bcf)} -f {qq(ref)} {qq(bam_sorted)} && "
            f"bcftools call --ploidy {ploidy} -m -v -o - {qq(raw_bcf)} | "
            f"bcftools filter -e 'QUAL<{min_qual} || DP<{min_dp}' - | "
            f"bgzip -c > {qq(vcf_all)}"
        )
        tabix_index(vcf_all)

        # Split SNVs / INDELs
        sh(f"bcftools view -v snps   {qq(vcf_all)} -Oz -o {qq(vcf_snvs)}")
        tabix_index(vcf_snvs)
        sh(f"bcftools view -v indels {qq(vcf_all)} -Oz -o {qq(vcf_indels)}")
        tabix_index(vcf_indels)

        # Extract TSVs (intermediate) for dNdScv
        snv_tsv   = outdir/"vcf"/"snvs"/f"{base}_snvs.tsv"
        indel_tsv = outdir/"vcf"/"indels"/f"{base}_indels.tsv"
        bcftools_query_to_tsv(vcf_snvs, snv_tsv)
        bcftools_query_to_tsv(vcf_indels, indel_tsv)

        # Load and label
        def load_tbl(p: Path):
            if not p.exists() or p.stat().st_size == 0:
                return pd.DataFrame(columns=["sampleID","chr","pos","ref","mut"])
            df = pd.read_csv(p, sep="\t", header=None, names=["chr","pos","ref","mut"])
            df.insert(0, "sampleID", base)
            return df

        snv_rows.append(load_tbl(snv_tsv))
        indel_rows.append(load_tbl(indel_tsv))

    # Combine outputs
    ddir = outdir/"dndscv"
    snv_all   = pd.concat(snv_rows,   ignore_index=True) if snv_rows   else pd.DataFrame(columns=["sampleID","chr","pos","ref","mut"])
    indel_all = pd.concat(indel_rows, ignore_index=True) if indel_rows else pd.DataFrame(columns=["sampleID","chr","pos","ref","mut"])

    dndscv_all = pd.concat([snv_all, indel_all], ignore_index=True)

    # Light schema validation for dNdScv
    required_cols = ["sampleID","chr","pos","ref","mut"]
    if list(dndscv_all.columns) != required_cols:
        sys.exit(f"ERROR: dNdScv schema mismatch. Got columns {list(dndscv_all.columns)} expected {required_cols}")
    if not dndscv_all["pos"].apply(lambda x: str(x).isdigit() and int(x) > 0).all():
        sys.exit("ERROR: Invalid POS values detected.")
    valid_nt = set(list("ACGTN") + ["-"])
    if not dndscv_all["ref"].astype(str).str.upper().apply(lambda s: set(s) <= valid_nt).all():
        sys.exit("ERROR: Invalid REF symbols detected.")
    if not dndscv_all["mut"].astype(str).str.upper().apply(lambda s: set(s) <= valid_nt).all():
        sys.exit("ERROR: Invalid MUT symbols detected.")

    # ------------------ Write outputs (CSV) ------------------
    ddir.mkdir(parents=True, exist_ok=True)

    # Optional per-type CSVs (for QC/inspection)
    snv_csv   = ddir / "combined_snv_variants.csv"
    indel_csv = ddir / "combined_indels_variants.csv"
    snv_all.to_csv(snv_csv, index=False)         # sampleID,chr,pos,ref,mut
    indel_all.to_csv(indel_csv, index=False)

    # REQUIRED: single combined CSV for dNdScv
    dndscv_csv = ddir / "dndscv_input.csv"
    dndscv_all.to_csv(dndscv_csv, index=False)

    # Simple per-sample summary (CSV)
    summary_rows = []
    if not snv_all.empty:
        c = snv_all.groupby("sampleID").size().rename("n").reset_index()
        c["type"] = "SNV"
        summary_rows.append(c)
    if not indel_all.empty:
        c = indel_all.groupby("sampleID").size().rename("n").reset_index()
        c["type"] = "INDEL"
        summary_rows.append(c)
    if summary_rows:
        summary = pd.concat(summary_rows, ignore_index=True)
        summary.to_csv(outdir/"reports_summary.csv", index=False)

    # ------------------ Optional: run dNdScv R engine ------------------
    if run_dndscv:
        rscript = shutil.which("Rscript")
        if rscript is None:
            print("[WARN] Rscript not found in PATH. Skipping dNdScv run.")
        else:
            cmd = [
                rscript, "R/run_dndscv.R",
                "--input", str(dndscv_csv),
                "--assembly", assembly,
                "--outdir", str(ddir)
            ]
            if assembly == "custom" and refdb_rds:
                cmd.extend(["--refdb_rds", refdb_rds])
            print("[info] Running dNdScv via Rscript...")
            sh(" ".join(shlex.quote(c) for c in cmd))

    print("\n Done.")
    print(f"  dNdScv input : {dndscv_csv}")
    print(f"  (optional) SNVs   : {snv_csv}")
    print(f"  (optional) INDELs : {indel_csv}")
    if (outdir/'reports_summary.csv').exists():
        print(f"  Summary      : {outdir/'reports_summary.csv'}")
    if run_dndscv:
        print(f"  dNdScv (R)   : outputs in {ddir} (e.g., globaldnds.csv, sel_cv.csv, genesmuts.csv, dndscv_plots.pdf)")

# --------------------------- entrypoint ---------------------------

def build_parser():
    p = argparse.ArgumentParser(
        description="Interactive or non-interactive FASTQ→BAM→VCF→dNdScv-table generator (+ optional dNdScv R analysis)"
    )
    # mode
    p.add_argument("--non-interactive", action="store_true",
                   help="Run with provided flags only (no prompts).")
    p.add_argument("--yes", "--assume-yes", dest="assume_yes", action="store_true",
                   help="Assume 'yes' to confirmations (useful for interactive run with auto indexing).")

    # inputs / outputs
    p.add_argument("--ref", help="Path to reference FASTA")
    p.add_argument("--r1-source", default="*_R1.fastq.gz",
                   help="Directory OR glob for R1 FASTQs (e.g., '/data/cohort' or '/data/*/*_R1.fastq.gz').")
    p.add_argument("--outdir", default="results", help="Output directory (default: results)")

    # params
    p.add_argument("--threads", type=int, default=4, help="Threads (default: 4)")
    p.add_argument("--ploidy",  type=int, default=2, help="bcftools call ploidy (default: 2)")
    p.add_argument("--min-qual", type=int, default=20, help="Min QUAL (default: 20)")
    p.add_argument("--min-dp",   type=int, default=10, help="Min DP (default: 10)")
    p.add_argument("--auto-index", action="store_true",
                   help="Auto-generate bwa/samtools indexes if missing, without prompt.")

    # dNdScv (R) optional step
    p.add_argument("--run-dndscv", action="store_true",
                   help="After building CSV, run the R dNdScv engine.")
    p.add_argument("--assembly", choices=["hg19","hg38","auto","custom"], default="auto",
                   help="Assembly for dNdScv refdb (default: auto). If 'custom', also pass --refdb-rds.")
    p.add_argument("--refdb-rds", default=None,
                   help="Path to custom refdb RDS (used only if --assembly=custom).")

    return p

def main():
    parser = build_parser()
    args = parser.parse_args()

    wd = Path.cwd()

    if args.non_interactive:
        # strict: require ref in non-interactive
        if not args.ref:
            sys.exit("ERROR: --non-interactive requires --ref <fasta>")
        run_pipeline(
            ref=Path(args.ref).expanduser().resolve(),
            r1_glob_or_dir=args.r1_source,
            outdir=Path(args.outdir),
            threads=args.threads,
            ploidy=args.ploidy,
            min_qual=args.min_qual,
            min_dp=args.min_dp,
            auto_index=args.auto_index,
            assume_yes=True,  # never prompt in non-interactive
            run_dndscv=args.run_dndscv,
            assembly=args.assembly,
            refdb_rds=args.refdb_rds
        )
        return

    # -------- INTERACTIVE FLOW (prompts) --------
    print("Accepted format for FASTQs is *_R1.fastq.gz and *_R2.fastq.gz")
    ok_fmt = prompt_yn("Is your input FASTQ in the accepted format?", default="y")
    if not ok_fmt:
        print("Please transform your FASTQs and re-run.")
        sys.exit(0)

    # Reference (allow any path)
    ref_default = None
    fasta_candidates = list(wd.glob("*.fa")) + list(wd.glob("*.fasta")) + list(wd.glob("*.fna"))
    if len(fasta_candidates) == 1:
        ref_default = str(fasta_candidates[0])
    ref_input = prompt_str("Enter path to reference FASTA", default=ref_default)
    ref_path = Path(ref_input).expanduser()
    if not ref_path.is_absolute():
        ref_path = (wd / ref_path)
    ref_path = ref_path.resolve()

    # FASTQ location (dir or glob)
    fq_source = prompt_str("FASTQ location (directory OR glob)", default=args.r1_source)

    outdir = Path(prompt_str("Choose output directory", default=args.outdir))
    threads = prompt_int("Threads to use", default=args.threads, minval=1)
    ploidy  = prompt_int("Ploidy for bcftools call", default=args.ploidy, minval=1)

    run_pipeline(
        ref=ref_path,
        r1_glob_or_dir=fq_source,
        outdir=outdir,
        threads=threads,
        ploidy=ploidy,
        min_qual=args.min_qual,
        min_dp=args.min_dp,
        auto_index=args.auto_index,
        assume_yes=args.assume_yes,
        run_dndscv=args.run_dndscv,
        assembly=args.assembly,
        refdb_rds=args.refdb_rds
    )

if __name__ == "__main__":
    main()
