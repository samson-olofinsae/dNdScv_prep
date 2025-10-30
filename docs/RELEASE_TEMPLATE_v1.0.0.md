# dNdScv_prep v1.0.0 - Release Notes (TEMPLATE)

> **Status:** Template for your first tagged release. Copy–paste into GitHub → Releases when publishing `v1.0.0`.

## Highlights
- **End-to-end prep** from paired-end FASTQs → BAM → VCF → **single combined CSV** (`dndscv_input.csv`) for the **dNdScv** R package.
- **Container/conda friendly**: minimal runtime deps (bwa, samtools, bcftools, bgzip, tabix, Python + pandas). Works on Linux, WSL2, and HPC.
- **Demo-ready**: includes toy `demo_ref/` and `demo_fastq/` to validate the full path.
- **Reproducible outputs**: deterministic folder structure, simple logs, and schema checks.

## What’s in this release
- `mutation-caller.py` (CLI) with **interactive** and **non-interactive** modes.
- Generates:
  - `results_*/dndscv/dndscv_input.csv` (**required**, merged SNVs+INDELs)
  - `results_*/dndscv/combined_snv_variants.csv` (optional QC)
  - `results_*/dndscv/combined_indels_variants.csv` (optional QC)
  - `results_*/reports_summary.csv` (per-sample variant counts)
- `environment.yml` for a pinned conda environment.
- `demo_ref/` (tiny toy FASTA) and `demo_fastq/` (A/B/C samples) for quick validation.

## Requirements
- Linux (or WSL2)
- Tools on PATH: `bwa`, `samtools`, `bcftools`, `bgzip`, `tabix`
- Python 3.8+ with `pandas` (or use `conda env create -f environment.yml`)

## Quickstart (Demo)
```bash
# 1) Clone
git clone https://github.com/samson-olofinsae/dNdScv_prep.git
cd dNdScv_prep

# 2) Create env (recommended)
conda env create -f environment.yml
conda activate dndscv-prep

# 3) Run (non-interactive demo)
python3 mutation-caller.py \
  --non-interactive \
  --ref demo_ref/demo.fa \
  --r1-source "demo_fastq/*_R1.fastq.gz" \
  --outdir results_demo \
  --threads 4 --ploidy 2 --min-qual 20 --min-dp 10 \
  --auto-index
```

**Outputs:** see `results_demo/dndscv/dndscv_input.csv`.

## Example artifacts (attach on the Release page)
- `results_demo/dndscv/dndscv_input.csv`
- `results_demo/reports_summary.csv`
- (optional) `results_demo/dndscv/combined_snv_variants.csv`
- (optional) `results_demo/dndscv/combined_indels_variants.csv`

> Tip: keeping demo outputs attached to the release helps reviewers validate quickly.

## Checksums (fill after attaching artifacts)
- `dndscv_input.csv` — `SHA256:48e6d6ed6bafc696fff1cd2d24695ab9128d2aadb5454d76b97579531aa2ee53`
- `reports_summary.csv` — `SHA256:3add603b662a85e7acc87ff5172b172df1631f32c5573b3af48b3664235a62d3`
## Reproducibility recipe
- OS: `<Ubuntu 22.04 / WSL2 Ubuntu>`
- Conda: `<23.x>`
- Tool versions: `bwa 0.7.x`, `samtools 1.1x`, `bcftools 1.1x`, `tabix htslib 1.1x`
- Command: (see Quickstart)
- Expected row counts (demo): provide counts from `reports_summary.csv`

## Known issues / Notes
- The **toy demo reference** uses `demo1` as a chromosome. It’s intentionally not human.
- For real human analyses, ensure **reference assembly, chromosome naming, and dNdScv refdb** are consistent (e.g., GRCh38 + `chr1` vs GRCh37 + `1`).
- This project **does not** re-implement the dNdScv engine; it **prepares** the input table.

## Acknowledgements
This tool **complements** the excellent **dNdScv** package by **Dr. Iñigo Martincorena** and colleagues (Wellcome Sanger Institute). Please cite their work when using this pipeline.

## Cite
If you use this pipeline, please cite:
- *dNdScv*: Martincorena I, et al. (provide canonical citation here).

And this repository:
```
Olofinsae S. dNdScv_prep: a reproducible pipeline to generate merged SNV+INDEL tables for dNdScv. v1.0.0. GitHub. https://github.com/samson-olofinsae/dNdScv_prep
```

---

### Changelog
- Initial public release `v1.0.0`.
- Stable CLI, demo data, and environment spec.
