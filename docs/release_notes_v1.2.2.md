# dNdScv_prep v1.2.2 - Beginner’s Note, Folder Setup & Troubleshooting

> **Release date:** 2025-11-04  
> **Type:** Documentation & UX Update  
> **Previous version:** v1.2.1  
> **Author:** Samson A. Olofinsae  
> **License:** MIT

---

## Highlights

- **Beginner-friendly onboarding:** clearer guidance for users new to command-line pipelines.  
- **Improved folder structure:** standardized directories (`user_ref/`, `user_fastq/`, `user_results/`) for real-data runs.  
- **Interactive clarity:** added precise examples on entering file paths interactively.  
- **Troubleshooting guide:** now includes common errors, file validation checks, and example bcftools filters.  
- **Refined documentation:** consistent markdown style, citation fixes, and clarified flag usage (`--auto-index`, `--yes`, etc.).  

---

## What’s in this release

### Documentation & UX Improvements
- Added a **Beginner’s Note** explaining how to prepare input directories and verify filenames.  
- Introduced standard directory setup:
  ```bash
  mkdir -p user_ref user_fastq user_results
  ```
  Place your reference FASTA and FASTQs there before running the pipeline.
- Updated **real-data command** example for clarity:
  ```bash
  python3 mutation-caller.py --ref user_ref/<your_reference_file>.fa --r1-source "user_fastq/*_R1.fastq.gz" --outdir user_results --threads 8 --ploidy 2 --min-qual 20 --min-dp 10 --auto-index --yes
  ```
  Replace `<your_reference_file>` with your actual FASTA filename (e.g. `hg19.fa`, `GRCh38.fa`).

### Troubleshooting section
Common beginner errors and their quick fixes:
- **Path errors** → use:
  ```bash
  ls user_ref/
  ls user_fastq/
  ```
  to confirm filenames.
- **Invalid MUT symbols detected** → caused by `ALT="*"` spanning deletions.  
  Automatically handled in v1.2.2, but manual filter if needed:
  ```bash
  bcftools view -e 'ALT="*"' input.vcf.gz -Oz -o filtered.vcf.gz
  ```
- **Indexing issues** → rerun with `--auto-index --yes`.

> **Purpose:** make the pipeline intuitive and forgiving for first-time users, ensuring reproducibility even on minimal setups.

---

## Quickstart (Demo)

```bash
# 1) Clone
git clone https://github.com/samson-olofinsae/dNdScv_prep.git
cd dNdScv_prep

# 2) Create env (recommended)
conda env create -f environment.yml
conda activate dndscv-prep

# 3) Run demo
python3 mutation-caller.py --non-interactive --ref demo_ref/demo.fa --r1-source "demo_fastq/*_R1.fastq.gz" --outdir results_demo --threads 4 --ploidy 2 --min-qual 20 --min-dp 10 --auto-index
```

**Outputs:**  
`results_demo/dndscv/dndscv_input.csv` (ready for `dNdScv` R package).  
Additional QC outputs include `reports_summary.csv`, `combined_snv_variants.csv`, and `combined_indels_variants.csv`.

---

## Reproducibility Recipe

| Component | Version / Notes |
|------------|----------------|
| OS | Ubuntu 22.04 / WSL2 Ubuntu |
| Conda | 23.x |
| Tools | bwa 0.7.x • samtools 1.1x • bcftools 1.1x • tabix (HTSlib 1.1x) |
| Command | See “Quickstart” section |
| Demo rows | Validate using `reports_summary.csv` |

---

## Known Issues / Notes

- The demo reference uses a toy chromosome `demo1`; for real analyses, ensure reference assembly and dNdScv database consistency (e.g. GRCh38 + `chr1`).  
- This tool **prepares** the dNdScv input table; it does **not** re-implement the dNdScv algorithm.  
- Generated index files (`.bwt`, `.fai`, `.pac`, `.ann`, `.amb`, `.sa`) remain **gitignored**.

---

## Acknowledgements

This tool complements the excellent **dNdScv** R package by  
**Dr Inigo Martincorena** and colleagues (Wellcome Sanger Institute).  
Please cite their work when using this pipeline.

---

## Citation

If you use this pipeline, please cite:

- *dNdScv*: Martincorena I et al. *Genome-wide determinants of somatic mutation rates in cancer.* Nature (2017).  
- *This repository*:
  ```
  Olofinsae S. dNdScv_prep: a reproducible pipeline to generate merged SNV + INDEL tables for dNdScv. v1.2.2. GitHub. https://github.com/samson-olofinsae/dNdScv_prep
  ```

---

## Changelog Summary

| Version | Date | Key Updates |
|----------|------|-------------|
| **v1.2.2** | 2025-11-04 | Beginner’s Note, folder setup, troubleshooting, improved clarity |
| **v1.2.1** | 2025-11-03 | Citation & visual fixes |
| **v1.2.0** | 2025-11-04 | Expanded README, auto-index explanation |
| **v1.1.0** | 2025-10-20 | New CLI flags, validation improvements |
| **v1.0.0** | 2025-09-15 | Initial public release |

---

**Maintainer:** Dr Samson A. Olofinsae  
**Repository:** [https://github.com/samson-olofinsae/dNdScv_prep](https://github.com/samson-olofinsae/dNdScv_prep)  
**License:** MIT  
**Contact:** [GitHub Issues](https://github.com/samson-olofinsae/dNdScv_prep/issues)

---

> _“Reproducibility begins with clarity — this release builds bridges for first-time users and reviewers alike.”_
