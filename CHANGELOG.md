# Changelog

All notable changes to this project will be documented in this file.

---

## [v1.2.2] - 2025-11-04

### Documentation & UX Improvements
- Added a **Beginner’s Note** with clear step-by-step guidance on how to enter file paths during interactive runs.
- Clarified that `user_ref/` and `user_fastq/` directories can (and should) be created manually before running the pipeline on real data.
- Updated **real data command** to use the placeholder `<your_reference_file>.fa`, reducing confusion for new users:
  ```bash
  python3 mutation-caller.py --ref user_ref/<your_reference_file>.fa --r1-source "user_fastq/*_R1.fastq.gz" --outdir user_results --threads 8 --ploidy 2 --min-qual 20 --min-dp 10 --auto-index --yes
  ```
  > Replace `<your_reference_file>` with your actual FASTA filename (e.g. `hg19.fa`, `GRCh38.fa`).
- Introduced a **Troubleshooting** section with tips for common beginner errors:
  - Using `ls user_ref/` and `ls user_fastq/` to confirm filenames.
  - Handling the “Invalid MUT symbols detected” error caused by `ALT="*"` (spanning deletions).
  - Optional bcftools filter fix:  
    ```bash
    bcftools view -e 'ALT="*"' input.vcf.gz -Oz -o filtered.vcf.gz
    ```
- Emphasized that generated index files (`.bwt`, `.fai`, `.pac`, `.ann`, `.amb`, `.sa`) remain **gitignored**.

> **Purpose:** Enhance accessibility and error resilience for new users running the tool on real datasets.

---

## [v1.2.1] - 2025-11-03

### Minor Documentation Fixes
- Corrected author name and citation formatting in **CITATION.cff**.
- Added “Latest” badge and link to v1.2.0 release in **README.md**.
- Improved alignment of markdown sections for GitHub rendering.
- No functional code changes; documentation-only release.

> **Purpose:** Visual and citation consistency improvements following the v1.2.0 release.

---

## [v1.2.0] - 2025-11-04

### Documentation Updates
- Expanded **README.md** with clearer guidance on running the pipeline interactively.
- Added detailed explanation for using `--auto-index --yes` when FASTA index files are missing.
- Clarified the difference between:
  - **Interactive mode** (`python3 mutation-caller.py`)
  - **Automated indexing mode** (`--auto-index --yes`)
  - **Non-interactive mode** (`--non-interactive`)
- Included explicit example paths and prompt responses for the demo dataset.
- Improved markdown structure for easier readability in GitHub.
- Added notes about `.bwt`, `.fai`, `.pac`, `.ann`, `.amb`, and `.sa` being **gitignored**.

> **Purpose:** Make onboarding simpler for new users and reproducibility reviewers, aligning documentation with the current behavior of the CLI.

---

## [v1.1.0] - 2025-10-20

### Features
- Introduced uniform prompt interface for CLI and interactive mode.
- Added `--yes` and `--auto-index` flags for simplified automation.
- Improved input validation for reference FASTA and FASTQ patterns.
- Enhanced summary output for reproducibility logs.

---

## [v1.0.0] - 2025-09-15

### Initial Release
- First public version of **dNdScv_prep**.
- Implements FASTQ → BAM → VCF → dNdScv input pipeline.
- Bundled demo FASTQs and reference FASTA for quick reproducibility.
- Released under MIT License.
