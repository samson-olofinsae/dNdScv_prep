# Changelog

All notable changes to this project will be documented in this file.

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
