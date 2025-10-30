# dNdScv_prep

### FASTQ → BAM → VCF → dNdScv Input (Container-Ready Workflow)

---

## Overview

**dNdScv_prep** is a lightweight, reproducible pipeline that automates generation of the `dNdScv_input.csv` file required by the [**dNdScv**](https://github.com/im3sanger/dndscv) R package developed by **Dr Inigo Martincorena** (Wellcome Sanger Institute).

Our tool does **not** re-implement dNdScv.  
Instead, it **complements** it - automating all upstream data-processing steps and simplifying the workflow for postdocs, early-career researchers, and clinicians who wish to run *dNdScv* reproducibly on their own datasets.

> **Core Workflow:**  
> FASTQ → BAM → VCF → Combined Variants → dNdScv Input CSV

### What does dNdScv do?
The **dNdScv** R package estimates the ratio of non-synonymous to synonymous substitutions (dN/dS) across genes, enabling the detection of positive selection in somatic mutations. It’s a powerful tool for identifying **driver genes** in cancer and somatic evolution studies.

This pipeline complements dNdScv by generating a clean, ready-to-use input file directly from raw FASTQ data - closing the gap between raw sequencing data and selection analysis.

---

## Clone the Repository

```bash
git clone https://github.com/samson-olofinsae/dNdScv_prep.git
cd dNdScv_prep
```

---

## Quick Start

### 1. Create environment
```bash
conda env create -f environment.yml
conda activate dndscv-prep
```

### 2 Run demo (tracked in repo)
```bash
python3 mutation-caller.py --ref demo_ref/demo.fa --r1-source "demo_fastq/*_R1.fastq.gz" --outdir results_demo --threads 4 --ploidy 2 --min-qual 20 --min-dp 10 --auto-index
```

This produces:
```
results_demo/
├── dndscv/
│   ├── dndscv_input.csv
│   ├── combined_snv_variants.csv
│   └── combined_indels_variants.csv
└── reports_summary.csv
```

### 3. Run on real data (ignored by Git)
```bash
python3 mutation-caller.py --ref user_ref/hg19.fa --r1-source "user_fastq/*_R1.fastq.gz" --outdir user_results --threads 8 --ploidy 2 --min-qual 20 --min-dp 10 --auto-index
```

---

## Repository Structure

| Folder | Purpose | Tracked |
|---------|----------|---------|
| `demo_ref/` | Tiny toy reference FASTA | Tracked |
| `demo_fastq/` | Small paired FASTQ demo files | Tracked |
| `results_demo/` | Example output for reviewers | Tracked |
| `R/demo/` | Example CSVs & demo outputs | Tracked |
| `user_ref/`, `user_fastq/`, `user_results/` | Real data (private) | Ignored |

This layout allows reproducible demo runs while keeping real genomic data out of version control.

---

## Parameters

| Flag | Description | Default |
|------|--------------|----------|
| `--ref` | Reference FASTA | *required* |
| `--r1-source` | Directory or glob for R1 FASTQs | `*_R1.fastq.gz` |
| `--outdir` | Output directory | `results` |
| `--threads` | Number of threads | `4` |
| `--ploidy` | Ploidy for variant calling | `2` |
| `--min-qual` | Minimum QUAL threshold | `20` |
| `--min-dp` | Minimum depth threshold | `10` |
| `--auto-index` | Auto-generate BWA / FASTA indexes | `false` |
| `--non-interactive` | Run without prompts | `false` |

---

## Containerisation

The `environment.yml` file specifies all dependencies (`bwa`, `samtools`, `bcftools`, `pandas`, `python >= 3.9`).  
This environment can be exported to **Singularity** or **Docker** images for cross-system reproducibility.

Example:
```bash
conda env export --no-builds > environment.yml
```

---

## Integration with dNdScv

After the pipeline generates `dndscv_input.csv`, it can be directly consumed by the R engine:

```r
library(dndscv)
m <- read.csv("results_demo/dndscv/dndscv_input.csv", stringsAsFactors = FALSE)
dndsout <- dndscv(m)
dndsout$sel_cv
```

For reproducible testing, a demo input is also provided at:
```
R/demo/dndscv_demo_input.csv
```

---

## Acknowledgements

Our collaborator and **dNdScv** creator  
**Dr Inigo Martincorena**, PhD - Group Leader, Somatic Evolution Group, Wellcome Sanger Institute  

Our mentors and collaborators  
- **Prof David Wedge**  -Cancer Research UK Manchester Centre  
- **Prof Rosalind Eeles** - Institute of Cancer Research, UK  
- **Prof Daniel Brewer** - University of East Anglia  
- **Prof Colin Cooper** - University of East Anglia  

> This pipeline is a collaborative complement to *dNdScv*, designed to promote openness, reproducibility, and accessibility in somatic mutation research.

---

## License

This project is released under the **MIT License** — see [LICENSE](./LICENSE) for details.
