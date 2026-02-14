# dNdScv_prep

[![Latest Release](https://img.shields.io/badge/release-v1.2.2-blue.svg)](https://github.com/samson-olofinsae/dNdScv_prep/releases/tag/v1.2.2)
[**Latest (v1.2.2) release notes →**](https://github.com/samson-olofinsae/dNdScv_prep/releases/tag/v1.2.2)

### FASTQ → BAM → VCF → dNdScv Input (Container‑Ready Workflow)

---

## Overview

**dNdScv_prep** is a lightweight, reproducible pipeline that automates generation of the `dNdScv_input.csv` file required by the [**dNdScv**](https://github.com/im3sanger/dndscv) R package developed by **Dr. Inigo Martincorena** (Wellcome Sanger Institute).

Our tool does **not** re‑implement dNdScv.
Instead, it **complements** it - automating upstream data‑processing steps and simplifying the workflow for postdocs, early‑career researchers, and clinicians who wish to run *dNdScv* reproducibly on their own datasets.

> **Core Workflow:**
> FASTQ → BAM → VCF → Combined Variants → **dNdScv Input CSV**

---

## Clone the Repository

```bash
git clone https://github.com/samson-olofinsae/dNdScv_prep.git
cd dNdScv_prep
```

---

## Quick Start

### 1) Create environment
```bash
conda env create -f environment.yml
conda activate dndscv-prep
```

### 2) Run demo (tracked in repo)
```bash
python3 mutation-caller.py --ref demo_ref/demo.fa --r1-source "demo_fastq/*_R1.fastq.gz" --outdir results_demo --threads 4 --ploidy 2 --min-qual 20 --min-dp 10 --auto-index --yes
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

### 3) Run on real data (ignored by Git)
```bash
python3 mutation-caller.py --ref user_ref/<your_reference_file>.fa --r1-source "user_fastq/*_R1.fastq.gz" --outdir user_results --threads 8 --ploidy 2 --min-qual 20 --min-dp 10 --auto-index --yes
```

> **Note:** Replace `<your_reference_file>` with your actual FASTA name (e.g. `hg19.fa`, `GRCh38.fa`, or `custom_ref.fa`).

---

## Preparing Your Own Data Folders

Before running the pipeline on your own data, make sure the following directories exist:

```
user_ref/       # store your reference FASTA file(s) here
user_fastq/     # store your paired FASTQ files here
```

If they don’t exist yet, simply create them manually:

```bash
mkdir -p user_ref user_fastq
```

Then copy your reference and FASTQ files into those folders.  
For example:

```bash
cp /path/to/your_reference.fa user_ref/
cp /path/to/sample*_R*.fastq.gz user_fastq/
```

Once your files are in place, you can run:

```bash
python3 mutation-caller.py --ref user_ref/<your_reference_file>.fa --r1-source "user_fastq/*_R1.fastq.gz" --outdir user_results  --threads 8 --ploidy 2 --min-qual 20 --min-dp 10 --auto-index --yes
```

> **Note:** Replace `<your_reference_file>` with your actual FASTA name  
> (e.g. `hg19.fa`, `GRCh38.fa`, or `custom_ref.fa`).  
> The `user_results/` folder will be created automatically when the pipeline runs.

---

## Interactive Mode (Recommended)

You can also run the pipeline **interactively** with no flags. This is helpful for first‑time users and for documentation screenshots.

```bash
python3 mutation-caller.py
```

When launched, the terminal will display prompts like:

```
Accepted format for FASTQs is *_R1.fastq.gz and *_R2.fastq.gz
TIP: Reference FASTA path example: ./demo_ref/demo.fa
TIP: FASTQ glob (R1) example:      ./demo_fastq/*_R1.fastq.gz

Enter path to reference FASTA [user_ref/<your_reference_file>.fa]: ./user_ref/<your_reference_file>.fa
FASTQ location (directory OR glob) [user_fastq/*_R1.fastq.gz]: ./user_fastq/*_R1.fastq.gz
Choose output directory [user_results]: user_results
Threads to use [4]: 4
Ploidy for bcftools call [2]: 2

Summary:
  Working dir : /path/to/dNdScv_prep
  Reference   : /path/to/dNdScv_prep/user_ref/<your_reference_file>.fa
  Outputs to  : /path/to/dNdScv_prep/user_results
  Threads     : 4
  Ploidy      : 2
  Filters     : QUAL>=20, DP>=10
  Samples     : 3
Proceed with processing? [Y/n]: y
```

> **Tip:** The pipeline expects paired files named like `sampleA_R1.fastq.gz` and `sampleA_R2.fastq.gz`.

---

## If your FASTA isn’t indexed (auto‑index)

If the required index files for your reference FASTA (e.g. `demo.fa.bwt`, `demo.fa.fai`, `demo.fa.pac`, `demo.fa.ann`, `demo.fa.amb`, `demo.fa.sa`) are **missing**, you can auto‑generate them while keeping the rest of the run interactive:

```bash
python3 mutation-caller.py --auto-index --yes
```

- `--auto-index` creates the **BWA** and **FASTA** indices on the fly, if absent.
- `--yes` auto‑accepts the yes/no confirmations so indexing proceeds, while still showing the normal prompts for paths and parameters.

> **Note:** Generated index files can be large and are intentionally **gitignored**.

---

## Parameters

| Flag | Description | Default |
|---|---|---|
| `--ref` | Path to reference FASTA | *required for `--non-interactive`* |
| `--r1-source` | Directory or glob for R1 FASTQs | `*_R1.fastq.gz` |
| `--outdir` | Output directory | `results` |
| `--threads` | Number of threads | `4` |
| `--ploidy` | Ploidy for variant calling | `2` |
| `--min-qual` | Minimum QUAL threshold | `20` |
| `--min-dp` | Minimum depth threshold | `10` |
| `--auto-index` | Auto-generate BWA / FASTA indexes if missing | `false` |
| `--yes` | Auto‑accept Y/N confirmations (still prints prompts) | `false` |
| `--non-interactive` | Fully silent mode (for CI/CD) | `false` |

---

## Containerisation

The `environment.yml` file specifies all dependencies (`bwa`, `samtools`, `bcftools`, `pandas`, `python >= 3.9`).
This environment can be exported to **Singularity** or **Docker** images for cross‑system reproducibility.

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
**Dr. Inigo Martincorena**, PhD - Group Leader, Somatic Evolution Group, Wellcome Sanger Institute

Our mentors and collaborators
- **Prof David Wedge** - Cancer Research UK Manchester Centre
- **Prof Rosalind Eeles** - Institute of Cancer Research, UK
- **Prof Daniel Brewer** - University of East Anglia
- **Prof Colin Cooper** - University of East Anglia

> This pipeline is a collaborative complement to *dNdScv*, designed to promote openness, reproducibility, and accessibility in somatic mutation research.

---

## License

This project is released under the **MIT License** - see [LICENSE](./LICENSE) for details.

---

_Last updated: 2025-11-04_




## Companion R analysis module (QC + visualisation + Shiny)

For downstream QC plots, interactive dashboards, and **VCF → `dndscv_input.csv`** conversion (without Bioconductor IO issues), see:

**dNdScv_R_analysis**  
(VCF → QC tables/plots → Shiny → dNdScv input)  
GitHub: https://github.com/samson-olofinsae/dNdScv_R_analysis

Together, these two repositories provide an end-to-end reproducible somatic-mutation workflow from FASTQs to selection analysis readiness.


_Last updated: 2026-12-14_