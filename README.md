This software aims to contribute to the ongoing scientific pursuit of detecting of the changes in the DNA sequence of genes that cause cells to become cancerous. 
Sequence analysis of tumour tissues for driver genes may not only help to plan treatment and cause cancer cells from growing, but can also help guide novel therapy (Ostroverkhova D, et al, 2023). 
Consequently, the development of suites of bioinformatics pipeline toward detecting novel cancer mutation drivers is urgent.

Variant calling is an essential downstream analysis in Next Generation Sequencing analysis. Filtering variants in sequenced DNA samples is an essential process towards determining variations in the samples relative to the reference sequence. While typically written in bash script, the pipeline script in this repo is coded in Python (filename: mutation-caller.py), which is the chosen programming language for this pipeline because of its powerful functionalities that enable developers to track bugs, unit-test codes, build control and automation. The pipeline employs an inital software-like user friendly interface which checks that the right file naming convention compatible with the pipeline is used by the user. Acceptable fastq file naming is : sample_(R1,R2)* where R1 and R2 denotes forward and reverse read, respectively. An unzipped reference genome file is also recommended.

The pipeline uses the following python modules for its processing: pandas, os, glob, os.path and subprocess ( importing call). The pipeline is able to process multiple samples (HPC server recommended).

The pipeline generates outputs in the user's working directory using the following folder and file structure (aslo attached as filename: Folder_structure.txt):

results/ ├── bam │   ├── sample.aligned.bam │   ├── sample.aligned.sorted.bam │   └── sample.aligned.sorted.bam.bai ├── bcf │   └── sample_raw.bcf ├── sam │   └── sample.aligned.sam └── vcf ├── sample_variants.vcf ├── indels │   └── samplle_indels.vcf.gz └── snvs └── sample_snvs.vcf.gz

From input fastq files, the pipeline generates VCF files, and also filters the files to seperately generate indels and snv VCFs. VCF file filtering is an essential step added for researh purpose. For instance, concatenation of indels and snvs variants from tumour samples may be used as inputs in the pipeline for generating mutation table used in somatic cancer driver analysis that employs Non Synonymous:Synonymous mutation ratio (Martincorena, et al., 2017)

Reference

Martincorena I, Raine KM, Gerstung M, Dawson KJ, Haase K, Van Loo P, Davies H, Stratton MR, Campbell PJ. Universal Patterns of Selection in Cancer and Somatic Tissues. Cell. 2017 Nov 16;171(5):1029-1041.e21. doi: 10.1016/j.cell.2017.09.042. Epub 2017 Oct 19. Erratum in: Cell. 2018 Jun 14;173(7):1823. PMID: 29056346; PMCID: PMC5720395.

Ostroverkhova D, Przytycka TM, Panchenko AR. Cancer driver mutations: predictions and reality. Trends Mol Med. 2023 Jul;29(7):554-566. doi: 10.1016/j.molmed.2023.03.007. Epub 2023 Apr 17. PMID: 37076339.
