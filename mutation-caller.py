########### This software works with unzipped reference genome file

import subprocess
import pandas as  pd
import os
import glob
import os.path
from subprocess import call
import rpy2
import rpy2.robjects as robjects
import matplotlib.pyplot as plt
import altair as alt

# Initial User interface

print ("Accepted format for fastq files is _(R1,R2).fastq.gz")
q = input ('Is your input fastq file in this format? enter y for Yes or n for No: ')
if q == 'y':
    ref_genome_file = input ("Accurately enter the name of your ref genome file: ")
elif q == 'n':
    print ('Please tranform your fastq file to the recommended format and restart the programme')
    quit()

wd = os.getcwd()

print (f'Your working directory: {wd}')

ref_genome_file_path = os.path.join(wd, ref_genome_file)

print (f'Your ref genome file path: {ref_genome_file_path}')


fastq_file_path = os.path.join(wd, '*.fastq.gz') 

print (f'Your fastq file path is {fastq_file_path}')


# Indexing the reference genome

indexing = subprocess.run (["bwa", "index", ref_genome_file_path], capture_output=True)


# Creating the outut folders and subfolders

os.mkdir ('results')

os.chdir('results')


output_folder = ['sam', 'bam', 'bcf', 'vcf']
for folder in output_folder:
    os.mkdir(folder)

os.chdir('vcf')
variant_folder = ['indels', 'snvs']
for variant in variant_folder:
    os.mkdir(variant)


os.chdir(wd)



for fq1 in os.path.join (wd, '*_R1.fastq.gz'):
    
    file = glob.glob("*R1.fastq.gz")
    for item in file:
        pattern = item.find("_")
   
        base = item[0:pattern]
        base_R1 = base + '_R1.fastq.gz'
        base_R2 = base + '_R2.fastq.gz'
        base_sam = base + '.aligned.sam'
        base_bam = base + '.aligned.bam'
        base_sorted_bam = base + '.aligned.sorted.bam'
        base_bcf = base + '_raw.bcf'
        base_variant = base + '_variants.vcf'
        base_final_variant = base + '_final_variants.vcf'
        base_indel = base + '_indels.vcf.gz'
        base_snv = base + '_snvs.vcf.gz'

        fq1 = os.path.join (wd, base_R1)
        fq2 = os.path.join (wd, base_R2)
        sam = os.path.join (wd, 'results', 'sam', base_sam)
        bam = os.path.join (wd,'results', 'bam', base_bam)
        sorted_bam = os.path.join (wd,'results', 'bam', base_sorted_bam)
        raw_bcf = os.path.join (wd,'results','bcf', base_bcf)
        variants = os.path.join (wd, 'results', 'vcf', base_variant)
        final_variants = os.path.join (wd, 'results', 'vcf', base_final_variant)
        indexed_variants = os.path.join (wd, 'results', 'vcf', base_final_variant)
        indels = os.path.join (wd, 'results', 'vcf', 'indels', base_indel)
        snvs = os.path.join (wd, 'results', 'vcf', 'snvs', base_snv)



        for command in (f"bwa mem {ref_genome_file_path} {fq1} {fq2} > {sam}", 
                        f"samtools view -S -b {sam} > {bam}",
                        f"samtools sort -o {sorted_bam} {bam}",
                        f"samtools index {sorted_bam}", 
                        f"bcftools mpileup -O b -o {raw_bcf} -f {ref_genome_file_path} {sorted_bam}",
                        f"bcftools call --ploidy 1 -m -v -o {variants} {raw_bcf}",
                        f"bgzip {final_variants}",
                        f"bcftools index {indexed_variants}",
                        f"bcftools view -v snps {indexed_variants} > {snvs}",
                        f"bcftools view -v indels {indexed_variants} > {indels}"):
              

              
            call(command, shell=True) 
