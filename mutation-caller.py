########### This software works with unzipped reference genome file

import subprocess
import pandas as  pd
import os
import glob
import os.path
from subprocess import call
# import rpy2
# import rpy2.robjects as robjects
# import matplotlib.pyplot as plt
# import altair as alt

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


output_folder = ['sam', 'bam', 'bcf', 'vcf', 'dndscv']
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
        print (f'base name is {base}')
        fq1= os.path.join (wd, f'{base}' + '_R1.fastq.gz')
        fq2= os.path.join (wd, f'{base}' + '_R2.fastq.gz')
        sam = os.path.join (wd, 'results', 'sam', f'{base}' + '.aligned.sam')
        bam = os.path.join (wd, 'results', 'bam', f'{base}' + '.aligned.bam')
        sorted_bam = os.path.join (wd, 'results', 'bam', f'{base}' + '.aligned.sorted.bam')
        raw_bcf = os.path.join (wd, 'results', 'bcf', f'{base}' + '_raw.bcf')
        variants = os.path.join (wd, 'results', 'vcf', f'{base}' + '_variants.vcf')
        final_variants = os.path.join (wd, 'results', 'vcf', f'{base}' + '_final_variants.vcf')
        indexed_variants = os.path.join (wd, 'results', 'vcf', f'{base}' + '_final_variants.vcf.gz')
        indels = os.path.join (wd, 'results', 'vcf', 'indels', f'{base}' + '_indels.vcf.gz')
        snvs = os.path.join (wd, 'results', 'vcf', 'snvs', f'{base}' + '_snvs.vcf.gz')
        


        




        
        for command in (f"bwa mem {ref_genome_file_path} {fq1} {fq2} > {sam}", 
                        f"samtools view -S -b {sam} > {bam}",
                        f"samtools sort -o {sorted_bam} {bam}",
                        f"samtools index {sorted_bam}", 
                        f"bcftools mpileup -O b -o {raw_bcf} -f {ref_genome_file_path} {sorted_bam}",
                        f"bcftools call --ploidy 1 -m -v -o {variants} {raw_bcf}",
                        f"vcfutils.pl varFilter {variants} > {final_variants}",
                        f"bgzip {final_variants}",
                        f"bcftools index {indexed_variants}",
                        f"bcftools view -v snps {indexed_variants} > -Oz -o {snvs}", # could we call variants from variants, not final_variants?
                        f"bcftools view -v indels {indexed_variants} > -Oz -o {indels}"):
              

              
            call(command, shell=True) 
        # Extracting indels variants

        print ("Contructing indels mutation table...") 

        
        os.chdir(os.path.join (wd, 'results', 'vcf', 'indels'))

        file = glob.glob("*vcf.gz")
        all_samples=[]
        for item in file:
                pattern = item.find("_")

                sampleID = item[0:pattern]
                df_variant=pd.read_csv(item,sep='\t',comment='#',usecols=[0,1,3,4],names=['CHR','POS','REF','ALT'])
                sampleids=[sampleID for i in range(df_variant.shape[0])]
                df_variant.insert(0,'SampleID',sampleids)
                all_samples.append(df_variant)
            
                print(f"Generating indels mutation table from {sampleID}")

            
       
        df_final=pd.concat(all_samples)
        df_final.to_csv('combined_indels_variants.csv',index=False)
        dndscv_path = os.path.join (wd, 'results', 'dndscv')

        status = subprocess.run (["cp", "combined_indels_variants.csv", dndscv_path], capture_output=True)

            
        # Extracting snvs variants

        print ("Contructing snvs mutation table...") 

        #wd = os.getcwd()
        os.chdir(os.path.join (wd, 'results', 'vcf', 'snvs'))

        file = glob.glob("*vcf.gz")
        all_samples=[]
        for item in file:
                pattern = item.find("_")

                sampleID = item[0:pattern]
                df_variant=pd.read_csv(item,sep='\t',comment='#',usecols=[0,1,3,4],names=['CHR','POS','REF','ALT'])
                sampleids=[sampleID for i in range(df_variant.shape[0])]
                df_variant.insert(0,'SampleID',sampleids)
                all_samples.append(df_variant)
            
                print(f"Generating snvs mutation table from {sampleID}")

            
        df_final=pd.concat(all_samples)
        df_final.to_csv('combined_snv_variants.csv',index=False, header=None)

        status = subprocess.run (["cp", "combined_snv_variants.csv", dndscv_path], capture_output=True)



        # Generate the dndsvc table
        
        os.chdir(os.path.join (wd, 'results', 'dndscv'))

              
        for command in ("cat combined_indels_variants.csv combined_snv_variants.csv > dndscv.csv",):
            call(command, shell=True) 

