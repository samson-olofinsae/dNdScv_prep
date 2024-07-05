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
