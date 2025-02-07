'''
Goal: 
This code aims to automatize bwa processing (bwa.py) and is calling a second file for analysing and extracting matching positions. 

Instructions:
The command line should be : python3 path/bwa.py path_to_data path_to_index/index.fasta Donor Recipient. 
For example: python3 bwa.py data/GalK/Rel index/K12_REL606.fasta Rel606 K12


(bwa.py) Juliette Bellengier
(ecoli.py) Abdelmajid Omarjee and Juliette Bellengier 

Contact  juliette.bellengier@inserm.fr if any question/comment
'''


import sys 
import os
import subprocess
import re #for regular expression in python 
import glob



#Associating variables with arguments
seq_path = sys.argv[1]
index_file_with_the_path = sys.argv[2]
recipient= sys.argv[3]
donor= sys.argv[4]
output_path = input("What is the output path", )


#Create a file 'fastaqsam' with the names of all the fastq files
fastaqsam= []
for file in os.listdir(seq_path):
    fastaqsam.append(file)
fastaqsam = sorted(fastaqsam)

#Filtering the files that don't have the fastq extension 
substring = 'fastq'
fastaqsam_fil = [s for s in fastaqsam if substring in s]

#Run bwa index reference genomes
subprocess.run(["bwa", "index", index_file_with_the_path], check=True)

#Split the R1 and R2 
odd = fastaqsam_fil[1 : : 2] #all the R1 in the name 
even = fastaqsam_fil[: : 2] #all the R2 in the name 

for i in range(len(odd)):
    R1 = list()
    R2= list()
    filename = odd[i]
    filename_sam = os.path.splitext(filename)[0]
    output_file = os.path.join(output_path, f"{os.path.splitext(filename)[0]}.sam")
    R1 = seq_path + '/'+ odd[i]
    R2 = seq_path + '/'+ even[i]
    
    #Run bwa mem
    bwa_run = ["bwa", "mem", index_file_with_the_path, R2, R1]
    print(bwa_run)
    with open(output_file, "w") as sam_file:
        subprocess.run(bwa_run, stdout=sam_file, check=True)

ext = 'sam'   
files = os.listdir(output_path)
names_sam = [o for o in files if ext in o]

#Calling ecoli.py to analysis, map and save the matching positions in a .txt file. 
for sam in names_sam:
    subprocess.run(["python", "ecoli.py", sam, recipient, donor])

    
    
    