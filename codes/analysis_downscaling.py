import os
import gzip
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#Parameters
genome_length = 4600000  
read_length = 128        
fastq_dir = "./RELHS"  

#Fonction for counting reads 
def count_reads(fastq_path):
    open_func = gzip.open if fastq_path.endswith(".gz") else open
    with open_func(fastq_path, "rt") as f:
        n_lines = sum(1 for _ in f)
    return n_lines // 4

#Step 1: READ ORIGINAL FILES:read R1 and the corresponding R2 file
original_coverages = {}
for fname in os.listdir(fastq_dir):
    if fname.endswith(".fastq") or fname.endswith(".fastq.gz"):
        if "_R1" not in fname:
            continue 

        base = fname.split("_R1")[0]
        parts = base.split("_")
        strain = parts[0] #the name of the file is split (it contains the strain name and the resolution)

        # Check if it's the original file (no fraction = no x) or a fragment downscaled in resolution (with an x)
        if any("x" in p for p in parts):
            continue
        path_R1 = os.path.join(fastq_dir, fname)
        path_R2 = os.path.join(fastq_dir, base + "_R2" + fname.split("_R1")[1])
        reads_R1 = count_reads(path_R1)
        reads_R2 = count_reads(path_R2)
        total_reads = reads_R1 + reads_R2
        #Compute the coverage of the sequenced sample
        coverage_mean = total_reads * read_length / genome_length
        original_coverages[strain] = coverage_mean

fractions_to_plot = [0.5, 0.4, 0.25, 0.1]  
data = []

#Step 2: READ DOWNSCALED FILES
for fname in os.listdir(fastq_dir):
    if fname.endswith(".fastq") or fname.endswith(".fastq.gz"):
        if "_R1" not in fname:
            continue  
        base = fname.split("_R1")[0]
        parts = base.split("_")
        strain = parts[0]
        fraction = None
        for p in parts:
            if p.endswith("x"):
                fraction = float(p.replace("x",""))
                break

        if fraction is None or fraction not in fractions_to_plot:
            fraction = "1" 
        path_R1 = os.path.join(fastq_dir, fname)
        path_R2 = os.path.join(fastq_dir, base + "_R2" + fname.split("_R1")[1])
        reads_R1 = count_reads(path_R1)
        reads_R2 = count_reads(path_R2)
        total_reads = reads_R1 + reads_R2

        coverage_mean = total_reads * read_length / genome_length
        diff_to_ref = coverage_mean - original_coverages[strain]

        data.append({
            "Strain": strain,
            "Fraction": fraction,
            "Coverage_diff": diff_to_ref, 
            "Coverage_estimated": coverage_mean
        })

Coverage_recap = pd.DataFrame(data) 
Coverage_recap_mean = Coverage_recap.groupby("Fraction")["Coverage_estimated"].mean().reset_index()
Coverage_recap.to_csv("coverage_tableRELHS.csv", index=False)