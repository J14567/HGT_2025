#!/bin/bash
set -euo pipefail

SEED=42 
r1_extension="_R1_001.fastq.gz"
GENOME_LEN=4600000
READ_LEN=128        

for R1 in *${r1_extension}; do
    sample_base=$(basename "$R1" "$r1_extension")
    sample="$sample_base"
    R2="${sample}_R2_001.fastq.gz"

    reads_R1=$(( $(zcat "$R1" |  wc -l) / 4 ))
    reads_R2=$(( $(zcat "$R2" | wc -l) / 4 ))
    total_reads=$((reads_R1+reads_R2))

    initial_cov=$(awk -v r=$total_reads -v l=$READ_LEN -v g=$GENOME_LEN 'BEGIN{print r*l/g}')

    echo "=== Downscaling ongoing for $sample ==="
    for frac in 0.5 0.4 0.25 0.1; do 
    #for frac in 0.175; do 
        target_cov=$(awk -v c=$initial_cov -v f=$frac 'BEGIN{print c*f}')
        cov_lab=$(awk -v v=$target_cov 'BEGIN{printf "%.0f", v}')
        out_R1="${sample}_${frac}x_R1_001.fastq.gz"
        out_R2="${sample}_${frac}x_R2_001.fastq.gz"
        echo "-> Target ${cov_lab}x (fraction=$frac)"
        seqtk sample -s$SEED "$R1" $frac | gzip > "$out_R1"
        seqtk sample -s$SEED "$R2" $frac | gzip > "$out_R2"
    done
done