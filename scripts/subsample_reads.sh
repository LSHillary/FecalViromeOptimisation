#!/bin/bash

forward_reads=$1
reverse_reads=$2
output_dir=$3

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Set the fixed seed
seed=100

# Generate a unique output file name
output_file="${output_dir}/subsampled_reads"

# Subsample the reads using seqtk with the fixed seed
seqtk sample -s "$seed" "$forward_reads" 10000 | gzip > "${output_file}_R1.fq.gz"
seqtk sample -s "$seed" "$reverse_reads" 10000 | gzip > "${output_file}_R2.fq.gz"

echo "Generated subsampled reads: ${output_file}_R1.fq.gz and ${output_file}_R2.fq.gz"
