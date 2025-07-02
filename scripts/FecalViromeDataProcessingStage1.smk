import os
import glob

# This script runs FastQC on the raw samples

# Define wildcards for samples
SAMPLES, = glob_wildcards('0-raw/{sample}_L002_R1_001.fastq.gz')

# Rule that defines the final targets
rule all_preprocessing:
    input:
        CheckRawFastQC = expand("Checks/1-RawFastQC_{sample}.done", sample = SAMPLES)

rule RawFastQC:
    input:
        ForRaw="0-raw/{sample}_L002_R1_001.fastq.gz",
        RevRaw="0-raw/{sample}_L002_R2_001.fastq.gz"
    output:
        CheckRawFastQC="Checks/1-RawFastQC_{sample}.done"
    params:
        PairedOutputDirectory="1-RawFastQC/FastQC/",
        tag="{sample}"
    threads: 8
    resources:
        mem_mb="8gb"
    shell:'''
        mkdir -p {params.PairedOutputDirectory} && \
        fastqc {input.ForRaw} -o {params.PairedOutputDirectory} -t {threads} && \
        fastqc {input.RevRaw} -o {params.PairedOutputDirectory} -t {threads} && \
        touch {output.CheckRawFastQC}
    '''