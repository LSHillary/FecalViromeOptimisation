import os
import yaml
from collections import defaultdict

configfile: "config/ViromeDataProcessing_config.yml"

with open(config['samples'], 'r') as f:
    samples = f.read().splitlines()

stool_samples = ['F3_S19','F4_S20', 'F5_S21', 'F6_S22', 'F7_S23', 'F8_S24', 'M1_S25', 'M2_S26', 'M3_S27', 'GP3_S28', 'GP4_S29', 'GP5_S30']

rule all:
    input:
        expand("Checks/11.1-UHGV_Index_UHGV.done"),
        expand("Checks/11.2-Biogeography_CoverM_UHGV_{sample}.done", sample=stool_samples),
        "Checks/Clustering_Biogeography_UHGV.done"

rule minimap2_index:
    input:
        fasta = "11-Biogeography/votus_full.fna"
    output:
        index = "11-Biogeography/Indexes/UHGV.mmi",
        Check = "Checks/11.1-UHGV_Index_UHGV.done"
    threads: 64
    params:
        tag = "UHGV_index"
    resources:
        mem_mb = "128gb",  # should be integer not "64gb" here
        time = "1-00:00:00",
        partition = "med2"
    shell: '''
        mkdir -p 11-Biogeography/Indexes Checks && \
        minimap2 -d {output.index} {input.fasta} && \
        touch {output.Check}
    '''

rule minimap2_and_coverm_UHGV:
    input:
        R1 = ancient("2-QC/2.2-ErrorCorrection/{sample}_EC_R1.fq.gz"),
        R2 = ancient("2-QC/2.2-ErrorCorrection/{sample}_EC_R2.fq.gz"),
        Index = "11-Biogeography/Indexes/UHGV.mmi"
    output:
        coverage = "11-Biogeography/CoverM/{sample}.tsv",
        coverage_mean = "11-Biogeography/CoverM/{sample}_mean.tsv",
        done = "Checks/11.2-Biogeography_CoverM_UHGV_{sample}.done"
    threads: 64
    params:
        min_coverage = 1,
        tag = "{sample}"
    benchmark: "Benchmarks/minimap2_coverm_{sample}.txt"
    resources:
        mem_mb = "128gb",
        partition = "med2",
        time = "2-00:00:00"
    shell:'''
        set -euo pipefail

        SCRATCH_WORKDIR=/scratch/$USER/$SLURM_JOBID/{wildcards.sample}
        mkdir -p $SCRATCH_WORKDIR

        echo "Running minimap2"
        minimap2 -ax sr -K 1G --sam-hit-only -2 -N 0 -t {threads} \
            {input.Index} \
            {input.R1} {input.R2} \
        | samtools view -u - \
        | samtools sort -@ {threads} -m 2G -o $SCRATCH_WORKDIR/aligned.bam

        echo "Running coverm"
        mkdir -p $(dirname {output.coverage})
        coverm contig \
            --bam-files $SCRATCH_WORKDIR/aligned.bam \
            --min-read-percent-identity 90 \
            --min-covered-fraction {params.min_coverage} \
            --methods covered_fraction \
            --output-format sparse \
            -t {threads} \
            -o {output.coverage}
        
        coverm contig \
            --bam-files $SCRATCH_WORKDIR/aligned.bam \
            --min-read-percent-identity 90 \
            --min-covered-fraction {params.min_coverage} \
            --methods mean trimmed_mean covered_bases variance rpkm tpm \
            --output-format sparse \
            -t {threads} \
            -o {output.coverage_mean}

        echo "Cleaning up scratch BAM"
        rm -rf $SCRATCH_WORKDIR

        touch {output.done}
    '''
