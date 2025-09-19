######################################################################################################
# SRA HUMAN READ SCRUBBING PIPELINE
# Author: Luke Hillary
# Date: 2025-09-19
#
# Description:
# This Snakemake pipeline performs human read removal on paired-end FASTQ files prior to SRA submission.
# It uses the NCBI SRA Human Scrubber tool, wrapped via an Apptainer (.sif) container for HPC use.
# Raw gzipped reads are first decompressed using `unpigz` for fast parallel decompression, and then
# passed to `scrub.sh` to remove human-derived reads.
#
# Input:
#   - Paired-end FASTQ files (.fastq.gz) located in: 0-raw/
#     e.g., 0-raw/F3_S19_L002_R1_001.fastq.gz, 0-raw/F3_S19_L002_R2_001.fastq.gz
#
# Output:
#   - Scrubbed, gzipped FASTQ files in: 13-Submission/
#     e.g., F3_S19.scrubbed.R1.fastq.gz, F3_S19.scrubbed.R2.fastq.gz
#   - Completion flag: Checks/13-HumanScrub_{sample}.done
#
# Requirements:
#   - Apptainer (module load apptainer)
#   - unpigz (module load pigz OR conda install -c conda-forge pigz)
#   - sra-human-scrubber.sif in /home/lhillary/apptainer/
#
# Notes:
#   - Each FASTQ file is processed independently by `scrub.sh`
#   - The paired-end reads are scrubbed separately, as required by the tool
#   - Temporary unzipped FASTQ files are cleaned automatically after the workflow completes
#   - Human k-mers are replaced with 'N' in the output files
#
# Usage:
#   snakemake -j 8 --use-singularity -s this_script.smk
######################################################################################################

rule unzip_fastq:
    input:
        ForRaw = "0-raw/{sample}_L002_R1_001.fastq.gz",
        RevRaw = "0-raw/{sample}_L002_R2_001.fastq.gz",
    output:
        ForUnzipped = temp("13-Submission/{sample}_R1.fastq"),
        RevUnzipped = temp("13-Submission/{sample}_R2.fastq"),
        Check = "Checks/13.1-Unzip_{sample}.done"
    threads: 8
    resources:
        mem_mb = 16000,
        partition = "med2",
        time = "1-00:00:00"
    params:
        tag = "{sample}"
    shell:
        """
        mkdir -p 13-Submission && \
        unpigz -p {threads} -c {input.ForRaw} > {output.ForUnzipped} && \
        unpigz -p {threads} -c {input.RevRaw} > {output.RevUnzipped} && \
        touch {output.Check}
        """

rule SRA_Human_Scrub_Pipe:
    input:
        ForUnzipped = "13-Submission/{sample}_R1.fastq",
        RevUnzipped = "13-Submission/{sample}_R2.fastq"
    output:
        ForScrubbed = "13-Submission/{sample}.scrubbed.R1.fastq.gz",
        RevScrubbed = "13-Submission/{sample}.scrubbed.R2.fastq.gz",
        done = "Checks/13.2-HumanScrub_Pipe_{sample}.done"
    threads: 16
    resources:
        mem_mb = 64000,
        time = "1-00:00:00",
        partition = "bmm"
    params:
        sif = "/home/lhillary/apptainer/sra-human-scrubber.sif",
        tag = "{sample}"
    shell:
        """
        apptainer exec {params.sif} scrub.sh -i {input.ForUnzipped} -o - | \
            pigz -p {threads} > {output.ForScrubbed}

        apptainer exec {params.sif} scrub.sh -i {input.RevUnzipped} -o - | \
            pigz -p {threads} > {output.RevScrubbed}

        touch {output.done}
        """