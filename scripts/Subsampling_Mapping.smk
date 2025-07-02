# Subsample mapping script

with open('subsampling/subsamples_16G.txt', 'r') as f:
    samples = f.read().splitlines()

rule all:
    input:
        Check_Mapping = expand("subsampling/Checks/6-Mapping_combined_{sample}.done", sample = samples),
        Check_coverm = "subsampling/Checks/6.2-CoverM_Combined_16G.done"

rule combined_minimap2_map:
    input:
        FwdReads="subsampling/2-QC/2.2-ErrorCorrection/{sample}_EC_R1.fq.gz",
        RevReads="subsampling/2-QC/2.2-ErrorCorrection/{sample}_EC_R2.fq.gz",
        CombinedContigs="subsampling/5-dereplication/combined_vOTUs.fna"
    output:
        CombinedBam=temp("subsampling/6-Mapping/{sample}_combined_mapped.bam"),
        CombinedCheckMinimap2Mapping="subsampling/Checks/6-Mapping_combined_{sample}.done"
    threads: 16
    params:
        tag="{sample}"
    resources:
        mem_mb=128000,
        partition="high2"
    shell:'''
        minimap2 -ax sr -t {threads} {input.CombinedContigs} {input.FwdReads} {input.RevReads} | \
        samtools view -u | \
        samtools sort -o {output.CombinedBam} && \
        touch {output.CombinedCheckMinimap2Mapping}
        '''

rule CombinedCoverM:
    wildcard_constraints:
        sample="|".join(samples)
    input:
        bam = expand("subsampling/6-Mapping/{sample}_combined_mapped.bam", sample = samples),
        CheckStool = expand("subsampling/Checks/6-Mapping_combined_{sample}.done", sample = samples),
    output:
        CheckCombinedCoverM = "subsampling/Checks/6.2-CoverM_Combined_16G.done",
        coverage = "subsampling/6-Mapping/Combined_coverage_tpm_16G.tsv",
        counts = "subsampling/6-Mapping/Combined_coverage_counts_16G.tsv",
    threads: 16
    params:
        tag = "CM_16G",
        CoverMPath = "6-Mapping/",
        min_coverage = 75
    resources:
        mem_mb = 16000,
        partition = "high2",
        time = "19:00:00"
    shell:'''
        coverm contig -t {threads} --methods mean trimmed_mean covered_bases variance rpkm tpm --min-read-percent-identity 90 --min-covered-fraction {params.min_coverage} \
        --bam-files {input.bam} --output-format dense -o {output.coverage} && \
        coverm contig -t {threads} --methods mean trimmed_mean count --min-read-percent-identity 90 \
        --bam-files {input.bam} --output-format dense -o {output.counts} && \
        touch {output.CheckCombinedCoverM}
        '''

