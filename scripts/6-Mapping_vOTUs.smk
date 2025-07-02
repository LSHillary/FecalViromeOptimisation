# Snakemake rule for mapping reads to contigs using minimap2, pipe to samtools to convert to sorted bam
rule stool_minimap2_map:
    input:
        FwdReads = "2-QC/2.2-ErrorCorrection/{sample}_EC_R1.fq.gz",
        RevReads = "2-QC/2.2-ErrorCorrection/{sample}_EC_R2.fq.gz",
        StoolContigs = "5-dereplication/combined_Genomad_vOTUs.fna",
    output:
        StoolBam = "6-Mapping/{sample}_stool_Genomad_mapped.bam",
        StoolCheckMinimap2Mapping = "Checks/6-Mapping_stool_Genomad_{sample}.done"
    threads: 16
    params:
        tag = "{sample}",
        
        CheckClusterStool = "Checks/LeidenClustering_Genomad_combined.done"
    resources:
        mem_mb = "64gb",
        partition = "high2"
    shell:'''
        minimap2 -axsr -t {threads} {input.StoolContigs} {input.FwdReads} {input.RevReads} | samtools view -u | samtools sort -o {output.StoolBam} \
        && \
        touch {output.StoolCheckMinimap2Mapping}
        '''

rule human_minimap2_map:
    input:
        FwdReads = "2-QC/2.2-ErrorCorrection/{sample}_EC_R1.fq.gz",
        RevReads = "2-QC/2.2-ErrorCorrection/{sample}_EC_R2.fq.gz",
        HumanContigs = "/home/lhillary/databases/HumanGenome/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
    output:
        HumanBam = temp("6-Mapping/{sample}_human_mapped.bam"),
        HumanCheckMinimap2Mapping = "Checks/6-Mapping_human_{sample}.done"
    threads: 16
    params:
        tag = "{sample}"
    resources:
        mem_mb = "64gb",
        partition = "med2"
    shell:'''
        minimap2 -axsr -t {threads} {input.HumanContigs} {input.FwdReads} {input.RevReads} | samtools view -u | samtools sort -o {output.HumanBam} \
        && \
        touch {output.HumanCheckMinimap2Mapping}
        '''

rule StoolCoverM:
    wildcard_constraints:
        stool_sample="|".join(stool_samples)
    input:
        bam = expand("6-Mapping/{stool_sample}_stool_Genomad_mapped.bam", stool_sample = stool_samples),
        CheckStool = expand("Checks/6-Mapping_stool_Genomad_{stool_sample}.done", stool_sample = stool_samples),
    output:
        CheckStoolCoverM = "Checks/6.1-CoverM_Stool_Genomad.done",
        coverage = "6-Mapping/Stool_Genomad_coverage_tpm.tsv",
    threads: 16
    params:
        tag = "Stool",
        CoverMPath = "6-Mapping/",
        min_coverage = 75
    resources:
        mem_mb = 16000,
        partition = "med2",
        time = "19:00:00"
    shell:'''
        coverm contig -t {threads} --methods mean trimmed_mean covered_bases variance rpkm tpm --min-read-percent-identity 90 --min-covered-fraction {params.min_coverage} \
        --bam-files {input.bam} --output-format dense -o {output.coverage} && \
        touch {output.CheckStoolCoverM}
        '''

rule HumanCoverM:
    input:
        bam = expand("6-Mapping/{sample}_human_mapped.bam", sample = stool_samples),
        CheckHuman = expand("Checks/6-Mapping_human_{sample}.done", sample = stool_samples),
    output:
        CheckHumanCoverM = "Checks/6.1-CoverM_Human.done",
        coverage = "6-Mapping/Human_coverage.tsv",
    threads: 16
    params:
        tag = "Human",
        CoverMPath = "6-Mapping/",
    resources:
        mem_mb = 16000,
        partition = "med2",
        time = "4:00:00"
    shell:'''
        coverm contig -t {threads} --methods covered_bases count --min-read-percent-identity 90\
        --bam-files {input.bam} --output-format dense -o {output.coverage} && \
        touch {output.CheckHumanCoverM}
        '''

rule StoolCoverMCounts:
    wildcard_constraints:
        stool_sample="|".join(stool_samples)
    input:
        bam = expand("6-Mapping/{stool_sample}_stool_Genomad_mapped.bam", stool_sample = stool_samples),
        CheckStool = expand("Checks/6-Mapping_stool_Genomad_{stool_sample}.done", stool_sample = stool_samples),
    output:
        CheckStoolCoverMCounts = "Checks/6.1-CoverM_Stool_Genomad_Counts.done",
        coverage = "6-Mapping/Stool_Genomad_coverage_counts.tsv",
    threads: 16
    params:
        tag = "Stool",
        CoverMPath = "6-Mapping/",
        min_coverage = 75
    resources:
        mem_mb = 16000,
        partition = "high2",
        time = "19:00:00"
    shell:'''
        coverm contig -t {threads} --methods mean trimmed_mean count --min-read-percent-identity 90 \
        --bam-files {input.bam} --output-format dense -o {output.coverage} && \
        touch {output.CheckStoolCoverMCounts}
        '''