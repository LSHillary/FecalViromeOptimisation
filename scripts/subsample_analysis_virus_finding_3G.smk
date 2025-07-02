with open('subsampling/subsamples_3G.txt', 'r') as f:
    samples = f.read().splitlines()

# ls Virome/subsampling/*_3G_R1.fastq.gz | sed 's/_R1.fastq.gz$//' | sed 's#.*/##' > Virome/subsampling/subsamples_3G.txt

rule all:
    input:
        CheckGenomad = expand("subsampling/Checks/Genomad_{sample}.done", sample = samples),
        CheckFilteredGenomadResults = expand("subsampling/Checks/4.2-FilterGenomadResults_{sample}.done", sample = samples),
        CheckDereplication = expand("subsampling/Checks/LeidenClustering_Genomad_combined.done")

rule GeNomad:
    input:
        ContigsIn = "subsampling/3-Assembly/Individual/{sample}_renamed_contigs.fa"
    output:
        CheckGenomad = "subsampling/Checks/Genomad_{sample}.done"
    threads: 16
    params:
        tag = "G{sample}",
        GenomadDB = "/group/jbemersogrp/databases/genomad/genomad_db",
        GenomadFolder = "subsampling/4-virus_identification/genomad/{sample}"
    resources:
        mem_mb = 64000,
        partition = "med2"
    shell:'''
        mkdir -p {params.GenomadFolder} && \
        micromamba run -n genomad_env genomad end-to-end --cleanup --enable-score-calibration -t {threads} {input.ContigsIn} \
        {params.GenomadFolder} {params.GenomadDB} && \
        touch {output.CheckGenomad}
        '''

rule FilterGenomad:
    input:
        CheckGenomad = "subsampling/Checks/Genomad_{sample}.done"
    output:
        FilteredGenomadResults = "subsampling/4-virus_identification/genomad/{sample}/{sample}_FilteredGenomadResults.tsv",
        CheckFilteredGenomadResults = "subsampling/Checks/4.2-FilterGenomadResults_{sample}.done",
        FilteredContigs = "subsampling/4-virus_identification/genomad/{sample}_FilteredContigs.fna"
    resources:
        partition = "med2"
    params:
        tag = "{sample}",
        Contigs = "subsampling/4-virus_identification/genomad/{sample}/{sample}_renamed_contigs_summary/{sample}_renamed_contigs_virus.fna",
        GenomadResults = "subsampling/4-virus_identification/genomad/{sample}/{sample}_renamed_contigs_summary/{sample}_renamed_contigs_virus_summary.tsv"
    shell:'''
        awk -F'\t' 'NR==1 || $2 >= 10000 || ($2 >= 1000 && $11 ~ /Monodnaviria/)' {params.GenomadResults} | \
        cut -f1 > {output.FilteredGenomadResults} && \
        seqtk subseq {params.Contigs} {output.FilteredGenomadResults} > {output.FilteredContigs} && \
        touch {output.CheckFilteredGenomadResults}
        '''

rule PreClusteringAggregation:
    input:
        combined_contigs = expand("subsampling/4-virus_identification/Genomad/{sample}_FilteredContigs.fna", sample=stool_samples),
    output:
        combined_contigs = "subsampling/4-virus_identification/combined_Genomad_FilteredContigs.fna",
    threads: 12
    params:
        tag="combined_Genomad"
    resources:
        mem_mb = 16000,
        partition = "med2",
        time = "5-00:00:00"
    shell:'''
        cat {input.combined_contigs} > {output.combined_contigs} && \
        touch {output.CheckAggregation}
        '''
 
rule Clustering:
    input:
        contigs = "subsampling/4-virus_identification/combined_Genomad_FilteredContigs.fna",
    output:
        Leiden_dereplicated = "subsampling/5-dereplication/combined_Genomad_3G_vOTUs.fna",
        CheckCluster = "subsampling/Checks/LeidenClustering_Genomad_combined.done"
    threads: 48
    params:
        tag="combined"
    resources:
        mem_mb = 64000,
        partition = "med2",
        time = "5-00:00:00"
    shell:'''
        snakemake --nolock --config input={input.contigs} \
                   output={output.Leiden_dereplicated} \
                   leiden_resolution=1.0 \
                   min_ani=0.95 \
                   min_cov=0.85 \
                   blast_max_evalue=1e-5 \
                   blast_threads={threads} \
                   --cores {threads} -s ../scripts/vOTU_clustering/contig-ani-leiden-clustering-pipeline.smk && \
        touch {output.CheckCluster}
        '''