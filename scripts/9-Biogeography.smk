rule BG_Clustering:
    input:
        contigs = expand("4-virus_identification/Genomad/{sample}_FilteredContigs.fna", sample = stool_samples),
        uhgvc_vOTUs = "11-Biogeography/votus_full.fna",
    output:
        Combined_Contigs = temp("11-Biogeography/AllContigs.fna"),
        Dereplicated_vOTUs = "9-Biogeography/BG_vOTUs.fna",
        CheckCluster = "Checks/Clustering_Biogeography.done"
    threads: 64
    params:
        tag="BG_Clustering"
    resources:
        mem_mb = "256gb",
        partition = "bmh",
        time = "3-00:00:00"
    shell:'''
        cat {input.contigs} {input.uhgvc_vOTUs} > {output.Combined_Contigs} && \
        snakemake --nolock --config input={output.Combined_Contigs} \
                   output={output.Dereplicated_vOTUs} \
                   leiden_resolution=1.0 \
                   min_ani=0.95 \
                   min_cov=0.85 \
                   blast_max_evalue=1e-5 \
                   blast_threads={threads} \
                   --cores {threads} -s ../scripts/vOTU_clustering/contig-ani-leiden-clustering-pipeline.smk && \
        touch {output.CheckCluster}
        '''
