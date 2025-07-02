#### 4 - Assembly ####

# Individual assembly of each sample using megahit
rule IndividualAssembly:
    input:
        CheckECFastQC="Checks/2.2-ErrorCorrectionFastQC_{sample}.done",
        R1 = "2-QC/2.3-Deduplication/{sample}_Dedup_R1.fq.gz",
        R2 = "2-QC/2.3-Deduplication/{sample}_Dedup_R2.fq.gz"
    output:
        CheckMegahitIndividual = "Checks/{sample}_IndividualAssembly.done",
        out_contig = "3-Assembly/Individual/{sample}.contigs.fa"
    threads: 16
    resources:
        mem_mb = 48000,
        cores = 16,
        partition = "bmh",
        time = "3-00:00:00"
    params:
        tag = "{sample}",
        output_folder = "3-Assembly/Individual",
        output_temp = "3-Assembly/megahit_temp"
    message: "paired end assembly on {params.tag}"
    shell:'''
    mkdir -p  3-Assembly/megahit_temp/

    # megahit does not allow force overwrite, so each assembly needs to take place in it's own directory. 
    megahit -1 {input.R1} -2 {input.R2} \
    -t 16 --continue --k-min 27 --min-contig-len 1000 --presets meta-large \
    --out-dir {params.output_temp}/{wildcards.sample} \
    --out-prefix {wildcards.sample} && \
    mv {params.output_temp}/{wildcards.sample}/{wildcards.sample}.contigs.fa \
    {params.output_folder} && touch {output.CheckMegahitIndividual}
    '''

# Create a rule to rename the contigs to {sample}_contig_<contig_number>
rule RenameContigs:
    input:
        contigs = "3-Assembly/Individual/{sample}.contigs.fa",
        CheckMegahitIndividual = "Checks/{sample}_IndividualAssembly.done"
    output:
        renamed_contigs = "3-Assembly/Individual/{sample}_renamed_contigs.fa",
        CheckRenameContigs = "Checks/RenameContigs_{sample}.done"
    params:
        tag = "{sample}"
    shell:'''
    awk '/^>/{{print ">" "{params.tag}" "_contig_" ++i; next}}{{print}}' < {input.contigs} > {output.renamed_contigs} && \
    touch {output.CheckRenameContigs}
    '''
