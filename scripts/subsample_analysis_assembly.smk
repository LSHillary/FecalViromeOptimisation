with open('subsampling/subsamples_13G.txt', 'r') as f:
    samples = f.read().splitlines()

# Creating a new list that includes only samples with specific substrings
metagenomes = [sample for sample in samples if "M1_S" in sample or "M2_S" in sample or "M3_S" in sample]

# Assuming 'samples' and 'metagenomes' are already defined
viromes = [sample for sample in samples if sample not in metagenomes]

# ls Virome/subsampling/*_1G_R1.fastq.gz | sed 's/_R1.fastq.gz$//' | sed 's#.*/##' > Virome/subsampling/subsamples_1G.txt


rule all:
    input:
        CheckEC = expand("subsampling/Checks/{sample}_ErrorCorrection.done", sample=samples),
        CheckDedup = expand("subsampling/Checks/{sample}_Deduplication.done", sample=samples),
        CheckRenameContigs = expand("subsampling/Checks/{sample}_RenameContigs.done", sample=samples),

rule QualityFiltering:
    input:
        ForRaw = "subsampling/{sample}_R1.fastq.gz",
        RevRaw = "subsampling/{sample}_R2.fastq.gz",
    output:
        ForQC=temp("subsampling/2-QC/2.1-Filtering/{sample}_QC_R1.fq.gz"),
        RevQC=temp("subsampling/2-QC/2.1-Filtering/{sample}_QC_R2.fq.gz"),
    resources:
        mem_mb=10000
    threads: 16
    params:
        tag = "QF{sample}"
    message: "QC filtering using bbduk"
    shell:'''
        mkdir -p 2-QC && \
        bbduk.sh in={input.ForRaw} in2={input.RevRaw} \
        ref=adapters,phix ktrim=r k=23 mink=11 hdist=1 tpe tbo \
        qtrim=r trimq=10 maxns=3 maq=3 minlen=50 mlf=0.333 \
        out={output.ForQC} out2={output.RevQC} \
        threads={threads}
    '''

# Run error correction on the QC filtered PE reads using tadpole
rule ErrorCorrection:
    input:
        ForQC="subsampling/2-QC/2.1-Filtering/{sample}_QC_R1.fq.gz",
        RevQC="subsampling/2-QC/2.1-Filtering/{sample}_QC_R2.fq.gz",
    output:
        ForEC="subsampling/2-QC/2.2-ErrorCorrection/{sample}_EC_R1.fq.gz",
        RevEC="subsampling/2-QC/2.2-ErrorCorrection/{sample}_EC_R2.fq.gz",
        CheckEC = "subsampling/Checks/{sample}_ErrorCorrection.done"
    resources:
        mem_mb=32000
    threads: 16
    params:
        tag = "EC{sample}"
    message: "Error correction using tadpole"
    shell:'''
    mkdir -p 2-QC/2.2-ErrorCorrection && \
    tadpole.sh \
        in={input.ForQC} \
        in2={input.RevQC} \
        out={output.ForEC} \
        out2={output.RevEC} \
        mode=correct ecc=t prefilter=1 \
        threads={threads} && \
    touch {output.CheckEC}
    '''

# Rule for removing PCR duplicates using clumpify
rule PcrDuplicateRemoval:
    input:
        ForEC="subsampling/2-QC/2.2-ErrorCorrection/{sample}_EC_R1.fq.gz",
        RevEC="subsampling/2-QC/2.2-ErrorCorrection/{sample}_EC_R2.fq.gz",
    output:
        ForDedup=temp("subsampling/2-QC/2.3-Deduplication/{sample}_Dedup_R1.fq.gz"),
        RevDedup=temp("subsampling/2-QC/2.3-Deduplication/{sample}_Dedup_R2.fq.gz"),
        CheckDedup="subsampling/Checks/{sample}_Deduplication.done"
    resources:
        mem_mb=32000
    threads: 16
    params:
        tag = "DD{sample}"
    message: "PCR duplicate removal using clumpify"
    shell:'''
        mkdir -p 2-QC/2.3-Deduplication && \
        clumpify.sh \
            in={input.ForEC} \
            in2={input.RevEC} \
            out={output.ForDedup} \
            out2={output.RevDedup} \
            dedupe subs=0 passes=2 \
            threads={threads} && \
        touch {output.CheckDedup}
    '''

# Individual assembly of each sample using megahit
rule IndividualAssembly:
    input:
        R1 = "subsampling/2-QC/2.3-Deduplication/{sample}_Dedup_R1.fq.gz",
        R2 = "subsampling/2-QC/2.3-Deduplication/{sample}_Dedup_R2.fq.gz"
    output:
        out_contig = temp("subsampling/3-Assembly/Individual/{sample}.contigs.fa"),
        CheckMegahitIndividual = "subsampling/Checks/{sample}_IndividualAssembly.done"
    threads: 16
    resources:
        mem_mb = 48000,
        cores = 16,
        partition = "high2",
        time = "3-00:00:00"
    params:
        tag = "AS{sample}",
        output_folder = "subsampling/3-Assembly/Individual",
        output_temp = "subsampling/3-Assembly/megahit_temp"
    shell:'''
    mkdir -p  subsampling/3-Assembly/megahit_temp/

    # megahit does not allow force overwrite, so each assembly needs to take place in it's own directory. 
    megahit -1 {input.R1} -2 {input.R2} \
    -t 16 --continue --k-min 27 --min-contig-len 1000 --presets meta-large \
    --out-dir {params.output_temp}/{wildcards.sample} \
    --out-prefix {wildcards.sample} && \
    cp {params.output_temp}/{wildcards.sample}/{wildcards.sample}.contigs.fa \
    {params.output_folder} && \
    touch {output.CheckMegahitIndividual}
    '''

# Create a rule to rename the contigs to {sample}_contig_<contig_number>
rule RenameContigs:
    input:
        contigs = "subsampling/3-Assembly/Individual/{sample}.contigs.fa"
    output:
        renamed_contigs = "subsampling/3-Assembly/Individual/{sample}_renamed_contigs.fa",
        CheckRenameContigs = "subsampling/Checks/{sample}_RenameContigs.done"
    params:
        tag = "RC{sample}"
    resources:
        partition = "high2"
    shell:'''
    awk '/^>/{{print ">" "{wildcards.sample}" "_contig_" ++i; next}}{{print}}' < {input.contigs} > {output.renamed_contigs} && \
    touch {output.CheckRenameContigs}
    '''
