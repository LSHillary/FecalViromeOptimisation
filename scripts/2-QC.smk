#### 2 - QC Filtering ####

# Rule 2.1 Adaptor removal and quality filtering on raw reads using bbduk
rule QualityFiltering:
    input:
        ForRaw = "0-raw/{sample}_L002_R1_001.fastq.gz",
        RevRaw = "0-raw/{sample}_L002_R2_001.fastq.gz",
        CheckRawFastQC = expand("Checks/1.1-RawFastQC_{sample}.done", sample = samples)
    output:
        ForQC="2-QC/2.1-Filtering/{sample}_QC_R1.fq.gz",
        RevQC="2-QC/2.1-Filtering/{sample}_QC_R2.fq.gz",
        SingleQC="2-QC/2.1-Filtering/{sample}_QC_U.fq.gz",
        CheckQC="Checks/2-QC_{sample}.done"
    resources:
        mem_mb=10000
    threads: 8
    params:
        tag="{sample}"
    message: "QC filtering using bbduk"
    shell:'''
        mkdir -p 2-QC && \
        bbduk.sh in={input.ForRaw} in2={input.RevRaw} \
        ref=adapters,phix ktrim=r k=23 mink=11 hdist=1 tpe tbo \
        qtrim=r trimq=10 maxns=3 maq=3 minlen=50 mlf=0.333 \
        out={output.ForQC} out2={output.RevQC} outs={output.SingleQC} \
        threads={threads} \
        && \
        touch {output.CheckQC}
    '''

# Rule 2.2 Run FastQC on processed reads
rule FilteringFastQC:
    input:
        ForProcessed="2-QC/2.1-Filtering/{sample}_QC_R1.fq.gz",
        RevProcessed="2-QC/2.1-Filtering/{sample}_QC_R2.fq.gz",
        SingleProcessed="2-QC/2.1-Filtering/{sample}_QC_U.fq.gz",
        CheckQC="Checks/2-QC_{sample}.done"
    output:
        CheckProcessedFastQC="Checks/FilteredFastQC_{sample}.done"
    params:
        PairedOutputDirectory="2-QC/FastQC/Filtered/Paired",
        SingletonsOutputDirectory="2-QC/FastQC/Filtered/Singletons",
        tag="{sample}"
    threads: 8
    resources:
        mem_gb=30
    shell:'''
        mkdir -p {params.PairedOutputDirectory} {params.SingletonsOutputDirectory} && \
        fastqc {input.ForProcessed} -o {params.PairedOutputDirectory} -t {threads} && \
        fastqc {input.RevProcessed} -o {params.PairedOutputDirectory} -t {threads} && \
        fastqc {input.SingleProcessed} -o {params.SingletonsOutputDirectory} -t {threads} && \
        touch {output.CheckProcessedFastQC}
    '''

