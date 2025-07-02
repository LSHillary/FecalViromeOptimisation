#### Note that this is an optional step and the pipeline has been designed to work without it.

# Run error correction on the QC filtered PE reads using tadpole
rule ErrorCorrection:
    input:
        ForQC="2-QC/2.1-Filtering/{sample}_QC_R1.fq.gz",
        RevQC="2-QC/2.1-Filtering/{sample}_QC_R2.fq.gz",
        CheckQC="Checks/2-QC_{sample}.done"
    output:
        ForEC="2-QC/2.2-ErrorCorrection/{sample}_EC_R1.fq.gz",
        RevEC="2-QC/2.2-ErrorCorrection/{sample}_EC_R2.fq.gz",
        CheckEC="Checks/2.2_EC_{sample}.done"
    resources:
        mem_mb=32000
    threads: 8
    params:
        tag="{sample}"
    message: "Error correction using tadpole"
    shell:'''
    mkdir -p 2-QC/2.2-ErrorCorrection && \
    tadpole.sh \
        in={input.ForQC} \
        in2={input.RevQC} \
        out={output.ForEC} \
        out2={output.RevEC} \
        mode=correct ecc=t prefilter=1 \
        threads={threads} \
    && \
    touch {output.CheckEC}
    '''

# Rule for removing PCR duplicates using clumpify
rule PcrDuplicateRemoval:
    input:
        ForEC="2-QC/2.2-ErrorCorrection/{sample}_EC_R1.fq.gz",
        RevEC="2-QC/2.2-ErrorCorrection/{sample}_EC_R2.fq.gz",
        CheckEC="Checks/2.2_EC_{sample}.done"
    output:
        ForDedup="2-QC/2.3-Deduplication/{sample}_Dedup_R1.fq.gz",
        RevDedup="2-QC/2.3-Deduplication/{sample}_Dedup_R2.fq.gz",
        CheckDedup="Checks/2.2_Dedup_{sample}.done"
    resources:
        mem_mb=10000
    threads: 8
    params:
        tag="{sample}"
    message: "PCR duplicate removal using clumpify"
    shell:'''
        mkdir -p 2-QC/2.3-Deduplication && \
        clumpify.sh \
            in={input.ForEC} \
            in2={input.RevEC} \
            out={output.ForDedup} \
            out2={output.RevDedup} \
            dedupe subs=0 passes=2 \
            threads={threads} \
        && \
        touch {output.CheckDedup}
    '''

