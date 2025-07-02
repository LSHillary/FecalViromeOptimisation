rule khmer_interleave:
    input:
        ForQC = "2-QC/2.2-ErrorCorrection/{sample}_EC_R1.fq.gz",
        RevQC = "2-QC/2.2-ErrorCorrection/{sample}_EC_R2.fq.gz",
        CheckQC = "Checks/2.2_EC_{sample}.done"
    output:
        InterleavedReads = "2-QC/2.4-Khmer/{sample}_EC_I.fq",
        Check = "Checks/2.4_Khmer-Interleave_{sample}.done"
    params:
        tag = "{sample}"
    threads: 8
    resources:
        mem_mb = 50000,
        partition = "bmh",
        time = "10:00:00"
    output:
    shell:'''
        interleave-reads.py {input.ForQC} {input.RevQC} --gzip -o {output.InterleavedReads} && \
        touch {output.Check}
    '''

rule khmer_trim:
    input:
        InterleavedReads = "2-QC/2.4-Khmer/{sample}_EC_I.fq",
        Check = "Checks/2.4_Khmer-Interleave_{sample}.done"
    output:
        KhmerTrimmedReads = "2-QC/2.4-Khmer/{sample}.trimmed.fq.gz",
        Check = "Checks/2.4_Khmer-Trimmed_{sample}.done"
    params:
        tag = "{sample}"
    threads: 64
    resources:
        mem_mb = 524288,
        time = "3-10:00:00",
        partition = "bmm"
    shell:'''
        trim-low-abund.py -V -C 2 -M {resources.mem_mb}M -k 31 --gzip -o {output.KhmerTrimmedReads} {input.InterleavedReads} && \
        touch {output.Check}
    '''

rule Khmer_count:
    input:
        KhmerTrimmedReads = "2-QC/2.4-Khmer/{sample}.trimmed.fq.gz",
        Check = "Checks/2.4_Khmer-Trimmed_{sample}.done"
    output:
        CountGraph = "2-QC/2.4-Khmer/{sample}.KmerCount.tsv",
        Check = "Checks/2.4_Khmer-Counted_{sample}.done"
    params:
        tag = "{sample}"
    threads: 64
    resources:
        mem_mb = 524288,
        time = "3-10:00:00",
        partition = "bmm"
    shell:'''
        load-into-counting.py -k 31 -M {resources.mem_mb}M {output.CountGraph} {input.KhmerTrimmedReads} && \
        touch {output.Check}
    '''

rule Khmer_abundance:
    input:
        KhmerTrimmedReads = "2-QC/2.4-Khmer/{sample}.trimmed.fq.gz",
        CountGraph = "2-QC/2.4-Khmer/{sample}.KmerCount.tsv",
        Check = "Checks/2.4_Khmer-Counted_{sample}.done"
    output:
        AbundanceTable = "2-QC/2.4-Khmer/{sample}.AbundanceTable.tsv",
        Check = "Checks/2.4_Khmer-Abundance_{sample}.done"
    params:
        tag = "{sample}"
    threads: 64
    resources:
        mem_mb = 524288,
        time = "3-10:00:00",
        partition = "bmm"
    shell:'''
        abundance-dist.py {input.CountGraph} {input.KhmerTrimmedReads} {output.AbundanceTable} && \
        touch {output.Check}
    '''