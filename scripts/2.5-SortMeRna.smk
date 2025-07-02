rule sortmerna:
    input:
        ForQC = "2-QC/2.2-ErrorCorrection/{sample}_EC_R1.fq.gz",
        RevQC = "2-QC/2.2-ErrorCorrection/{sample}_EC_R2.fq.gz",
        CheckQC = "Checks/2.2_EC_{sample}.done"
    output:
        CheckSortmeRNA = "Checks/2.5-SortMeRna_{sample}.done"
    threads: 16
    resources:
        mem_mb = 131072,
        tag = "{sample}",
        partition = "bmh",
        time = "2-00:00:00"
    params:
        tag = "{sample}"
    message:
        "Running sortmerna on {input.ForQC} and {input.RevQC}"
    shell:'''
    sortmerna --ref /group/jbemersogrp/databases/sortmerna/smr_v4.3_sensitive_db.fasta \
    --reads {input.ForQC} \
    --reads {input.RevQC} \
    --workdir 2-QC/2.5-SortMeRna/{params.tag}_wd \
    --idx-dir /home/lhillary/sortmerna/run/idx/ \
    --aligned 2-QC/2.5-SortMeRna/{params.tag}_paired_rrna \
    --other 2-QC/2.5-SortMeRna/{params.tag}_paired_other \
    --paired_in \
    --num_alignments 1 \
    --out2 \
    --fastx \
    -v --threads 40 && \
    touch {output.CheckSortmeRNA}
    '''

rule SortMeRnaFastQC:
    input:
        ForRRNA = "2-QC/2.5-SortMeRna/{sample}_paired_rrna_fwd.fq.gz",
        RevRRNA = "2-QC/2.5-SortMeRna/{sample}_paired_rrna_rev.fq.gz",
        ForOther = "2-QC/2.5-SortMeRna/{sample}_paired_other_fwd.fq.gz",
        RevOther = "2-QC/2.5-SortMeRna/{sample}_paired_other_rev.fq.gz",
        CheckSortmeRNA = "Checks/2.5-SortMeRna_{sample}.done"
    output:
        CheckSortMeRnaFastQC = "Checks/2.5-SortMeRnaFastQC_{sample}.done"
    params:
        OutputDirectory = "2-QC/2.5-SortMeRna/FastQC/",
        tag = "{sample}"
    threads: 16
    message:
        "Running FastQC on {wildcards.sample}"
    shell:'''
        mkdir -p {params.OutputDirectory} && \
        fastqc {input.ForRRNA} -o {params.OutputDirectory} -t {threads} && \
        fastqc {input.RevRRNA} -o {params.OutputDirectory} -t {threads} && \
        fastqc {input.ForOther} -o {params.OutputDirectory} -t {threads} && \
        fastqc {input.RevOther} -o {params.OutputDirectory} -t {threads} && \
        touch {output.CheckSortMeRnaFastQC}
    '''