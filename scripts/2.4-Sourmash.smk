rule SourmashSketchReads:
    input:
        R1 = "2-QC/2.2-ErrorCorrection/{sample}_EC_R1.fq.gz",
        R2 = "2-QC/2.2-ErrorCorrection/{sample}_EC_R2.fq.gz"
    output:
        signature="2-QC/2.4-Sourmash/{sample}.sig",
        Check = "Checks/2.4-SourmashSketchReads_{sample}.done"
    threads: 4
    resources:
        memory = "32G",
        time = "8:00:00",
        partition = "high2"
    params:
        tag="{sample}"
    priority: 50
    shell:'''
        micromamba run -n SOURMASH_ENV sourmash sketch dna -p scaled=1000,k=31,abund --merge {params.tag} -o {output.signature} {input.R1} {input.R2} && \
        touch {output.Check}
    '''

rule SourmashSketchContigs:
    input:
        contigs = "5-dereplication/combined_Genomad_vOTUs.fna"
    output:
        signature = "2-QC/2.4-Sourmash/Genomad_vOTUs.sig",
        check = "Checks/2.4-SourmashSketchContigs.done"
    threads: 4
    resources:
        memory = "32G",
        time = "8:00:00",
        partition = "high2"
    params:
        tag="Contigs"
    priority: 100
    shell:'''
        micromamba run -n SOURMASH_ENV sourmash sketch dna -p scaled=1000,k=31,abund --singleton -o {output.signature} {input.contigs} && \
        touch {output.check}
    '''

rule sourmash_abundance_histogram:
    input:
        sig = "2-QC/2.4-Sourmash/{sample}.sig"
    output:
        tsv = "2-QC/2.4-Sourmash/{sample}_abund_hist.tsv",
        check = "Checks/2.4-Sourmash_abund_hist_{sample}.done"
    threads: 2
    resources:
        memory = "8G",
        time = "2:00:00",
        partition = "high2"
    params:
        tag = "{sample}"
    shell:'''bash ../scripts/run_sourmash_abund_hist.sh {wildcards.sample}'''

rule SourmashGather:
    input:
        sample_sig = "2-QC/2.4-Sourmash/{sample}.sig",
        contigs_sig = "2-QC/2.4-Sourmash/Genomad_vOTUs.sig",
        CheckContigs = "Checks/2.4-SourmashSketchContigs.done",
        CheckReads = "Checks/2.4-SourmashSketchReads_{sample}.done"
    output:
        results = "2-QC/2.4-Sourmash/{sample}_gather.csv",
        check = "Checks/2.4-Sourmash_Gather_{sample}.done"
    threads: 4
    resources:
        memory = "32G",
        time = "8:00:00",
        partition = "high2"
    params:
        tag = "{sample}"
    shell:'''
        sourmash gather --scaled 1000 -k 31 \
        -o {output.results} {input.sample_sig} {input.contigs_sig} && \
        touch {output.check}
    '''