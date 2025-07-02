SAMPLES = ['F3_S19','F4_S20', 'F5_S21', 'F6_S22', 'F7_S23', 'F8_S24', 'M1_S25', 'M2_S26', 'M3_S27', 'GP3_S28', 'GP4_S29', 'GP5_S30']

rule all:
    input:
        expand("subsampling/{sample}_subsampled_depth_{depth}G_R1.fastq.gz", sample=SAMPLES, depth=range(1, 18)),
        expand("subsampling/{sample}_subsampled_depth_{depth}G_R2.fastq.gz", sample=SAMPLES, depth=range(1, 18))

rule subsample_reads:
    input:
        R1 = "0-raw/{sample}_raw_R1.fq.gz",
        R2 = "0-raw/{sample}_raw_R2.fq.gz"
    output:
        R1 = "subsampling/{sample}_subsampled_depth_{depth}G_R1.fastq.gz",
        R2 = "subsampling/{sample}_subsampled_depth_{depth}G_R2.fastq.gz"
    params:
        num_reads = lambda wildcards: str(int(wildcards.depth) * 1000000000 / (150 * 2)),  # Adjust calculation as necessary,
        tag = "{sample}"
    threads: 8
    resources:
        mem_mb = 64000,
        partition = "bmh"
    shell: '''        # Use seqtk to subsample reads
        seqtk sample -s100 {input.R1} {params.num_reads} | gzip > {output.R1}
        seqtk sample -s100 {input.R2} {params.num_reads} | gzip > {output.R2}
        '''
