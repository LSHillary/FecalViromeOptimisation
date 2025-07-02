rule BacPhlip:
    input:
        vOTUs = "12-Lifestyle/{sample}_FilteredContigs.fna"
    output:
        CheckBacPhlip = "Checks/12-BacPhlip_{sample}.done"
    threads: 32
    params:
        tag = "BacPhlip",
    resources:
        mem_mb = "64gb",
        partition = "high2",
        time = "2-24:00:00"
    shell:'''
        mkdir -p 12-Lifestyle && \
        micromamba run -n BACPHLIP_ENV bacphlip \
        -i {input.vOTUs} --multi_fasta -f && \
        touch {output.CheckBacPhlip}
        '''
