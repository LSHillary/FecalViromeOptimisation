rule BacPhlip_vOTU:
    input:
        vOTUs = "5-dereplication/combined_Genomad_vOTUs.fna"
    output:
        CheckBacPhlipVOTU = "Checks/12-BacPhlip_vOTUs.done"
    threads: 64
    params:
        tag = "BacPhlip",
    resources:
        mem_mb = "128gb",
        partition = "med2",
        time = "2-24:00:00"
    shell:'''
        mkdir -p 12-Lifestyle_vOTUs && \
        micromamba run -n BACPHLIP_ENV bacphlip \
        -i {input.vOTUs} --multi_fasta -f && \
        touch {output.CheckBacPhlipVOTU}
        '''