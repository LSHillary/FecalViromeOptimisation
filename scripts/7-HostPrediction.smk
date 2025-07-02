# Snakefile for predicting virus hosts using iPhop

rule iPhop_predict:
    input:
        ContigsIn = "5-dereplication/combined_Genomad_vOTUs.fna"
    output:
        Check_iPhop = "Checks/iphop_combined_Genomad.done",
        OutDirectory = directory("7-HostPrediction/combined")
    threads: 12
    resources:
        mem_mb = "64gb",
        partition = "bmh",
        time = "24:00:00"
    params:
        tag = "SC_combined",
        Database = "/group/jbemersogrp/databases/iphop/latest/Aug_2023_pub_rw/",
    shell:'''
        micromamba run -n iphop_env iphop predict -f {input.ContigsIn} -o {output.OutDirectory} -t {threads} -d {params.Database} && \
        touch {output.Check_iPhop}
        '''