rule Pharokka_all:
    input:
        ProteinsIn = "5-dereplication/{sample}_FilteredProteins.faa"
    output:
        OutDirectory = directory("8-FunctionalAnnotation/{sample}_proteins"),
        CheckPharokkaAll = "Checks/8.1-Pharokka_{sample}_proteins.done"
    params:
        tag = "{sample}",
        PharokkaDatabase = "/home/lhillary/Software/Pharokka/pharokka_v1.4.0_databases/."
    threads: 32
    resources:
        mem_mb = "64gb",
        time = "2-4:00:00",
        partition = "bmm"
    shell:'''
        mkdir -p {output.OutDirectory} && \
        micromamba run -n PHAROKKA_env pharokka_proteins.py -i {input.ProteinsIn} -o {output.OutDirectory} -d {params.PharokkaDatabase} -t {threads} \
        -f && \
        touch {output.CheckPharokkaAll}
    '''

rule DefenseFinder_all:
    input:
        ProteinsIn = "5-dereplication/{sample}_FilteredProteins.faa"
    output:
        OutDirectory = directory("8-FunctionalAnnotation/DefenseFinder/{sample}"),
        CheckDefenseFinderAll = "Checks/8.4-DefenseFinder_{sample}.done"
    params:
        tag = "{sample}",
        DF_ModelsDir = "/group/jbemersogrp/databases/defensefinder"
    threads: 32
    resources:
        mem_mb = "64gb",
        time = "2-4:00:00",
        partition = "bmm"
    shell:'''
        micromamba run -n DEFENSEFINDER_ENV defense-finder run \
        -o {output.OutDirectory} \
        --models-dir {params.DF_ModelsDir} \
        --antidefensefinder \
        --db-type gembase \
        {input.ProteinsIn} && \
        touch {output.CheckDefenseFinderAll}
        '''