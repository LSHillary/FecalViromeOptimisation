rule vContact3:
    input:
        ProteinsIn = "5-Dereplication/COBRA/COBRA_vOTUs_proteins.faa",
        vOTUs = "5-Dereplication/COBRA_vOTUs.fna",
        G2G_file = "10-Taxonomy/COBRA_vOTUs_g2g.tsv",
        Length_file = "10-Taxonomy/vOTUs_length.tsv"
    output:
        Check = "Checks/10.2-vContact3_COBRA.done"
    params:
        Output_Directory = "10-Taxonomy/COBRA",
        tag = "vContact3",
        vContact3_db = "/group/jbemersogrp/databases/vContact3/db_223"
    threads: 64
    resources:
        mem_mb = "500gb",
        partition = "bmh",
        time = "3-00:00:00"
    shell:'''
        mkdir -p {params.Output_Directory} && \
        micromamba run -n VCONTACT3_CONDA_ENV vcontact3 run --proteins {input.ProteinsIn} --gene2genome {input.G2G_file} --output {params.Output_Directory} --db-domain "prokaryotes" --db-version 223 --db-path {params.vContact3_db} \
        --gene2genome {input.G2G_file} --len-nucleotide {input.Length_file} --threads {threads} && \
        touch {output.Check}
    '''

rule vContact3_nucleotide:
    input:
        ContigsIn = "5-dereplication/{group}_Genomad_vOTUs.fna",
    output:
        CheckVC3 = "Checks/10.1-vContact3_nucleotide_{group}_Genomad.done"
    params:
        Output_Directory = "10-Taxonomy/{group}",
        tag = "{group}",
        vContact3_db = "/group/jbemersogrp/databases/vContact3/db_223"
    threads: 64
    resources:
        mem_mb = "64gb",
        partition = "med2",
        time = "3-00:00:00"
    shell:'''
        mkdir -p {params.Output_Directory} && \
        micromamba run -n VCONTACT3_CONDA_ENV vcontact3 run --nucleotide {input.ContigsIn} --output {params.Output_Directory} --db-domain "prokaryotes" --db-version 223 --db-path {params.vContact3_db} \
        --threads {threads} --pyrodigal-gv && \
        touch {output.CheckVC3}
    '''
