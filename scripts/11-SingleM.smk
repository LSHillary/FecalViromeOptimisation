rule SingleM:
    input:
        FwdReads = "0-raw/{sample}_raw_R1.fq.gz",
        RevReads = "0-raw/{sample}_raw_R2.fq.gz",
    output:
        OtuTable = "11-MicrobialProfiling/{sample}_otu_table.tsv",
        TaxProfile = "11-MicrobialProfiling/{sample}_tax_profile.tsv",
        CheckSingleM = "Checks/11-SingleM_{sample}.done"
    threads: 32
    params:
        tag = "{sample}",
    resources:
        mem_mb = "64gb",
        partition = "high2",
        time = "2-24:00:00"
    shell:'''
        mkdir -p 11-MicrobialProfiling && \
        micromamba run -n SINGLEM_ENV singlem pipe -1 {input.FwdReads} -2 {input.RevReads} \
        -p {output.TaxProfile} --otu-table {output.OtuTable} --threads {threads} && \
        touch {output.CheckSingleM}
        '''

rule SingleM_Summarise:
    input:
        TaxProfile = "11-MicrobialProfiling/{sample}_tax_profile.tsv"
    output:
        SummarisedTaxProfile = "11-MicrobialProfiling/{sample}_summary.tsv",
        CheckSingleM_Summarise = "Checks/11-SingleM_{sample}_summarised.done"
    threads: 32
    params:
        tag = "{sample}",
    resources:
        partition = "high2",
        time = "2-24:00:00"
    shell:'''
        micromamba run -n SINGLEM_ENV singlem summarise --input-taxonomic-profile {input.TaxProfile} \
        --output-taxonomic-profile-with-extras {output.SummarisedTaxProfile} && \
        touch {output.CheckSingleM_Summarise}
    '''