# Define the rule for running GeNomad
rule GeNomad:
    input:
        CheckMegahitIndividual = "Checks/{sample}_IndividualAssembly.done",
        CheckRenameContigs = "Checks/{sample}_RenameContigs.done",
        ContigsIn = "3-Assembly/Individual/{sample}_renamed_contigs.fa"
    output:
        CheckGenomad = "Checks/4.1_Genomad_{sample}.done"
    threads: 16
    params:
        tag = "{sample}",
        GenomadDB = "/group/jbemersogrp/databases/genomad/genomad_db",
        GenomadFolder = "4-virus_identification/Genomad/{sample}"
    resources:
        mem_mb = 64000,
        partition = "high2"
    shell:'''
        mkdir -p {params.GenomadFolder} && \
        micromamba run -n genomad_env genomad end-to-end --cleanup --enable-score-calibration -t {threads} {input.ContigsIn} \
        {params.GenomadFolder} {params.GenomadDB} &&
        touch {output.CheckGenomad}
        '''

rule FilterGenomad:
    input:
        CheckGenomad = "Checks/4.1_Genomad_{sample}.done"
    output:
        FilteredGenomadResults = "4-virus_identification/Genomad/{sample}/{sample}_FilteredGenomadResults.tsv",
        CheckFilteredGenomadResults = "Checks/4.2-FilterGenomadResults_{sample}.done",
        FilteredContigs = "4-virus_identification/Genomad/{sample}_FilteredContigs.fna"
    resources:
        partition = "med2"
    params:
        tag = "{sample}",
        Contigs = "4-virus_identification/Genomad/{sample}/{sample}_renamed_contigs_summary/{sample}_renamed_contigs_virus.fna",
        GenomadResults = "4-virus_identification/Genomad/{sample}/{sample}_renamed_contigs_summary/{sample}_renamed_contigs_virus_summary.tsv"
    shell:'''
        awk -F'\t' 'NR==1 || $2 >= 10000 || ($2 >= 1000 && $11 ~ /Monodnaviria/)' {params.GenomadResults} | \
        cut -f1 > {output.FilteredGenomadResults} && \
        seqtk subseq {params.Contigs} {output.FilteredGenomadResults} > {output.FilteredContigs} && \
        touch {output.CheckFilteredGenomadResults}
        '''