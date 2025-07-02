#### 3.4 - Assembly cleanup ####

# Detect which assembly strategy was used. Options include
# 1 - individual assemblies for each sample - this is the default: identified by folder 3-Assembly/Individual
# 2 - a single CoAssembly for all samples: identified by folder 3-Assembly/CoAssembly
# 3 - a single CoAssembly of digitally normalized reads for all samples: identified by folder 3-Assembly/DigiNormCoAssembly

# Detect which assembly strategy was used and list all contig files
assembly_list = []
folders = ['Individual']
existing_folders = [folder for folder in folders if os.path.exists(f'3-Assembly/{folder}')]
for folder in existing_folders:
    assembly_list.extend(glob.glob(f'3-Assembly/{folder}/*_renamed_.contigs.fa'))

# Run metaQuast on all identified assemblies with contigs >1000 bp
rule Quast:
    input:
        assembly <- assembly_list
    output:
        report="3-Assembly/MetaQuast/" # You can specify the desired output location here
    shell:'''
        metaquast {input} -o {output.report} --min-contig 1000
        '''

# Rule to rename individual sample contigs
rule rename_individual_contigs:
    input:
        contigsIn="3-Assembly/Individual/{sample}.contigs.fa"
    output:
        contigsOut="3-Assembly/Individual/{sample}_renamed.contigs.fa",
        contig_names="3-Assembly/Individual/{sample}_contig_names.csv",
        CheckRenamed="Checks/{sample}_renamed_done.txt"
    params:
        sname="{sample}"
    threads: 8
    resources:
        mem_mb=24000,
        partition="high2"
    shell:'''
        awk '/^>/{{original=$0; print ">" "{params.sname}" "_contig_" ++i; newname=$0; print newname >> "{output.contig_names}"; print original","newname >> "{output.contig_names}"; next}} {{print}}' {input.contigsIn} > {output.contigsOut} && \
        touch {output.CheckRenamed}
        '''

# Rule to concatenate all the renamed contigs
rule concatenate_contigs:
    input:
        individual_contigs=expand("3-Assembly/Individual/{sample}_renamed.contigs.fa", sample=samples) if 'Individual' in existing_folders else [],
        coassembly_contigs=expand("3-Assembly/{folder}/coassembly_renamed.contigs.fa", folder=[folder for folder in existing_folders if folder != 'Individual']),
        csv_files=expand("3-Assembly/{folder}/{type}_contig_names.csv", folder=existing_folders, type=['{sample}', 'coassembly']) if 'Individual' in existing_folders else expand("3-Assembly/{folder}/coassembly_contig_names.csv", folder=[folder for folder in existing_folders if folder != 'Individual'])

    output:
        contigsOut="3-Assembly/{folder}/coassembly_renamed.contigs.fa",
        contig_names="3-Assembly/{folder}/coassembly_contig_names.csv",
        CheckRenamed="Checks/coassembly_renamed_done.txt"
    params:
        prefix=lambda wildcards: 'CoAssembly_' if wildcards.folder == 'CoAssembly' else 'DigiNorm_CoAssembly_'
    threads: 8
    resources:
        mem_mb=24000,
        partition="high2"
    shell:'''
        awk '/^>/{{original=$0; print ">" "{params.prefix}" "contig_" ++i; newname=$0; print newname >> "{output.contig_names}"; print original","newname >> "{output.contig_names}"; next}} {{print}}' {input.contigsIn} > {output.contigsOut} && \
        touch {output.CheckRenamed}
        '''

# Rule to concatenate all the renamed contigs
rule concatenate_contigs:
    input:
        individual_contigs=expand("3-Assembly/Individual/{sample}_renamed.contigs.fa", sample=samples),
        coassembly_contigs=expand("3-Assembly/{folder}/coassembly_renamed.contigs.fa", folder=['CoAssembly', 'DigiNormCoAssembly']),
        csv_files=expand("3-Assembly/{folder}/{type}_contig_names.csv", folder=folders, type=['{sample}', 'coassembly']),
    output:
        contigsOut="3-Assembly/individual/all_individual_assembly_contigs.fa",
        final_csv="3-Assembly/all_contig_names.csv",
        CheckCat="Checks/individual_concatenate_contigs_done.txt"
    shell:'''
        cat {input.individual_contigs} {input.coassembly_contigs} > {output.contigsOut} && \
        echo "Original_Contig_Name,New_Contig_Name" > {output.final_csv} && \
        cat {input.csv_files} >> {output.final_csv} && \
        rm {input.individual_contigs} {input.coassembly_contigs} && \
        touch {output.CheckCat}
        '''

