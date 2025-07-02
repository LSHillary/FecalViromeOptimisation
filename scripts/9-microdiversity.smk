rule update_gene_headers:
    input:
        genes_fna="5-dereplication/combined_vOTUs_genes.fna",  # Input gene sequences
        proteins_faa="5-dereplication/combined_vOTUs_proteins.faa"  # Input protein sequences
    output:
        updated_genes_fna="5-dereplication/combined_vOTUs_genes_with_headers.fna",  # Output updated gene sequences
        CheckUpdateGeneHeaders="Checks/5.3-UpdateGeneHeaders.done"  # Output check file
    params:
        tag="combined"
    run:
        # Python code to update gene headers
        genes_fna = input.genes_fna
        proteins_faa = input.proteins_faa
        updated_genes_fna = output.updated_genes_fna

        # Function to parse the protein headers and store them in a dictionary
        def parse_protein_headers(proteins_faa):
            protein_headers = {}
            with open(proteins_faa, 'r') as prot_file:
                for line in prot_file:
                    if line.startswith('>'):
                        # Extract the gene identifier (e.g., "GP5_S30_contig_5330_13") and the rest of the header
                        gene_id = line.split()[0][1:]  # Remove ">" and extract the first part
                        additional_info = " ".join(line.split()[1:])  # Extract everything after the first space
                        protein_headers[gene_id] = additional_info
            return protein_headers

        # Function to update the gene headers
        def update_gene_headers(genes_fna, protein_headers, output_file):
            with open(genes_fna, 'r') as genes_file, open(output_file, 'w') as out_file:
                for line in genes_file:
                    if line.startswith('>'):
                        gene_id = line.strip()[1:]  # Extract the gene identifier
                        if gene_id in protein_headers:
                            # Write the new header with additional information from the protein file
                            out_file.write(f">{gene_id} {protein_headers[gene_id]}\n")
                        else:
                            # If no match is found, just write the original header
                            out_file.write(line)
                    else:
                        # Write the sequence lines as they are
                        out_file.write(line)

        # Parse protein headers
        protein_headers = parse_protein_headers(proteins_faa)

        # Update gene headers
        update_gene_headers(genes_fna, protein_headers, updated_genes_fna)

        # Create the check file
        with open(output.CheckUpdateGeneHeaders, 'w') as check_file:
            check_file.write("")

rule instrain_profile:
    input:
        bam = "6-Mapping/{sample}_stool_Genomad_mapped.bam",
        CheckUpdateGeneHeaders="Checks/5.3-UpdateGeneHeaders.done"
    output:
        Check = "Checks/instrain_profile.{sample}.done"
    params:
        tag = "{sample}",
        vOTU_fasta = "5-dereplication/combined_Genomad_vOTUs.fna",
        genes_fasta = "5-dereplication/combined_vOTUs_genes_with_headers.fna",
    threads: 12
    resources:
        mem_mb = "128gb",
        partition = "bmh",
        time = "96:00:00"
    shell: '''
        mkdir -p 9-microdiversity && \
        micromamba run -n inStrain_env inStrain profile {input.bam} {params.vOTU_fasta} -o {params.tag} -p {threads} \
        --min_read_ani 0.95 --gene_file {params.genes_fasta} && \
        touch {output.Check}
    '''

rule instrain_compare:
    input:
        bams = expand("6-Mapping/{sample}_combined_mapped.bam", sample = stool_samples),
        ProfileCheck = expand("Checks/instrain_profile.{sample}.done", sample = stool_samples),
    output:
        Check = "Checks/instrain_compare.done"
    params:
        tag = "combined",
        profiles = expand("9-microdiversity/{sample}_instrain_profile", sample = stool_samples)
    threads: 16
    resources:
        mem_mb = 64000,
        partition = "bmh",
        time = "96:00:00"
    shell: '''
        micromamba run -n instrain_env inStrain compare -i {params.profiles} -o {sample} -p {threads} \
        --breadth 0.75 -bams  && \
        mv {sample}_compare 9-microdiversity/{sample}_instrain_compare && \
        touch {output.Check}
    '''