rule PreClusteringAggregation:
    input:
        combined_contigs = expand("4-virus_identification/Genomad/{sample}_FilteredContigs.fna", sample=stool_samples),
    output:
        combined_contigs = "4-virus_identification/combined_Genomad_FilteredContigs.fna",
    threads: 12
    params:
        tag="combined_Genomad"
    resources:
        mem_mb = 16000,
        partition = "med2",
        time = "5-00:00:00"
    shell:'''
        cat {input.combined_contigs} > {output.combined_contigs} && \
        touch {output.CheckAggregation}
        '''
 
rule Clustering:
    input:
        contigs = "4-virus_identification/combined_Genomad_FilteredContigs.fna",
    output:
        Leiden_dereplicated = "5-dereplication/combined_Genomad_vOTUs.fna",
        CheckCluster = "Checks/LeidenClustering_Genomad_combined.done"
    threads: 48
    params:
        tag="combined"
    resources:
        mem_mb = 64000,
        partition = "med2",
        time = "5-00:00:00"
    shell:'''
        snakemake --nolock --config input={input.contigs} \
                   output={output.Leiden_dereplicated} \
                   leiden_resolution=1.0 \
                   min_ani=0.95 \
                   min_cov=0.85 \
                   blast_max_evalue=1e-5 \
                   blast_threads={threads} \
                   --cores {threads} -s ../scripts/vOTU_clustering/contig-ani-leiden-clustering-pipeline.smk && \
        touch {output.CheckCluster}
        '''

rule ExtractGenesAndProteins:
    input:
        GeNomadFilteredContigs = "4-virus_identification/Genomad/{sample}_FilteredContigs.fna",
        GenomadGenes = "4-virus_identification/Genomad/{sample}/{sample}_renamed_contigs_summary/{sample}_renamed_contigs_virus_genes.tsv",
        GenomadProteins = "4-virus_identification/Genomad/{sample}/{sample}_renamed_contigs_summary/{sample}_renamed_contigs_virus_proteins.faa"
    output:
        Filtered_proteins = "5-dereplication/{sample}_FilteredProteins.faa",
        Filtered_genes_tsv = "5-dereplication/{sample}_FilteredGenes.tsv",
        Filtered_genes_fasta = "5-dereplication/{sample}_FilteredGenes.fna",
        CheckProteins = "Checks/5.3-ExtractProteins_{sample}.done",
    params:
        tag = "{sample}"
    threads: 8
    resources:
        mem_mb = "32gb",
        time = "2-00:00:00",
        partition = "med2"
    run:
        import os
        from Bio import SeqIO
        import pandas as pd

        # Step 1: Extract all contig names from filtered contigs fasta file
        filtered_contigs = set()
        with open(input.GeNomadFilteredContigs, "r") as contigs_file:
            for record in SeqIO.parse(contigs_file, "fasta"):
                filtered_contigs.add(record.id)

        # Step 2: Add contig name column to genes file by removing the last underscore and following number from the 'gene' column
        genes_df = pd.read_csv(input.GenomadGenes, sep="\t")
        genes_df['contig_name'] = genes_df['gene'].str.rsplit('_', n=1).str[0]  # Corrected line

        # Step 3: Filter genes file for contigs present in the filtered contigs file
        filtered_genes_df = genes_df[genes_df['contig_name'].isin(filtered_contigs)]
        
        # Write filtered genes to output TSV file
        filtered_genes_df.to_csv(output.Filtered_genes_tsv, sep="\t", index=False)

        # Step 4: Filter proteins file for genes in the filtered genes set and write to FASTA
        filtered_genes = set(filtered_genes_df['gene'])
        with open(input.GenomadProteins, "r") as proteins_file, open(output.Filtered_proteins, "w") as filtered_proteins_file:
            for record in SeqIO.parse(proteins_file, "fasta"):
                if record.id in filtered_genes:
                    SeqIO.write(record, filtered_proteins_file, "fasta")

        # Step 5: Extract the nucleotide sequence for each gene from the filtered contig file
        contig_records = SeqIO.to_dict(SeqIO.parse(input.GeNomadFilteredContigs, "fasta"))

        with open(output.Filtered_genes_fasta, "w") as gene_fasta:
            for _, row in filtered_genes_df.iterrows():
                contig_id = row['contig_name']
                start_pos = int(row['start']) - 1  # Convert to 0-based indexing
                end_pos = int(row['end'])
                strand = row['strand']

                if contig_id in contig_records:
                    contig_seq = contig_records[contig_id].seq
                    gene_seq = contig_seq[start_pos:end_pos]

                    # Reverse complement if strand is -1
                    if strand == -1:
                        gene_seq = gene_seq.reverse_complement()

                    # Write gene sequence to output FASTA
                    gene_fasta.write(f">{row['gene']}\n{gene_seq}\n")

        # Step 6: Create a check file to signal successful completion
        with open(output.CheckProteins, "w") as check_file:
            check_file.write("ExtractGenesAndProteins completed successfully.\n")

import os
from Bio import SeqIO
import pandas as pd
from multiprocessing import Pool

def extract_votu_contigs(votu_fasta):
    """Extract contig names from the vOTUs FASTA file."""
    votu_contigs = set()
    with open(votu_fasta, "r") as votu_file:
        for record in SeqIO.parse(votu_file, "fasta"):
            votu_contigs.add(record.id)
    return votu_contigs

def process_sample(sample, votu_contigs):
    """Filter genes and proteins per sample."""
    # Read gene TSV and filter for vOTU contigs
    genes_df = pd.read_csv(f"5-dereplication/{sample}_FilteredGenes.tsv", sep="\t")
    filtered_genes_df = genes_df[genes_df['contig_name'].isin(votu_contigs)]

    # Filter proteins FASTA for genes in the filtered genes set
    proteins_fasta_in = f"5-dereplication/{sample}_FilteredProteins.faa"
    filtered_genes = set(filtered_genes_df['gene'])
    filtered_proteins = []
    with open(proteins_fasta_in, "r") as proteins_file:
        for record in SeqIO.parse(proteins_file, "fasta"):
            if record.id in filtered_genes:
                filtered_proteins.append(record)

    # Filter gene FASTA file for genes in filtered genes set
    genes_fasta_in = f"5-dereplication/{sample}_FilteredGenes.fna"
    filtered_gene_sequences = []
    with open(genes_fasta_in, "r") as genes_fasta_file:
        for record in SeqIO.parse(genes_fasta_file, "fasta"):
            if record.id in filtered_genes:
                filtered_gene_sequences.append(record)

    return filtered_genes_df, filtered_proteins, filtered_gene_sequences

def write_merged_output(output, filtered_genes_dfs, filtered_proteins, filtered_gene_sequences):
    """Write the final merged genes and proteins files."""
    # Merge all filtered genes TSV files and write to output
    merged_genes_df = pd.concat(filtered_genes_dfs)
    merged_genes_df.to_csv(output.vOTU_Genes, sep="\t", index=False)

    # Merge and write filtered proteins to output FASTA
    with open(output.vOTU_Proteins, "w") as proteins_out:
        for record in filtered_proteins:
            SeqIO.write(record, proteins_out, "fasta")

    # Merge and write filtered gene sequences to output FASTA
    with open(output.vOTU_GenesFasta, "w") as genes_fasta_out:
        for record in filtered_gene_sequences:
            SeqIO.write(record, genes_fasta_out, "fasta")

rule Filter_vOTU_Proteins:
    input:
        ProteinsIn = expand("5-dereplication/{sample}_FilteredProteins.faa", sample=stool_samples),
        GenesIn = expand("5-dereplication/{sample}_FilteredGenes.tsv", sample=stool_samples),
        FilteredGenesFasta = expand("5-dereplication/{sample}_FilteredGenes.fna", sample=stool_samples),
        vOTUsIn = "5-dereplication/combined_Genomad_vOTUs.fna"
    output:
        vOTU_Proteins = "5-dereplication/combined_vOTUs_proteins.faa",
        vOTU_Genes = "5-dereplication/combined_vOTUs_genes.tsv",
        vOTU_GenesFasta = "5-dereplication/combined_vOTUs_genes.fna",
        Check_vOTU_Filtering = "Checks/5.3-ExtractProteins_combined_vOTUs.done"
    params:
        tag = "vOTUs"
    threads: 8
    resources:
        mem_mb = "32gb",
        time = "2-00:00:00",
        partition = "med2"
    run:
        # Step 1: Extract all contig names from vOTUs file
        votu_contigs = extract_votu_contigs(input.vOTUsIn)

        # Step 2: Process each sample in parallel
        pool = Pool(threads)
        results = pool.starmap(process_sample, [(sample, votu_contigs) for sample in stool_samples])
        pool.close()
        pool.join()

        # Unpack the results
        filtered_genes_dfs = [result[0] for result in results]
        filtered_proteins = [protein for result in results for protein in result[1]]
        filtered_gene_sequences = [gene_seq for result in results for gene_seq in result[2]]

        # Step 3: Write merged results to output files
        write_merged_output(output, filtered_genes_dfs, filtered_proteins, filtered_gene_sequences)

        # Step 4: Create the check file indicating success
        with open(output.Check_vOTU_Filtering, "w") as check_file:
            check_file.write("Filter_vOTU_Proteins completed successfully.\n")