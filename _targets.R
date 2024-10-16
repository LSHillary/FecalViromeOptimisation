# Load packages
library(targets)

# Source R functions
source("scripts/Functions_FVO.R")

# Set target-specific options (optional)
tar_option_set(
  packages = c("tidyverse", "scales","ggVennDiagram", "ggpubr", "vegan"),
  format = "qs"
)

# Define the pipeline (targets)
list(
  # 1 - Data Import ----
  tar_target(
  ## Metadata
    metadata,
    import_metadata("data/MetaData/FecalViromeOptimisationMetaData.csv")
  ),
  ## Read data
  tar_target(
    raw_multiqc,
    import_multiqc("data/RawReads/multiqc_general_stats.txt", metadata)
  ),
  tar_target(
    kmers,
    read_and_merge_tsvs("data/Kmers", ".AbundanceTable.tsv")
  ),
  tar_target(
    human_read_counts,
    import_read_counts("data/ReadMapping/Human_coverage.tsv", "human")
  ),
  tar_target(
    viral_read_counts,
    import_read_counts("data/ReadMapping/Stool_Genomad_coverage_counts.tsv", "stool_Genomad")
  ),
  tar_target(
    vOTU_tpm,
    import_tpm("data/ReadMapping/Stool_Genomad_coverage_tpm.tsv", "stool_Genomad")
  ),
  ## vOTU Cluster data
  tar_target(
    ClusteringData,
    import_vOTU_clusters("data/ViralContigs/combined_Genomad_FilteredContigs_clusters.tsv")
  ),
  tar_target(
    GenomadData,
    read_and_merge_tsvs("data/ViralContigs", "_renamed_contigs_virus_summary.tsv")
  ),
  tar_target(
    BacPhlipData,
    import_bacphlip("data/ViralContigs/combined_Genomad_vOTUs.fna.bacphlip")
  ),
  tar_target(
    DefenseFinderData,
    import_defensefinder("data/DefenseFinder/combined_vOTUs_proteins_defense_finder_systems.tsv")
  ),
  tar_target(
    iPhopData,
    import_iphop_data("data/HostPrediction/Host_prediction_to_genus_m90.csv")
  ),
  # 2 - Data Processing----
  tar_target(
    processed_kmers,
    process_kmers(kmers, metadata)  
  ),
  tar_target(
    all_read_counts,
    merge_read_counts(viral_read_counts, human_read_counts, raw_multiqc)
  ),
  tar_target(
    vira_contig_data,
    merge_viral_contig_data(GenomadData, ClusteringData, metadata)
  ),
  tar_target(
    alpha_diversity,
    calculate_alpha_diversity(vOTU_tpm, metadata)
  ),
  tar_target(
    PA_table,
    create_PA_table(vOTU_tpm, metadata, ClusteringData)
  ),
  # 3 - Data Visualisation ----
  tar_target(
    plt_raw_reads,
    plot_raw_reads(raw_multiqc)
  ),
  tar_target(
    plt_kmers,
    plot_Kmers(processed_kmers)
  ),
  tar_target(
    plt_read_percentages,
    plot_read_percentages(all_read_counts)
  ),
  tar_target(
    plt_length,
    plot_contig_length(vira_contig_data)
  ),
  tar_target(
    plt_alpha_diversity,
    plot_alpha_diversity(alpha_diversity)
  ),
  tar_target(
    plt_Venns,
    plot_venn_diagrams(PA_table)
  ),
  tar_target(
    plt_stacked,
    plot_vOTU_stacked_barplot(vOTU_tpm, metadata, GenomadData)
  )
)

