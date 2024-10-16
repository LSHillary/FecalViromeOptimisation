# Load packages
library(targets)

# Source R functions
source("scripts/Functions_FVO.R")

# Set target-specific options (optional)
tar_option_set(
  packages = c("tidyverse", "scales"),
  format = "qs"
)

# Define the pipeline (targets)
list(
  # 1 - Data Import ----
  tar_target(
    metadata,
    import_metadata("data/MetaData/FecalViromeOptimisationMetaData.csv")
  ),
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
  # 2 - Data Processing----
  tar_target(
    processed_kmers,
    process_kmers(kmers, metadata)  
  ),
  tar_target(
    all_read_counts,
    merge_read_counts(viral_read_counts, human_read_counts, raw_multiqc)
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
  )
)