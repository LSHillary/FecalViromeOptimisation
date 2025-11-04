###############################################################################
# Fecal Virome Optimization (FVO) Manuscript – Targets Pipeline
# Author: Luke Hillary
# Institution: University of California, Davis – Emerson Lab
# Repository: https://github.com/LSHillary/FecalViromeOptimisation
# License: [MIT / GPL-3 / CC-BY-4.0]
# Version: v1.0 (Version of Record for Zenodo DOI)
# Date: 2025-11-04
###############################################################################
# Description:
# This `_targets.R` file defines the reproducible analysis workflow used in the
# Fecal Virome Optimization (FVO) manuscript. It imports raw and processed data,
# executes preprocessing and diversity calculations, and produces key summary
# figures (K-mer profiles, Venn diagrams, stacked barplots, etc.).
#
# The workflow is divided into three main stages:
#   1. Data Import – reads metadata, coverage tables, taxonomy, and annotations
#   2. Data Processing – merges, normalizes, and summarizes analysis tables
#   3. Data Visualization – generates figures and summary plots
#
# All objects are cached by {targets} to ensure reproducible and incremental runs.
# Outputs are serialized as `.rds` for stability and downstream re-use.
# The objects are then used by Run_FVO_manuscript.R to produce the final figures
# and tables for the manuscript.
###############################################################################

# 0. Setup

# Load core package
library(targets)

# Source shared function definitions
source("scripts/Functions_FVO.R")

# Set pipeline options (package dependencies, output format)
tar_option_set(
  packages = c(
    "tidyverse", "vegan", "scales", "ggpubr",
    "ggVennDiagram", "pairwiseAdonis", "multcompView", "rcompanion"
  ),
  format = "rds"  # safer for ggplot and list outputs
)

# Log versions for reproducibility
message("Loaded {targets} version: ", packageVersion("targets"))

# Target pipeline definition

list(
  # 1. Data Import -------------------------------------------------------
  
  # Reads all required inputs (metadata, mapping summaries, taxonomy, etc.)
  # from the `/data` directory into tidy, analysis-ready data frames.
  
  # --- Metadata & QC
  tar_target(
    metadata,
    import_metadata("data/MetaData/FecalViromeOptimisationMetaData.csv"),
    description = "Sample metadata (IDs, treatments, and descriptive variables)."
  ),
  tar_target(
    raw_multiqc,
    import_multiqc("data/RawReads/multiqc_general_stats.txt", metadata),
    description = "MultiQC summary statistics merged with metadata."
  ),
  tar_target(
    kmers,
    read_and_merge_tsvs("data/Kmers", ".AbundanceTable.tsv"),
    description = "Merged per-sample K-mer abundance tables."
  ),
  
  # --- Read Mapping Summaries
  tar_target(
    human_read_counts,
    import_read_counts("data/ReadMapping/Human_coverage.tsv", "human"),
    description = "Read-pair counts mapping to the human genome."
  ),
  tar_target(
    viral_read_counts,
    import_read_counts("data/ReadMapping/Stool_Genomad_coverage_counts.tsv", "stool_Genomad"),
    description = "Read-pair counts mapping to viral contigs."
  ),
  tar_target(
    vOTU_tpm,
    import_tpm("data/ReadMapping/Stool_Genomad_coverage_tpm.tsv", "stool_Genomad"),
    description = "TPM (transcripts-per-million) values for viral contigs."
  ),
  
  # --- Viral Contig & Annotation Data
  tar_target(
    ClusteringData,
    import_vOTU_clusters("data/ViralContigs/combined_Genomad_FilteredContigs_clusters.tsv"),
    description = "Mapping of viral contigs to vOTU clusters."
  ),
  tar_target(
    GenomadData,
    read_and_merge_tsvs("data/ViralContigs", "_renamed_contigs_virus_summary.tsv"),
    description = "Merged geNomad viral summary tables."
  ),
  tar_target(
    PharokkaData,
    import_pharokka("data/Pharokka"),
    description = "Merged Pharokka protein annotations."
  ),
  
  # --- Host & Lifestyle Predictions
  tar_target(
    BacPhlipData_vOTU,
    import_bacphlip_vOTU("data/Lifestyle/combined_Genomad_vOTUs.fna.bacphlip"),
    description = "BACPHLIP lifestyle predictions (temperate vs virulent)."
  ),
  tar_target(
    iPhopData,
    import_iphop_data("data/HostPrediction/Host_prediction_to_genus_m90.csv"),
    description = "iPHoP host predictions at the genus level."
  ),
  
  # --- Microbial Profiling
  tar_target(
    SingleMData,
    import_singlem("data/MicrobialProfiling", metadata),
    description = "Phylum-level bacterial relative abundances from SingleM."
  ),
  
  # 2. Data Processing ---------------------------------------------------

  # Combines, cleans, and summarizes imported data into analysis-ready tables:
  #   - merges read counts
  #   - computes alpha/beta diversity
  #   - creates presence/absence tables
  
  # Process kmer data
  tar_target(
    processed_kmers,
    process_kmers(kmers, metadata),
    description = "Cleaned and metadata-merged K-mer tables."
  ),
  # Merge read counts
  tar_target(
    all_read_counts,
    merge_read_counts(viral_read_counts, human_read_counts, raw_multiqc),
    description = "Unified summary of viral, human, and total reads per sample."
  ),
  # Merge viral contig data
  tar_target(
    viral_contig_data,
    merge_viral_contig_data(GenomadData, ClusteringData, metadata),
    description = "Filtered viral contigs merged with clustering and metadata."
  ),
  # Calculate alpha diversity
  tar_target(
    alpha_diversity,
    calculate_alpha_diversity(vOTU_tpm, metadata),
    description = "Sample-level alpha diversity (vOTU richness)."
  ),
  # Produce a Presence/Absence table
  tar_target(
    PA_table,
    create_PA_table(vOTU_tpm, metadata, ClusteringData),
    description = "Presence/absence table for vOTUs across samples."
  ),
  
  # 3. Data Visualization -------------------------------------------------
  
  # Generates main manuscript figures and summary plots.
  # All ggplot objects are saved as `.rds` for portability and figure export.
  
  # Plot kmer data
  tar_target(
    plt_kmers,
    plot_Kmers(processed_kmers),
    description = "Log–log K-mer abundance curves by replicate."
  ),
  # Plot Venn diagrams - note that this was tweaked in Affinity to align labels
  # for the publication figure
  tar_target(
    plt_Venns,
    plot_venn_diagrams(PA_table),
    description = "Venn diagrams of assembled, mapped, and mapped-only vOTUs."
  ),
  # Plot stacked barplots of viral community composition
  tar_target(
    plt_stacked,
    plot_vOTU_stacked_barplot(vOTU_tpm, metadata, GenomadData),
    description = "Stacked barplots of viral community composition."
  )
)