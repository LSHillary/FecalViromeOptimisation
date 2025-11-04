###############################################################################
# Fecal Virome Optimization (FVO) Manuscript – Core Functions
# Author: Luke Hillary
# Institution: University of California, Davis – Emerson Lab
# Repository: https://github.com/LSHillary/FecalViromeOptimisation
# License: MIT
# Version: v1.0 (Version of Record for Zenodo DOI)
# Date: 2025-11-04
###############################################################################
# Description:
# This script defines the core R functions used in the Fecal Virome Optimization
# (FVO) manuscript pipeline. It supports all targets and figure-generation steps
# executed by the `_targets.R` and `Run_FVO_manuscript.R` scripts.
#
# The functions are grouped into the following sections:
#   0. Parameter definitions (labels, levels, color palettes)
#   1. Data import utilities for metadata, coverage, taxonomy, and annotation
#   2. Data processing functions for read merging, diversity calculations,
#      and table generation (alpha/beta diversity, presence–absence)
#   3. Visualization utilities for figure generation (e.g., stacked barplots,
#      PCoA ordinations, Venn diagrams, K-mer profiles)
#
# Key design principles:
#   • Modular and self-contained – each function is callable independently
#   • Targets-compatible – all outputs serialize cleanly within {targets}
#   • Reproducible – parameterized filepaths, consistent factor levels, and
#     fixed palettes to ensure reproducible figures
#   • Minimal dependencies – tidyverse core + vegan + ggpubr + ggVennDiagram
#
# Usage:
#   Source this file in `_targets.R`:
#     tar_source("Functions_FVO.R")
#
#   or interactively:
#     source("Functions_FVO.R")
#
# Notes:
#   - Designed for R ≥ 4.3 and tidyverse ≥ 2.0
#   - All ggplot objects should be saved with `format = "rds"` in {targets}
#     to avoid serialization errors.
###############################################################################


# 0. Parameter settings ----

# Define consistent labels and factor levels for sample descriptions used 
# throughout the Fecal Virome Optimization (FVO) analyses and figures.
# These are used for axis labels, legend text, and ordered facets.

# Mapping of short sample codes (used in filenames/metadata) to human-readable labels
DescriptionModifiers <- c(
  "UTV" = "Untreated\nVirome",
  "DTV" = "DNase\nTreated\nVirome",
  "MDA" = "MDA\nVirome",
  "MtG" = "Meta\ngenome"
)

# Define the display order of descriptions in all plots and summaries
DescriptionLevels <- c(
  "Untreated\nVirome",
  "DNase\nTreated\nVirome",
  "MDA\nVirome",
  "Meta\ngenome"
)



# 1. Data Import Functions ----

# These functions read and preprocess external data files (metadata, MultiQC,
# mapping summaries, TPM tables, taxonomy, and host/lifestyle predictions)
# into tidy data frames for downstream analysis.
# All assume UTF-8 CSV/TSV files and return tidyverse-style tibbles.

# Read metadata table
# Arguments:
#   filepath — path to metadata CSV
# Returns:
#   Data frame of sample metadata.
import_metadata <- function(filepath) {
  read.csv(filepath, sep = ",", header = TRUE)
}

# Read MultiQC FastQC summary
# Cleans column names, removes R2 reads, and merges with metadata.
# Arguments:
#   filepath — path to MultiQC general stats TSV
#   metadata — data frame of sample metadata
# Returns:
#   Data frame of cleaned MultiQC stats merged with metadata.
import_multiqc <- function(filepath, metadata) {
  read.csv(filepath, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    rename_with(~ gsub("FastQC_mqc\\.generalstats\\.fastqc\\.", "", .x)) %>%
    filter(!grepl("R2", Sample)) %>%
    mutate(Sample = gsub("_raw_R1", "", Sample)) %>%
    select(Sample, percent_gc, total_sequences) %>%
    merge(metadata, by = "Sample", all.x = TRUE)
}

# Read and merge multiple TSV file
# Reads all files matching a pattern within a folder and combines them.
# Arguments:
#   folder — directory path
#   file_extension — regex pattern (e.g. "_tax_profile.tsv$")
# Returns:
#   Data frame combining all TSVs with sample identifiers.
read_and_merge_tsvs <- function(folder, file_extension) {
  files <- list.files(path = folder, pattern = file_extension)
  df_list <- lapply(files, function(f) {
    df <- read_tsv(file.path(folder, f), show_col_types = FALSE)
    df$sample <- gsub(file_extension, "", f)
    df
  })
  bind_rows(df_list)
}

# Import read-count summaries
# Sums mapped read counts across contigs and converts to read pairs.
# Arguments:
#   filename — path to read-count table
#   TEXT — descriptor (e.g. "human" or "viral")
# Returns:
#   Data frame of total read pairs per sample.
import_read_counts <- function(filename, TEXT) {
  read.csv(filename, sep = "\t") %>%
    select(contains("Read.Count")) %>%
    rename_with(~ gsub(paste0("_", TEXT, "_mapped\\.Read\\.Count"), "", .x)) %>%
    summarise(across(everything(), sum)) %>%
    t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>%
    rename(Reads = V1) %>%
    mutate(Reads = Reads / 2)
}

# Import TPM tables
# Reads TPM values for each contig/sample.
# Arguments:
#   filename — path to TPM table
#   TEXT — descriptor (e.g. "human" or "viral")
# Returns:
#   Data frame of TPM values per contig/sample.
import_tpm <- function(filename, TEXT) {
  read.csv(filename, sep = "\t") %>%
    select(contains("TPM"), contains("Contig")) %>%
    rename_with(~ gsub(paste0("_", TEXT, "_mapped\\.TPM"), "", .x))
}

# Import vOTU cluster assignments
# Reads vOTU-to-cluster mapping from clustering output.
# Arguments:
#   filename — path to vOTU cluster file
# Returns:
#   Data frame of vOTU and corresponding cluster.
import_vOTU_clusters <- function(filename) {
  read.csv(filename, sep = "\t", header = FALSE) %>%
    rename(vOTU = V1, Cluster = V2)
}

# Import SingleM outputs
# Aggregates coverage to phylum level and merges with metadata.
# Arguments:
#   folder — directory path containing SingleM TSVs
#   metadata — data frame of sample metadata
# Returns:
#   Data frame of phylum-level relative abundances per sample.
import_singlem <- function(folder, metadata) {
  read_and_merge_tsvs(folder, "_tax_profile.tsv") %>%
    separate(taxonomy, into = c("Root","Domain","Phylum","Class","Order",
                                "Family","Genus","Species"), sep = ";", fill = "right") %>%
    mutate(
      Phylum = ifelse(is.na(Phylum), paste("Unclassified", Domain), Phylum),
      Phylum = gsub("^d__|^p__|\\s", "", Phylum),
      Phylum = ifelse(grepl("Unclassified", Phylum), "Unclassified", Phylum),
      Phylum = ifelse(grepl("Bacillota", Phylum), "Bacillota", Phylum)
    ) %>%
    group_by(sample, Phylum) %>%
    summarise(coverage = sum(coverage), .groups = "drop_last") %>%
    mutate(RelAbund = coverage / sum(coverage) * 100) %>%
    ungroup() %>%
    rename(Sample = sample) %>%
    merge(metadata, by = "Sample")
}

# Import BACPHLIP lifestyle predictions
# Adds Temperate/Virulent/Unclassified status and marks proviruses.
# Arguments:
#   filepath — path to BACPHLIP output file
# Returns:
#   Data frame of vOTU lifestyle predictions
import_bacphlip_vOTU <- function(filepath) {
  read.csv(filepath, sep = "\t", header = TRUE) %>%
    rename(vOTU = X) %>%
    mutate(
      status = case_when(
        Temperate >= 0.95 ~ "Temperate",
        Virulent >= 0.95  ~ "Virulent",
        TRUE              ~ "Unclassified"
      ),
      Provirus = ifelse(grepl("provirus", vOTU, ignore.case = TRUE), "Provirus", "Virus")
    )
}

# Import Pharokka annotation tables
# Reads and merges Pharokka protein annotations.
# Arguments:
#   folder — directory path containing Pharokka TSVs
# Returns:
#   Data frame of Pharokka protein annotations with contig IDs.
import_pharokka <- function(folder) {
  read_and_merge_tsvs(folder, "_pharokka_proteins_full_merged_output.tsv") %>%
    mutate(Contig = gsub("_[0-9]*$", "", ID))
}

# Import iPHoP host predictions
# Cleans and resolves taxonomy, keeping top-confidence hits only.
# Arguments:
#   filename — path to iPHoP host prediction CSV
# Returns:
#   Data frame of vOTU host predictions with cleaned taxonomy.
import_iphop_data <- function(filename) {
  read.csv(filename) %>%
    rename(vOTU = Virus) %>%
    group_by(vOTU) %>%
    filter(Confidence.score == max(Confidence.score)) %>%
    separate(Host.genus,
             into = c("HostDomain","HostPhylum","HostClass","HostOrder",
                      "HostFamily","HostGenus"),
             sep = ";", fill = "right") %>%
    mutate(across(starts_with("Host"), ~ gsub("^[a-z]__", "", .x))) %>%
    mutate(
      HostPhylum = ifelse(grepl("Bacillota", HostPhylum), "Bacillota", HostPhylum),
      HostGenus  = ifelse(grepl("Clostridium", HostGenus), "Clostridium", HostGenus)
    ) %>%
    group_by(vOTU, Confidence.score) %>%
    mutate(across(starts_with("Host"),
                  ~ ifelse(n_distinct(.x) > 1, NA, .x))) %>%
    ungroup() %>%
    select(-List.of.methods) %>%
    distinct()
}


# 2. Data Processing Functions ----

# These functions perform internal data transformations used in the FVO
# analysis pipeline, including read-merging, normalization, and calculation of
# diversity and presence–absence metrics.

# Process K-mer data
# Cleans and merges K-mer abundance data with metadata.
# Arguments:
#   df_Kmers     — data frame of K-mer abundance and counts
#   df_metadata  — data frame of sample metadata
# Returns:
#   Data frame of merged K-mer abundances with metadata annotations.
process_kmers <- function(df_Kmers, df_metadata) {
  df_Kmers %>%
    separate(
      col = colnames(df_Kmers)[1],
      into = strsplit(colnames(df_Kmers)[1], ",")[[1]],
      sep = ","
    ) %>%
    mutate(
      abundance = as.numeric(abundance),
      count = as.numeric(count)
    ) %>%
    merge(df_metadata, by.x = "sample", by.y = "Sample", all.x = TRUE) %>%
    filter(Description != "Soil Virome") %>%
    rename(Sample = sample)
}

# Merge read-count summaries
# Combines viral, human, and raw MultiQC data; computes read percentages.
# Arguments:
#   viral_read_counts — data frame of viral read counts
#   human_read_counts — data frame of human read counts
#   raw_multiqc       — data frame of MultiQC read statistics
# Returns:
#   Long-format data frame of read-type percentages per sample.
merge_read_counts <- function(viral_read_counts, human_read_counts, raw_multiqc) {
  viral_read_counts <- viral_read_counts %>% rename(ViralReads = Reads)
  human_read_counts <- human_read_counts %>% rename(HumanReads = Reads)
  
  raw_multiqc %>%
    merge(human_read_counts, by = "Sample", all.x = TRUE) %>%
    merge(viral_read_counts, by = "Sample", all.x = TRUE) %>%
    mutate(
      HumanPercent = HumanReads / total_sequences * 100,
      ViralPercent = ViralReads / total_sequences * 100
    ) %>%
    pivot_longer(
      cols = c("HumanPercent", "rRNApercentage", "ViralPercent"),
      names_to = "Type", values_to = "ReadPercentage"
    ) %>%
    mutate(
      Type = case_when(
        Type == "HumanPercent" ~ "Human",
        Type == "rRNApercentage" ~ "rRNA",
        Type == "ViralPercent" ~ "Viral"
      )
    ) %>%
    filter(Description != "Soil Virome")
}

# Merge viral-contig data
# Filters GeNomad annotations for viral contigs and merges metadata.
# Arguments:
#   GenomadData   — data frame of annotated contigs
#   ClusteringData— vOTU cluster mapping
#   metadata       — data frame of sample metadata
# Returns:
#   Annotated data frame of viral contigs with metadata and taxonomy.
merge_viral_contig_data <- function(GenomadData, ClusteringData, metadata) {
  
  # Extract all viral contig IDs
  viral_contigs <- ClusteringData$Cluster %>%
    strsplit(",") %>% unlist()
  
  GenomadData %>%
    rename(Contig = seq_name, Sample = sample) %>%
    filter(Contig %in% viral_contigs) %>%
    merge(metadata, by = "Sample", all.x = TRUE) %>%
    mutate(
      DescriptionLong = DescriptionModifiers[Description],
      DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels)
    ) %>%
    separate(
      taxonomy,
      into = c("Viruses","Realm","Kingdom","Phylum","Class","Order","Family"),
      sep = ";", fill = "right"
    ) %>%
    filter(Contig %in% viral_contigs)
}

# Calculate alpha diversity
# Computes viral richness (vOTU count per sample).
# Arguments:
#   df_mapping — TPM matrix (rows = contigs, columns = samples)
#   metadata   — data frame of sample metadata
# Returns:
#   Data frame of richness values merged with metadata.
calculate_alpha_diversity <- function(df_mapping, metadata) {
  df_mapping %>%
    pivot_longer(-Contig, names_to = "Sample", values_to = "TPM") %>%
    filter(TPM > 0) %>%
    group_by(Sample) %>%
    summarise(Richness = n_distinct(Contig), .groups = "drop") %>%
    merge(metadata, by = "Sample", all.x = TRUE)
}

# Create presence–absence table
# Generates binary (0/1) presence–absence calls per vOTU and sample.
# Adds an assembled-presence flag based on vOTU cluster membership.
# Arguments:
#   df           — TPM matrix (rows = contigs, columns = samples)
#   metadata     — data frame of sample metadata
#   df_Clustering— vOTU cluster mapping
# Returns:
#   Data frame with vOTU, sample, PA, and assembled-presence indicators.
create_PA_table <- function(df, metadata, df_Clustering) {
  df %>%
    pivot_longer(-Contig, names_to = "Sample", values_to = "TPM") %>%
    mutate(PA = ifelse(TPM > 0, 1, 0)) %>%
    select(-TPM) %>%
    merge(metadata, by = "Sample", all.x = TRUE) %>%
    rename(vOTU = Contig) %>%
    merge(df_Clustering, by = "vOTU", all.x = TRUE) %>%
    mutate(
      AssembledPresence = mapply(function(sample, cluster)
        grepl(sample, cluster, fixed = TRUE), Sample, Cluster),
      AssembledPresence = as.integer(AssembledPresence)
    )
}

# Calculate beta diversity
# Removes singleton contigs (present in only one sample) and renormalizes TPM.
# Arguments:
#   df — TPM matrix (rows = contigs, columns = samples)
# Returns:
#   Normalized TPM matrix excluding singleton contigs.
calculate_beta_diversity <- function(df) {
  
  # Remove contigs present in only one sample
  df_filtered <- df %>%
    rowwise() %>%
    filter(sum(c_across(-Contig) > 0) > 1) %>%
    ungroup()
  
  # Renormalize TPM values to sum to 1 × 10⁶ per sample
  df_filtered %>%
    mutate(across(-Contig, ~ .x / sum(.x) * 1e6))
}

# 3. Data Visualisation Functions ----

# These functions generate all publication-quality figures for the FVO manuscript,
# including K-mer profiles, barplots with compact-letter displays, PCoA ordinations,
# Venn diagrams, and stacked taxonomic barplots.
# All outputs are ggplot objects that can be serialized in {targets}.

# 3.1  Color Palettes                                                          #

# Replicate-level palette (used in K-mer plots)
Replicates_palette <- c(
  "-DV1" = "#56B4E9", "-DV2" = "#56B4E9", "-DV3" = "#56B4E9",
  "+DV1" = "#0072B2", "+DV2" = "#0072B2", "+DV3" = "#0072B2",
  "MDA1" = "#009E73", "MDA2" = "#009E73", "MDA3" = "#009E73",
  "MtG1" = "#FF7F0E", "MtG2" = "#FF7F0E", "MtG3" = "#FF7F0E"
)

# Sample-type palette (used across most figures)
Description_palette <- c(
  "UTV" = "#56B4E9",
  "DTV" = "#0072B2",
  "MDA" = "#009E73",
  "MtG" = "#FF7F0E"
)

# Baltimore-class / lifestyle palette (stacked barplots)
Stacked_palette <- c(
  "dsDNA" = "#0072B2", "ssDNA" = "#56B4E9", "Unassigned" = "#808080",
  "Other Caudoviricetes" = "#56B4E9", "Crassvirales" = "#0072B2",
  "Microviridae" = "#D55E00", "Other Monodnaviria" = "#E69F00",
  "Herpesviridae" = "#009E73",
  "Temperate" = "#56B4E9", "Virulent" = "#0072B2",
  "Provirus" = "#56B4E9", "Virus" = "#0072B2"
)

# Bacterial phyla palette (SingleM outputs)
Phyla_palette <- c(
  "Actinomycetota"    = "#E69F00",
  "Bacillota"         = "#56B4E9",
  "Bacteroidota"      = "#009E73",
  "Desulfobacterota"  = "#F0E442",
  "Pseudomonadota"    = "#0072B2",
  "Verrucomicrobiota" = "#D55E00",
  "Unclassified"      = "#808080"
)

# 3.2  Plot Functions                                                          #

# Plot K-mer abundance profiles
# Arguments:
#   df — K-mer abundance data frame (abundance, count, ShortSamples)
# Returns:
#   ggplot object showing log–log K-mer abundance curves by replicate.
plot_Kmers <- function(df) {
  df %>%
    filter(abundance > 0, count > 0) %>%
    mutate(ShortSamples = factor(
      ShortSamples,
      levels = c("+DV1","+DV2","+DV3","-DV1","-DV2","-DV3",
                 "MDA1","MDA2","MDA3","MtG1","MtG2","MtG3")
    )) %>%
    ggplot(aes(x = abundance, y = count, colour = ShortSamples)) +
    geom_line(alpha = 0.5, linewidth = 0.3) +
    scale_x_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_colour_manual(values = Replicates_palette) +
    labs(x = "K-mer abundance", y = "K-mers") +
    theme_bw() +
    theme(
      axis.text.x = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

# Plot group means with Compact Letter Display
# Performs ANOVA + Tukey HSD, then plots grouped barplots with CLD letters.
# Arguments:
#   df         — data frame containing values and groups
#   value_col  — name of column with numeric values (string)
#   group_col  — name of column with grouping factor (string)
#   PALETTE    — named vector of colors for groups
# Returns
#   a list: $plot, $tukey, $cld.
plot_means_with_cld <- function(df, value_col, group_col, PALETTE) {
  aov_model <- aov(as.formula(paste(value_col, group_col, sep = " ~ ")), data = df)
  tukey_result <- TukeyHSD(aov_model)
  cld <- multcompLetters(tukey_result[[group_col]][, "p adj"])$Letters
  
  df_summary <- df %>%
    group_by(.data[[group_col]]) %>%
    summarise(mean = mean(.data[[value_col]]),
              sd = sd(.data[[value_col]]), .groups = "drop") %>%
    mutate(cld = cld[as.character(.data[[group_col]])])
  
  y_pos <- max(df_summary$mean + df_summary$sd) * 1.1
  
  plt <- ggplot(df_summary, aes(x = .data[[group_col]], y = mean, fill = .data[[group_col]])) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    geom_text(aes(label = cld, y = y_pos), vjust = 0, size = 3) +
    scale_fill_manual(values = PALETTE) +
    labs(x = "Processing method", y = "Mean ± SD") +
    theme_bw()
  
  list(plot = plt, tukey = tukey_result, cld = cld)
}


# Plot vOTU Venn diagrams
# Produces three Venn diagrams: Assembled, Mapped, and Mapped-only vOTUs.
# Arguments:
#   df — presence–absence data frame with vOTU, Description, PA, AssembledPresence
# Returns:
#   a list of ggplot objects and underlying vOTU sets.
plot_venn_diagrams <- function(df) {
  library(ggVennDiagram); library(ggpubr)
  
  order_vec <- c("UTV","DTV","MDA","MtG")
  
  make_sets <- function(d, flag, label_map) {
    d %>% filter(.data[[flag]] == 1) %>%
      select(vOTU, Description) %>%
      split(.$Description) %>%
      lapply(\(x) unique(x$vOTU)) %>%
      .[order_vec] %>%
      setNames(label_map)
  }
  
  label_map <- c("DNase\nTreated\nVirome","Untreated\nVirome",
                 "MDA\nVirome","Meta\ngenome")
  
  mapped_sets    <- make_sets(df, "PA", label_map)
  assembled_sets <- make_sets(df, "AssembledPresence", label_map)
  
  mapped_only <- setNames(lapply(order_vec, function(desc) {
    d <- df %>% filter(Description == desc)
    setdiff(
      d$vOTU[d$PA == 1],
      d$vOTU[d$AssembledPresence == 1]
    )
  }), label_map)
  
  venn_plot <- function(sets, title) {
    ggVennDiagram(sets, label = "count", edge_size = 0, set_size = 0) +
      scale_fill_distiller(limits = c(0, 200), palette = "Blues", direction = 1) +
      labs(title = title, fill = "Count") +
      theme_bw(base_size = 9) +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5))
  }
  
  list(
    plt_Assembled   = venn_plot(assembled_sets, "Assembled vOTUs"),
    plt_Mapped      = venn_plot(mapped_sets,    "Mapped vOTUs"),
    plt_MappedOnly  = venn_plot(mapped_only,    "Mapped-only vOTUs"),
    sets_Assembled  = assembled_sets,
    sets_Mapped     = mapped_sets,
    sets_MappedOnly = mapped_only
  )
}

# Plot PCoA ordination and PERMANOVA results
# Performs Bray–Curtis PCoA and PERMANOVA.
# Arguments:
#   df       — TPM matrix (rows = contigs, columns = samples)
#   metadata — data frame of sample metadata
#   Description_palette — named vector of colors for sample descriptions
# Returns:
#   a list of: $pcoa_plot, $permanova, $pairwise_permanova.
plot_pcoa <- function(df, metadata, Description_palette = Description_palette) {
  metadata <- metadata %>%
    mutate(DescriptionLong = factor(
      DescriptionModifiers[Description],
      levels = DescriptionLevels
    ))
  
  dist_mat <- df %>%
    column_to_rownames("Contig") %>%
    t() %>%
    vegdist(method = "bray")
  
  pcoa_res <- cmdscale(dist_mat, eig = TRUE, k = 2)
  axes <- as.data.frame(pcoa_res$points)
  colnames(axes) <- c("PCoA1","PCoA2")
  axes$Sample <- rownames(axes)
  
  axes <- merge(axes, metadata, by = "Sample") %>%
    mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels))
  
  var_exp <- 100 * pcoa_res$eig[1:2] / sum(pcoa_res$eig)
  
  p <- ggplot(axes, aes(PCoA1, PCoA2, color = DescriptionLong, shape = DescriptionLong)) +
    geom_point(size = 4, alpha = 0.6) +
    scale_colour_manual(values = Description_palette) +
    scale_shape_manual(values = c(16,17,15,18)) +
    coord_fixed() +
    theme_bw() +
    labs(
      x = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y = sprintf("PC2 (%.1f%%)", var_exp[2]),
      color = "Processing method"
    ) +
    theme(panel.grid = element_blank())
  
  perm <- adonis2(dist_mat ~ Description, data = axes, permutations = 999)
  pair <- pairwise.adonis2(dist_mat ~ Description, data = axes,
                           p.adjust.method = "BH", permutations = 999)
  
  list(pcoa_plot = p, permanova = perm, pairwise_permanova = pair)
}

# Generate multiple PCoAs (all, no ssDNA, no Monodnaviria)
# Wrapper to produce standard and filtered ordinations.
# Arguments:
#   df_mapping     — TPM matrix (rows = contigs, columns = samples)
#   stacked        — list containing $df_Baltimore with DNA types
#   GenomadData    — data frame of annotated contigs
#   metadata       — data frame of sample metadata
# Returns:
#   a list of PCoA plots and PERMANOVA results.
produce_pcoas <- function(df_mapping, stacked, GenomadData, metadata,
                          Description_palette = Description_palette) {
  df_all  <- calculate_beta_diversity(df_mapping)
  res_all <- plot_pcoa(df_all, metadata, Description_palette)
  
  monodna <- GenomadData %>% filter(grepl("Monodnaviria", taxonomy)) %>% pull(seq_name)
  ssDNA   <- stacked$df_Baltimore %>% filter(DNA == "ssDNA") %>% pull(Contig)
  
  df_no_mono <- df_mapping %>% filter(!Contig %in% monodna) %>% calculate_beta_diversity()
  df_no_ss   <- df_mapping %>% filter(!Contig %in% ssDNA)   %>% calculate_beta_diversity()
  
  res_no_ss <- plot_pcoa(df_no_ss, metadata, Description_palette)
  res_no_m  <- plot_pcoa(df_no_mono, metadata, Description_palette)
  
  list(
    pcoa_all          = res_all$pcoa_plot,
    pcoa_no_ssDNA     = res_no_ss$pcoa_plot,
    pcoa_no_monodna   = res_no_m$pcoa_plot,
    permanova_all     = res_all$permanova,
    permanova_no_ssDNA= res_no_ss$permanova,
    permanova_no_mono = res_no_m$permanova,
    pairwise_all      = res_all$pairwise_permanova
  )
}

# Plot stacked barplots of viral composition
# Generates Baltimore-class and family-level stacked barplots.
# Arguments:
#   df_mapping  — TPM matrix (rows = contigs, columns = samples)
#   metadata    — data frame of sample metadata
#   df_Genomad  — data frame of annotated contigs
# Returns:
#   A list of: $plt_stacked, $plt_DNA_stacked, $df_Baltimore, $dsDNA_Monodnaviria, $ssDNA_extra.
plot_vOTU_stacked_barplot <- function(df_mapping, metadata, df_Genomad) {
  
  df_tax <- df_Genomad %>%
    rename(Contig = seq_name) %>%
    separate(taxonomy,
             into = c("Viruses","Realm","Kingdom","Phylum",
                      "Class","Order","Family"), sep = ";", fill = "right")
  
  dsDNA_Monodnaviria <- df_tax %>% filter(Class == "Papovaviricetes") %>% pull(Contig)
  ssDNA_extra <- df_tax %>%
    filter(Family %in% c("Alphasatellitidae","Spiraviridae",
                         "Anelloviridae","Tolecusatellitidae")) %>% pull(Contig)
  
  df_tax <- df_tax %>%
    mutate(DNA = case_when(
      Realm == "Monodnaviria" & Contig %in% dsDNA_Monodnaviria ~ "dsDNA",
      Realm == "Monodnaviria"                                  ~ "ssDNA",
      Realm %in% c("Duplodnaviria","Adnaviria","Varidnaviria") ~ "dsDNA",
      Realm == "Unassigned" | is.na(Realm)                     ~ "Unassigned",
      TRUE                                                     ~ "Unassigned"
    ))
  
  df_Baltimore <- df_mapping %>%
    pivot_longer(-Contig, names_to = "Sample", values_to = "TPM") %>%
    filter(TPM > 0) %>%
    merge(metadata, by = "Sample") %>%
    mutate(
      DescriptionLong = factor(DescriptionModifiers[Description],
                               levels = DescriptionLevels),
      rep = substr(ShortSamples, nchar(ShortSamples), nchar(ShortSamples))
    ) %>%
    merge(df_tax, by = "Contig", all.x = TRUE) %>%
    mutate(DNA = replace_na(DNA, "Unassigned"))
  
  # Baltimore-class barplot
  df_DNA_stacked <- df_Baltimore %>%
    group_by(DNA, rep, DescriptionLong) %>%
    summarise(TPM = sum(TPM), .groups = "drop") %>%
    mutate(RelAbund = TPM / 10000)
  
  plt_DNA_stacked <- ggplot(df_DNA_stacked, aes(rep, RelAbund, fill = DNA)) +
    geom_bar(stat = "identity") +
    facet_wrap(~DescriptionLong, nrow = 1) +
    scale_fill_manual(values = Stacked_palette) +
    labs(x = "Processing method", y = "Relative abundance (%)") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "bottom",
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Family-level barplot
  df_stacked <- df_Baltimore %>%
    mutate(
      Family = coalesce(Family, Order, Class, Phylum, Kingdom, Realm, Viruses),
      Order  = coalesce(Order,  Class, Phylum, Kingdom, Realm, Viruses),
      Family = case_when(
        Class == "Caudoviricetes" & Family != "Crassvirales" ~ "Other Caudoviricetes",
        Realm == "Monodnaviria" & Family != "Microviridae"   ~ "Other Monodnaviria",
        is.na(Family)                                        ~ "Unassigned",
        TRUE                                                 ~ Family
      )
    ) %>%
    group_by(Family, ShortSamples, rep, DescriptionLong) %>%
    summarise(TPM = sum(TPM), .groups = "drop") %>%
    mutate(RelAbund = TPM / 10000)
  
  plt_stacked <- ggplot(df_stacked, aes(rep, RelAbund, fill = Family)) +
    geom_bar(stat = "identity") +
    facet_wrap(~DescriptionLong, nrow = 1) +
    scale_fill_manual(values = Stacked_palette) +
    labs(y = "Relative abundance (%)") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(colour = "black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "bottom") +
    guides(fill = guide_legend(title = "Lifestyle", nrow = 3))
  
  list(
    plt_stacked        = plt_stacked,
    plt_DNA_stacked    = plt_DNA_stacked,
    df_Baltimore       = df_Baltimore,
    dsDNA_Monodnaviria = dsDNA_Monodnaviria,
    ssDNA_extra        = ssDNA_extra
  )
}

# Helper function to streamline repetitive calls for QUAST metric plots
# Arguments:
#   df       — data frame with Metric, Value, DescriptionLong columns
#   metric   — specific metric to plot (string)
#   y_label  — y-axis label (string)
#   title    — plot title (string)
# Returns:
#   ggplot object of the specified QUAST metric plot
make_quast_plot <- function(df, metric, y_label, title) {
  plot_means_with_cld(df %>% filter(Metric == metric), "Value", "DescriptionLong", Description_palette)$plot +
    labs(y = y_label, title = title) +
    theme_bw(base_size = 10) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

# Helper: one-way ANOVA + Tukey extractor
# Arguments:
#   df            — data frame with value and group columns
#   value_col     — name of column with numeric values (string)
#   group_col     — name of column with grouping factor (string)
#   panel_label   — label for the panel/analysis (string)
#   measure_label — optional label for the measure (string)
# Returns:
#  a list: $anova (data frame of ANOVA results), $tukey (data frame of Tukey results)
one_way_stats <- function(df, value_col, group_col, panel_label, measure_label = NULL) {
  stopifnot(value_col %in% names(df), group_col %in% names(df))
  
  d <- df %>%
    dplyr::select(all_of(c(value_col, group_col))) %>%
    dplyr::rename(value = !!value_col, group = !!group_col) %>%
    dplyr::mutate(
      value = as.numeric(value),
      group = factor(group)
    ) %>%
    dplyr::filter(stats::complete.cases(value, group))
  
  # Skip if fewer than two groups have data
  if (dplyr::n_distinct(d$group) < 2) {
    return(list(
      anova = tibble::tibble(
        panel = panel_label,
        measure = measure_label %||% value_col,
        k_groups = dplyr::n_distinct(d$group),
        df1 = NA_integer_, df2 = NA_integer_,
        F = NA_real_, p_value = NA_real_
      ),
      tukey = tibble::tibble(
        panel = panel_label,
        measure = measure_label %||% value_col,
        contrast = character(),
        diff = numeric(), lwr = numeric(), upr = numeric(), p_adj = numeric()
      )
    ))
  }
  
  fit <- stats::aov(value ~ group, data = d)
  sm  <- summary(fit)[[1]]
  
  df1 <- sm[["Df"]][1]
  df2 <- sm[["Df"]][2]
  Fv  <- sm[["F value"]][1]
  pv  <- sm[["Pr(>F)"]][1]
  
  anova_line <- tibble::tibble(
    panel   = panel_label,
    measure = measure_label %||% value_col,
    k_groups = dplyr::n_distinct(d$group),
    df1 = df1, df2 = df2,
    F = as.numeric(Fv),
    p_value = as.numeric(pv)
  )
  
  tk <- TukeyHSD(fit)[["group"]] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("contrast") %>%
    dplyr::as_tibble() %>%
    dplyr::transmute(
      panel   = panel_label,
      measure = measure_label %||% value_col,
      contrast,
      diff = diff, lwr = lwr, upr = upr, p_adj = `p adj`
    )
  
  list(anova = anova_line, tukey = tk)
}
