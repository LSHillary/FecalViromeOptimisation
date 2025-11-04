###############################################################################
# Fecal Virome Optimisation Analysis – Parent Script
# Author: Luke Hillary
# Institution: University of California, Davis – Emerson Lab
# Project Repository: https://github.com/LSHillary/FecalViromeOptimisation
# Version: v1.0 (Manuscript Release)
# Date: 2025-11-04
###############################################################################
# Description:
# This parent script reproduces all analyses, figures, and tables associated
# with the *Fecal Virome Optimisation* (FVO) study.
#
# The pipeline integrates viral, microbial, and quality-control data to compare
# four processing methods for human gut virome characterisation:
#   - Untreated viromes (UTV)
#   - DNase-treated viromes (DTV)
#   - MDA-amplified viromes (MDA)
#   - Total metagenomes (MtG)
#
# The workflow performs the following:
#   0. Setup and environment configuration
#   1. Execution of the {targets} pipeline for reproducible data import and preprocessing
#   2. Loading and wrangling of processed datasets for figure generation
#   3. Main Figure Production (Figs 1–4)
#   4. Supplementary Figure Production (Figs S1–S7)
#   5. Statistical summary export (ANOVA and Tukey HSD tables)
#   6. Graphical Abstract assembly
#
# Each figure and table corresponds directly to those in the published manuscript.
# The script depends on intermediate objects produced by `_targets.R`
# and functions defined in `scripts/Functions_FVO.R`.
#
# Output:
#   - Publication-quality figures in /figures/
#   - Statistical summaries in /stats/
#   - Reproducible DAG visualisation (targets_dag.html)
#
# Reproducibility:
#   - All steps are deterministic and version-locked via {targets}.
#   - Figures and statistics regenerate automatically when input data change.
#
# Running the code:
#   - Typically I will run this code interactively line-by-line in RStudio.
#   - Alternatively, the entire script can be sourced in an R session, however
#     this is not recommended for debugging or iterative development, and has
#     not been extensively tested.
#
# Recommended citation:
#   <INSERT CITATION DETAILS HERE>
#
# License: MIT
###############################################################################

# 0 - Setup ----

## Clear environment
rm(list = ls())

# Set seed
set.seed(42)

## Core tidyverse (includes dplyr, ggplot2, tidyr, readr, etc.)
library(conflicted)
library(tidyverse)

# Ensure dplyr functions take precedence (avoids conflicts with plyr, reshape2, etc.)
lapply(c("filter", "lag", "select", "rename", "summarise", "arrange", "mutate"),
       conflict_prefer, winner = "dplyr")

## Visualisation and figure assembly
library(ggh4x)          # Extended ggplot2 functionality
library(ggpubr)         # Figure layout and publication tools
library(ggrastr)        # Rasterisation for large ggplots

## Statistical modelling and comparisons
library(lme4)           # Linear mixed-effects models
library(lmerTest)       # P-values for lme4 models
library(multcompView)   # Compact letter displays for pairwise tests
library(pairwiseAdonis) # Pairwise PERMANOVA tests
library(vegan)          # Ordination and ecological distance metrics

## Data wrangling and pipeline management
library(reshape2)       # Legacy reshape support
library(targets)        # Reproducible pipeline management

## Load custom functions
source("scripts/Functions_FVO.R")

# These definitions override the defaults in Functions_FVO.R
# and ensure consistent formatting of processing method labels and colors
# This is necessary for the figures that are generated for the first time
# in this script.

DescriptionModifiers <- c(
  "UTV" = "Untreated\nVirome",
  "DTV" = "DNase\nTreated\nVirome",
  "MDA" = "MDA\nVirome",
  "MtG" = "Meta\ngenome"
)

DescriptionLevels <- c(
  "Untreated\nVirome",
  "DNase\nTreated\nVirome",
  "MDA\nVirome",
  "Meta\ngenome"
)

Description_palette <- c(
  "Untreated\nVirome" = "#56B4E9",
  "DNase\nTreated\nVirome" = "#0072B2",
  "MDA\nVirome" = "#009E73",
  "Meta\ngenome" = "#FF7F0E"
)

# Visualize the {targets} dependency graph interactively in RStudio
tar_visnetwork()

message("Setup complete. Targets version: ", packageVersion("targets"))

# 1 - Run targets pipeline for Fecal Virome Optimisation Analysis ----

# Run the full analysis pipeline
# This executes all targets defined in _targets.R, rebuilding only those out of date
tar_make()

# Save the dependency graph as a standalone HTML for reproducibility documentation
htmlwidgets::saveWidget(
  tar_visnetwork(targets_only = TRUE),
  file = "figures/targets_dag.html",
  selfcontained = TRUE
)

message("Pipeline completed and DAG exported to figures/targets_dag.html")

# 2 - Load and wrangle data ----

## Load required objects for figure generation
metadata        <- tar_read(metadata)      # Sample metadata
df_mapping      <- tar_read(vOTU_tpm)      # vOTU relative abundance data (TPM)
df_reads        <- tar_read(all_read_counts) # Read count summaries
df_host         <- tar_read(iPhopData)     # Host predictions
df_Bacphlip   <- tar_read(BacPhlipData_vOTU)
df_Genomad    <- tar_read(GenomadData)
df_clustering <- tar_read(ClusteringData)
df_SingleM    <- tar_read(SingleMData)
ls_stacked     <- tar_read(plt_stacked)
ls_Venns   <- tar_read(plt_Venns)
df_clustering <- tar_read(ClusteringData)

# Prepare vOTU mapping table for relative abundance plots
df_mapping_cleaned <- df_mapping %>%
  rename(vOTU = Contig) %>%
  pivot_longer(-vOTU, names_to = "Sample", values_to = "TPM") %>%
  mutate(RelAbund = TPM / 10000)

# Clean and relabel read-count data for figure summaries
df_reads <- df_reads %>%
  mutate(
    DescriptionLong = factor(DescriptionModifiers[Description],
                             levels = DescriptionLevels),
    total_sequences = total_sequences / 1e6
  )

# Import and tidy QUAST assembly summary (includes viral contig counts)
df_quast <- read.csv("data/Quast/Quast_report.tsv", sep = "\t") %>%
  rename_with(~ str_remove(., "_renamed_contigs"), contains("_renamed_contigs")) %>%
  t() %>% as.data.frame() %>%
  setNames(.[1,]) %>% .[-1,] %>%
  rownames_to_column("Sample") %>%
  merge(metadata, by = "Sample") %>%
  select(Sample, Description, `# contigs`, `Total length`, `Largest contig`,
         N50, L50, auN, FilteredViralContigs) %>%
  mutate(across(-c(Sample, Description), as.numeric)) %>%
  pivot_longer(-c(Sample, Description), names_to = "Metric", values_to = "Value") %>%
  mutate(
    DescriptionLong = factor(DescriptionModifiers[Description],
                             levels = DescriptionLevels)
  )

## Compute summary statistics for key read metrics

# Mean total reads per processing method (in millions)
df_read_means <- df_reads %>%
  distinct(Description, Sample, total_sequences) %>%
  summarise(
    Mean_total_reads = mean(total_sequences, na.rm = TRUE),
    SD_total_reads   = sd(total_sequences, na.rm = TRUE)
  )

# Mean viral reads and richness per processing method
df_viral_means <- df_reads %>%
  distinct(Description, Sample, total_sequences, FilteredViralContigs, ViralReads) %>%
  group_by(Description) %>%
  summarise(
    mean_viral_reads       = mean(ViralReads / 1e6, na.rm = TRUE),
    mean_viral_richness    = mean(FilteredViralContigs, na.rm = TRUE),
    mean_total_sequences   = mean(total_sequences, na.rm = TRUE),
    SD_viral_reads         = sd(ViralReads / 1e6, na.rm = TRUE),
    SD_viral_richness      = sd(FilteredViralContigs, na.rm = TRUE),
    SD_total_sequences     = sd(total_sequences, na.rm = TRUE)
  )

# Mean and SD of read-type proportions (Human, rRNA, Viral)
df_mean_read_types <- df_reads %>%
  distinct(Description, Sample, total_sequences, Type, ReadPercentage) %>%
  group_by(Description, Type) %>%
  summarise(
    mean_reads = mean(ReadPercentage, na.rm = TRUE),
    SD_reads   = sd(ReadPercentage, na.rm = TRUE)
  )

## Prepare Genomad-filtered and read-depth data for Figs 3–4

# Filter GeNomad contigs to those present in viral clusters
all_filtered_contigs <- df_clustering %>%
  summarise(all_clusters = paste(Cluster, collapse = ",")) %>%
  pull(all_clusters) %>%
  strsplit(",") %>%
  unlist()

df_Genomad_filtered <- df_Genomad %>%
  filter(seq_name %in% all_filtered_contigs) %>%
  rename(Contig = seq_name, Sample = sample) %>%
  left_join(metadata, by = "Sample") %>%
  mutate(
    DescriptionLong = DescriptionModifiers[Description],
    DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels)
  ) %>%
  filter(!is.na(DescriptionLong))

# Import and process read-depth subsampling data
df_read_depth_100 <- read_and_merge_tsvs("data/Subsampling/", "G.tsv") %>%
  select(Contig, contains("mapped Mean")) %>%
  pivot_longer(cols = -Contig, names_to = "Sample", values_to = "Coverage") %>%
  mutate(
    Depth = as.numeric(str_extract(Sample, "(?<=depth_).*?(?=G)")),
    Sample = str_extract(Sample, ".*?(?=_subsampled)")
  ) %>%
  left_join(metadata, by = "Sample") %>%
  filter(!is.na(Coverage)) %>%
  group_by(Description, Depth, Sample) %>%
  summarise(
    TotalCount = n(),
    CountOver100 = sum(Coverage > 100, na.rm = TRUE),
    Percentage = (CountOver100 / TotalCount) * 100,
    .groups = "drop"
  ) %>%
  complete(Sample, Description,
           fill = list(TotalCount = 0, CountOver100 = 0, Percentage = 0)) %>%
  mutate(
    Description = replace_na(Description, "MtG"),
    DescriptionLong = DescriptionModifiers[Description]
  )

# 3 - Figure Generation ----
## Figure 1 — How the four methods compare in terms of viral recovery ----

# 3.1  Venn diagrams (top panel)

# Extract Venn plot
plt_Venns  <- ls_Venns$plt_Venn
Fig1_top   <- plt_Venns

# 3.2  Mapped vs. assembled virus counts (middle panels)

# Count mapped vOTUs per sample
df_mapped_viruses <- df_mapping_cleaned %>%
  filter(RelAbund != 0) %>%
  group_by(Sample) %>%
  summarise(Count = n()) %>%
  merge(metadata, by = "Sample") %>%
  mutate(
    Description = DescriptionModifiers[Description],
    Description = factor(Description, levels = DescriptionLevels)
  ) %>%
  select(Sample, Description, Count)

# Extract assembled vOTU counts per sample
df_assembled_viruses <- df_quast %>%
  pivot_wider(names_from = Metric, values_from = Value) %>%
  select(Sample, DescriptionLong, FilteredViralContigs) %>%
  rename(Description = DescriptionLong)

# Summarise mapped virus counts (mean ± SD)
df_mapped_viruses_summary <- df_mapped_viruses %>%
  group_by(Description) %>%
  summarise(
    mean_mapped_viruses = mean(Count, na.rm = TRUE),
    SD_mapped_viruses   = sd(Count, na.rm = TRUE)
  )

# Combine viral richness and mapped-virus summaries
df_viral_means_2 <- df_viral_means %>%
  mutate(Description = DescriptionModifiers[Description]) %>%
  merge(df_mapped_viruses_summary, by = "Description")

# Prepare long + wide summary tables (used for ratios and manuscript text)
df_long <- df_viral_means_2 %>%
  pivot_longer(cols = c(mean_viral_richness, mean_mapped_viruses),
               names_to = "Metric", values_to = "Mean") %>%
  pivot_longer(cols = c(SD_viral_richness, SD_mapped_viruses),
               names_to = "Metric_SD", values_to = "SD") %>%
  filter(str_remove(Metric_SD, "SD_") == str_remove(Metric, "mean_")) %>%
  mutate(
    Metric = recode(Metric,
                    "mean_viral_richness" = "Viral Richness",
                    "mean_mapped_viruses" = "Mapped Viruses"),
    Description = factor(Description, levels = DescriptionLevels)
  )

# Create wide summary table for ratio calculations
df_wide <- df_long %>%
  filter(Metric %in% c("Viral Richness", "Mapped Viruses") == TRUE) %>%
  select(Description, Metric, Mean) %>%
  pivot_wider(names_from = Metric, values_from = Mean) %>%
  mutate(Ratio = `Mapped Viruses` / `Viral Richness`)

# Plot mapped vOTU counts
plt_list_mapped <- plot_means_with_cld(df_mapped_viruses, "Count", "Description", Description_palette)

# Prepare mapped vOTU plot
plt_mapped <- plt_list_mapped$plot +
  labs(y = "Mapped vOTUs", fill = "Processing Method") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  lims(y = c(0, 525))

# Plot assembled vOTU counts
plt_list_assembled <- plot_means_with_cld(df_assembled_viruses, "FilteredViralContigs",
                                          "Description", Description_palette)

# Prepare assembled vOTU plot
plt_assembled <- plt_list_assembled$plot +
  labs(y = "Assembled vOTUs", fill = "Processing Method") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  lims(y = c(0, 525))

# Combine mapped + assembled panels
Fig1_middle <- ggarrange(
  plt_assembled, plt_mapped,
  ncol = 2, nrow = 1,
  labels = c("c", "d"),
  common.legend = TRUE,
  legend = "none"
)

# 3.3  Beta-diversity (PCoA) and viral reads (bottom panels)

# Generate PCoA plots and PERMANOVA tests
plt_pcoas <- produce_pcoas(df_mapping, ls_stacked, df_Genomad, metadata,
                           Description_palette = Description_palette)

# Bray–Curtis dissimilarity and dispersion tests
df_bray <- calculate_beta_diversity(df_mapping) %>%
  column_to_rownames(var = "Contig") %>%
  t() %>%
  vegdist(method = "bray")

groups <- metadata %>%
  filter(Sample %in% labels(df_bray)) %>%
  arrange(match(Sample, labels(df_bray))) %>%
  pull(Description)

bd <- betadisper(df_bray, group = groups)
anova(bd)
TukeyHSD(bd)

# Prepare main PCoA plot
pcoa_all <- plt_pcoas$pcoa_all +
  labs(colour = "Processing Method", fill = "Processing Method",
       shape = "Processing Method") +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  coord_fixed() +
  theme_bw() +
  # Remove grid lines
  theme(panel.grid = element_blank())

# Viral-read percentage plots
plt_list_viral_reads <- plot_means_with_cld(
  df_reads %>% filter(Type == "Viral"),
  "ReadPercentage", "DescriptionLong", Description_palette
)

plt_viral_reads <- plt_list_viral_reads$plot +
  labs(y = "Viral reads (%)", fill = "Processing Method") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

# Combine bottom panels
Fig1_bottom <- ggarrange(
  plt_viral_reads, pcoa_all,
  ncol = 2, nrow = 1,
  labels = c("e", "f")
)

# 3.4  Assemble and save Figure 1

Fig1 <- ggarrange(
  Fig1_top, Fig1_middle, Fig1_bottom,
  ncol = 1, nrow = 3,
  heights = c(1, 1, 1)
)

# Export in both PDF and SVG for publication
for (ext in c("pdf", "svg")) {
  ggsave(
    filename = file.path("figures", paste0("Fig1.", ext)),
    plot = Fig1,
    width = 170, height = 190, units = "mm"
  )
}

# Note: The manuscript version was lightly edited in Affinity Designer
#       for label legibility—no data or layout changes.

## Figure 2 - The MDA-viromes are super different, but we can correct for the differences ----

# Import stacked bar plot of vOTU relative abundance
plt_stacked <- tar_read(plt_stacked)$plt_stacked +
  theme(legend.title = element_blank())

# Assemble top panel
Fig2_top <- plt_stacked +
  theme(legend.position = "right")

# Import and style k-mer composition plot
plt_kmers <- tar_read(plt_kmers) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Create PCoA plot with ssDNA viruses removed
plt_pcoa_ssdna <- plt_pcoas$pcoa_no_ssDNA +
  labs(
    colour = "Sample\nType",
    fill   = "Sample\nType",
    shape  = "Sample\nType",
    y = "PC2\n(3.2%)"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.margin = ggplot2::margin(t = 0, r = 40, b = 0, l = 0)
  ) +
  guides(
    colour = guide_legend(nrow = 1),
    shape  = guide_legend(nrow = 1),
    fill   = guide_legend(nrow = 1)
  )

# Combine bottom panels (k-mers + PCoA)
Fig2_bottom <- ggarrange(
  plt_kmers, plt_pcoa_ssdna,
  ncol = 2, nrow = 1,
  labels = c("b", "c"),
  widths = c(1, 2)
)

# Assemble full Figure 2
Fig2 <- ggarrange(
  Fig2_top, Fig2_bottom,
  nrow = 2,
  labels = c("a", "")
)

# Save both PDF and SVG versions
for (ext in c("pdf", "svg")) {
  ggsave(
    filename = file.path("figures", paste0("Fig2.", ext)),
    plot = Fig2,
    width = 170, height = 100, units = "mm"
  )
}


## Figure 3 - The metagenomes recover different viral community members ----

# Add host predictions
df_Genomad_filtered_host <- df_Genomad_filtered %>%
  mutate(Provirus = ifelse(grepl("provirus", Contig), "Provirus", "Virus")) %>%
  merge(df_host, by.x = "Contig", by.y = "vOTU", all.x = TRUE)

# Summarise viral host distribution
df_Genomad_filtered_host_summary <- df_Genomad_filtered_host %>%
  select(Contig, Provirus, HostPhylum, Sample) %>%
  merge(metadata, by = "Sample", all.x = TRUE) %>%
  filter(HostPhylum %in% c("Bacteroidota", "Bacillota")) %>%
  rename(vOTU = Contig) %>%
  group_by(Provirus, HostPhylum, Description) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(
    Description = DescriptionModifiers[Description],
    Description = factor(Description, levels = DescriptionLevels)
  )

# Merge BACPHLIP with host predictions
df_Bacphlip_host <- df_Bacphlip %>%
  merge(df_host, by = "vOTU", all.x = TRUE)

# Merge lifestyle data with mapping
df_Lifestyle_mapped <- df_Bacphlip %>%
  merge(df_mapping_cleaned, by = "vOTU", all.y = TRUE) %>%
  merge(metadata, by = "Sample", all.x = TRUE) %>%
  select(vOTU, ShortSamples, RelAbund, status, Provirus, Description) %>%
  mutate(RelAbund = replace_na(RelAbund, 0))

# Merge lifestyle + host data
df_Lifestyle_mapped_host <- df_Bacphlip_host %>%
  merge(df_mapping_cleaned, by = "vOTU", all.y = TRUE) %>%
  merge(metadata, by = "Sample", all.x = TRUE) %>%
  select(vOTU, ShortSamples, RelAbund, status, Provirus, Description, HostPhylum) %>%
  mutate(
    RelAbund = replace_na(RelAbund, 0),
    DescriptionLong = DescriptionModifiers[Description]
  )

# Summarise lifestyle data
df_Lifestyle_mapped_summary <- df_Lifestyle_mapped %>%
  group_by(status, ShortSamples, Description) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop") %>%
  mutate(
    status = factor(status, c("Temperate", "Virulent", "Unclassified")),
    Value = RelAbund,
    rep = substr(ShortSamples, nchar(ShortSamples), nchar(ShortSamples)),
    DescriptionLong = DescriptionModifiers[Description]
  )

# Summarise lifestyle + host data
df_Lifestyle_mapped_host_summary <- df_Lifestyle_mapped_host %>%
  group_by(status, ShortSamples, Description, HostPhylum) %>%
  summarise(RelAbund = sum(RelAbund), Count = n(), .groups = "drop") %>%
  mutate(
    status = factor(status, c("Temperate", "Virulent", "Unclassified")),
    Value = RelAbund,
    rep = substr(ShortSamples, nchar(ShortSamples), nchar(ShortSamples)),
    DescriptionLong = DescriptionModifiers[Description]
  )

# Calculate temperate:virulent ratios
df_Lifestyle_mapped_host_ratio_summary <- df_Lifestyle_mapped_host %>%
  filter(
    HostPhylum %in% c("Bacillota", "Bacteroidota"),
    status %in% c("Temperate", "Virulent"),
    RelAbund > 0
  ) %>%
  group_by(ShortSamples, DescriptionLong, HostPhylum, status) %>%
  summarise(
    TotalAbundance = sum(RelAbund, na.rm = TRUE),
    UniqueVOTUs = n_distinct(vOTU),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = status,
    values_from = c(TotalAbundance, UniqueVOTUs),
    values_fill = 0
  ) %>%
  mutate(
    AbundanceRatio = TotalAbundance_Temperate / (TotalAbundance_Temperate + TotalAbundance_Virulent),
    CountRatio = UniqueVOTUs_Temperate / (UniqueVOTUs_Temperate + UniqueVOTUs_Virulent)
  ) %>%
  group_by(DescriptionLong, HostPhylum) %>%
  summarise(
    mean_abundance_ratio = mean(AbundanceRatio, na.rm = TRUE),
    sd_abundance_ratio   = sd(AbundanceRatio, na.rm = TRUE),
    se_abundance_ratio   = sd_abundance_ratio / sqrt(n()),
    mean_count_ratio     = mean(CountRatio, na.rm = TRUE),
    sd_count_ratio       = sd(CountRatio, na.rm = TRUE),
    se_count_ratio       = sd_count_ratio / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# Reshape for plotting
df_ratio_long <- df_Lifestyle_mapped_host_ratio_summary %>%
  pivot_longer(
    cols = c(mean_abundance_ratio, se_abundance_ratio, mean_count_ratio, se_count_ratio),
    names_to = c("stat", "metric"),
    names_pattern = "(.+?)_(abundance|count)_ratio",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    metric = recode(metric, "abundance" = "Relative\nAbundance", "count" = "Unique\nvOTUs"),
    DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels)
  )

# Lifestyle stacked barplot
plt_Lifestyle <- df_Lifestyle_mapped_summary %>%
  mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels)) %>%
  ggplot(aes(x = rep, y = Value, fill = status)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  labs(x = "Sample", y = "vOTU\nRelative Abundance (%)", fill = "Lifestyle") +
  scale_fill_manual(values = Stacked_palette) +
  guides(fill = guide_legend(title = "Lifestyle", nrow = 2)) +
  facet_wrap(~DescriptionLong, nrow = 1, drop = TRUE)

# Subset for proviruses (temperate only)
df_Provirus_mapped_summary <- df_Bacphlip %>%
  merge(df_mapping_cleaned, by = "vOTU", all.y = TRUE) %>%
  merge(metadata, by = "Sample", all.x = TRUE) %>%
  filter(status == "Temperate") %>%
  select(vOTU, ShortSamples, RelAbund, Provirus, Description) %>%
  mutate(RelAbund = replace_na(RelAbund, 0)) %>%
  group_by(Provirus, ShortSamples, Description) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop") %>%
  group_by(ShortSamples) %>%
  mutate(RelAbund = RelAbund / sum(RelAbund) * 100) %>%
  mutate(
    Provirus = factor(Provirus, c("Virus", "Provirus")),
    Value = RelAbund,
    rep = substr(ShortSamples, nchar(ShortSamples), nchar(ShortSamples)),
    DescriptionLong = DescriptionModifiers[Description]
  )

# Provirus plot
plt_list_Proviruses <- df_Provirus_mapped_summary %>%
  filter(Provirus == "Provirus") %>%
  mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels)) %>%
  plot_means_with_cld("Value", "DescriptionLong", Description_palette)

plt_Proviruses <- plt_list_Proviruses$plot +
  theme_bw(base_size = 10) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(y = "Provirus\nProportion (%)") +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(title = "Sample\nType", nrow = 2)) +
  scale_x_discrete(labels = function(x) str_wrap(x, 10))

# Assemble top panels
Fig3_top <- ggarrange(plt_Lifestyle, plt_Proviruses, ncol = 2, labels = c("a", "b"))

# Host prediction stacked barplot
df_host_stacked <- df_host %>%
  merge(df_mapping_cleaned, by = "vOTU", all.y = TRUE) %>%
  select(vOTU, Sample, HostPhylum, RelAbund) %>%
  mutate(
    HostPhylum = replace_na(HostPhylum, "Unclassified"),
    HostPhylum = factor(HostPhylum,
                        c("Actinomycetota","Bacillota","Bacteroidota",
                          "Desulfobacterota","Pseudomonadota","Verrucomicrobiota","Unclassified"))
  ) %>%
  filter(RelAbund > 0) %>%
  group_by(Sample, HostPhylum) %>%
  summarise(RelAbund = sum(RelAbund), Count = n(), .groups = "drop") %>%
  merge(metadata, by = "Sample") %>%
  mutate(
    rep = substr(ShortSamples, nchar(ShortSamples), nchar(ShortSamples)),
    DescriptionLong = DescriptionModifiers[Description],
    Data = "Viruses (by predicted host)"
  )

# Summary of host data
df_host_stacked_summary <- df_host_stacked %>%
  group_by(DescriptionLong, Data, rep) %>%
  mutate(Count = (Count / sum(Count)) * 100) %>%
  group_by(DescriptionLong, HostPhylum, Data) %>%
  summarise(
    Mean_RelAbund = mean(RelAbund, na.rm = TRUE),
    SD_RelAbund = sd(RelAbund, na.rm = TRUE),
    Mean_Count = mean(Count),
    SD_Count = sd(Count)
  )

# Load and process SingleM data for metagenomes
df_Microbiome <- df_SingleM %>%
  mutate(
    rep = substr(ShortSamples, nchar(ShortSamples), nchar(ShortSamples)),
    DescriptionLong = DescriptionModifiers[Description]
  ) %>%
  filter(DescriptionLong == "Meta\ngenome") %>%
  rename(HostPhylum = Phylum) %>%
  mutate(HostPhylum = str_remove(HostPhylum, "p__")) %>%
  select(DescriptionLong, HostPhylum, RelAbund, rep) %>%
  mutate(Data = "Bacteria")

# Merge bacterial and viral host data
df_Microbiome2 <- df_Microbiome %>%
  merge(df_host_stacked, by = c("DescriptionLong","HostPhylum","Data","rep"), all = TRUE) %>%
  mutate(RelAbund = coalesce(RelAbund.x, RelAbund.y)) %>%
  select(-RelAbund.x, -RelAbund.y) %>%
  mutate(
    DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels),
    Data = factor(Data, levels = c("Bacteria", "Viruses (by predicted host)"))
  )

# Calculate Bacillota:Bacteroidota ratio
df_Bacillota_ratio <- df_Microbiome2 %>%
  filter(HostPhylum %in% c("Bacillota", "Bacteroidota") & Data == "Bacteria") %>%
  pivot_wider(names_from = HostPhylum, values_from = RelAbund) %>%
  mutate(Ratio = Bacillota / Bacteroidota)

mean(df_Bacillota_ratio$Ratio)
sd(df_Bacillota_ratio$Ratio)

# Microbiome and host stacked barplot
plt_microbiome <- df_Microbiome2 %>%
  ggplot(aes(x = rep, y = RelAbund, fill = HostPhylum)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "white")
  ) +
  labs(x = "Sample", y = "vOTU\nRelative Abundance (%)", fill = "Predicted\nHost Phylum") +
  scale_fill_manual(values = Phyla_palette) +
  facet_nested(~ Data + DescriptionLong, drop = TRUE) +
  theme(
    legend.position = "right",
    legend.justification = c("center", "top"),
    legend.box.margin = margin(t = -62, r = 0, b = 0, l = 0)
  )

# Assemble full Figure 3
Fig3 <- ggarrange(Fig3_top, plt_microbiome, ncol = 1, labels = c("", "c"), heights = c(1, 1))

# Save Figure 3
for (ext in c("pdf", "svg")) {
  ggsave(
    filename = file.path("figures", paste0("Fig3.", ext)),
    plot = Fig3,
    width = 170, height = 150, units = "mm"
  )
}

## Figure 4 - Technical Recovery Differences ----

# Create a plot of the Genomad confidence scores for just the contigs that were
plt_confidence_filtered <- df_Genomad_filtered %>%
  ggplot(aes(x = virus_score, color = DescriptionLong, group = Sample)) +  # use group = Sample to show 3 lines per condition
  stat_ecdf(geom = "step", alpha = 0.7, linewidth = 0.5) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right") +
  labs(
    x = "Genomad Prediction Confidence",
    y = "Cumulative\nFrequency",
    color = "Processing Method"
  ) +
  scale_color_manual(values = Description_palette)

# Create Fig4 as a plot of the read depth data
Fig4b <- df_read_depth_100 %>%
  ggplot(aes(x = Depth, y = Percentage, colour = DescriptionLong, fill = DescriptionLong)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_smooth() +
  geom_vline(xintercept = c(3.2, 10), linetype = "dashed", colour = "gray50") +
  labs(x = "Sequencing Depth (Gbp)", 
       y = "vOTUs With\n>= 100x Coverage (%)",
       fill = "Processing Method",
       colour = "Processing Method") +
  scale_color_manual(values = Description_palette) +
  scale_fill_manual(values = Description_palette) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2), colour = guide_legend(nrow = 2))

# Combine both panels into Fig4
Fig4 <- ggarrange(plt_confidence_filtered, Fig4b, labels = c("a","b"),
                  common.legend = TRUE, legend = "bottom", ncol = 2, nrow = 1)

# Save Fig4 as a PDF (note, this was lightly edited in Affinity Designer for
# label legibility)
ggsave("figures/Fig4.pdf", plot = Fig4,
       width = 170, height = 60, units = "mm", device = "pdf")

# Save Fig4 as a SVG
ggsave("figures/Fig4.svg", plot = Fig4,
       width = 170, height = 60, units = "mm", device = "svg")

# 4 - Supplementary Figure Generation ----
## Figure S1 - Experimental design ----
### Note - this was made in Biorender ###
### This section is here as a placeholder ###

## Figure S2 – Contig data ----

# Individual panels
plt_contigs        <- make_quast_plot(df_quast, "# contigs",           "Contigs",           "Number of contigs")
plt_viral_contigs  <- make_quast_plot(df_quast, "FilteredViralContigs", "Viral contigs",     "Number of viral contigs")
plt_largest        <- make_quast_plot(df_quast, "Largest contig",       "Largest contig (bp)", "Largest contig")
plt_auN            <- make_quast_plot(df_quast, "auN",                  "auN",               "auN")
plt_N50            <- make_quast_plot(df_quast, "N50",                  "N50",               "N50")
plt_L50            <- make_quast_plot(df_quast, "L50",                  "L50",               "L50")
plt_total_length   <- make_quast_plot(df_quast, "Total length",         "Total length (bp)", "Total length")

# Assemble into Fig S2
FigS2 <- ggarrange(
  plt_contigs, plt_largest, plt_auN,
  plt_N50, plt_L50, plt_total_length,
  ncol = 2, nrow = 3,
  labels = c("a", "b", "c", "d", "e", "f"),
  common.legend = TRUE, legend = "bottom"
)

# Save Figure S2
for (ext in c("pdf", "svg")) {
  ggsave(
    filename = file.path("figures", paste0("FigS2.", ext)),
    plot = FigS2,
    width = 170, height = 170, units = "mm"
  )
}

## Figure S3 – Ratio of mapped to assembled ----

# Venn diagram of mapped-only vOTUs
plt_mapped_only <- ls_Venns$plt_MappedOnly

# Merge and calculate ratio of mapped to assembled vOTUs
df_all_viral <- df_mapped_viruses %>%
  merge(df_assembled_viruses, by = c("Sample", "Description")) %>%
  rename(
    AssembledViruses = FilteredViralContigs,
    MappedViruses = Count
  ) %>%
  mutate(Ratio = MappedViruses / AssembledViruses)

# Create plot list for ratio of mapped to assembled vOTUs
plt_list_ratio <- plot_means_with_cld(
  df_all_viral,
  "Ratio",
  "Description",
  Description_palette
)

# Plot ratio
plt_ratio <- plt_list_ratio$plot +
  labs(
    y = "Mapped / Assembled Viruses",
    title = "Mapped vOTUs / Assembled viral contigs"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Assemble Figure S3
FigS3 <- ggarrange(
  plt_mapped_only, plt_ratio,
  ncol = 1, labels = c("a", "b"),
  common.legend = TRUE, legend = "bottom"
)

# Save Figure S3
for (ext in c("pdf", "svg")) {
  ggsave(
    filename = file.path("figures", paste0("FigS3.", ext)),
    plot = FigS3,
    width = 160, height = 160, units = "mm"
  )
}

## Figure S4 – Read-based data ----

# Plot total reads per sample
plt_list_raw_reads <- plot_means_with_cld(
  df_reads %>%
    select(DescriptionLong, total_sequences) %>%
    distinct(),
  "total_sequences", "DescriptionLong", Description_palette
)

plt_raw_reads <- plt_list_raw_reads$plot +
  labs(y = "Total reads (millions)", title = "Total reads", fill = "Processing Method") +
  theme_bw(base_size = 10) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Plot human reads
plt_list_human_reads <- plot_means_with_cld(
  df_reads %>% filter(Type == "Human"),
  "ReadPercentage", "DescriptionLong", Description_palette
)

plt_human_reads <- plt_list_human_reads$plot +
  labs(y = "Percentage (%)", title = "Human reads", fill = "Processing Method") +
  theme_bw(base_size = 10) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Plot rRNA reads
plt_list_rRNA_reads <- plot_means_with_cld(
  df_reads %>% filter(Type == "rRNA"),
  "ReadPercentage", "DescriptionLong", Description_palette
)

plt_rRNA_reads <- plt_list_rRNA_reads$plot +
  labs(y = "Percentage (%)", title = "rRNA gene reads", fill = "Processing Method") +
  theme_bw(base_size = 10) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Unique k-mers
df_unique_kmers <- tar_read(processed_kmers) %>%
  filter(abundance == 1) %>%
  mutate(
    DescriptionLong = DescriptionModifiers[Description],
    DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels)
  )

plt_list_unique_kmers <- plot_means_with_cld(df_unique_kmers, "count", "DescriptionLong", Description_palette)

plt_unique_kmers <- plt_list_unique_kmers$plot +
  labs(y = "Unique kmers", title = "Unique kmers", fill = "Processing Method") +
  theme_bw(base_size = 10) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# Assemble Figure S4
FigS4 <- ggarrange(
  plt_raw_reads, plt_rRNA_reads,
  plt_human_reads, plt_unique_kmers,
  ncol = 2, nrow = 2,
  labels = c("a", "b", "c", "d"),
  common.legend = TRUE, legend = "bottom"
)

# Save Figure S4
for (ext in c("pdf", "svg")) {
  ggsave(
    filename = file.path("figures", paste0("FigS4.", ext)),
    plot = FigS4,
    width = 160, height = 160, units = "mm"
  )
}

## Figure S5 – Mean Bray–Curtis dissimilarities ----

# Convert Bray–Curtis distance matrix to long format
df_bray_long <- as.matrix(df_bray) %>%
  melt(varnames = c("Sample1", "Sample2"), value.name = "BrayCurtis") %>%
  filter(Sample1 != Sample2)  # Remove diagonal

# Add Processing Method metadata
df_bray_long <- df_bray_long %>%
  left_join(metadata %>% select(Sample1 = Sample, SampleType1 = Description), by = "Sample1") %>%
  left_join(metadata %>% select(Sample2 = Sample, SampleType2 = Description), by = "Sample2") %>%
  mutate(
    SampleType1 = DescriptionModifiers[SampleType1],
    SampleType2 = DescriptionModifiers[SampleType2]
  )

# Remove duplicate pairs (A–B vs B–A)
df_bray_long <- df_bray_long %>%
  rowwise() %>%
  mutate(pair_id = paste(sort(c(Sample1, Sample2)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair_id, .keep_all = TRUE)

# Create within-group comparisons
within_group_long <- df_bray_long %>%
  filter(SampleType1 == SampleType2) %>%
  mutate(
    Comparison = SampleType1,
    ComparisonType = "Within"
  )

# Create between-group comparisons
between_group_long <- df_bray_long %>%
  filter(SampleType1 != SampleType2) %>%
  rowwise() %>%
  mutate(
    Comparison = paste(sort(c(SampleType1, SampleType2)), collapse = " vs "),
    ComparisonType = "Between"
  ) %>%
  ungroup()

# Combine within- and between-group tables
df_bray_comparisons <- bind_rows(within_group_long, between_group_long)

# Summarise Bray–Curtis dissimilarities
df_bray_summary <- df_bray_comparisons %>%
  group_by(Comparison, ComparisonType) %>%
  summarise(
    Mean_BrayCurtis = mean(BrayCurtis, na.rm = TRUE),
    SD_BrayCurtis   = sd(BrayCurtis, na.rm = TRUE),
    Count           = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Comparison = str_replace_all(as.character(Comparison), "\n", " "),
    Comparison = str_replace_all(Comparison, " vs ", " vs\n"),
    Comparison = str_replace_all(Comparison, "Meta genome", "Metagenome"),
    Comparison = factor(Comparison, levels = unique(Comparison)),
    ComparisonType = factor(ComparisonType, levels = c("Within", "Between"))
  )

# Plot Bray–Curtis comparisons
plt_bray_comparisons <- ggplot(df_bray_summary, aes(x = Comparison, y = Mean_BrayCurtis, fill = ComparisonType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = Mean_BrayCurtis - SD_BrayCurtis, ymax = Mean_BrayCurtis + SD_BrayCurtis),
    width = 0.2, position = position_dodge(width = 0.9)
  ) +
  labs(x = "", y = "Mean Bray–Curtis\nDissimilarity", fill = "Comparison Type") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = c("Within" = "#56B4E9", "Between" = "#D55E00"))

# Save Figure S5
for (ext in c("pdf", "svg")) {
  ggsave(
    filename = file.path("figures", paste0("FigS5.", ext)),
    plot = plt_bray_comparisons,
    width = 160, height = 100, units = "mm"
  )
}
## Figure S6 – Annotation ----

# Prepare GeNomad data (unfiltered and filtered)
df_Genomad_unfiltered <- df_Genomad %>%
  rename(Contig = seq_name, Sample = sample) %>%
  left_join(metadata, by = "Sample") %>%
  mutate(
    DescriptionLong = DescriptionModifiers[Description],
    DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels)
  ) %>%
  filter(!is.na(DescriptionLong))

# Identify removed contigs (not retained after filtering)
df_Genomad_removed <- df_Genomad_unfiltered %>%
  filter(!Contig %in% df_Genomad_filtered$Contig)

# ECDF plots for Genomad prediction confidence
plt_confidence_unfiltered <- df_Genomad_unfiltered %>%
  ggplot(aes(x = virus_score, colour = DescriptionLong, group = Sample)) +
  stat_ecdf(geom = "step", alpha = 0.7, linewidth = 0.5) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(
    x = "GeNomad Prediction\nConfidence",
    y = "Cumulative Frequency",
    colour = "Processing Method",
    title = "All viral\ncontigs"
  ) +
  scale_colour_manual(values = Description_palette)

plt_confidence_removed <- df_Genomad_removed %>%
  ggplot(aes(x = virus_score, colour = DescriptionLong, group = Sample)) +
  stat_ecdf(geom = "step", alpha = 0.7, linewidth = 0.5) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(
    x = "GeNomad Prediction\nConfidence",
    y = "Cumulative Frequency",
    colour = "Processing Method",
    title = "Contigs removed\nby filtering"
  ) +
  scale_colour_manual(values = Description_palette)

plt_confidence_vOTU <- df_Genomad_filtered %>%
  filter(Contig %in% df_clustering$vOTU) %>%
  ggplot(aes(x = virus_score, colour = DescriptionLong, group = Sample)) +
  stat_ecdf(geom = "step", alpha = 0.7, linewidth = 0.5) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(
    x = "GeNomad Prediction\nConfidence",
    y = "Cumulative Frequency",
    colour = "Processing Method",
    title = "vOTU representative\ncontigs"
  ) +
  scale_colour_manual(values = Description_palette)

# Combine unfiltered and removed ECDF plots
plt_confidence_combined <- ggarrange(
  plt_confidence_unfiltered + theme(legend.position = "none"),
  plt_confidence_removed + theme(legend.position = "none"),
  ncol = 2, labels = c("a", "b"),
  common.legend = TRUE, legend = "bottom"
)

# Summarise mean confidence per sample and test
df_summary <- df_Genomad_filtered %>%
  group_by(DescriptionLong, Sample) %>%
  summarise(mean_confidence = mean(virus_score, na.rm = TRUE), .groups = "drop")

aov_model <- aov(mean_confidence ~ DescriptionLong, data = df_summary)
summary(aov_model)
TukeyHSD(aov_model)

# Load AMR data
df_AMR <- read.csv("data/ViralContigs/merged_amr_filtered.tsv", sep = "\t") %>%
  mutate(Contig = sub("_[^_]+$", "", gene)) %>%
  filter(Contig %in% all_filtered_contigs)

# AMR subset for vOTU representatives
df_AMR_vOTU <- df_AMR %>% filter(Contig %in% df_mapping$Contig)

# Identify AMR-associated vOTUs and clusters
AMR_vOTUs <- unique(df_AMR_vOTU$Contig)
AMR_clusters <- df_clustering %>% filter(vOTU %in% AMR_vOTUs)

# Load Pharokka data
df_Pharokka <- tar_read(PharokkaData)

# Subset for moron / AMG / host takeover genes
df_morons <- df_Pharokka %>%
  filter(category == "moron, auxiliary metabolic gene and host takeover") %>%
  rename(vOTU = Contig)

# Aggregate Pharokka annotations by category
df_Pharokka_category <- df_Pharokka %>%
  rename(Sample = sample) %>%
  count(Sample, category, name = "Count") %>%
  left_join(metadata, by = "Sample") %>%
  mutate(
    DescriptionLong = DescriptionModifiers[Description],
    DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels)
  )

# Compute per-sample percentage of moron/AMG/host takeover genes
df_Pharokka_moron <- df_Pharokka_category %>%
  group_by(Sample, DescriptionLong) %>%
  summarise(
    TotalGenes = sum(Count),
    TotalMoron = sum(Count[category == "moron, auxiliary metabolic gene and host takeover"]),
    Percentage = (TotalMoron / TotalGenes) * 100,
    .groups = "drop"
  )

# Summarise across samples
df_Pharokka_moron_summary <- df_Pharokka_moron %>%
  group_by(DescriptionLong) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    .groups = "drop"
  )

# Plot moron/AMG/host takeover gene percentage
plt_list_morons <- plot_means_with_cld(df_Pharokka_moron, "Percentage", "DescriptionLong", Description_palette)

plt_morons <- plt_list_morons$plot +
  labs(
    y = "Viral contigs carrying a\nmoron, auxiliary metabolic gene\nand host takeover gene (%)",
    fill = "Processing Method"
  ) +
  scale_fill_manual(values = Description_palette) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

# Final plot
FigS6 <- plt_morons

# Save Figure S6
for (ext in c("pdf", "svg")) {
  ggsave(
    filename = file.path("figures", paste0("FigS6.", ext)),
    plot = FigS6,
    width = 170, height = 100, units = "mm"
  )
}
## Figure S7 – Unified Human Gut Virome (UHGV) Catalog comparison ----

# Load and process clustering data
df_UHGV <- read.csv("data/Biogeography/AllContigs_UHGV_all_clusters.tsv", header = FALSE, sep = "\t") %>%
  rename(Contig = V1, Cluster = V2) %>%
  filter(grepl("contig", Cluster)) %>%
  mutate(
    UHGV = grepl("UHGV", Cluster),
    NumSeqs = str_count(Cluster, ",") + 1,
    NumUHGVs = str_count(Cluster, "UHGV"),
    NumStudy = NumSeqs - NumUHGVs
  )

# Plot number of clusters containing UHGV vs not
plt_UHGV_Clusters <- df_UHGV %>%
  group_by(UHGV) %>%
  summarise(NumClusters = n_distinct(Cluster), .groups = "drop") %>%
  ggplot(aes(x = UHGV, y = NumClusters, fill = UHGV)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = NumClusters), vjust = 1.5, size = 3) +
  labs(x = "Cluster contains a UHGV vOTU", y = "Clusters") +
  scale_fill_manual(values = c("TRUE" = "#F5C875", "FALSE" = "#A6D5F5")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", panel.grid = element_blank())

# Identify contig sets
UHGV_filtered_contigs <- df_UHGV %>%
  filter(UHGV) %>%
  summarise(all_clusters = paste(Cluster, collapse = ",")) %>%
  pull(all_clusters) %>%
  strsplit(",") %>%
  unlist() %>%
  intersect(all_filtered_contigs)

NonUHGV_filtered_contigs <- df_UHGV %>%
  filter(!UHGV) %>%
  summarise(all_clusters = paste(Cluster, collapse = ",")) %>%
  pull(all_clusters) %>%
  strsplit(",") %>%
  unlist() %>%
  intersect(all_filtered_contigs)

UHGV_contigs <- df_UHGV %>%
  filter(UHGV) %>%
  summarise(all_clusters = paste(Cluster, collapse = ",")) %>%
  pull(all_clusters) %>%
  strsplit(",") %>%
  unlist() %>%
  .[grepl("UHGV", .)]

# Calculate percentage of UHGV-linked contigs per processing method
df_UHGVs <- df_Genomad_filtered %>%
  mutate(UHGV = Contig %in% UHGV_filtered_contigs) %>%
  left_join(metadata, by = "Sample") %>%
  group_by(Sample, DescriptionLong) %>%
  summarise(
    UHGV = sum(UHGV, na.rm = TRUE),
    Total = n(),
    Percentage = (UHGV / Total) * 100,
    .groups = "drop"
  )

# ANOVA
anova_result <- aov(Percentage ~ DescriptionLong, data = df_UHGVs)
summary(anova_result)

# Plot percentage UHGV data
plt_list_UHGVs <- plot_means_with_cld(df_UHGVs, "Percentage", "DescriptionLong", Description_palette)
plt_UHGVs <- plt_list_UHGVs$plot +
  labs(x = "Processing Method", y = "Percentage of\nviral contigs (%)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

# Combine cluster and percentage plots
FigS7_Top <- ggarrange(plt_UHGV_Clusters, plt_UHGVs, ncol = 2, labels = c("a", "b"))

# Load mapping data
folder <- "data/Biogeography/Mapping_UHGV/"
files <- list.files(folder, full.names = TRUE)
mean_files <- files[grepl("_mean", basename(files), ignore.case = TRUE)]
non_mean_files <- files[!grepl("_mean", basename(files), ignore.case = TRUE)]

process_files <- function(file_list) {
  file_list %>%
    set_names() %>%
    map_dfr(~ {
      df <- read_tsv(.x, col_types = cols(.default = "c"))
      last_col <- names(df)[ncol(df)]
      df <- df %>% filter(.data[[last_col]] != "0")  # use tidy eval to select last column
      df$SourceFile <- basename(.x)
      df
    }) %>%
    select(-Sample) %>%
    rename(Sample = SourceFile) %>%
    mutate(Sample = gsub(".tsv", "", Sample))
}

df_UHGV_mean <- process_files(mean_files) %>% mutate(Sample = gsub("_mean", "", Sample))
df_UHGV_nonmean <- process_files(non_mean_files)

df_UHGV_mapping <- inner_join(df_UHGV_mean, df_UHGV_nonmean, by = c("Contig", "Sample")) %>%
  filter(`Covered Fraction` >= 0.75)

MappedOnlyUHGVs <- unique(df_UHGV_mapping$Contig[!df_UHGV_mapping$Contig %in% UHGV_contigs])
df_UHGV_mapping <- df_UHGV_mapping %>% mutate(MappedOnlyUHGV = Contig %in% MappedOnlyUHGVs)

df_UHGV_mapping_only <- df_UHGV_mapping %>% filter(MappedOnlyUHGV)

# Metadata for UHGVs
# Note that this file is not hosted on GitHub due to size constraints but has been added
# to the Zenodo archive - see README for details
df_UHGV_metadata <- read.csv("data/Biogeography/votus_metadata_extended.tsv", sep = "\t") %>%
  filter(uhgv_genome %in% c(MappedOnlyUHGVs, UHGV_contigs)) %>%
  mutate(MappedOnlyUHGV = uhgv_genome %in% MappedOnlyUHGVs)

# Summary for mapped UHGVs
df_UHGV_mapping_summary <- df_UHGV_mapping %>%
  filter(`Covered Fraction` >= 0.75) %>%
  count(Sample, name = "Count") %>%
  left_join(metadata, by = "Sample") %>%
  mutate(
    DescriptionLong = DescriptionModifiers[Description],
    DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels)
  )

df_UHGV_mapping_only_summary <- df_UHGV_mapping_only %>%
  count(Sample, name = "Count") %>%
  left_join(metadata, by = "Sample") %>%
  mutate(
    DescriptionLong = DescriptionModifiers[Description],
    DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels)
  )

# Plot UHGV mapping summaries
plt_UHGV_mapping <- plot_means_with_cld(df_UHGV_mapping_summary, "Count", "DescriptionLong", Description_palette)$plot +
  labs(x = "Processing Method", y = "Count") +
  scale_y_continuous(limits = c(0, 265), expand = expansion(mult = c(0, 0.1))) +
  theme(panel.grid = element_blank(), legend.position = "bottom")

plt_UHGV_mapping_only <- plot_means_with_cld(df_UHGV_mapping_only_summary, "Count", "DescriptionLong", Description_palette)$plot +
  labs(x = "Processing Method", y = "Count") +
  scale_y_continuous(limits = c(0, 265), expand = expansion(mult = c(0, 0.1))) +
  theme(panel.grid = element_blank(), legend.position = "bottom")

FigS7_Middle <- ggarrange(plt_UHGV_mapping, plt_UHGV_mapping_only, ncol = 2, labels = c("c", "d"), common.legend = TRUE)

# Coverage comparison
df_UHGV_coverage <- df_UHGV_mapping %>%
  group_by(Contig, MappedOnlyUHGV) %>%
  summarise(MaxMean = max(as.numeric(`Trimmed Mean`)), .groups = "drop")

wilcox_result <- wilcox.test(MaxMean ~ MappedOnlyUHGV, data = df_UHGV_coverage)
formatted_p <- if (wilcox_result$p.value < 1e-5) "p < 0.00001" else paste0("p = ", signif(wilcox_result$p.value, 2))

plt_UHGV_coverage <- ggplot(df_UHGV_coverage, aes(x = MappedOnlyUHGV, y = MaxMean, fill = MappedOnlyUHGV)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_jitter(aes(color = MappedOnlyUHGV), width = 0.2, alpha = 0.5, size = 0.5) +
  scale_fill_manual(values = c("TRUE" = "#F5C875", "FALSE" = "#A6D5F5")) +
  scale_color_manual(values = c("TRUE" = "#F5C875", "FALSE" = "#A6D5F5")) +
  scale_y_log10(labels = scales::label_log()) +
  annotate("text", x = 2, y = max(df_UHGV_coverage$MaxMean, na.rm = TRUE) * 1.1, label = paste("Wilcox test,", formatted_p), size = 4) +
  labs(x = "Only detected by mapping", y = "Coverage") +
  theme_bw(base_size = 10) + theme(panel.grid = element_blank(), legend.position = "bottom")

# Genome length comparison
wilcox_result_length <- wilcox.test(genome_length ~ MappedOnlyUHGV, data = df_UHGV_metadata)
formatted_p_length <- if (wilcox_result_length$p.value < 1e-5) "p < 0.00001" else paste0("p = ", signif(wilcox_result_length$p.value, 2))

plt_UHGV_length <- df_UHGV_metadata %>%
  ggplot(aes(x = MappedOnlyUHGV, y = genome_length, fill = MappedOnlyUHGV)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_jitter(aes(color = MappedOnlyUHGV), size = 0.5, alpha = 0.5, width = 0.2) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_y_log10(labels = scales::label_log()) +
  annotate("text", x = 1.85, y = max(df_UHGV_metadata$genome_length, na.rm = TRUE) * 1.2,
           label = paste("Wilcox test,", formatted_p_length), size = 4) +
  scale_fill_manual(values = c("TRUE" = "#F5C875", "FALSE" = "#A6D5F5")) +
  scale_color_manual(values = c("TRUE" = "#F5C875", "FALSE" = "#A6D5F5")) +
  labs(y = "Genome Length", x = "Only detected by mapping") +
  theme_bw(base_size = 10) + theme(panel.grid = element_blank(), legend.position = "bottom")

FigS7_Bottom <- ggarrange(plt_UHGV_coverage, plt_UHGV_length, ncol = 2, labels = c("e", "f"), common.legend = TRUE)

# Final assembly
FigS7 <- ggarrange(FigS7_Top, FigS7_Middle, FigS7_Bottom, nrow = 3)

# Save Figure S7
for (ext in c("pdf", "svg")) {
  ggsave(
    filename = file.path("figures", paste0("FigS7.", ext)),
    plot = FigS7,
    width = 170, height = 170, units = "mm"
  )
}

# Summary tables
df_UHGVs_summary_table <- df_UHGVs %>%
  group_by(DescriptionLong) %>%
  summarise(
    Mean_Percentage = mean(Percentage, na.rm = TRUE),
    SD_Percentage = sd(Percentage, na.rm = TRUE),
    Count = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Percentage))

df_UHGV_mapping_summary_table <- df_UHGV_mapping_summary %>%
  group_by(DescriptionLong) %>%
  summarise(
    Mean_Count = mean(Count, na.rm = TRUE),
    SD_Count = sd(Count, na.rm = TRUE),
    Count = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Count))
# Stats tables for SI ----

## Define panels to include in summary tables ----
panels <- list(
  # Fig 1c – Assembled vOTUs
  list(df = df_assembled_viruses, value = "FilteredViralContigs", group = "Description",
       panel = "Fig1c_Assembled_vOTUs", measure = "Assembled vOTUs"),
  
  # Fig 1d – Mapped vOTUs
  list(df = df_mapped_viruses, value = "Count", group = "Description",
       panel = "Fig1d_Mapped_vOTUs", measure = "Mapped vOTUs"),
  
  # Fig 1e – Viral reads (%)
  list(df = dplyr::filter(df_reads, Type == "Viral"), value = "ReadPercentage", group = "DescriptionLong",
       panel = "Fig1e_Viral_Reads", measure = "Viral reads (%)"),
  
  # Fig 3b – Provirus proportions
  list(df = df_Provirus_mapped_summary %>%
         filter(Provirus == "Provirus") %>%
         mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels)),
       value = "Value", group = "DescriptionLong",
       panel = "Fig3b_Provirus_Proportion", measure = "Provirus proportion (%)"),
  
  # Fig S2 – Assembly metrics
  list(df = filter(df_quast, Metric == "# contigs"),
       value = "Value", group = "DescriptionLong", panel = "FigS2a_NumContigs", measure = "# contigs"),
  list(df = filter(df_quast, Metric == "FilteredViralContigs"),
       value = "Value", group = "DescriptionLong", panel = "FigS2b_ViralContigs", measure = "Viral contigs"),
  list(df = filter(df_quast, Metric == "Largest contig"),
       value = "Value", group = "DescriptionLong", panel = "FigS2c_LargestContig", measure = "Largest contig (bp)"),
  list(df = filter(df_quast, Metric == "auN"),
       value = "Value", group = "DescriptionLong", panel = "FigS2d_auN", measure = "auN"),
  list(df = filter(df_quast, Metric == "N50"),
       value = "Value", group = "DescriptionLong", panel = "FigS2e_N50", measure = "N50"),
  list(df = filter(df_quast, Metric == "L50"),
       value = "Value", group = "DescriptionLong", panel = "FigS2f_L50", measure = "L50"),
  list(df = filter(df_quast, Metric == "Total length"),
       value = "Value", group = "DescriptionLong", panel = "FigS2g_TotalLength", measure = "Total length (bp)"),
  
  # Fig S3 – Ratio of mapped to assembled
  list(df = df_all_viral, value = "Ratio", group = "Description",
       panel = "FigS3b_Mapped_to_Assembled_Ratio", measure = "Mapped/Assembled"),
  
  # Fig S4 – Read-based metrics
  list(df = df_reads %>% select(DescriptionLong, total_sequences) %>% distinct(),
       value = "total_sequences", group = "DescriptionLong",
       panel = "FigS4a_Total_Reads", measure = "Total reads (millions)"),
  list(df = filter(df_reads, Type == "rRNA"),
       value = "ReadPercentage", group = "DescriptionLong",
       panel = "FigS4b_rRNA_Reads", measure = "rRNA reads (%)"),
  list(df = filter(df_reads, Type == "Human"),
       value = "ReadPercentage", group = "DescriptionLong",
       panel = "FigS4c_Human_Reads", measure = "Human reads (%)"),
  list(df = df_unique_kmers, value = "count", group = "DescriptionLong",
       panel = "FigS4d_UniqueKmers", measure = "Unique kmers"),
  
  # Fig S6 – Moron/AMG/host-takeover genes
  list(df = df_Pharokka_moron, value = "Percentage", group = "DescriptionLong",
       panel = "FigS6a_Moron_AMG_HostTakeover", measure = "Moron/AMG/host takeover (%)"),
  
  # Fig S7 – UHGV comparison
  list(df = df_UHGVs, value = "Percentage", group = "DescriptionLong",
       panel = "FigS7b_UHGV_Percentage", measure = "% contigs in UHGV clusters"),
  list(df = df_UHGV_mapping_summary, value = "Count", group = "DescriptionLong",
       panel = "FigS7c_UHGV_Mapping_Counts", measure = "UHGV mapping count"),
  list(df = df_UHGV_mapping_only_summary, value = "Count", group = "DescriptionLong",
       panel = "FigS7d_UHGV_MappingOnly_Counts", measure = "UHGV mapping-only count")
)

## Run ANOVA + Tukey tests for all panels ----
stats_list <- lapply(panels, function(p) {
  one_way_stats(p$df, p$value, p$group, p$panel, p$measure)
})

anova_tbl <- bind_rows(lapply(stats_list, `[[`, "anova")) %>%
  mutate(
    F_formatted = ifelse(is.na(F), NA_character_,
                         sprintf("F[%d,%d] = %.3f", df1, df2, F)),
    p_formatted = case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 1e-5 ~ "p < 0.00001",
      TRUE ~ paste0("p = ", signif(p_value, 3))
    )
  ) %>%
  select(panel, measure, k_groups, df1, df2, F, p_value, F_formatted, p_formatted)

tukey_tbl <- bind_rows(lapply(stats_list, `[[`, "tukey")) %>%
  mutate(
    p_formatted = case_when(
      is.na(p_adj) ~ NA_character_,
      p_adj < 1e-5 ~ "p < 0.00001",
      TRUE ~ paste0("p = ", signif(p_adj, 3))
    ),
    contrast = stringr::str_squish(stringr::str_replace_all(contrast, "[\r\n]+", " "))
  ) %>%
  tidyr::separate(
    contrast,
    into = c("group1", "group2"),
    sep = " - ", remove = FALSE,
    extra = "merge", fill = "right"
  ) %>%
  arrange(panel, measure, contrast)

## Export statistical summaries ----
if (!dir.exists("stats")) dir.create("stats", recursive = TRUE)
readr::write_tsv(anova_tbl, "stats/anova_summary.tsv")
readr::write_tsv(tukey_tbl, "stats/tukey_summary.tsv")

# Graphical abstract ----
plt_stacked_GA <- plt_stacked +
  theme_bw(base_size = 15) +
  labs(fill = "Viral\ntaxonomy") +
  guides(fill = guide_legend(title = "Viral\ntaxonomy", ncol = 1)) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )

plt_microbiome_GA <- plt_microbiome +
  theme_bw(base_size = 15) +
  labs(y = "Relative\nAbundance (%)") +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    legend.justification = c("center", "top"),
    legend.box.margin = margin(t = -80, r = 0, b = 0, l = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )

GraphicalAbstract <- ggarrange(
  plt_stacked_GA, plt_microbiome_GA,
  ncol = 1, nrow = 2
)

ggsave(
  "figures/GraphicalAbstractFigure.pdf",
  plot = GraphicalAbstract,
  width = 200, height = 150, units = "mm", device = "pdf"
)

