# Fecal Viromes Optimisation Analysis Parent Script
# This pipeline will run all the necessary steps to generate the figures and
# tables for the Fecal Viromes Optimisation Analysis

# 0 - Setup ----

## Clear environment
rm(list=ls())

## Load libraries
library(ggh4x)          # For extended ggplot2 functionality
library(ggpubr)         # For figure generation
library(multcompView)   # For pairwise comparisons
library(pairwiseAdonis) # For pairwise adonis tests
library(reshape2)       # For reshaping data
library(targets)        # For pipeline management
library(tidyverse)      # For data manipulation and visualisation
library(vegan)          # For beta diversity calculations

## Load custom functions
source("scripts/Functions_FVO.R") 

# 1 - Run targets pipeline for Fecal Viromes Optimisation Analysis ----

# Set new names for mutated Sample Types
DescriptionModifiers <- c("UTV" = "Untreated\nVirome",
                          "DTV" = "DNase\nTreated\nVirome",
                          "MDA" = "MDA\nVirome",
                          "MtG" = "Meta\ngenome")

# Set levels for mutated Sample Types
DescriptionLevels <- c("Untreated\nVirome","DNase\nTreated\nVirome", "MDA\nVirome", "Meta\ngenome")

# Set color palette for mutated Sample Types
Description_palette <- c("Untreated\nVirome" = "#56B4E9",
                         "DNase\nTreated\nVirome" = "#0072B2",
                         "MDA\nVirome" = "#009E73",
                         "Meta\ngenome" = "#FF7F0E")

## Visualise pipeline
tar_visnetwork()

## Run pipeline
tar_make()

# 2 - load and wrangle data ----
## Load required objects for figure production
# Import metadata
metadata <- tar_read(metadata) # Metadata

# Import read counts
df_mapping <- tar_read(vOTU_tpm) # vOTU relative abundance data

# Import read counts for viral reads
df_mapping_cleaned <- df_mapping %>%
  rename(vOTU = Contig) %>%
  # Pivot longer
  pivot_longer(cols = -vOTU, names_to = "Sample", values_to = "TPM") %>%
  mutate(RelAbund = TPM / 10000)

# Import read counts
df_reads <- tar_read(all_read_counts) %>% # Raw read counts
  mutate(DescriptionLong = DescriptionModifiers[Description],
         total_sequences = total_sequences / 1e6) %>%
  mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels))

# Import quast results that also contains viral contig counts
df_quast <- read.csv("data/Quast/Quast_report.tsv", sep = "\t") %>%
  # Remove "_renamed_contigs" from column names
  rename_with(~str_remove(., "_renamed_contigs"), contains("_renamed_contigs")) %>%
  # transpose the df and make column 1 the new column headers
  t() %>%
  as.data.frame() %>%
  # First row to column names
  setNames(.[1,]) %>%
  # Remove first row
  .[-1,] %>%
  # Convert rownames to a column called Sample
  rownames_to_column(var = "Sample") %>%
  # Merge with metadata
  merge(metadata, by = "Sample") %>%
  select(Sample, Description, `# contigs`, `Total length`, `Largest contig`, N50, L50, auN, FilteredViralContigs) %>%
  # Make all columns to numeric except columns Description and Sample
  mutate(across(-c(Sample, Description), as.numeric)) %>%
  pivot_longer(cols = -c(Sample, Description), names_to = "Metric", values_to = "Value") %>%
  mutate(DescriptionLong = DescriptionModifiers[Description]) %>%
  mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels))


# Calculate mean reads for in-text summary 

# Calculate mean raw reads for each sample type
df_read_means <- df_reads %>%
  select(Description, Sample, total_sequences) %>%
  unique() %>%
  summarise(
    Mean = mean(total_sequences / 1e6, na.rm = TRUE),
    SD = sd(total_sequences / 1e6, na.rm = TRUE)
  )

# Calculate mean viral reads and richness for each sample type
df_viral_means <- df_reads %>%
  select(Description, Sample, total_sequences, FilteredViralContigs, ViralReads) %>%
  unique() %>%
  group_by(Description) %>%
  summarise(
    mean_viral_reads = mean(ViralReads/ 1e6, na.rm = TRUE),
    mean_viral_richness = mean(FilteredViralContigs, na.rm = TRUE),
    mean_total_sequences = mean(total_sequences / 1e6, na.rm = TRUE),
    SD_viral_reads = sd(ViralReads / 1e6, na.rm = TRUE),
    SD_viral_richness = sd(FilteredViralContigs, na.rm = TRUE),
    SD_total_sequences = sd(total_sequences / 1e6, na.rm = TRUE)
  )

# Calculate mean and SD for read percentages for each sample type
df_mean_read_types <- df_reads %>%
  select(Description, Sample, total_sequences, Type, ReadPercentage) %>%
  unique() %>%
  group_by(Description, Type) %>%
  summarise(
    mean_reads = mean(ReadPercentage, na.rm = TRUE),
    SD_reads = sd(ReadPercentage, na.rm = TRUE)
  )

# 3 - Figure Generation ----
## Figure 1 - How the four methods compare in terms of viral recovery ----

# Plot Venn diagram
ls_Venns <- tar_read(plt_Venns)

# Set the Venn diagram as an object
plt_Venns <- ls_Venns$plt_Venn

# Set the top of Fig1 as plt_Venns
Fig1_top <- plt_Venns

# Create a df of just the mapped virus data
df_mapped_viruses <- df_mapping_cleaned %>%
  filter(RelAbund!=0) %>%
  # Count the number of viruses in each sample
  group_by(Sample) %>%
  summarise(
    # Count the number of viruses in each sample
    Count = n()) %>%
  # Merge with metadata
  merge(metadata, by = "Sample") %>%
  mutate(Description = DescriptionModifiers[Description]) %>%
  mutate(Description = factor(Description, levels = DescriptionLevels)) %>%
  select(Sample, Description, Count)

# Create a df of just the assembled virus data
df_assembled_viruses <- df_quast %>%
  pivot_wider(
    names_from = Metric,
    values_from = Value
  ) %>%
  select(Sample, DescriptionLong, FilteredViralContigs) %>%
  rename(Description= DescriptionLong)

# Create a summary table        
df_mapped_viruses_summary <- df_mapped_viruses%>%
  select(Sample, Description, Count) %>%
  group_by(Description) %>%
  summarise(
    mean_mapped_viruses = mean(Count, na.rm = TRUE),
    SD_mapped_viruses = sd(Count, na.rm = TRUE))

# Create a copy of df_viral_means and merge with df_mapped_viruses_summary
df_viral_means_2 <- df_viral_means %>% 
  # Change Description
  mutate(Description = DescriptionModifiers[Description]) %>%
  merge(df_mapped_viruses_summary, by = "Description")

# # Create a long version
df_long <- df_viral_means_2 %>%
  pivot_longer(cols = c(mean_viral_richness, mean_mapped_viruses),
    names_to = "Metric",values_to = "Mean") %>%
  pivot_longer(
    cols = c(SD_viral_richness, SD_mapped_viruses),
    names_to = "Metric_SD",
    values_to = "SD") %>%
  filter(str_remove(Metric_SD, "SD_") == str_remove(Metric, "mean_")) %>%
  mutate(Metric = recode(Metric,
                         "mean_viral_richness" = "Viral Richness",
                         "mean_mapped_viruses" = "Mapped Viruses")) %>%
  # Set levels
  mutate(Description = factor(Description, levels = DescriptionLevels)) %>%
  # Change Viral Richness to Assembled Viruses
  mutate(Description = recode(Description,
                              "Viral Richness" = "Assembled Viruses",
                              "Mapped Viruses" = "Mapped Viruses"))

df_wide <- df_long %>% filter(Metric %in% c("Viral Richness", "Mapped Viruses") == TRUE) %>%
  select(Description, Metric, Mean) %>%
  pivot_wider(
    names_from = Metric,
    values_from = Mean
  ) %>% 
  mutate(Ratio = `Mapped Viruses`/ `Viral Richness`)

# Plot the mapped vOTU data
plt_mapped <- plot_means_with_cld(df_mapped_viruses,"Count", "Description", Description_palette) +
  labs(y = "Mapped vOTUs", title = "", fill = "Sample Type") +
  theme_bw() +
  # Remove major and minor grid lines
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # add y limits 0-550
  lims(y = c(0, 525))

# Plot the assembled vOTU data
plt_assembled <- plot_means_with_cld(df_assembled_viruses,"FilteredViralContigs", "Description", Description_palette) +
  labs(y = "Assembled vOTUs", title = "", fill = "Sample Type") +
  theme_bw() +
  # Remove major and minor grid lines
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  lims(y = c(0, 525))

# Assemble Fig1 middle
Fig1_middle <- ggarrange(plt_assembled, plt_mapped, ncol = 2, nrow = 1,
                      labels = c("c", "d"), common.legend = TRUE, legend = "none")

# Import PCoA plots
GenomadData <- tar_read(GenomadData)
stacked <- tar_read(plt_stacked)

# Import PCoA plots
plt_pcoas <- produce_pcoas(df_mapping, stacked,
                           GenomadData, metadata,
                           Description_palette = Description_palette)

# Create a Bray-Curtis dissimilarity matrix
df_bray <- calculate_beta_diversity(df_mapping) %>%
  column_to_rownames(var = "Contig") %>%
  t() %>% vegdist(method = "bray")

# Match group labels to samples in distance matrix
groups <- metadata %>%
  filter(Sample %in% labels(df_bray)) %>%
  arrange(match(Sample, labels(df_bray))) %>%
  pull(Description)  # or DescriptionModifiers if you've modified them


# Calculate distances to group centroids
bd <- betadisper(df_bray, group = groups)

# Test for significant differences in dispersion
anova(bd)
TukeyHSD(bd)

# Plot PCoA with group centroids
pcoa_all <- plt_pcoas$pcoa_all +
  labs(title = "", colour = "Sample Type", fill = "Sample Type", shape = "Sample Type")+
  # make x and y axis 10% bigger
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  coord_fixed()

# Plot viral reads data
plt_viral_reads <- plot_means_with_cld(df_reads %>% filter(Type == "Viral"),"ReadPercentage", "DescriptionLong", Description_palette) +
  labs(y = "Percentage (%)", title = "Viral reads", fill = "Sample Type") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(title = "", y = "Viral reads (%)") +
  theme(legend.position = "none")

# Arrange the plots for Fig1
Fig1_bottom <- ggarrange(plt_viral_reads, pcoa_all, ncol = 2, nrow = 1,
                      labels = c("e", "f"))

# Assemble Fig1
Fig1 <- ggarrange(Fig1_top, Fig1_middle, Fig1_bottom,
                      ncol = 1, nrow = 3,
                      heights = c(1,1,1))

# Save Fig1 as a PDF
ggsave("figures/Fig1.pdf", plot = Fig1,
       width = 170, height = 190, units = "mm", device = "pdf")

# Save Fig1 as a SVG
ggsave("figures/Fig1.svg", plot = Fig1,
       width = 170, height = 190, units = "mm", device = "svg")

## Figure 2 - The MDA-viromes are super different, but we can correct for the differences ----

# Import stacked bar plot of vOTU relative abundance
plt_stacked <- tar_read(plt_stacked)$plt_stacked +
  # Remove legend title
  theme(legend.title = element_blank())

# Assemble Figure 2 top
Fig2_top <-plt_stacked + 
  # put legend on right
  theme(legend.position = "right")

# Plot k-mers data
plt_kmers <- tar_read(plt_kmers) +
  theme_bw() + 
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plt_pcoa_ssdna <- plt_pcoas$pcoa_no_ssDNA +
  labs(colour = "Sample\nType", fill = "Sample\nType", shape = "Sample\nType",
       y = "PC2\n(3.2%)") +
  # Make x and y axes 20% bigger
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.margin = ggplot2::margin(t = 0, r = 40, b = 0, l = 0)) +
  guides(colour = guide_legend(nrow = 1),
         shape = guide_legend(nrow = 1),
         fill = guide_legend(nrow = 1))

# Assemble Figure 2 bottom
Fig2_bottom <- ggarrange(plt_kmers, plt_pcoa_ssdna, ncol = 2, nrow = 1,
                      labels = c("b", "c"),
                      # increase ratio
                      widths = c(1, 2))

Fig2 <- ggarrange(Fig2_top, Fig2_bottom,nrow = 2, labels = c("a",""))

ggsave("figures/Fig2.pdf", plot = Fig2,
       width = 170, height = 100, units = "mm", device = "pdf")

# Save SVG version
ggsave("figures/Fig2.svg", plot = Fig2,
       width = 170, height = 100, units = "mm", device = "svg")

## Figure 3 - The metagenomes recover different viral community members ----
# Import BacPhlip data
df_Bacphlip <- tar_read(BacPhlipData)

# First merge df_BacPhlip and df_mapping by both "vOTU" and "Sample" (left join to keep all BacPhlip data)
df_Lifestyle_mapped_summary <- df_Bacphlip %>% 
  select(-Sample) %>%
  merge(df_mapping_cleaned, by = c("vOTU"), all.y = TRUE) %>%
  merge(metadata, by = "Sample", all.x = TRUE) %>%
  select(vOTU, ShortSamples, RelAbund, status, Provirus, ShortSamples, Description) %>%
  mutate(RelAbund = ifelse(is.na(RelAbund), 0, RelAbund)) %>%
  group_by(status, ShortSamples, Description) %>%
  summarise(RelAbund = sum(RelAbund)) %>% ungroup() %>%
  mutate(status = factor(status, levels = c("Temperate", "Virulent", "Unclassified"))) %>%
  rename(Value = RelAbund) %>%
  # Create a column "rep" from last character from ShortSamples
  mutate(rep = substr(ShortSamples, nchar(ShortSamples), nchar(ShortSamples))) %>%
  # Create a new column DescriptionLong using values in Description and DescriptionModifiers
  mutate(DescriptionLong = DescriptionModifiers[Description])

# Create Lifestyle stacked barplot
plt_Lifestyle <- df_Lifestyle_mapped_summary %>% mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels)) %>%
  ggplot(aes(x = rep, y = Value, fill = status)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Sample", y = "vOTU\nRelative Abundance (%)", fill = "Lifestyle") +
  theme(strip.background = element_rect(fill = "white")) +
  scale_fill_manual(values = Stacked_palette) +
  guides(fill = guide_legend(title = "Lifestyle", nrow = 2)) +
  facet_wrap(~DescriptionLong, drop = TRUE, nrow = 1) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")

# Subset lifestyle data for just proviruses that are predicted to be temperate
df_Provirus_mapped_summary <- df_Bacphlip %>%
  select(-Sample) %>%
  merge(df_mapping_cleaned, by = c("vOTU"), all.y = TRUE) %>%
  merge(metadata, by = "Sample", all.x = TRUE) %>%
  filter(status == "Temperate") %>%
  select(vOTU, ShortSamples, RelAbund, Provirus, ShortSamples, Description) %>%
  mutate(RelAbund = ifelse(is.na(RelAbund), 0, RelAbund)) %>%
  group_by(Provirus, ShortSamples, Description) %>%
  summarise(RelAbund = sum(RelAbund)) %>% ungroup() %>%
  group_by(ShortSamples) %>%
  mutate(RelAbund = RelAbund / sum(RelAbund) * 100) %>%
  mutate(Provirus = factor(Provirus, levels = c("Virus", "Provirus"))) %>%
  rename(Value = RelAbund) %>%
  mutate(Variable = "vOTU Mapped Relative Abundance (%)") %>%
  mutate(rep = substr(ShortSamples, nchar(ShortSamples), nchar(ShortSamples))) %>%
  mutate(DescriptionLong = DescriptionModifiers[Description])

# Create provirus plot
plt_Proviruses <- df_Provirus_mapped_summary %>%
  filter(Provirus == "Provirus",
         Variable == "vOTU Mapped Relative Abundance (%)") %>%
  mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels)) %>%
  plot_means_with_cld("Value", "DescriptionLong", Description_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_bw(base_size = 10) +
  labs(y = "Provirus\nProportion (%)") +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(title = "Sample\nType", nrow = 2)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))

# Assemble Figure 3
Fig3_top <- ggarrange (plt_Lifestyle, plt_Proviruses, ncol = 2, labels = c("a", "b"))

# Import SingleM data
df_SingleM <- tar_read(SingleMData) %>%
  mutate(rep = substr(ShortSamples, nchar(ShortSamples), nchar(ShortSamples))) %>%
  mutate(DescriptionLong = DescriptionModifiers[Description])

# Import host prediction data
df_host <- tar_read(iPhopData)

# Create host prediction stacked barplot
df_host_stacked <- df_host %>%
  merge(df_mapping_cleaned, by = "vOTU", all.y = TRUE) %>%
  select(c(vOTU, Sample, HostPhylum, RelAbund)) %>%
  mutate(HostPhylum = ifelse(is.na(HostPhylum), "Unclassified", HostPhylum)) %>%
  mutate(HostPhylum = factor(HostPhylum,
                             levels = c("Actinomycetota","Bacillota", "Bacteroidota",
                                        "Desulfobacterota", "Pseudomonadota", "Verrucomicrobiota",
                                        "Unclassified"))) %>%
  group_by(Sample, HostPhylum) %>% filter(RelAbund > 0) %>%
  summarise(RelAbund = sum(RelAbund),
            Count = n()) %>%
  merge(metadata, by = "Sample") %>%
  mutate(rep = substr(ShortSamples, nchar(ShortSamples), nchar(ShortSamples))) %>%
  mutate(DescriptionLong = DescriptionModifiers[Description]) %>%
  select(DescriptionLong, HostPhylum, Count, RelAbund, rep) %>%
  mutate(Data = "Viruses (by predicted host)")

# Create a summary df of the host data
df_host_stacked_summary <- df_host_stacked %>%
  group_by(DescriptionLong, Data, rep) %>%
  mutate(Count = (Count/sum(Count))*100) %>%
  group_by(DescriptionLong, HostPhylum, Data) %>%
  summarise(Mean_RelAbund = mean(RelAbund, na.rm = TRUE),
            SD_RelAbund = sd(RelAbund, na.rm = TRUE),
            Mean_Count = mean(Count),
            SD_Count = sd(Count))

# Load the SingleM data and filter for metagenomes
df_Microbiome <- df_Microbiome <- tar_read(SingleMData) %>%
  mutate(rep = substr(ShortSamples, nchar(ShortSamples), nchar(ShortSamples))) %>%
  mutate(DescriptionLong = DescriptionModifiers[Description]) %>%
  filter(DescriptionLong == "Meta\ngenome") %>%
  rename(HostPhylum = Phylum) %>%
  select(DescriptionLong, HostPhylum, RelAbund, rep) %>%
  mutate(Data = "Bacteria")

# Create a copy of df_Microbiome and merge with df_host_stacked
df_Microbiome2 <- df_Microbiome %>%
  merge(df_host_stacked, by = c("DescriptionLong","HostPhylum", "Data", "rep"), all = TRUE) %>%
  mutate(
    RelAbund = coalesce(RelAbund.x, RelAbund.y)
  ) %>%
  select(-RelAbund.x, -RelAbund.y) %>% # Drop the original split columns
  mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels),
         Data = factor(Data, levels = c("Bacteria", "Viruses (by predicted host)")))

# Calculate the Bacillota to Bacteroidota ratio
df_Bacillota_ratio <- df_Microbiome2 %>%
  filter(HostPhylum %in% c("Bacillota", "Bacteroidota") &
           Data == "Bacteria") %>%
  pivot_wider(names_from = HostPhylum, values_from = RelAbund) %>%
  mutate(Ratio = Bacillota / Bacteroidota)

# Calculate mean and SD of the Bacillota to Bacteroidota ratio
mean(df_Bacillota_ratio$Ratio)  
sd(df_Bacillota_ratio$Ratio)

# Create a plot of the microbiome and host data
plt_microbiome <- df_Microbiome2 %>%
  ggplot(aes(x = rep, y = RelAbund, fill = HostPhylum)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "white")) +
  labs(x = "Sample", 
       y = "vOTU\nRelative Abundance (%)", 
       fill = "Predicted\nHost Phylum") +
  scale_fill_manual(values = Phyla_palette) +
  facet_nested(~ Data + DescriptionLong, drop = TRUE)

# Import read depth data
df_read_depth_100 <- read_and_merge_tsvs("data/Subsampling/", "G.tsv") %>%
  select(Contig, contains("mapped Mean")) %>%
  pivot_longer(cols = -Contig, names_to = "Sample", values_to = "Coverage") %>%
  mutate(Depth = as.numeric(str_extract(Sample, "(?<=depth_).*?(?=G)")),
         Sample = str_extract(Sample, ".*?(?=_subsampled)")) %>%
  left_join(metadata, by = "Sample") %>%
  filter(!is.na(Coverage)) %>%
  group_by(Description, Depth, Sample) %>%
  summarise(TotalCount = n(),
            CountOver100 = sum(Coverage > 100, na.rm = TRUE),
            Percentage = (CountOver100 / TotalCount) * 100,
            .groups = "drop") %>%
  complete(Sample, Description, fill = list(TotalCount = 0, CountOver100 = 0, Percentage = 0)) %>%
  mutate(Description = replace_na(Description, "MtG"),
         DescriptionLong = DescriptionModifiers[Description])

# Create Fig3_bottom as a plot of the read depth data
Fig3_bottom <- df_read_depth_100 %>%
  ggplot(aes(x = Depth, y = Percentage, colour = DescriptionLong, fill = DescriptionLong)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_smooth() +
  geom_vline(xintercept = c(3.2, 10), linetype = "dashed", colour = "gray50") +
  labs(x = "Sequencing Depth (Gbp)", 
       y = "vOTUs With\n>= 100x Coverage (%)",
       fill = "Sample Type",
       colour = "Sample Type") +
  scale_color_manual(values = Description_palette) +
  scale_fill_manual(values = Description_palette) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2), colour = guide_legend(nrow = 2))

# Assemble Fig 3
Fig3 <- ggarrange(Fig3_top, plt_microbiome, Fig3_bottom, ncol = 1, nrow = 3,
                      labels = c("", "c","d"), legend = "right",heights = c(1.5,1.5,1))

# Save Fig3 as a PDF
ggsave("figures/Fig3_new.pdf", plot = Fig3,
       width = 170, height = 190, units = "mm", device = "pdf")

# Save Fig3 as a SVG
ggsave("figures/Fig3_new.svg", plot = Fig3,
       width = 170, height = 190, units = "mm", device = "svg")

# 4 - Supplementary Figure Generation ----
## Figure S1 - Experimental design ----
### Note - this was made in Biorender ###
### This section is here as a placeholder ###

## Figure S2 - Contig data ----
# Create plot of contigs
plt_contigs <- plot_means_with_cld(df_quast %>% filter(Metric == "# contigs"), "Value", "DescriptionLong", Description_palette) +
  labs(y = "Contigs", title = "Number of contigs") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  # Remove major and minor grid lines
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom")

# Create plot of filtered viral contigs
plt_viral_contigs <- plot_means_with_cld(df_quast %>% filter(Metric == "FilteredViralContigs"), "Value", "DescriptionLong", Description_palette) +
  labs(y = "Viral contigs", title = "Number of viral contigs") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  # Remove major and minor grid lines
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom")

# Create plot of largest contig
plt_Largest_contig <- plot_means_with_cld(df_quast %>% filter(Metric == "Largest contig"), "Value", "DescriptionLong", Description_palette) +
  labs(y = "Largest contig (bp)", title = "Largest contig") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  # Remove major and minor grid lines
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom")

# Create plot of auN
plt_auN <- plot_means_with_cld(df_quast %>% filter(Metric == "auN"), "Value", "DescriptionLong", Description_palette) +
  labs(y = "auN", title = "auN") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  # Remove major and minor grid lines
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom")

# Create plot of N50
plt_N50 <- plot_means_with_cld(df_quast %>% filter(Metric == "N50"), "Value", "DescriptionLong", Description_palette) +
  labs(y = "N50", title = "N50") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  # Remove major and minor grid lines
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom")

# Create plot of L50
plt_L50 <- plot_means_with_cld(df_quast %>% filter(Metric == "L50"), "Value", "DescriptionLong", Description_palette) +
  labs(y = "L50", title = "L50") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  # Remove major and minor grid lines
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom")

# Create plot of Total length
plt_Total_Length <- plot_means_with_cld(df_quast %>% filter(Metric == "Total length"), "Value", "DescriptionLong", Description_palette) +
  labs(y = "Total length (bp)", title = "Total length") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  # Remove major and minor grid lines
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom")

# Arrange plots into FigS2
FigS2 <- ggarrange(plt_contigs, plt_Largest_contig, plt_auN, plt_N50, plt_L50, plt_Total_Length, ncol = 2, nrow = 3,
                   labels = c("a", "b", "c", "d", "e", "f"), common.legend = TRUE, legend = "bottom")

# Save plot as pdf
ggsave("figures/FigS2.pdf", plot = FigS2,
       width = 160, height = 160, units = "mm", device = "pdf")

# Save svg version
ggsave("figures/FigS2.svg", plot = FigS2,
       width = 170, height = 170, units = "mm", device = "svg")

## Figure S3 - Ratio of mapped to assembled ----

# Set plt_mapped_only as the venn diagram for mapped only vOTUs
plt_mapped_only <- ls_Venns$plt_MappedOnly

# Merge the dfs and calculate the mapped to assembled ratio
df_all_viral<- df_mapped_viruses %>%
  merge(df_assembled_viruses, by = c("Sample", "Description")) %>%
  rename(AssembledViruses = FilteredViralContigs, MappedViruses = Count) %>%
  mutate(Ratio = MappedViruses / AssembledViruses)

# Plot the ratio of mapped to assembled vOTUs
plt_ratio <- plot_means_with_cld(df_all_viral, "Ratio", "Description", Description_palette) +
  labs(y = "Mapped / Assembled Viruses", title = "Mapped vOTUs/ Assembled viral contigs") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_bw() +
  # Remove major and minor grid lines
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Assemble Figure S3
FigS3 <- ggarrange(plt_ratio, plt_mapped_only,
                   ncol = 1, labels = c("a", "b"))

# Save plot as pdf
ggsave("figures/FigS3.pdf", plot = FigS3,
       width = 160, height = 160, units = "mm", device = "pdf")

# Save svg version
ggsave("figures/FigS3.svg", plot = FigS3,
       width = 160, height = 160, units = "mm", device = "svg")

## Figure S4 - Read-based data ----

# Plot raw reads data
plt_raw_reads <- plot_means_with_cld(df_reads %>% select(DescriptionLong, total_sequences) %>% unique(),"total_sequences", "DescriptionLong", Description_palette) +
  labs(y = "Total reads (millions)", title = "Total reads", fill = "Sample Type") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  # Remove major and minor grid lines
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

# Plot human reads data
plt_human_reads <- plot_means_with_cld(df_reads %>% filter(Type == "Human"),"ReadPercentage", "DescriptionLong", Description_palette) +
  labs(y = "Percentage (%)", title = "Human reads", fill = "Sample Type") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  # Remove major and minor grid lines
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Plot rRNA reads data
plt_rRNA_reads <- plot_means_with_cld(df_reads %>% filter(Type == "rRNA"),"ReadPercentage", "DescriptionLong", Description_palette) +
  labs(y = "Percentage (%)", title = "rRNA gene reads", fill = "Sample Type") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  # Remove major and minor grid lines
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Unique K-mer counts
# Import and subset unique k-mer data
df_unique_kmers <- tar_read(processed_kmers) %>%
  filter(abundance == 1) %>%
  mutate(DescriptionLong = DescriptionModifiers[Description]) %>%
  mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels))

# Plot unique k-mer data
plt_unique_kmers <- plot_means_with_cld(df_unique_kmers, "count", "DescriptionLong", Description_palette) +
  labs(y = "Unique kmers", title = "Unique kmers", fill = "Sample Type") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  # Remove major and minor grid lines
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom")

# Assemble Figure S4
FigS4 <- ggarrange(plt_raw_reads, plt_rRNA_reads,
                   plt_human_reads, plt_unique_kmers,
                   ncol = 2, nrow = 2, labels = c("a", "b", "c","d"),
                   common.legend = TRUE, legend = "bottom")

# Save plot as pdf
ggsave("figures/FigS4.pdf", plot = FigS4,
       width = 160, height = 160, units = "mm", device = "pdf")

# Save svg version
ggsave("figures/FigS4.svg", plot = FigS4,
       width = 160, height = 160, units = "mm", device = "svg")

## Figure S5 - Mean Bray-Curtis dissimilarities ----

# Convert distance matrix to long format
df_bray_long <- as.matrix(df_bray) %>%
  melt(varnames = c("Sample1", "Sample2"), value.name = "BrayCurtis") %>%
  filter(Sample1 != Sample2)  # Remove diagonal

# Add sample type metadata
df_bray_long <- df_bray_long %>%
  left_join(metadata %>% select(Sample1 = Sample, SampleType1 = Description), by = "Sample1") %>%
  left_join(metadata %>% select(Sample2 = Sample, SampleType2 = Description), by = "Sample2") %>%
  mutate(SampleType1 = DescriptionModifiers[SampleType1],
         SampleType2 = DescriptionModifiers[SampleType2])

# Remove duplicate pairs (e.g., A-B and B-A)
df_bray_long <- df_bray_long %>%
  rowwise() %>%
  mutate(pair_id = paste(sort(c(Sample1, Sample2)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair_id, .keep_all = TRUE)

# Add a Comparison column for within-group
within_group_long <- df_bray_long %>%
  filter(SampleType1 == SampleType2) %>%
  mutate(Comparison = SampleType1,ComparisonType = "Within")

# Add a Comparison column for between-group (sorted to avoid duplicates)
between_group_long <- df_bray_long %>%
  filter(SampleType1 != SampleType2) %>%
  rowwise() %>%
  mutate(Comparison = paste(sort(c(SampleType1, SampleType2)), collapse = " vs "),
         ComparisonType = "Between") %>%
  ungroup()

# Combine into a single tidy table
df_bray_comparisons <- bind_rows(within_group_long, between_group_long)

# Summarise comparisons
df_bray_summary <- df_bray_comparisons %>%
  group_by(Comparison, ComparisonType) %>%
  summarise(Mean_BrayCurtis = mean(BrayCurtis, na.rm = TRUE),
            SD_BrayCurtis = sd(BrayCurtis, na.rm = TRUE),
            Count = n(),
            .groups = "drop") %>%
  mutate(Comparison = factor(Comparison, levels = unique(Comparison)),
         ComparisonType = factor(ComparisonType, levels = c("Within", "Between"))) %>%
  mutate(Comparison = str_replace_all(as.character(Comparison), "\n", " "),
         Comparison = str_replace_all(Comparison, " vs ", " vs\n"),
         Comparison = str_replace_all(Comparison, "Meta genome", "Metagenome")) %>%
  mutate(Comparison = factor(Comparison, levels = unique(Comparison)))

# Plot Comparisons as a boxplot
plt_bray_comparisons <- ggplot(df_bray_summary, aes(x = Comparison, y = Mean_BrayCurtis, fill = ComparisonType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Mean_BrayCurtis - SD_BrayCurtis, ymax = Mean_BrayCurtis + SD_BrayCurtis),
                width = 0.2, position = position_dodge(width = 0.9)) +
  labs(x = "", y = "Mean Bray-Curtis\nDissimilarity", fill = "Comparison Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = c("Within" = "#56B4E9", "Between" = "#D55E00"))

# Save the Bray-Curtis comparisons plot as FigS5
ggsave("figures/FigS5.pdf", plot = plt_bray_comparisons,
       width = 160, height = 100, units = "mm", device = "pdf")

# Save the Bray-Curtis comparisons plot as SVG
ggsave("figures/FigS5.svg", plot = plt_bray_comparisons,
       width = 160, height = 100, units = "mm", device = "svg")
## Figure S6 - Annotation ----

# Load Genomad data output
df_Genomad <- tar_read(GenomadData)

# Import Clustering data
df_clustering <- tar_read(ClusteringData)

# Get a list of all viral sequences (not just vOTU representative sequences)
all_filtered_contigs <- df_clustering %>%
  # Concatenate all Cluster rows with a comma separator
  summarise(all_clusters = paste(Cluster, collapse = ",")) %>%
  # Split the concatenated string into a vector
  pull(all_clusters) %>%
  strsplit(split = ",") %>%
  unlist()

# Filter Genomad data for only contigs in all_filtered_contigs
df_Genomad_filtered <- df_Genomad %>%
  filter(seq_name %in% all_filtered_contigs) %>%
  rename(Contig = seq_name, Sample = sample) 

# Create a tsv of AMR data
df_AMR <- read.csv("data/ViralContigs/merged_amr_filtered.tsv", sep = "\t") %>%
  mutate(Contig = sub("_[^_]+$", "", gene)) %>%
  filter(Contig %in% all_filtered_contigs)

# Filter the AMR df for just vOTU representative sequences
df_AMR_vOTU <- df_AMR %>%
  filter(Contig %in% df_mapping$Contig)

# Extract vOTU names from the AMR data
AMR_vOTUs <- df_AMR_vOTU %>% pull(Contig) %>% unique()

# Identify vOTU clusters
AMR_clusters <- df_clustering %>%
  filter(vOTU %in% AMR_vOTUs)

# Join with metadata
df_Genomad_meta <- df_Genomad %>%
  rename(Contig = seq_name, Sample = sample) %>%
  left_join(metadata, by = "Sample") %>%
  mutate(DescriptionLong = DescriptionModifiers[Description]) %>%
  mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels)) %>%
  filter (is.na(SampleType) == FALSE)

# Import DefenseFinder data
df_Defense <- tar_read(DefenseFinderData)

# Merge with df_Genomad_filtered, keeping all rows in df_Defense
df_Defense_merged <- df_Defense %>%
  left_join(df_Genomad_filtered, by = c("vOTU" = "Contig")) %>%
  # Merge with metadata by Sample
  left_join(metadata, by = "Sample") %>%
  select(Sample, Description, type, subtype, activity) %>%
  mutate(DescriptionLong = DescriptionModifiers[Description]) %>%
  mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels))

# Create a df that contains the counts and percentages of sequences carrying defense genes
df_defense_activity <- df_Defense_merged %>%
  group_by(Sample, DescriptionLong, activity) %>%
  summarise(
    Count = n(),
    Percentage = (Count / nrow(df_Defense_merged)) * 100
  )

# Create a plot of defense systems
plt_defense <- df_defense_activity %>% filter(activity == "Defense") %>%
  plot_means_with_cld("Percentage", "DescriptionLong", Description_palette) +
  labs(y = "Viral contigs carrying\n bacteriophage defense systems (%)",
       fill = "Sample Type") +
  # Remove major and minor grid lines
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))  # Alternate labels into two rows

# Create a plot of antidefense systems
plt_antidefense <- df_defense_activity %>% filter(activity == "Antidefense") %>%
  plot_means_with_cld("Percentage", "DescriptionLong", Description_palette) +
  labs(y = "Viral contigs carrying\nantidefense systems (%)",
       fill = "Sample Type") +
  # Remove major and minor grid lines
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))  # Alternate labels into two rows

# Read in Pharokka annotation data
df_Pharokka <- tar_read(PharokkaData)

# Filter Pharokka data for moron, auxiliary metabolic gene and host takeover genes
df_morons <- df_Pharokka %>% filter(category == "moron, auxiliary metabolic gene and host takeover") %>%
  rename(vOTU = Contig)

# Aggregate Pharokka data into gene categories
df_Pharokka_category <- df_Pharokka %>%
  rename(Sample = sample) %>%
  select(category, Sample) %>%
  group_by(Sample, category) %>%
  summarise(Count = n()) %>%
  merge(metadata, by = "Sample") %>%
  mutate(DescriptionLong = DescriptionModifiers[Description]) %>%
  mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels))

# Create a plot of moron, auxiliary metabolic gene and host takeover genes
df_Pharokka_moron <- df_Pharokka_category %>%
  group_by(Sample, DescriptionLong) %>%
  summarise(TotalGenes = sum(Count),
            TotalMoron = sum(Count[category == "moron, auxiliary metabolic gene and host takeover"]),
            Percentage = (TotalMoron / TotalGenes) * 100)

# Create a summary df of moron, auxiliary metabolic gene and host takeover genes
df_Pharokka_moron_summary <- df_Pharokka_moron %>% ungroup() %>%
  group_by(DescriptionLong) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE)
  )

# Create a plot of moron, auxiliary metabolic gene and host takeover genes
plt_morons <- df_Pharokka_moron %>% 
  plot_means_with_cld("Percentage", "DescriptionLong", Description_palette) +
  labs(y = "Viral contigs carrying a\nmoron, auxiliary metabolic gene\nand host takeover gene (%)",
       fill = "Sample Type") +
  scale_fill_manual(values = Description_palette) +
  # Remove major and minor grid lines
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))  # Alternate labels into two rows

# Arrange plots into FigS6
FigS6 <- ggarrange(plt_morons, plt_defense, plt_antidefense, nrow = 1, heights = c(1, 1),
                   labels = c("a", "b", "c"), common.legend = TRUE, legend = "bottom")

# Save plot as pdf
ggsave("figures/FigS6.pdf", plot = FigS6,
       width = 160, height = 110, units = "mm", device = "pdf")

# Save svg version
ggsave("figures/FigS6.svg", plot = FigS6,
       width = 160, height = 110, units = "mm", device = "svg")

## Figure S7 - Unified Human Gut Virome Catalog comparison ----

# Load in contig clustering data from clustering with the Gut Phage Database
df_UHGV <- read.csv("data/Biogeography/AllContigs_UHGV_all_clusters.tsv", header = FALSE, sep = "\t") %>%
  rename(Contig = V1, Cluster = V2) %>%
  # filter for rows that contain text "contig" in Cluster
  filter(grepl("contig", Cluster)) %>%
  # Create column that is true/ false on if "UHGV" in Cluster
  mutate(UHGV = grepl("UHGV", Cluster),
         NumSeqs = str_count(Cluster, ",") + 1,
         NumUHGVs = str_count(Cluster, "UHGV"),
         NumStudy = NumSeqs-NumUHGVs)

# Summarise the number of clusters with a UHGV vOTU in them
table(df_UHGV$UHGV)

# Plot the number of clusters with a UHGV vOTU in them and the number without
plt_UHGV_Clusters <- df_UHGV %>%
  group_by(UHGV) %>%
  summarise(NumClusters = n_distinct(Cluster)) %>%
  ggplot(aes(x = UHGV, y = NumClusters, fill = UHGV)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = NumClusters), vjust = 1.5, size = 3) +
  labs(x = "Cluster Contains a UHGV vOTU", y = "Clusters") +
  scale_fill_manual(values = c("TRUE" = "#F5C875",
                               "FALSE" = "#A6D5F5")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        # Remove grid lines
        panel.grid = element_blank())

# Create a list of contigs from this study that occur within a cluster that also
# contains a UHGV from the Gut Phage Database
UHGV_filtered_contigs <- df_UHGV %>%
  filter(UHGV == TRUE) %>%
  summarise(all_clusters = paste(Cluster, collapse = ",")) %>%
  pull(all_clusters) %>%
  strsplit(split = ",") %>%
  unlist() %>%
  .[. %in% all_filtered_contigs]

# Create a list of contigs from this study that do not occur within a cluster that
# also contains a UHGV vOTU
NonUHGV_filtered_contigs <- df_UHGV %>%
  filter(UHGV == FALSE) %>%
  summarise(all_clusters = paste(Cluster, collapse = ",")) %>%
  pull(all_clusters) %>%
  strsplit(split = ",") %>%
  unlist() %>%
  .[. %in% all_filtered_contigs]

# Create a list of UHGV contigs from this study that occur within a cluster that
# also contains a UHGV vOTU
UHGV_contigs <- df_UHGV %>%
  filter(UHGV == TRUE) %>%
  summarise(all_clusters = paste(Cluster, collapse = ",")) %>%
  pull(all_clusters) %>%
  strsplit(split = ",") %>%
  unlist() %>%
  .[grepl("UHGV", .)]

# Filter Genomad data for viral contigs from this study also clustering with
# the UHGV vOTUs
df_UHGVs <- df_Genomad_filtered %>%
  mutate(UHGV = Contig %in% UHGV_filtered_contigs) %>%
  left_join(metadata, by = "Sample") %>%
  mutate(DescriptionLong = DescriptionModifiers[Description]) %>%
  mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels)) %>%
  select(Sample, DescriptionLong, ShortSamples, UHGV) %>%
  group_by(Sample, DescriptionLong) %>%
  summarise(UHGV = sum(UHGV, na.rm = TRUE),
            Total = n(),
            Percentage = (UHGV/ Total)*100)

# Summarise the percentage GPD data
df_UHGVs_summary <- df_UHGVs %>% ungroup() %>%
  summarise(Mean = mean(Percentage, na.rm = TRUE),
            SD = sd(Percentage, na.rm = TRUE))

# Plot Percentage GPD data
plt_UHGVs <- df_UHGVs %>% plot_means_with_cld("Percentage", "DescriptionLong", Description_palette) +
  labs(x = "Sample Type", y = "Percentage of\nviral contigs (%)", fill = "Sample Type") +
  # Expand y-axis limits
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  # Remove grid lines
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # Move legend to the bottom
  theme(legend.position = "none")

# Assemble the top of Fig S7
FigS7_Top <- ggarrange(plt_UHGV_Clusters, plt_UHGVs, labels = c("a","b"), ncol = 2, nrow = 1)

# Read all files in folder, filter out rows with values of 0 in the last column and merge them together, with a column for filename
folder <- "data/Biogeography/Mapping_UHGV/"

# Get list of files in the folder
files <- list.files(folder, full.names = TRUE)

# Split files based on "_mean" presence
mean_files <- files[grepl("_mean", basename(files), ignore.case = TRUE)]

# Filter out files that do not contain "_mean"
non_mean_files <- files[!grepl("_mean", basename(files), ignore.case = TRUE)]

# Function to read, filter, and process each file
process_files <- function(file_list) {
  file_list %>%
    set_names() %>%
    map_dfr(~ {
      df <- read_tsv(.x, col_types = cols(.default = "c"))
      last_col <- ncol(df)
      df <- df %>% filter(df[[last_col]] != "0")
      df$SourceFile <- basename(.x)
      return(df)
    }) %>%
    select(-Sample) %>%
    rename(Sample = SourceFile) %>%
    mutate(Sample = gsub(".tsv", "", Sample))}

# Process mean files
df_UHGV_mean     <- process_files(mean_files) %>%
  mutate(Sample = gsub("_mean", "", Sample))

# Process non-mean files
df_UHGV_nonmean  <- process_files(non_mean_files)

# Merge the mean and non-mean dataframes and filter for covered fraction >= 0.75
df_UHGV_mapping <- inner_join(
  df_UHGV_mean,
  df_UHGV_nonmean,
  by = c("Contig", "Sample")) %>%
  filter(`Covered Fraction` >= 0.75)

# Create a new column to indicate if the contig is a UHGV contig
MappedOnlyUHGVs <- df_UHGV_mapping  %>%
  pull(Contig) %>%
  unique() %>%
  # Filter for those not in UHGV_contigs
  .[! . %in% UHGV_contigs]

# Create a new column to indicate if the contig is a MappedOnlyUHGV
df_UHGV_mapping <- df_UHGV_mapping %>%
  mutate(MappedOnlyUHGV = Contig %in% MappedOnlyUHGVs)

# Filter the df_UHGV_mapping to only include MappedOnlyUHGVs
df_UHGV_mapping_only <- df_UHGV_mapping %>%
  filter(MappedOnlyUHGV == TRUE)

# Read in the metadata for UHGVs
df_UHGV_metadata <- read.csv("data/Biogeography/votus_metadata_extended.tsv", sep = "\t") %>%
  filter(uhgv_genome %in% MappedOnlyUHGVs |
           uhgv_genome %in% UHGV_contigs) %>%
  mutate(MappedOnlyUHGV = uhgv_genome %in% MappedOnlyUHGVs)

# Create a summary table for df_UHGV_mapping and df_UHGV_mapping_only
df_UHGV_mapping_summary <- df_UHGV_mapping %>%
  filter(`Covered Fraction` >= 0.75) %>%
  group_by(Sample) %>%
  summarise(Count = n()) %>%
  merge(metadata, by = "Sample") %>%
  select(Description, Sample, Count) %>%
  mutate(DescriptionLong = DescriptionModifiers[Description]) %>%
  mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels))

# Create a summary table for df_UHGV_mapping_only
df_UHGV_mapping_only_summary <- df_UHGV_mapping_only %>%
  group_by(Sample) %>%
  summarise(Count = n()) %>%
  merge(metadata, by = "Sample") %>%
  select(Description, Sample, Count) %>%
  mutate(DescriptionLong = DescriptionModifiers[Description]) %>%
  mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels))

# Create a plot of the UHGV mapping summary
plt_UHGV_mapping <- plot_means_with_cld(df_UHGV_mapping_summary, "Count", "DescriptionLong", Description_palette) +
  labs(x = "Sample Type", y = "Count", fill = "Sample Type") +
  scale_y_continuous(
    limits = c(0, 265),
    expand = expansion(mult = c(0, 0.1))
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = Description_palette)

plt_UHGV_mapping_only <- plot_means_with_cld(df_UHGV_mapping_only_summary, "Count", "DescriptionLong", Description_palette) +
  labs(x = "Sample Type", y = "Count", fill = "Sample Type") +
  # Expand y-axis limits
  scale_y_continuous(
    limits = c(0, 265),
    expand = expansion(mult = c(0, 0.1))
  ) +
  # Remove grid lines
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # Move legend to the bottom
  theme(legend.position = "bottom") +
  scale_fill_manual(values = Description_palette)

# Assemble the middle of Fig S7
FigS7_Middle <- ggarrange(plt_UHGV_mapping, plt_UHGV_mapping_only, 
                          ncol = 2, nrow = 1, labels = c("c", "d"),
                          common.legend = TRUE, legend = "bottom")


# Calculate the maximum mean coverage for each contig and whether it is a MappedOnlyUHGV
df_UHGV_coverage <- df_UHGV_mapping %>%
  group_by(Contig, MappedOnlyUHGV) %>%
  summarise(MaxMean = max(`Trimmed Mean`)) %>%
  mutate(MaxMean = as.numeric(MaxMean)) %>%
  ungroup()

# Run the Wilcoxon test
wilcox_result <- wilcox.test(MaxMean ~ MappedOnlyUHGV, data = df_UHGV_coverage)

# Extract and format the p-value
formatted_p <- if (wilcox_result$p.value < 1e-5) {"p < 0.00001"} else {
  paste0("p = ", signif(wilcox_result$p.value, digits = 2))}

# Create the plot
plt_UHGV_coverage <- ggplot(df_UHGV_coverage, aes(x = MappedOnlyUHGV, y = MaxMean, fill = MappedOnlyUHGV)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
  geom_jitter(aes(color = MappedOnlyUHGV), width = 0.2, alpha = 0.5, size = 0.5) +
  scale_fill_manual(values = c("TRUE" = "#F5C875",
                               "FALSE" = "#A6D5F5")) +
  scale_color_manual(values = c("TRUE" = "#F5C875",
                                "FALSE" = "#A6D5F5")) +
  scale_y_log10(labels = scales::label_log()) +
  annotate("text", x = 2, y = max(df_UHGV_coverage$MaxMean, na.rm = TRUE)+1.5e5, 
           label = paste("Wilcox test,", formatted_p), size = 4) +
  labs(x = "Only detected by mapping", y = "Coverage", fill = "Only detected by mapping",
       colour = "Only detected by mapping") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        panel.grid = element_blank())

# Perform a Wilcoxon test for genome length
wilcox_result_length <- wilcox.test(genome_length ~ MappedOnlyUHGV, data = df_UHGV_metadata)

# Extract and format the p-value for genome length
formatted_p_length <- if (wilcox_result_length$p.value < 1e-5) {
  "p < 0.00001"
} else {
  paste0("p = ", signif(wilcox_result_length$p.value, digits = 2))
}

# Build the plot
plt_UHGV_length <- df_UHGV_metadata %>%
  mutate(
    MappedOnlyUHGV = as.factor(MappedOnlyUHGV),
    phylum_name = fct_infreq(phylum_name)
  ) %>%
  ggplot(aes(x = MappedOnlyUHGV, y = genome_length, fill = MappedOnlyUHGV)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.6, color = NA) +
  geom_jitter(aes(color = MappedOnlyUHGV), size = 0.5, alpha = 0.5, width = 0.2, height = 0) +
  geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
  scale_y_log10(labels = scales::label_log()) +
  annotate("text", 
           x = 1.85, 
           y = max(df_UHGV_metadata$genome_length, na.rm = TRUE) * 1.2, 
           label = paste("Wilcox test,", formatted_p_length), 
           size = 4) +
  scale_fill_manual(values = c("TRUE" = "#F5C875",
                               "FALSE" = "#A6D5F5")) +
  scale_color_manual(values = c("TRUE" = "#F5C875",
                                "FALSE" = "#A6D5F5")) +
  labs(y = "Genome Length",
       x = "Only detected by mapping",
       fill = "Mapped Only UHGV",
       colour = "Mapped Only UHGV") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        panel.grid = element_blank())

# Arrange the plots into a single figure
FigS7_Bottom <- ggarrange(plt_UHGV_coverage, plt_UHGV_length,
                          ncol = 2, nrow = 1, labels = c("e", "f"),
                          common.legend = TRUE, legend = "bottom")

# Assemble the final figure
FigS7 <- ggarrange(FigS7_Top, FigS7_Middle, FigS7_Bottom, nrow = 3)

# Save plot as pdf
ggsave("figures/FigS7.pdf", plot = FigS7,
       width = 170, height = 170, units = "mm", device = "pdf")

# Save svg version
ggsave("figures/FigS7.svg", plot = FigS7,
       width = 170, height = 170, units = "mm", device = "svg")

# Produce a summary table of df_UHGVs by DescriptionLong
df_UHGVs_summary_table <- df_UHGVs %>%
  group_by(DescriptionLong) %>%
  summarise(
    Mean_Percentage = mean(Percentage, na.rm = TRUE),
    SD_Percentage = sd(Percentage, na.rm = TRUE),
    Count = n()
  ) %>%
  arrange(desc(Mean_Percentage))

# Print the summary statistics for df_UHGVs
mean(df_UHGVs$Percentage)
sd(df_UHGVs$Percentage)

# Produce a summary table of df_UHGV_mapping_summary by DescriptionLong
df_UHGV_mapping_summary_table <- df_UHGV_mapping_summary %>%
  group_by(DescriptionLong) %>%
  summarise(Mean_Count = mean(Count, na.rm = TRUE),
            SD_Count = sd(Count, na.rm = TRUE),
            Count = n()) %>%
  arrange(desc(Mean_Count))