# Clear environment
rm(list=ls())
source("scripts/Functions_FVO.R")
library(tidyverse)

metadata <- tar_read(metadata)

df_mapping <- tar_read(vOTU_tpm)

df_Quast <- read.csv("data/ViralContigs/report.tsv", sep = "\t")

df_Clustering <- read.csv("data/ViralContigs/combined_Genomad_FilteredContigs_clusters.tsv", sep = "\t", header = FALSE) %>%
  rename("vOTU" = "V1", "Cluster" = "V2")

viral_contigs <- df_Clustering$Contigs %>%
  # Merge all contigs into a single vector, split individual entries by ","
  strsplit(",") %>%
  # Unlist the vector
  unlist()

df_Genomad <- read_and_merge_tsvs("data/ViralContigs/", "_renamed_contigs_virus_summary.tsv") %>%
  rename(Contig = seq_name,
         Sample = sample) %>%
  filter(Contig %in% viral_contigs) %>%
  merge(metadata, by = "Sample", all.x = TRUE) %>%
  mutate(Description = factor(Description, levels = c("DNase Treated Virome", "Untreated Virome", "MDA Amplified Virome", "Metagenome"))) %>%
  # Split taxonomy into individual columns
  separate(col = "taxonomy", into = c("Viruses", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family"), sep = ";", fill = "right")

df_Genomad %>% ggplot(aes(x = ShortSamples, y = length, fill = Description)) +
  geom_violin() +
  labs(title = "Contig Length",
       x = "Sample",
       y = "Number of contigs") +
  theme_bw() +
  theme(legend.position = "none",
        # remove Grid lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10() +
  scale_fill_manual(values = Description_palette)

df_alpha_diversity <- df_mapping %>%
  # Pivot columns longer except for Contig and convert names to Sample and values to TPM
  pivot_longer(cols = -c(Contig), names_to = "Sample", values_to = "TPM") %>%
  # Summarise number of vOTUs with TPM > 0
  filter(TPM > 0) %>%
  group_by(Sample) %>%
  summarise(Richness = n_distinct(Contig)) %>% merge(metadata, by = "Sample", all.x = TRUE)

# Plot alpha diversity
plt_alpha <- df_alpha_diversity %>% ggplot(aes(x = ShortSamples, y = Richness, fill = Sample)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample",
       y = "Number of assembled\nviral contigs") +
  theme_bw() +
  theme(legend.position = "none",
        # remove Grid lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = Replicates_palette)

# Create Presence/ Absence table
df_PA <- df_mapping %>%
  # Pivot columns longer except for Contig and convert names to Sample and values to TPM
  pivot_longer(cols = -c(Contig), names_to = "Sample", values_to = "TPM") %>%
  # Create a binary presence/absence table
  mutate(PA = ifelse(TPM > 0, 1, 0)) %>% 
  select(-TPM) %>%
  merge(metadata, by = "Sample", all.x = TRUE) %>%
  rename(vOTU = Contig) %>%
  merge(df_Clustering, by = "vOTU", all.x = TRUE) %>%
  # Apply grepl row-wise using mapply
  mutate(AssembledPresence = mapply(function(sample, cluster) grepl(sample, cluster), Sample, Cluster)) %>%
  # Convert logical to 1/0
  mutate(AssembledPresence = ifelse(AssembledPresence, 1, 0))

library(ggVennDiagram)

# Filter by MappedPresence == 1, then split vOTU based on Description
df_MappedVenn <- df_PA %>% 
  filter(PA == 1) %>%
  select(vOTU, Description) %>%
  split(.$Description)

df_AssembledVenn <- df_PA %>% 
  filter(AssembledPresence == 1) %>%
  select(vOTU, Description) %>%
  split(.$Description)

DescriptionLevelOrder <- c("DNase Treated Virome", "Untreated Virome", "MDA Amplified Virome", "Metagenome")

Mapped_sets <- lapply(df_MappedVenn, function(x) unique(x$vOTU)) %>%
  # Change order of list elements to match DescriptionLevelOrder
  .[DescriptionLevelOrder] %>%
  setNames(c("+DNase", "-DNase", "MDA", "MG"))

Assembled_sets <- lapply(df_AssembledVenn, function(x) unique(x$vOTU)) %>%
  # Change order of list elements to match DescriptionLevelOrder
  .[DescriptionLevelOrder] %>%
  setNames(c("+DNase", "-DNase", "MDA", "MG"))

Venn_Mapped <- ggVennDiagram(Mapped_sets, label = "count", edge_size = 0) +
  scale_fill_distiller(palette = "Blues", direction = 1) + 
  theme(legend.position = "bottom") +
  # Set fill scale to 0-200 and going from white to blue
  scale_fill_distiller(limits = c(0, 200), palette = "Blues", direction = 1)

Venn_Assembled <- ggVennDiagram(Assembled_sets, label = "count", edge_size = 0) +
  scale_fill_distiller(palette = "Blues", direction = 1) + 
  theme(legend.position = "bottom") +
  # Set fill scale to 0-200 and going from white to blue
  scale_fill_distiller(limits = c(0, 200), palette = "Blues", direction = 1)

ggarrange(Venn_Mapped, Venn_Assembled, ncol = 2, labels = c("a", "b"),
          common.legend = TRUE, legend = "bottom")

library(vegan)
library(dplyr)

library(dplyr)

create_df_beta <- function(df_mapping) {
  
  # Step 1: Remove contigs (singletons) present in only one sample
  df_filtered <- df_mapping %>%
    rowwise() %>%
    filter(sum(c_across(-Contig) > 0) > 1) %>%  # Keep contigs present in more than one sample
    ungroup()
  
  # Step 2: Renormalize TPM values so that they sum to 1 million in each sample
  df_renormalized <- df_filtered %>%
    mutate(across(-Contig, ~ . / sum(.) * 1e6))  # Renormalize TPM values to sum to 1 million per sample
  
  return(df_renormalized)
}

df_beta <- create_df_beta(df_mapping)

# Beta diversity
plot_pcoa <- function(df, metadata, Description_palette) {
  
  # Step 1: Calculate Bray-Curtis Dissimilarity
  df_beta_diversity <- df %>%
    column_to_rownames(var = "Contig") %>%   # Convert contigs to row names
    t() %>%                                  # Transpose the data so rows are samples, columns are contigs
    vegdist(method = "bray")                 # Calculate Bray-Curtis dissimilarity on samples
  
  # Step 2: Perform PCoA on Bray-Curtis dissimilarity matrix
  pcoa_result <- cmdscale(df_beta_diversity, eig = TRUE, k = 2)
  
  # Step 3: Extract PCoA axes
  pcoa_axes <- as.data.frame(pcoa_result$points)
  colnames(pcoa_axes) <- c("PCoA1", "PCoA2")
  
  # Add Sample names as a column for plotting
  pcoa_axes$Sample <- rownames(pcoa_axes)
  
  # Step 4: Merge with metadata to add 'Description' column
  pcoa_axes <- merge(pcoa_axes, metadata, by = "Sample") %>%
    mutate(Description = factor(Description, levels = c("DNase Treated Virome", "Untreated Virome", "MDA Amplified Virome", "Metagenome")))
  
  # Step 5: Calculate the percentage of variance explained by each axis
  variance_explained <- 100 * pcoa_result$eig[1:2] / sum(pcoa_result$eig)
  
  # Step 6: Plot the PCoA results, colored by 'Description'
  pcoa_plot <- ggplot(pcoa_axes, aes(x = PCoA1, y = PCoA2, color = Description, shape = Description)) +
    geom_point(size = 4, alpha = 0.5) +
    labs(x = paste0("PCoA1\n(", round(variance_explained[1], 1), "% variance)"),
         y = paste0("PCoA2\n(", round(variance_explained[2], 1), "% variance)"),) +
    scale_shape_manual(values = c(16, 17, 15, 18)) +  # Assign different shapes
    theme_bw() +
    # Remove grid lines
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_colour_manual(values = Description_palette) +
    # Fix co-oridnate axis
    coord_fixed()
  
  # Return the plot
  return(pcoa_plot)
}

# Plot PCoA
plt_pcoa <- plot_pcoa(df_beta, metadata, Description_palette)

plt_pcoa

monodnaviria <- df_Genomad %>%
  filter(Realm == "Monodnaviria")

df_beta_no_mondna <- df_mapping %>% filter(!Contig %in% monodnaviria$Contig) %>%
  create_df_beta()

plt_pcoa_no_mondna <- plot_pcoa(df_beta_no_mondna, metadata, Description_palette) +
  lims(x = c(-0.255, 0.69), y = c(-0.1, 0.15))

plt_pcoa_no_mondna
plt_pcoas <- ggarrange(plt_pcoa, plt_pcoa_no_mondna, ncol = 1, labels = c("a", "b"),
          common.legend = TRUE, legend = "right", heights = c(2.5,1))

df_tax <- df_Genomad %>%
  select(Contig, Viruses, Realm, Kingdom, Phylum, Class, Order, Family)

# Stacked barplots
df_stacked <- df_mapping %>%
  pivot_longer(cols = -c(Contig), names_to = "Sample", values_to = "TPM") %>%
  filter(TPM > 0) %>%
  merge(metadata, by = "Sample", all.x = TRUE) %>%
  mutate(Description = factor(Description, levels = c("DNase Treated Virome", "Untreated Virome", "MDA Amplified Virome", "Metagenome"))) %>%
  # Merge with df_Genomad
  merge(df_tax, by = "Contig", all.x = TRUE) %>%
  # Modify NA in Realm to "Unassigned"
  mutate(Realm = ifelse(is.na(Realm), "Unassigned", Realm)) %>%
  # Summarise TPM values in each sample to Realm level
  group_by(Realm, ShortSamples, Description) %>%
  summarise(TPM = sum(TPM)) %>%
  mutate(RelAbund = TPM / 10000)

plt_stacked <- df_stacked %>%
  ggplot(aes(x = ShortSamples, y = RelAbund, fill = Realm)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample",
       y = "Relative Abundance (x 10^4)") +
  theme_bw() +
  theme(legend.position = "bottom",
        # remove Grid lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = Description_palette)

plt_stacked
