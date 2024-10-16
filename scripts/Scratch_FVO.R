# Clear environment
rm(list=ls())
source("scripts/Functions_FVO.R")
library(tidyverse)

metadata <- tar_read(metadata)

df_mapping <- tar_read(vOTU_tpm)

df_Quast <- read.csv("data/ViralContigs/report.tsv", sep = "\t")

df_SingleM <- read_and_merge_tsvs("data/MicrobialProfiling", "_tax_profile.tsv") %>%
  # Split taxonomy
  separate(taxonomy, into = c("Root","Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") %>%
  # If Phylum = NA, then take value from Domain and fill as paste("Unclassified", Domain)
  mutate(Phylum = ifelse(is.na(Phylum), paste("Unclassified", Domain), Phylum)) %>%
  # Remove "d_" and "p_" from Phylum
  mutate(Phylum = gsub("d__", "", Phylum)) %>%
  mutate(Phylum = gsub("p__", "", Phylum)) %>%
  # Rename all Phyla containing text "Bacillota" to Bacillota
  mutate(Phylum = ifelse(grepl("Bacillota", Phylum), "Bacillota", Phylum)) %>%
  select(c(sample, Phylum, coverage)) %>%
  group_by(sample, Phylum) %>%
  summarise(coverage = sum(coverage)) %>%
  group_by(sample) %>%
  mutate(RelAbund = coverage/sum(coverage) *100) %>%
  rename(Sample = sample) %>%
  # Merge with metadata
  merge(metadata, by = "Sample")

plt_SingleM_coverage <- df_SingleM %>% filter(Description != "Soil Virome") %>%
  ggplot(aes(x = ShortSamples, y = coverage, fill = Phylum)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Sample", y = "Relative Abundance (%)", fill = "Phylum") +
  ggtitle("Microbial Profiling")

plt_SingleM_relative_abundance <- df_SingleM %>% filter(Description != "Soil Virome") %>%
  ggplot(aes(x = ShortSamples, y = RelAbund, fill = Phylum)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Sample", y = "Relative Abundance (%)", fill = "Phylum")

plt_SingleM_relative_abundance
