# Fecal Viromes Optimisation Analysis Parent Script

# 0 - Setup ----

## Clear environment
rm(list=ls())

## Load libraries
library(targets)
library(tidyverse)

## Load functions
source("scripts/Functions_FVO.R")
#tar_destroy()
# 1 - Run targets pipeline
## Visualise pipeline
tar_visnetwork()

## Run pipeline
tar_make()

# Figure 1 - raw reads ----
plt_kmers <- tar_read(plt_kmers) + 
  # Remove legend
  theme(legend.position = "none")

plt_raw_reads <- tar_read(plt_raw_reads) +
  # Remove legend
  theme(legend.position = "none") + 
  # Remove x-axis ticks, labels and title
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

plt_raw_reads

plt_read_percentages <- tar_read(plt_read_percentages)

plt_top <- ggarrange(plt_raw_reads, plt_kmers, ncol = 2, nrow = 1,
                     labels = c("a", "b"))

plt_fig_1 <- ggarrange(plt_top, plt_read_percentages, ncol = 1, nrow = 2,
                       labels = c("", "c"))

#plt_fig_1

# Figure 2 - Viral abundance comparisons


# Figure 2 - contigs ----
plt_length <- tar_read(plt_length)

