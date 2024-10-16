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
  theme(legend.position = "none")

plt_read_percentages <- tar_make(plt_read_percentages)

plt_top <- ggarrange(plt_raw_reads, plt_kmers, ncol = 2, nrow = 1,
                     labels = c("a", "b"))

plt_fig_1 <- ggarrange(plt_top, plt_read_percentages, ncol = 1, nrow = 2,
                       labels = c("", "c"))



# Figure 2 - Viral abundance comparisons
