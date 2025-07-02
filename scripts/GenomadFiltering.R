#### Genomad Data Processing ####

# 0 - Setup ----
## Clear Environment
rm(list = ls())

## Load Libraries
library(tidyverse)

## Functions ----

## Function to load and merge all genomad tsv files in <folder> with an additional
## column for the sample name, taken from the filename with the basename and file
## extension removed
read_and_merge_tsvs <- function(folder, file_extension){
  # Get list of files in folder
  files <- list.files(path = folder, pattern = file_extension)
  # Create empty list to store dataframes
  df_list <- list()
  # Loop through files
  for (i in 1:length(files)){
    # Read file
    df <- read_tsv(paste(folder, "/", files[i], sep = ""))
    # Add sample name column
    df$sample <- gsub(file_extension, "", files[i])
    # Add dataframe to list
    df_list[[i]] <- df
  }
  # Merge all dataframes in list
  df_merged <- bind_rows(df_list)
  # Return merged dataframe
  return(df_merged)
}

## Function to write filtered viral contig names to a file, or create an empty
## file if the df is empty
write_seq_files <- function(input_df, virus_type, sample_names = sample_names) {
  library(tidyverse)
  library(readr)
  
  if (!dir.exists("data")) {
    dir.create("data")
  }
  
  # Debug: Print the list of unique samples
  print(paste("All samples:", toString(sample_names)))
  
  for (sample_name in sample_names) {
    sample_filtered_df <- input_df %>% filter(sample == sample_name)
    
    output_filename <- file.path("data", paste0(sample_name, "_", virus_type, ".tsv"))
    empty_filename <- file.path("data", paste0(sample_name, "_NO_HITS_", virus_type, ".tsv"))
    
    if (!is.null(sample_filtered_df) && nrow(sample_filtered_df) > 0) {
      write_tsv(sample_filtered_df %>% select(seq_name), output_filename, col_names = FALSE)
    } else {
      # Debug: Print the empty_filename
      print(paste("Attempting to create empty file:", empty_filename))
      
      system(paste("touch", shQuote(empty_filename)))
    }
  }
}

# 1 - Data import and filtering ----
setwd("~/Projects/FecalViromeOptimisation/Virome/4-genomad/Individual/filtered")
## Import Genomad Data
df_gmad_raw <- read_and_merge_tsvs("data", ".contigs_virus_summary.tsv")

sample_names <- unique(df_gmad_raw$sample)

## Filter quality
df_gmad_Q <- df_gmad_raw[df_gmad_raw$fdr < 0.05,]

## Filter ssDNA and non-ssDNA viruses into separate lists
ssDNA_filter <- c("Monodnaviria","Ainoaviricetes","Tolecusatellitidae", "Anelloviridae", "Alphasatellitidae")
dsDnaMonovirales <- c("Papovaviricetes", "Trapavirae")
# Select rows where the column taxonomy contains at least one of the terms in ssDNA_filter
df_gmad_ssDNA <- df_gmad_Q[grepl(paste(ssDNA_filter, collapse = "|"), df_gmad_Q$taxonomy),]
 
# Remove entries that are from the dsDNA Monovirales
df_gmad_ssDNA <- df_gmad_ssDNA[!grepl(paste(dsDnaMonovirales, collapse = "|"), df_gmad_ssDNA$taxonomy),]

# Select Trapavirae
df_gmad_Trapavirae <- df_gmad_Q[grepl("Trapavirae", df_gmad_Q$taxonomy),]

# Select all viruses from df_gmad_Q not in df_gmad_ssDNA or df_gmad_Trapavirae
df_gmad_filtered <- df_gmad_Q[!df_gmad_Q$seq_name %in% df_gmad_ssDNA$seq_name & !df_gmad_Q$seq_name %in% df_gmad_Trapavirae$seq_name,]

# Also filter df_gmad_filtered to length >10000 kb
df_gmad_filtered <- df_gmad_filtered[df_gmad_filtered$length > 10000,]

# 2 - Writing list of seq_names to file ----
## Write ssDNA viral contig names
write_seq_files(df_gmad_ssDNA, "ssDNA", sample_names)

## Write Trapavirae viral contig names
write_seq_files(df_gmad_Trapavirae, "Trapavirae", sample_names)

## Write all other filtered viral contig names
write_seq_files(df_gmad_filtered, "filtered", sample_names)