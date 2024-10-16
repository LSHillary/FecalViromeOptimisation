# 1 - Importing data ----
# Import Metadata
import_metadata <- function(filepath){
  df <- read.csv(filepath, sep = ",", header = TRUE)
  return(df)
}

# Import Raw MultiQC
import_multiqc <- function(filepath, metadata){
  df <- read.csv(filepath, sep = "\t", header = TRUE, stringsAsFactors = TRUE) %>%
    # Remove "FastQC_mqc.generalstats.fastqc." from column names
    rename_all(~gsub("FastQC_mqc.generalstats.fastqc.", "", .)) %>%
    # Remove R2
    filter(!grepl("R2", Sample)) %>%
    # Remove "_raw_R1" from sample name %>%
    mutate(Sample = gsub("_raw_R1", "", Sample)) %>%
    # Select relevant columns
    select(Sample, percent_gc, total_sequences) %>%
    merge(metadata, by = "Sample", all.x = TRUE)
  return(df)
}

# Import and merge multiple tsvs
read_and_merge_tsvs <- function(folder, file_extension){
  # Get list of files in folder
  files <- list.files(path = folder, pattern = file_extension)
  # Create empty list to store dataframes
  df_list <- list()
  # Loop through files
  for (i in 1:length(files)){
    # Read file
    df <- read_tsv(paste(folder, "/", files[i], sep = ""), show_col_types = FALSE)
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

# Import human read coverage data
import_read_counts <- function(filename, TEXT){
  df <- read.csv(filename, sep = "\t") %>%
    # Filter columns for those containing "Read.Count"
    select(contains("Read.Count")) %>%
    # Remove "_human_mapped.Read.Count" from column names
    rename_all(~gsub(paste0("_", TEXT, "_mapped.Read.Count"), "", .)) %>%
    # Summarise total counts for each column
    summarise_all(sum) %>%
    # Transpose and turn column names into a new column called Sample
    t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>%
    rename("Reads" = "V1") %>%
    # Divide HumanReads by 2 to convert to read pairs
    mutate(Reads = Reads / 2)
  return(df)
}

# Import TPM values
filename <- "data/ReadMapping/Stool_Genomad_coverage_tpm.tsv"
TEXT <- "stool_Genomad"
import_tpm <- function(filename, TEXT){
  df <- read.csv(filename, sep = "\t") %>%
    # Filter columns that contain "TPM" or "Contig"
    select(contains("TPM"), contains("Contig")) %>%
    # Remove "_human_mapped.Read.Count" from column names
    rename_all(~gsub(paste0("_", TEXT, "_mapped.TPM"), "", .))
  return(df)
}
# 2 - Data Processing ----
## Process Kmer data
process_kmers <- function(df_Kmers, df_metadata){
  df <- df_Kmers %>%
    separate(col = colnames(df_Kmers)[1], into = strsplit(colnames(df_Kmers)[1], ",")[[1]], sep = ",") %>%
    mutate(abundance = as.numeric(abundance),
           count = as.numeric(count)) %>%
    merge(df_metadata, by = "sample", by.y = "Sample", all.x = TRUE) %>%
    filter(Description != "Soil Virome") %>%
    rename("Sample" = "sample")
  return(df)
}

# Merge read count data
merge_read_counts <- function(viral_read_counts, human_read_counts, raw_multiqc){
  viral_read_counts <- viral_read_counts %>%
    rename("ViralReads" = "Reads")
  human_read_counts <- human_read_counts %>%
    rename("HumanReads" = "Reads")
  df_all <- merge(raw_multiqc, human_read_counts, by = "Sample", all.x = TRUE) %>%
    merge(viral_read_counts, by = "Sample", all.x = TRUE) %>%
    #Calculate Human %
    mutate(HumanPercent = HumanReads / total_sequences * 100) %>%
    # Calculate Viral %
    mutate(ViralPercent = ViralReads / total_sequences * 100) %>%
    # Pivot longer for HumanPercent, rRNApercentage
    pivot_longer(cols = c("HumanPercent", "rRNApercentage", "ViralPercent"),
                 names_to = "Type", values_to = "ReadPercentage") %>%
    # Change HumanPercent to "Human", rRNApercentage to "rRNA" and ViralPercent to "Viral"
    mutate(Type = case_when(
      Type == "HumanPercent" ~ "Human",
      Type == "rRNApercentage" ~ "rRNA",
      Type == "ViralPercent" ~ "Viral"
    )) %>%
    filter(Description != "Soil Virome") %>%
    mutate(Description = factor(Description, levels = c("Untreated Virome",
                                                        "DNase Treated Virome",
                                                        "MDA Amplified Virome",
                                                        "Metagenome")))
  return(df_all)
}

# 3 - Data Visualisation ----
# Palettes
Replicates_palette <- c("#0077BB","#0077BB","#0077BB",
                        "#33BBEE", "#33BBEE", "#33BBEE",
                        "#009988", "#009988", "#009988",
                        "#CC3311", "#CC3311", "#CC3311")

Description_palette <- c("#0077BB", "#33BBEE", "#009988", "#CC3311")

# Theme for the manuscript
theme_paper <- function(){
  theme_bw() +
    theme(
      plot.title = element_text(size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)) +
    # Remove grid lines
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

## Plot Kmer data
library(ggplot2)
library(scales)

plot_Kmers <- function(df){
  plt <- df %>%
    filter(abundance > 0 & count > 0) %>%
    mutate(ShortSamples = factor(ShortSamples, levels = c("+DNase 1", "+DNase 2", "+DNase 3",
                                                          "-DNase 1", "-DNase 2", "-DNase 3",
                                                          "MDA 1", "MDA 2", "MDA 3",
                                                          "MetaG 1", "MetaG 2", "MetaG 3"))) %>%
    ggplot(aes(x = abundance, y = count)) +
    geom_line(aes(colour = ShortSamples), alpha = 0.5, linewidth = 1) +
    labs(x = "Kmer Abundance", y = "Number of Kmers") +
    scale_x_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_colour_manual(values = Replicates_palette) +
    theme_bw() +
    theme(
      axis.text.x = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  
  return(plt)
}

plot_raw_reads <- function(df){
  plt <- df %>%
    filter(Description != "Soil Virome") %>%
    mutate(Description = factor(Description, levels = c("Untreated Virome",
                                                        "DNase Treated Virome",
                                                        "MDA Amplified Virome",
                                                        "Metagenome"))) %>%
    ggplot(aes(x = ShortSamples, y = total_sequences/1e6, fill = Description)) +
    geom_bar(stat = "identity") +
    labs(y = expression("Number of Reads (x"~10^6~")")) +
    theme_bw() +
    # Remove grid lines
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_fill_manual(values = Description_palette)
  return(plt)
}

plot_read_percentages <- function(df){plt <-ggplot(df %>% filter(Description != "Soil Virome"), aes(x = ShortSamples, y = ReadPercentage, fill = Description)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    labs(y = "Percentage (%)",
         fill = "Type") +
    scale_fill_manual(values = Description_palette) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    facet_wrap(~Type, scales = "free_y") +
    # Make facet label background white
    theme(strip.background = element_rect(fill = "white"))
  return(plt)
}