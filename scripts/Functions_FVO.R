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
import_tpm <- function(filename, TEXT){
  df <- read.csv(filename, sep = "\t") %>%
    # Filter columns that contain "TPM" or "Contig"
    select(contains("TPM"), contains("Contig")) %>%
    # Remove "_human_mapped.Read.Count" from column names
    rename_all(~gsub(paste0("_", TEXT, "_mapped.TPM"), "", .))
  return(df)
}

# Import vOTU clusters
import_vOTU_clusters <- function(filename){
  df <- read.csv(filename, sep = "\t", header = FALSE) %>%
  rename("vOTU" = "V1", "Cluster" = "V2")
  return(df)
}

# Import Lifestyle data
import_bacphlip <- function(filepath){df <- read.csv(filepath, sep = "\t") %>%
  # Create column status with Temperate, Virulent or Unknown
  # Temperate = Temperate < 0.05, Virulent. = Virulent < 0.05, Unknown = other
  mutate(status = ifelse(Temperate >= 0.95, "Temperate", ifelse(Virulent >= 0.95, "Virulent", "Unknown"))) %>%
  rename(vOTU = X)
return(df)
}

# Import Defense genes data
import_defensefinder <- function(filepath){read.csv(filepath, sep = "\t") %>%
    #Create vOTU from "sys_beg" and remove the last underscore and the text after
    mutate(vOTU = gsub("_[0-9]*$", "", sys_beg))
  return(df)
}

# Import Host Data
import_iphop_data <- function(filename){
  df <- read.csv("data/HostPrediction/Host_prediction_to_genus_m90.csv") %>%
    rename(vOTU = Virus) %>%
    
    # Step 1: Group by vOTU and filter by the highest Confidence.score
    group_by(vOTU) %>%
    filter(Confidence.score == max(Confidence.score)) %>%
    
    # Step 2: Split taxonomy information into individual levels
    separate(Host.genus, into = c("HostDomain", "HostPhylum", "HostClass", "HostOrder",
                                  "HostFamily", "HostGenus"), sep = ";", fill = "right") %>%
    
    # Step 3: Clean up taxonomy by removing prefixes like "d__", "p__", etc.
    mutate(across(starts_with("Host"), ~ gsub("^[a-z]__", "", .))) %>%
    
    # Step 4: Collapse certain host taxa for consistency
    mutate(HostPhylum = ifelse(grepl("Bacillota", HostPhylum), "Bacillota", HostPhylum),
           HostGenus = ifelse(grepl("Clostridium", HostGenus), "Clostridium", HostGenus)) %>%
    
    # Step 5: Handle multiple hits with the same Confidence.score
    # For vOTUs with multiple hits and same Confidence.score, keep the first hit
    # and remove the lowest taxonomic levels that don't match
    group_by(vOTU, Confidence.score) %>%
    arrange(vOTU, Confidence.score) %>%
    
    # Identify duplicates and remove mismatching lower taxonomic levels
    mutate(HostGenus = ifelse(n_distinct(HostGenus) > 1, NA, HostGenus),
           HostFamily = ifelse(n_distinct(HostFamily) > 1, NA, HostFamily),
           HostOrder  = ifelse(n_distinct(HostOrder) > 1, NA, HostOrder),
           HostClass  = ifelse(n_distinct(HostClass) > 1, NA, HostClass),
           HostPhylum = ifelse(n_distinct(HostPhylum) > 1, NA, HostPhylum)) %>%
    
    # Step 6: Ungroup after processing
    ungroup() %>% select(-List.of.methods) %>% unique()
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

# Merge ViralContig data
merge_viral_contig_data <- function(GenomadData, ClusteringData, metadata){
  # Identify viral contigs
  viral_contigs <- ClusteringData$Cluster %>%
    # Merge all contigs into a single vector, split individual entries by ","
    strsplit(",") %>%
    # Unlist the vector
    unlist()
  
  # Filter GenomadData for viral contigs
  df <- GenomadData %>%
    rename(Contig = seq_name,
           Sample = sample) %>%
    filter(Contig %in% viral_contigs) %>%
    merge(metadata, by = "Sample", all.x = TRUE) %>%
    mutate(Description = factor(Description, levels = c("DNase Treated Virome", "Untreated Virome", "MDA Amplified Virome", "Metagenome"))) %>%
    # Split taxonomy into individual columns
    separate(col = "taxonomy", into = c("Viruses", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family"), sep = ";", fill = "right") %>%
    filter(Contig %in% viral_contigs)
  return(df)
}

#df_mapping <- tar_read(vOTU_tpm)
#df_Genomad <- tar_read(GenomadData)
#df_BacPhlip <- tar_read(BacPhlipData)
#df_DefenseFinder <- tar_read(DefenseFinderData)
# Merge vOTU data
merge_vOTU_data <- function(df_mapping, df_Genomad, df_BacPhlip,
                            df_DefenseFinder, metadata){
  df <- df_mapping %>%
    merge(df_Genomad, by = "Contig", all.x = TRUE) %>%
    merge(metadata, by = "Sample", all.x = TRUE) %>%
    mutate(Description = factor(Description, levels = c("DNase Treated Virome", "Untreated Virome", "MDA Amplified Virome", "Metagenome"))) %>%
    rename(vOTU = Contig) %>%
    merge(df_BacPhlip, by = "vOTU", all.x = TRUE) %>%
    merge(df_DefenseFinder, by = "vOTU", all.x = TRUE)
  return(df)
}

# Process Alpha Diversity
calculate_alpha_diversity <- function(df_mapping, metadata){
  df_alpha_diversity <- df_mapping %>%
    # Pivot columns longer except for Contig and convert names to Sample and values to TPM
    pivot_longer(cols = -c(Contig), names_to = "Sample", values_to = "TPM") %>%
    # Summarise number of vOTUs with TPM > 0
    filter(TPM > 0) %>%
    group_by(Sample) %>%
    summarise(Richness = n_distinct(Contig)) %>% merge(metadata, by = "Sample", all.x = TRUE)
  return(df_alpha_diversity)
}

# Create Presence/ Absence table
create_PA_table <- function(df, metadata, df_Clustering){
  df_PA <- df %>%
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
  return(df_PA)
}

# Create Beta Diversity table
calculate_beta_diversity <- function(df){
  
  # Step 1: Remove contigs (singletons) present in only one sample
  df_filtered <- df %>%
    rowwise() %>%
    filter(sum(c_across(-Contig) > 0) > 1) %>%  # Keep contigs present in more than one sample
    ungroup()
  
  # Step 2: Renormalize TPM values so that they sum to 1 million in each sample
  df_renormalized <- df_filtered %>%
    mutate(across(-Contig, ~ . / sum(.) * 1e6))  # Renormalize TPM values to sum to 1 million per sample
  
  return(df_renormalized)
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

plot_contig_length <- function(df){
  plt <- df %>% ggplot(aes(x = ShortSamples, y = length, fill = Description)) +
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
  
  return(plt)
}

plot_alpha_diversity <- function(df){
  plt_alpha <- df %>% ggplot(aes(x = ShortSamples, y = Richness, fill = Sample)) +
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
  return(plt_alpha)
}

plot_venn_diagrams <- function(df){
  df_MappedVenn <- df %>% 
    filter(PA == 1) %>%
    select(vOTU, Description) %>%
    split(.$Description)
  
  df_AssembledVenn <- df %>% 
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
  
  plt <- ggarrange(Venn_Mapped, Venn_Assembled, ncol = 2, labels = c("d", "e"),
            common.legend = TRUE, legend = "bottom")
  return(plt)
}

# Plot PCoA
plot_pcoa <- function(df, metadata, Description_palette = Description_palette) {
  
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
#df_mapping <- tar_read(vOTU_tpm)
#GenomadData <- tar_read(GenomadData)
#metadata <- tar_read(metadata)

produce_pcoas <-  function(df_mapping, GenomadData, metadata,
                           Description_palette = Description_palette){
  df_beta <- calculate_beta_diversity(df_mapping)
  plt_pcoa <- plot_pcoa(df_beta, metadata, Description_palette)
  monodnaviria <- GenomadData %>%
    filter(grepl("Monodnaviria", taxonomy))
  df_beta_no_mondna <- df_mapping %>% filter(!Contig %in% monodnaviria$seq_name) %>%
    calculate_beta_diversity()
  
  plt_pcoa_no_mondna <- plot_pcoa(df_beta_no_mondna, metadata, Description_palette) +
    lims(x = c(-0.255, 0.69), y = c(-0.1, 0.15))
  
  ls_pcoas <- ggarrange(plt_pcoa, plt_pcoa_no_mondna, ncol = 1, labels = c("b", "c"),
                        common.legend = TRUE, legend = "right", heights = c(2.5,1))
  return(ls_pcoas)
}

plot_vOTU_stacked_barplot <- function(df_mapping, metadata, df_Genomad){
  df_tax <- df_Genomad %>%
    rename(Contig = seq_name) %>%
    # Split taxonomy
    separate(col = "taxonomy", into = c("Viruses", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family"), sep = ";", fill = "right") %>%
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
         y = "Relative Abundance (%)") +
    theme_bw() +
    theme(legend.position = "bottom",
          # remove Grid lines
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = Description_palette)
  
  return(plt_stacked)
}



