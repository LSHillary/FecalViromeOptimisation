# 0 Setting parameters ----

# Changes Short Sample names to be more readable
DescriptionModifiers <- c("UTV" = "Untreated\nVirome",
                          "DTV" = "DNase\nTreated\nVirome",
                          "MDA" = "MDA\nVirome",
                          "MtG" = "Meta\ngenome")

# Sets the order of the Description levels
DescriptionLevels <- c("Untreated\nVirome","DNase\nTreated\nVirome", "MDA\nVirome", "Meta\ngenome")


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

# Import SingleM data
import_singlem <- function(filename, metadata = metadata){
  df <- read_and_merge_tsvs(filename, "_tax_profile.tsv") %>%
    # Split taxonomy
    separate(taxonomy, into = c("Root","Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") %>%
    # If Phylum = NA, then take value from Domain and fill as paste("Unclassified", Domain)
    mutate(Phylum = ifelse(is.na(Phylum), paste("Unclassified", Domain), Phylum)) %>%
    # Remove "d_" and "p_" from Phylum
    mutate(Phylum = gsub("d__", "", Phylum)) %>%
    mutate(Phylum = gsub("p__", "", Phylum)) %>%
    mutate(Phylum = gsub(" ","", Phylum)) %>%
    mutate(Phylum = ifelse(grepl("Unclassified", Phylum), "Unclassified", Phylum)) %>%
    # Rename all Phyla containing text "Bacillota" to Bacillota
    mutate(Phylum = ifelse(grepl("Bacillota", Phylum), "Bacillota", Phylum)) %>%
    select(c(sample, Phylum, coverage)) %>%
    group_by(sample, Phylum) %>%
    summarise(coverage = sum(coverage)) %>%
    group_by(sample) %>%
    mutate(RelAbund = coverage/sum(coverage) *100,
           Phylum = factor(Phylum, levels = c("Actinomycetota","Bacillota", "Bacteroidota",
                                              "Desulfobacterota", "Pseudomonadota", "Verrucomicrobiota",
                                              "Unclassified"))) %>%
    rename(Sample = sample) %>%
    # Merge with metadata
    merge(metadata, by = "Sample")
  return(df)}

# Import Lifestyle data
import_bacphlip <- function(filepath){
  df <- read_and_merge_tsvs(filepath, "_FilteredContigs.fna.bacphlip") %>%
    rename(vOTU = `...1`) %>%
    # Create column status with Temperate, Virulent or Unknown
    # Temperate = Temperate < 0.05, Virulent. = Virulent < 0.05, Unknown = other
    mutate(status = ifelse(Temperate >= 0.95, "Temperate", ifelse(Virulent >= 0.95, "Virulent", "Unclassified"))) %>%
    # Create column Provirus if Provirus is in the vOTU name
    mutate(Provirus = ifelse(grepl("provirus", vOTU), "Provirus", "Virus")) %>%
    rename(Sample = sample)
  return(df)
}

import_bacphlip_vOTU <- function(filepath){
  df <- read.csv(filepath, sep = "\t", header = TRUE) %>%
    rename(vOTU = X) %>%
    # Create column status with Temperate, Virulent or Unknown
    # Temperate = Temperate < 0.05, Virulent. = Virulent < 0.05, Unknown = other
    mutate(status = ifelse(Temperate >= 0.95, "Temperate", ifelse(Virulent >= 0.95, "Virulent", "Unclassified"))) %>%
    # Create column Provirus if Provirus is in the vOTU name
    mutate(Provirus = ifelse(grepl("provirus", vOTU), "Provirus", "Virus"))
  return(df)
  }

# Import Pharokka data
import_pharokka <- function(filepath){
  df <- read_and_merge_tsvs(filepath, "_pharokka_proteins_full_merged_output.tsv") %>%
    mutate(Contig = gsub("_[0-9]*$", "", ID))
  return(df)
}

import_pharokka_vOTU <- function(filepath){
  df <- read.csv(filepath, sep = "\t") %>%
    mutate(Contig = gsub("_[0-9]*$", "", ID))
  return(df)
}

# Import Defense genes data
import_defensefinder <- function(filepath){
  df <- read_and_merge_tsvs(filepath, "_FilteredProteins_defense_finder_systems.tsv") %>%
    #Create ViralSequence from "sys_beg" and remove the last underscore and the text after
    mutate(ViralSequence = gsub("_[0-9]*$", "", sys_beg))
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
    filter(Description != "Soil Virome")
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
    mutate(DescriptionLong = DescriptionModifiers[Description]) %>%
    mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels)) %>%
    # Split taxonomy into individual columns
    separate(col = "taxonomy", into = c("Viruses", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family"), sep = ";", fill = "right") %>%
    filter(Contig %in% viral_contigs)
  return(df)
}

# Merge vOTU data
merge_vOTU_data <- function(df_mapping, df_Genomad, df_BacPhlip,
                            df_DefenseFinder, metadata){
  df <- df_mapping %>%
    merge(df_Genomad, by = "Contig", all.x = TRUE) %>%
    merge(metadata, by = "Sample", all.x = TRUE) %>%
    mutate(Description = factor(Description, levels = DescriptionLevels)) %>%
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
Replicates_palette <- c(
  "-DV1" = "#56B4E9", "-DV2" = "#56B4E9", "-DV3" = "#56B4E9",
  "+DV1" = "#0072B2", "+DV2" = "#0072B2", "+DV3" = "#0072B2",
  "MDA1" = "#009E73", "MDA2" = "#009E73", "MDA3" = "#009E73",
  "MtG1" = "#FF7F0E", "MtG2" = "#FF7F0E", "MtG3" = "#FF7F0E")

Description_palette <- c("UTV" = "#56B4E9",
                         "DTV" = "#0072B2",
                         "MDA" = "#009E73",
                         "MtG" = "#FF7F0E")

Stacked_palette <- c("dsDNA" = "#0072B2", "ssDNA" = "#56B4E9",
                     "Unassigned" = "#808080",
                     "Other Caudoviricetes" = "#56B4E9",
                     "Crassvirales" = "#0072B2",
                     "Microviridae" = "#D55E00",
                     "Other Monodnaviria" = "#E69F00",
                     "Herpesviridae" = "#009E73",
                     "Temperate" = "#56B4E9",
                     "Virulent" = "#0072B2",
                     "Provirus" = "#56B4E9",
                     "Virus" = "#0072B2",
                     "Unclassified" = "#808080")



Phyla_palette <- c(
  "Actinomycetota"    = "#E69F00", 
  "Bacillota"         = "#56B4E9",
  "Bacteroidota"      = "#009E73",
  "Desulfobacterota"  = "#F0E442",
  "Pseudomonadota"    = "#0072B2",
  "Verrucomicrobiota" = "#D55E00",
  "Unclassified"      = "#808080"
)

# Theme for the manuscript

## Plot Kmer data
plot_Kmers <- function(df){
  plt <- df %>%
    filter(abundance > 0 & count > 0) %>%
    mutate(ShortSamples = factor(ShortSamples, levels = c("+DV1", "+DV2", "+DV3",
                                                          "-DV1", "-DV2", "-DV3",
                                                          "MDA1", "MDA2", "MDA3",
                                                          "MtG1", "MtG2", "MtG3"))) %>%
    ggplot(aes(x = abundance, y = count)) +
    geom_line(aes(colour = ShortSamples), alpha = 0.5, linewidth = 0.3) +
    labs(x = "Kmer Abundance", y = "Kmers") +
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

# Plots barplots with CLD letter grouping
plot_means_with_cld <- function(df, value_col, group_col, PALETTE) {
  # Perform ANOVA
  aov_model <- aov(as.formula(paste(value_col, group_col, sep = " ~ ")), data = df)
  
  # Perform Tukey HSD post hoc test
  tukey_result <- TukeyHSD(aov_model)
  
  # Generate Compact Letter Display
  cld <- multcompLetters(tukey_result[[group_col]][, "p adj"])$Letters
  
  # Summarize data to calculate mean and standard deviation
  df_summary <- df %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      mean = mean(.data[[value_col]]),
      sd = sd(.data[[value_col]])
    )
  
  # Add the CLD to the summary dataframe
  df_summary$cld <- cld[as.character(df_summary[[group_col]])]
  
  # Calculate dynamic y-limit for text placement
  y_min <- max(df_summary$mean + df_summary$sd) # Minimum y-value
  text_y_position <- y_min + 0.1 * abs(y_min)   # 10% above y-limit
  
  # Generate barplot
  plt <- ggplot(df_summary, aes(x = .data[[group_col]], y = mean, fill = .data[[group_col]])) +
    geom_bar(stat = "identity", show.legend = TRUE) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    geom_text(aes(label = cld, y = text_y_position), vjust = 0, size = 3) + # Position labels
    labs(x = group_col, y = "Mean ± SD") +
    theme_bw() +
    labs(x = "Processing Method") +
    scale_fill_manual(values = PALETTE)
  
  return(list(
    plot = plt,
    tukey = tukey_result,
    cld = cld
  ))
}

# Plot vOTU Venn diagrams
plot_venn_diagrams <- function(df) {
  library(dplyr)
  library(tidyr)
  library(ggVennDiagram)
  library(ggpubr)
  
  DescriptionLevelOrder <- c("UTV", "DTV", "MDA", "MtG")
  
  # Mapped vOTUs
  df_MappedVenn <- df %>%
    filter(PA == 1) %>%
    select(vOTU, Description) %>%
    split(.$Description)
  
  Mapped_sets <- lapply(df_MappedVenn, function(x) unique(x$vOTU)) %>%
    .[DescriptionLevelOrder] %>%
    setNames(c("DNase\nTreated\nVirome", "Untreated\nVirome", "MDA\nVirome", "Meta\ngenome"))
  
  # Assembled vOTUs
  df_AssembledVenn <- df %>%
    filter(AssembledPresence == 1) %>%
    select(vOTU, Description) %>%
    split(.$Description)
  
  Assembled_sets <- lapply(df_AssembledVenn, function(x) unique(x$vOTU)) %>%
    .[DescriptionLevelOrder] %>%
    setNames(c("DNase\nTreated\nVirome", "Untreated\nVirome", "MDA\nVirome", "Meta\ngenome"))
  
  # Correct Mapped-only logic
  get_mapped_only_votus <- function(df, description_label) {
    df_group <- df %>% filter(Description == description_label)
    mapped_any <- df_group %>%
      filter(PA == 1) %>%
      pull(vOTU) %>%
      unique()
    assembled_any <- df_group %>%
      filter(AssembledPresence == 1) %>%
      pull(vOTU) %>%
      unique()
    setdiff(mapped_any, assembled_any)
  }
  
  MappedOnly_sets <- setNames(lapply(DescriptionLevelOrder, function(desc) {
    get_mapped_only_votus(df, desc)
  }), c("DNase\nTreated\nVirome", "Untreated\nVirome", "MDA\nVirome", "Meta\ngenome"))
  
  # Venn Diagrams ------------------------------------------------------------
  
  Venn_Mapped <- ggVennDiagram(Mapped_sets, label = "count", edge_size = 0, set_size = 0) +
    theme(legend.position = "bottom") +
    annotate("text", x = 0.14, y = 0.81, label = "DNase\nTreated\nVirome", size = 2.5) +
    annotate("text", x = 0.32, y = 0.84, label = "Untreated\nVirome", size = 2.5) +
    annotate("text", x = 0.65, y = 0.84, label = "MDA\nVirome", size = 2.5) +
    annotate("text", x = 0.85, y = 0.81, label = "Meta\ngenome", size = 2.5) +
    scale_fill_distiller(limits = c(0, 200), palette = "Blues", direction = 1) +
    theme(plot.margin = margin(0.5, 0, 0, 0, "cm")) +
    labs(title = "Mapped vOTUs", fill = "Count") +
    theme(plot.title = element_text(hjust = 0.5, margin = margin(0, 0, 0.4, 0, "cm")))
  
  Venn_Assembled <- ggVennDiagram(Assembled_sets, label = "count", edge_size = 0, set_size = 0) +
    scale_fill_distiller(palette = "Blues", direction = 1) + 
    annotate("text", x = 0.14, y = 0.81, label = "DNase\nTreated\nVirome", size = 2.5) +
    annotate("text", x = 0.32, y = 0.85, label = "Untreated\nVirome", size = 2.5) +
    annotate("text", x = 0.65, y = 0.85, label = "MDA\nVirome", size = 2.5) +
    annotate("text", x = 0.85, y = 0.81, label = "Meta\ngenome", size = 2.5) +
    theme(legend.position = "bottom") +
    scale_fill_distiller(limits = c(0, 200), palette = "Blues", direction = 1) +
    theme(plot.margin = margin(0.5, 0, 0, 0, "cm")) +
    labs(title = "Assembled vOTUs", fill = "Count") +
    theme(plot.title = element_text(hjust = 0.5, margin = margin(0, 0, 0.4, 0, "cm")))
  
  Venn_MappedOnly <- ggVennDiagram(MappedOnly_sets, label = "count", edge_size = 0, set_size = 0) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    annotate("text", x = 0.14, y = 0.81, label = "DNase\nTreated\nVirome", size = 2.5) +
    annotate("text", x = 0.32, y = 0.85, label = "Untreated\nVirome", size = 2.5) +
    annotate("text", x = 0.65, y = 0.85, label = "MDA\nVirome", size = 2.5) +
    annotate("text", x = 0.85, y = 0.81, label = "Meta\ngenome", size = 2.5) +
    theme(legend.position = "bottom") +
    scale_fill_distiller(limits = c(0, 200), palette = "Blues", direction = 1) +
    theme(plot.margin = margin(0.5, 0, 0, 0, "cm")) +
    labs(title = "Mapped-only vOTUs (failed assembly filters)", fill = "Count") +
    theme(plot.title = element_text(hjust = 0.5, margin = margin(0, 0, 0.4, 0, "cm")))
  
  # Combine main two into one panel
  plt_Venn <- ggarrange(Venn_Assembled, Venn_Mapped, ncol = 2, labels = c("a", "b"),
                        common.legend = TRUE, legend = "bottom")
  
  return(list(
    Venn_Assembled   = Assembled_sets,
    Venn_Mapped      = Mapped_sets,
    Venn_MappedOnly  = MappedOnly_sets,
    plt_Venn         = plt_Venn,
    plt_MappedOnly   = Venn_MappedOnly
  ))
}

# Plot PCoA
plot_pcoa <- function(df, metadata, Description_palette = Description_palette) {
  metadata <- metadata %>% 
    mutate(DescriptionLong = DescriptionModifiers[Description]) %>%
    mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels))
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
    mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels))
  
  # Step 5: Calculate the percentage of variance explained by each axis
  variance_explained <- 100 * pcoa_result$eig[1:2] / sum(pcoa_result$eig)
  
  # Step 6: Plot the PCoA results, colored by 'Description'
  pcoa_plot <- ggplot(pcoa_axes, aes(x = PCoA1, y = PCoA2, color = DescriptionLong, shape = DescriptionLong)) +
    geom_point(size = 4, alpha = 0.5) +
    labs(x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
         y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
         fill = "Processing Method") +
    scale_shape_manual(values = c(16, 17, 15, 18)) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    scale_colour_manual(values = Description_palette) +
    guides(fill = guide_legend(title = "Sample\nType", nrow = 2))  +
    scale_x_continuous(expand = expansion(mult = 0.05)) +  # Expand x-axis by 5%
    scale_y_continuous(expand = expansion(mult = 0.05)) +    # Expand y-axis by 5%
    coord_fixed()
  # Calculate PERMANOVA
  pcoa_permanova <- adonis2(df_beta_diversity ~ Description, data = pcoa_axes, permutations = 999)
  pcoa_permanova

  # Pairwise PERMANOVA
  pairwise_permanova <- pairwise.adonis2(df_beta_diversity ~ Description, data = pcoa_axes, p.adjust.method = "BH", permutations = 999)
  pairwise_permanova
  
  # Create a list of results
  results_list <- list(
    pcoa_plot = pcoa_plot,
    permanova = pcoa_permanova,
    pairwise_permanova = pairwise_permanova
  )
  # Return the plot
  return(results_list)
}

# Produces all PCoA plots and permanova results
produce_pcoas <-  function(df_mapping, stacked, GenomadData, metadata,
                           Description_palette = Description_palette){
  df_beta <- calculate_beta_diversity(df_mapping)
  plt_pcoa <- plot_pcoa(df_beta, metadata, Description_palette)
  monodnaviria <- GenomadData %>%
    filter(grepl("Monodnaviria", taxonomy))
  df_beta_no_mondna <- df_mapping %>% filter(!Contig %in% monodnaviria$seq_name) %>%
    calculate_beta_diversity()
  
  ssDNA <- stacked$df_Baltimore %>%
    filter(DNA == "ssDNA") %>% pull(Contig)
  df_beta_no_ssDNA <- df_mapping %>% filter(!Contig %in% ssDNA) %>%
    calculate_beta_diversity()
    
  plt_pcoa_no_ssDNA <- plot_pcoa(df_beta_no_ssDNA, metadata, Description_palette)
  
  plt_pcoas <- list(
    pcoa_all = plt_pcoa$pcoa_plot,
    pcoa_no_ssDNA = plt_pcoa_no_ssDNA$pcoa_plot,
    permanova_all = plt_pcoa$permanova,
    permanova_no_ssDNA = plt_pcoa_no_ssDNA$permanova,
    pairwise_permanova_all = plt_pcoa$pairwise_permanova,
    pairwise_permanova_no_ssDNA = plt_pcoa_no_ssDNA$pairwise_permanova
  )
  return(plt_pcoas)
}

# Plot stacked barplot for vOTU data
plot_vOTU_stacked_barplot <- function(df_mapping, metadata, df_Genomad){
  df_tax <- df_Genomad %>%
    rename(Contig = seq_name) %>%
    # Split taxonomy
    separate(col = "taxonomy", into = c("Viruses", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family"), sep = ";", fill = "right") %>%
    select(Contig, Viruses, Realm, Kingdom, Phylum, Class, Order, Family)
  
  dsDNA_Monodnaviria <- df_tax %>%
    filter(Class == "Papovaviricetes") %>%
    pull(Contig)
  
  ssDNA_extra <- df_tax %>%
    filter(Family %in% c("Alphasatellitidae","Spiraviridae","Anelloviridae","Tolecusatellitidae")) %>%
    pull(Contig)
  
  df_tax <- df_tax %>%
    mutate(DNA =
             # If realm == monodnaviria and class != Papovaviricetes, DNA = ssDNA, otherwise dsDNA
             ifelse(Realm == "Monodnaviria" & !(Contig %in% dsDNA_Monodnaviria), "ssDNA",
                    ifelse(Realm == "Monodnaviria" & (Contig %in% dsDNA_Monodnaviria), "dsDNA",
                           ifelse(Realm == "Unassigned","Unassigned",
                                  ifelse(Realm %in% c("Duplodnaviria", "Adnaviria", "Varidnaviria"), "dsDNA", "Unassigned")))))
  
  # Stacked barplots
  df_Baltimore <- df_mapping %>%
    pivot_longer(cols = -c(Contig), names_to = "Sample", values_to = "TPM") %>%
    filter(TPM > 0) %>%
    merge(metadata, by = "Sample", all.x = TRUE) %>%
    mutate(DescriptionLong = DescriptionModifiers[Description],
           rep = substr(ShortSamples, nchar(ShortSamples), nchar(ShortSamples))) %>%
    mutate(DescriptionLong = factor(DescriptionLong, levels = DescriptionLevels)) %>%
    # Merge with df_Genomad
    merge(df_tax, by = "Contig", all.x = TRUE) %>%
    # Modify NA in Realm to "Unassigned"
    mutate(DNA = ifelse(is.na(DNA), "Unassigned", DNA))
  
  df_DNA_stacked <- df_Baltimore %>%
    # Summarise TPM values in each sample to Realm level
    group_by(DNA, rep, DescriptionLong) %>%
    summarise(TPM = sum(TPM)) %>%
    mutate(RelAbund = TPM / 10000)
  
  plt_DNA_stacked <- df_DNA_stacked %>%
    ggplot(aes(x = rep, y = RelAbund, fill = DNA)) +
    geom_bar(stat = "identity") +
    labs(x = "Processing Method",
         y = "Relative Abundance (%)") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          # remove Grid lines
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = Stacked_palette) +
    facet_wrap(~DescriptionLong, nrow = 1)

  df_stacked <- df_Baltimore %>%
    mutate(
      Family = coalesce(Family, Order, Class, Phylum, Kingdom, Realm, Viruses),
      Order = coalesce(Order, Class, Phylum, Kingdom, Realm, Viruses)) %>%
    # Change Family where Class = "Caudoviricetes" and Family != "Crassvirales" to "Other Caudoviricetes"
    mutate(Family = ifelse(Class == "Caudoviricetes" & Family != "Crassvirales", "Other Caudoviricetes", Family)) %>%
    # Change Family where Realm = "Monodnaviria" and Family != "Microviridae" to "Other Monodnaviria"
    mutate(Family = ifelse(Realm == "Monodnaviria" & Family != "Microviridae", "Other Monodnaviria", Family)) %>%
    # Change Family = NA to Unassigned
    mutate(Family = ifelse(is.na(Family), "Unassigned", Family)) %>%
    group_by(Family, ShortSamples, rep, DescriptionLong) %>%
    summarise(TPM = sum(TPM)) %>%
    mutate(RelAbund = TPM / 10000) %>%
    mutate(Family = factor(Family, levels = c("dsDNA", "ssDNA","Crassvirales",
                                              "Other Caudoviricetes",
                                              "Herpesviridae",
                                              "Microviridae", "Other Monodnaviria","Unassigned")))
  
  df_stacked_summary <- df_stacked %>%
    group_by(DescriptionLong, Family) %>%
    summarise(mean_RelAbund = mean(RelAbund),
              SD_RelAbund = sd(RelAbund))
  
  plt_stacked <- df_stacked %>% ggplot(aes(x = rep, y = RelAbund, fill = Family)) +
    geom_bar(stat = "identity") +
    labs(x = "Processing Method",
         y = "Relative\nAbundance (%)") +
    theme_bw(base_size = 10) +
    theme(legend.position = "right",
          # remove Grid lines
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = Stacked_palette) +
    facet_wrap(~DescriptionLong, nrow = 1) +
    # Remove x axis ticks, axis label and tick labels
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
    # Make facet label backgrounds white
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(colour = "black")) +
    guides(fill = guide_legend(title = "Lifestyle", nrow = 3)) +  # Ensure colors are in one row
    theme(legend.position = "bottom")
  
  ls_stacked <- list(plt_stacked = plt_stacked,
                      plt_DNA_stacked = plt_DNA_stacked,
                      df_Baltimore = df_Baltimore,
                      dsDNA_Monodnaviria = dsDNA_Monodnaviria,
                      ssDNA_extra = ssDNA_extra)
  
  return(ls_stacked)
}



