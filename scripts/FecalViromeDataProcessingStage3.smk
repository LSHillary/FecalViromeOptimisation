import os
import yaml
import glob

#### SET UP ####

stool_samples = ['F3_S19','F4_S20', 'F5_S21', 'F6_S22', 'F7_S23', 'F8_S24', 'M1_S25', 'M2_S26', 'M3_S27', 'GP3_S28', 'GP4_S29', 'GP5_S30']

#### Check Rule ####

rule all:
    input:
        # Error correction
        #CheckEC = expand("Checks/2.2_EC_{sample}.done", sample = stool_samples),
        # PCR read duplicate removal
        #CheckDedup = expand("Checks/2.2_Dedup_{sample}.done", sample = stool_samples),
        # K-mer counting and analysis
        #Check_Khmer_Interleave = expand("Checks/2.4_Khmer-Interleave_{sample}.done", sample = stool_samples),
        #Check_Khmer_Trim = expand("Checks/2.4_Khmer-Trimmed_{sample}.done", sample = stool_samples),
        #Check_Khmer_Counted = expand("Checks/2.4_Khmer-Counted_{sample}.done", sample = stool_samples),
        #Check_Kmer_Abundance = expand("Checks/2.4_Khmer-Abundance_{sample}.done", sample = stool_samples),
        # RNA read quantification
        #Check_SortMeRna = expand("Checks/2.5-SortMeRna_{sample}.done", sample=samples),
        #CheckSortMeRnaFastQC = expand("Checks/2.5-SortMeRnaFastQC_{sample}.done", sample = stool_samples),
        # Assembly
        #CheckMegahitIndividual = expand("Checks/{sample}_IndividualAssembly.done", sample = stool_samples),
        #CheckRenameContigs = expand("Checks/RenameContigs_{sample}.done", sample = stool_samples),
        # Virus Identification
        #CheckGenomad = expand("Checks/4.1_Genomad_{sample}.done", sample = stool_samples),
        #CheckFilteredGenomadResults = expand("Checks/4.2-FilterGenomadResults_{sample}.done", sample = stool_samples),
        # Viral contig clustering
        #CheckClustering = expand("Checks/LeidenClustering_Genomad_{group}.done", group = cluster_groups),
        # Identification of viral protein sequences
        #CheckProteins = expand("Checks/5.3-ExtractProteins_{sample}.done", sample = stool_samples),
        #Check_vOTU_Filtering = "Checks/5.3-ExtractProteins_combined_vOTUs.done",
        # Read Mapping and quantification
        #CombinedCheckMinimap2Mapping = expand("Checks/6-Mapping_combined_{sample}.done", sample = stool_samples),
        #CheckCoverM = "Checks/6.1-CoverM_Stool_Genomad.done",
        #CheckCoverMCounts = "Checks/6.1-CoverM_Stool_Genomad_Counts.done",
        #CheckHumanCoverM = "Checks/6.1-CoverM_Human.done",
        # Host Prediction
        #Check_iPhop = expand("Checks/iphop_Genomad.done", group = iphop_groups),
        # Viral Contig Annotation
        #CheckPharokka = "Checks/8.1-Pharokka_combined_vOTU_proteins.done",
        #CheckDefenseFinder = "Checks/8.4-DefenseFinder_vOTUs.done",
        #CheckPharokkaAll = expand("Checks/8.1-Pharokka_{sample}_proteins.done", sample = stool_samples),
        #CheckDefenseFinder = expand("Checks/8.4-DefenseFinder_{sample}.done", sample = stool_samples),
        # Prokaryote read quantification
        #CheckSingleM = expand("Checks/11-SingleM_{sample}.done", sample = stool_samples),
        # Virus lifestyle prediction
        #CheckBacPhlip = expand("Checks/12-BacPhlip_{sample}.done", sample = stool_samples),

#### PIPELINE ####

#include: "2.2-EC_Dedup.smk"
#include: "2.4-Khmer.smk"
#include: "2.5-SortMeRna.smk" # Indentifies rRNA reads
#include: "3-IndividualAssembly.smk"
#include: "4-VirusIdentificationGenomad.smk"
#include: "5-DereplicationGenomad.smk" # Dereplicates contigs
#include: "6-Mapping_vOTUs.smk" # Maps reads back to the assemblies
#include: "7-HostPrediction.smk"
#include: "8-Annotation.smk" # Annotates proteins
#include: "9-Biogeography.smk" # Clusters vOTUs
#include: "11-SingleM.smk" # Runs SingleM
#include: "12-BacPhlip.smk" # Runs BacPhlip