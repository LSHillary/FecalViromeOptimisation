######################################################################################################
# HUMAN FECAL VIROME OPTIMISATION PIPELINE (PARENT WORKFLOW)
# Author: Luke Hillary
# Date: 2025-09-19
#
# Description:
# This Snakemake pipeline orchestrates a modular, end-to-end analysis of human fecal viromes
# prepared with different methods. It includes modules for read preprocessing, assembly, 
# viral contig identification, dereplication, functional annotation, host prediction,
# and preparation for SRA submission.
#
# The pipeline is designed for **flexible activation** of specific workflow modules via
# commenting/uncommenting the `include:` statements under the "PIPELINE" section.
#
# Stage 3 contains all data processing steps after initial read QC
#
# Sample metadata (stool_samples) is defined in the script and should be updated as needed.
#
# Output:
#   - Each step writes its own `.done` flag in the `Checks/` directory
#
# Requirements:
#   - Apptainer/Singularity and required containers (e.g., sra-human-scrubber.sif)
#   - Conda environment or modules providing required tools (e.g., unpigz, sortmerna, genomad, etc.)
#
# Notes:
#   - All steps are sample-parallelized via `stool_samples` defined at the top.
#   - This PARENT Snakefile assumes modular structure and follows a “comment-to-disable” approach.
#   - Ensure any required parameters, groupings (e.g. `cluster_groups`, `iphop_groups`) are defined.
######################################################################################################

import os
import yaml
import glob

#### SET UP ####
#'F3_S19',
stool_samples = ['F4_S20', 'F5_S21', 'F6_S22', 'F7_S23', 'F8_S24', 'M1_S25', 'M2_S26', 'M3_S27', 'GP3_S28', 'GP4_S29', 'GP5_S30']

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
        #CheckPharokkaAll = expand("Checks/8.1-Pharokka_{sample}_proteins.done", sample = stool_samples),
        # Prokaryote read quantification
        #CheckSingleM = expand("Checks/11-SingleM_{sample}.done", sample = stool_samples),
        # Virus lifestyle prediction
        #CheckBacPhlipVOTU = "Checks/12-BacPhlip_vOTUs.done",
        #CheckHumanReadScrubbing = expand("Checks/13.2-HumanScrub_{sample}.done", sample = stool_samples),
        #CheckCompressedScrubbedReads = expand("Checks/13.3-Compress_{sample}.done", sample = stool_samples),
        CheckHumanReadScrubbingPipe = expand("Checks/13.2-HumanScrub_Pipe_{sample}.done", sample = stool_samples),

#### PIPELINE ####
# Turn on/off different parts of the pipeline by commenting/uncommenting the relevant include lines
#include: "2.2-EC_Dedup.smk" # Error correction and PCR duplicate removal
#include: "2.4-Khmer.smk" # K-mer profiling
#include: "2.5-SortMeRna.smk" # Indentifies rRNA reads
#include: "3-IndividualAssembly.smk" # Assembles reads into contigs
#include: "4-VirusIdentificationGenomad.smk" # Identifies viral contigs with Genomad
#include: "5-DereplicationGenomad.smk" # Dereplicates contigs
#include: "6-Mapping_vOTUs.smk" # Maps reads back to the assemblies
#include: "7-HostPrediction.smk" # Predicts hosts for viral contigs using iPhop
#include: "8-Annotation.smk" # Annotates proteins
#include: "9-Biogeography.smk" # Clusters vOTUs
#include: "11-SingleM.smk" # Runs SingleM to profile prokaryotic communities
#include: "12-BacPhlip.smk" # Runs BacPhlip for virus lifestyle prediction
include: "13-NCBI_Submission.smk" # Scrubs human reads and gathers annotations on all viral sequences for NCBI submission