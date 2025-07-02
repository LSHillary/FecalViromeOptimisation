import os
import yaml
import glob

#### SET UP ####

stool_samples = ['F3_S19','F4_S20', 'F5_S21', 'F6_S22', 'F7_S23', 'F8_S24', 'M1_S25', 'M2_S26', 'M3_S27', 'GP3_S28', 'GP4_S29', 'GP5_S30']

#### Check Rule ####

rule all:
    input:
        CheckProcessedFastQC = expand("Checks/FilteredFastQC_{sample}.done", sample = stool_samples)

#### PIPELINE ####

#include: 2-QC.smk