# ViromeDataProcessing
Central codebase for virome data processing pipeline

## Step 0 - Set up and fetch raw data

### Installing Conda

Personally, I'm going to try micromamba:

https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html

`curl micro.mamba.pm/install.sh | bash`

### Setting up the ViromeProcessing conda environment

`micromamba create -n ViromeDataProcessing -f ViromeDataProcessing.yml`

### Non-conda resources

Follow these instructions to install bbtools
`https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/`

## Step 0 - Get data and verify

Run the following command:

`conda activate ViromeDataProcessing`

### Download raw data
For people wanting to reproduce this analysis post publication, download raw reads from the SRA

For me, I copied the data from the UC Davis Genome Center server into

### Run 0-file_renaming.smk
This file will reorganise the data and should be customised for each project.

## Step 1 - Preprocessing
### Run 1-preprocessing.smk
This is done by hashing out all steps but 1-preprocessing in the parent Snakefile ViromeDataProcessing.smk
This step should be run on its own, as the results need to be inspected in order to customise the QC steps.

### Run multiqc on the results from 1-preprocessing.smk
I haven't yet been able to get multiqc to run within a conda environment so I'm relying on the installed version on the UC Davis Farm cluster. This is a bit of a hack, but it works. To load this run:
`module load deprecated/multiqc`

Feel free to adjust to the specifics of your cluster.

Then this command from your working directory:
`cd Virome/FastQC/Raw && multiqc . && mv multiqc_report.html multiqc_raw.html && cd ../../..`

## Step 2 - QC
### Adjust the settings in 2-QC.smk and run
Settings in 2-QC.smk should be adjusted based on the results of 1-preprocessing.smk but the supplied settings are a good starting point.
Run 2-QC.smk by removing the hash in the parent Snakefile ViromeDataProcessing.smk.

### Run multiqc on the results from 2-QC.smk
Repeat the multiqc step as above, but this time on the QC results (don't forget to set up multiqc first).:
`cd Virome/FastQC/Processed/Paired && multiqc . && mv multiqc_report.html multiqc_QC.html && cd ../../..`
You can also run multiqc on the unpaired reads:
`cd Virome/FastQC/Processed/Single && multiqc . && mv multiqc_report.html multiqc_QC.html && cd ../../..`

Inspect the MultiQC report before progressing any further. You want to verify that your filtering settings have had the desired effect.
A word of caution here. If one read is removed, the other read can often be poor quality anyway. If the quality of the singleton reads is poor, it's worth considering not using them further. The rest of the pipeline will run without them but feel free to adapt the code to your needs.

### Run error correction and deduplication
There is some evidence that error correction and deduplication can improve assembly. This step is optional and can be skipped if you want to save time or have other concerns. If you want to run it, remove the hash in the parent Snakefile ViromeDataProcessing.smk for file "2.2-EC_Dedup.smk". Something to note though, as error rates can vary between Illumina runs, if you have resequenced a library, run error correction on each run separately or skip it. If you choose to include error correction and deduplication, read mapping should be performed using the error corrected reads pre-deduplication. This is because the deduplication step will remove reads that are identical, which will cause problems for getting accurate read counts. Duplicate reads can come from both biological and PCR duplicates, so I prefer to keep duplicates in before mapping. To run this step, remove the hash in the parent Snakefile ViromeDataProcessing.smk.

## Step 3 - Assembly