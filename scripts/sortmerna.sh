#!/bin/bash
#SBATCH --job-name=sortmerna
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=lhillary@ucdavis.edu
#SBATCH --job-name=sortmerna
#SBATCH --nodes=1
#SBATCH -t 10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --partition=high2


# for calculating the amount of time the job takes
begin=`date +%s`
echo $HOSTNAME

source ~/.bashrc
conda activate sortmerna

path=${1}
sample=${2}

cd ${path}

sortmerna --ref /group/jbemersogrp/databases/sortmerna/smr_v4.3_sensitive_db.fasta \
--reads ${path}/${sample} \
--reads ${path}/${sample} \
--idx-dir /home/lhillary/sortmerna/run/idx/ \
--workdir ./sortmerna/${sample}_paired \
--aligned ./sortmerna/${sample}_paired_rrna \
--other ./sortmerna/${sample}_paired_other \
--paired_in \
--out2 \
--fastx \
-v --threads 40
# getting end time to calculate time elapsed
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed

