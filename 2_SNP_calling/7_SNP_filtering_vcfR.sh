#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=20G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:20G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 24:00:00         # time limit in format dd-hh:mm:ss

# Specify the tasks to run:
#SBATCH --array=1-6  # Create 6 tasks, numbers 1 to 6
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here

#Modules
module load r/4.2.1
module load bioinformatics
module load bcftools

#Name of output SNP library

line_num=$(expr $SLURM_ARRAY_TASK_ID)

echo "$line_num"

# Get library and genome names

Library_name=$(sed -n "${line_num}p" /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/2_SNP_calling/library_combinations/library_name)
genome=$(sed -n "${line_num}p" /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/2_SNP_calling/library_combinations/genome)

library_version=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_v3/)

echo "Processing database $Library_name using $genome"

Rscript /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/2_SNP_calling/7_SNP_filtering_vcfR.R $Library_name $library_version

## Get the most highest quality samples
cd $library_version/$Library_name/

bcftools stats -s - $Library_name.all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.vcf.gz > ${Library_name}_bcfstats_snps.txt

plot-vcfstats -p bcfstats ${Library_name}_bcfstats_snps.txt

