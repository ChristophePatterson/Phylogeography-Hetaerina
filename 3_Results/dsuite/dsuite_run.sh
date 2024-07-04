#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=2G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:2G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 02:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out
#SBATCH --array=6   # Create 32 tasks, numbers 1 to 32

# Commands to execute start here
module purge
module load gcc/native
module load gdal
module load geos/3.10.1
module load r/4.2.1

# Get library name from slurm array
# line_num=(6)
line_num=$(expr $SLURM_ARRAY_TASK_ID)
echo "$line_num"
# Library name
Library_name=$(sed -n "${line_num}p" /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/library_combinations/library_name)
# Version of library
library_version=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/)

# Run R script
# Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/dsuite/dsuite_input.R $Library_name $library_version

cd $library_version/$Library_name/dsuite/
# run dsuite across all snps
/nobackup/tmjj24/apps/Dsuite/Build/Dsuite Dquartets hybrid_samples.vcf.gz \
pop_file.txt \
-o ${Library_name}_filtered_allChr.txt 