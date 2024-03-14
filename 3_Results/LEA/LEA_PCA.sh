#!/bin/bash

#SBATCH -c 3 
#SBATCH --mem=10G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:10G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 72:00:00         # time limit in format dd-hh:mm:ss

#SBATCH --array=1-6   # Create 32 tasks, numbers 1 to 32
#SBATCH --output=/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_v2/slurm-%x.%j.out

# Commands to execute start here
module purge
module load gcc/native
module load gdal
module load geos/3.10.1
module load r/4.2.1

# Get library name from slurm array
line_num=$(expr $SLURM_ARRAY_TASK_ID)
echo "$line_num"
# Library name
Library_name=$(sed -n "${line_num}p" /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/2_SNP_calling/library_combinations/library_name)
# Version of library
library_version=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_v2/)

# Run R script
Rscript /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/3_Results/LEA/LEA_PCA.R $Library_name $library_version