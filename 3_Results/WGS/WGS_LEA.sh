#!/bin/bash

#SBATCH -c 5 
#SBATCH --mem=20G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:20G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 24:00:00         # time limit in format dd-hh:mm:ss

# Specify the tasks to run:
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here

#Modules
module load r/4.2.1

# Library
Library_name="WGS_titia_PAC"
link_filt="10000"

## Input directory
input_dir=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/${Library_name}_r$link_filt/)
cd $input_dir

Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/WGS/WGS_LEA.R $Library_name $input_dir $SLURM_CPUS_PER_TASK
