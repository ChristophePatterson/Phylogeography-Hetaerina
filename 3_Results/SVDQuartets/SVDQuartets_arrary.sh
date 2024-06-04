#!/bin/bash

#SBATCH -c 20 
#SBATCH --mem=1G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:1G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 24:00:00         # time limit in format dd-hh:mm:ss

#SBATCH --array=1-6   # Create 32 tasks, numbers 1 to 32
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here
module purge
module load gcc/native
module load gdal
module load geos/3.10.1
module load r/4.2.1
module load bioinformatics
module load paup/4.0a

#Name of output SNP library
line_num=$(expr $SLURM_ARRAY_TASK_ID)

# Get library name
SNP_library=$(sed -n "${line_num}p" /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/library_combinations/library_name)
input_dir=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/)

# Run SVDquartet
paup $input_dir/$SNP_library/SVDQuartets/Paup_excute_file.nex -n

cd ~

Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/SVDQuartets/SVDQuartets_plot.R $SNP_library $input_dir






