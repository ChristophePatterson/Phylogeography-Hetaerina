#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=50G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:50G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 24:00:00         # time limit in format dd-hh:mm:ss

#SBATCH --array=6   # Create 32 tasks, numbers 1 to 32
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here
module purge
module load gcc/native
module load gdal
module load geos/3.10.1
module load r/4.2.1
module load bioinformatics
module load bcftools/1.15

#Name of output SNP library
line_num=$(expr $SLURM_ARRAY_TASK_ID)
# line_num=(3)
# Get library name
SNP_library=$(sed -n "${line_num}p" /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/library_combinations/library_name)
#Output directory
input_dir=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/)
output_dir=($input_dir/$SNP_library/BPP/model_runs/)

rm -r $input_dir/$SNP_library/BPP/
mkdir -p $input_dir/$SNP_library/BPP
mkdir -p $output_dir

# Number of loci to subsample
select_N=(2000)

# Create g-phocs input file
Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/BPP/BPP_seq_config_generation.R $SNP_library $input_dir $select_N
