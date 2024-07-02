#!/bin/bash

#SBATCH -c 24 
#SBATCH --mem=50G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:50G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 24:00:00         # time limit in format dd-hh:mm:ss

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
# Get library and genome names
Library_name=(VCF_chrom_r10000)
library_version=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/)
output_dir=($library_version/$Library_name/G-Phocs/model_runs/)
mkdir -p $output_dir

/nobackup/tmjj24/apps/G-PhoCS/bin/G-PhoCS /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/WGS/G-Phocs_WGS/G-Phocs-WGS-Iso.config -n $SLURM_CPUS_PER_TASK \
> $output_dir/G-Phocs-WGS-Iso.log




