#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=50G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:50G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 24:00:00         # time limit in format dd-hh:mm:ss

#SBATCH --array=4,6   # Create 32 tasks, numbers 1 to 32
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

# Get library name
SNP_library=$(sed -n "${line_num}p" /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/2_SNP_calling/library_combinations/library_name)
#Output directory
input_dir=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_v2/)
mkdir -p $input_dir/$SNP_library/G-Phocs
mkdir -p $input_dir/$SNP_library/G-Phocs/gphocs-loci
# Convert bcf into vcf.gz
bcftools view -O z $input_dir/$SNP_library/$SNP_library.all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.bcf > $input_dir/$SNP_library/$SNP_library.all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.vcf.gz

# Number of samples to select from
select_N=(2)

# Remove any previous gphocs files
## rm $input_dir/$SNP_library/G-Phocs/gphocs-loci/*
# Create g-phocs input file
## Rscript /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/3_Results/G-Phocs/vcfR2g-phocs.R $SNP_library $input_dir $select_N

cd $input_dir/$SNP_library/G-Phocs/
# Concat all files into one
cat gphocs-loci/Gphocs_header.txt gphocs-loci/*.gphocs > $SNP_library.gphocs

## Gerenation input files for gphocs
Rscript /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/3_Results/G-Phocs/g-phocs-config-generation.R $SNP_library $input_dir $select_N



