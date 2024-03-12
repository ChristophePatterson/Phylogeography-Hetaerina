#!/bin/bash

#SBATCH -c 15 
#SBATCH --mem=5G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:5G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 72:00:00         # time limit in format dd-hh:mm:ss

#SBATCH --array=1-2   # Create 32 tasks, numbers 1 to 32
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here
module purge
module load gcc/native
module load gdal
module load geos/3.10.1
module load r/4.2.1
module load bioinformatics
module load raxml/8.2.12 

# Convert vcf2phy
Rscript /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/ddRAD_Durham_and_Sheffield/9_vcf2phy_ASC_array.R

#Name of output SNP library
line_num=$(expr $SLURM_ARRAY_TASK_ID)

# Get library name

SNP_library=$(sed -n "${line_num}p" /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/ddRAD_Durham_and_Sheffield/library_combinations/library_name)
phy_file=$SNP_library

#Output directory
output_dir=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_Pool_X_HetTit1.0_HetAmer1.0/$SNP_library/RAxML)

mkdir -p $output_dir

cd $output_dir

rm *$phy_file.raxmlASC_GTRGAM_HPCaT15_noCUAJa02_v3*

# Run RAxML
raxmlHPC-PTHREADS -s $phy_file.ASC.noCAUJa02.phy -f a -m ASC_GTRGAMMA --asc-corr=lewis -n $phy_file.raxmlASC_GTRGAM_HPCaT15_noCUAJa02_v3 -p 12345 -x 12345 -N 100 -T 15

cd ~
Rscript /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/ddRAD_Durham_and_Sheffield/9_RAxML_plot_tree_array.R







