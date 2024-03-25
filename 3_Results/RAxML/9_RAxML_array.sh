#!/bin/bash

#SBATCH -c 24 
#SBATCH --mem=5G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:5G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 72:00:00         # time limit in format dd-hh:mm:ss

#SBATCH --array=1-6   # Create 32 tasks, numbers 1 to 32
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here
module purge
module load gcc/native
module load gdal
module load geos/3.10.1
module load r/4.2.1
module load bioinformatics
module load raxml/8.2.12 

#Name of output SNP library
line_num=$(expr $SLURM_ARRAY_TASK_ID)

# Get library name
SNP_library=$(sed -n "${line_num}p" /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/2_SNP_calling/library_combinations/library_name)
phy_file=$SNP_library

#Output directory
input_dir=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_v3/)

mkdir -p $input_dir/$SNP_library/RAxML/

cd $input_dir/$SNP_library/RAxML/

rm *$phy_file.raxmlASC_GTRGAM_HPCaT15_noCUAJ*

# Run RAxML
pwd
echo $phy_file.noX.noCAUJ.phy
ls

raxmlHPC-PTHREADS -s $phy_file.noX.noCUAJ.phy -f a -m ASC_GTRGAMMA --asc-corr=lewis -n $phy_file.raxmlASC_GTRGAM_HPCaT15_noCUAJ -p 12345 -x 12345 -N 100 -T 24

cd ~
Rscript /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/3_Results/RAxML/9_RAxML_plot_tree_array.R $SNP_library $input_dir

