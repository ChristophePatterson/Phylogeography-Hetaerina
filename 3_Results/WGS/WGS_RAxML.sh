#!/bin/bash

#SBATCH -c 15 
#SBATCH --mem=10G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:10G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 00:15:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here
module purge
module load gcc/native
module load gdal
module load geos/3.10.1
module load r/4.2.1
module load bioinformatics
module load raxml/8.2.12 


#Output directory
## Input directory
input_dir=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/VCF_chrom_r10000)
phy_file=(WGS_titia_chr1-12)
cd $input_dir
output_dir=($input_dir/results/RAxML)

mkdir -p $output_dir

Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/WGS/vcf_2_phy.R $input_dir $output_dir $phy_file

cd $output_dir

# Run RAxML
rm *$phy_file.raxmlASC_GTRGAM_HPCaT15*
raxmlHPC-PTHREADS -s $phy_file.phy -f a -m ASC_GTRGAMMA --asc-corr=lewis -n $phy_file.raxmlASC_GTRGAM_HPCaT15_v1 -p 12345 -x 12345 -N 100 -T 15

cd ~
# Rscript /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/ddRAD_Durham_and_Sheffield/9_RAxML_plot_tree_array.R







