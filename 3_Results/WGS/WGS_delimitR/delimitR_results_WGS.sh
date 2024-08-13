#!/bin/bash

#SBATCH -c 12 
#SBATCH --mem=100G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:100G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 24:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

module load r/4.2.1


# Path to working directory
# Get library name
library=(VCF_chrom_r10000)
VCF=(WGS_titia_chr1-12_filtered)
dir_path=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/$library/)
species="titia"

cd $dir_path

# Assign values to sp1, sp2, and sp3 based on the value of library
sp_0=2
sp_1=2
sp_2=4

output_dir=(delimitR/${sp_0}_${sp_1}_${sp_2}_nreps100_MIG0_trees_3_sbatch)

# Select directory that delimitR has inputed files into
cd $output_dir

# convert .obs name that delimitR can interpret.
cp ${VCF}_MSFS.obs ${library}_MSFS.obs

echo "Running Rscript"
# Rscript that produces reduced priors
### Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/delimitR/delimR_run_delimR_priors.R $dir_path $output_dir $SLURM_CPUS_PER_TASK $library
#  Once this has been created results plotted using this script
Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/delimitR/delimitR_plot_results.R $dir_path $output_dir $SLURM_CPUS_PER_TASK $library

cd ~