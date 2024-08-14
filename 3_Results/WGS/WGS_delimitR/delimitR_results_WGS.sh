#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=10G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:10G       # $TMPDIR space required on each compute node, up to 400G.
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

output_dir=(delimitR/${sp_0}_${sp_1}_${sp_2}_nreps1000_MIG1_trees_3_sbatch)

# Select directory that delimitR has inputed files into
cd $output_dir

# convert .obs name that delimitR can interpret.
cp ${VCF}_MSFS.obs ${library}_MSFS.obs

# Remove previous prior runs
rm -r Prior/
rm Binned_Processed_${library}_*_MSFS.obs
rm ${library}_*_MSFS.obs

echo "Running Rscript"
# Rscript that produces reduced priors
Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/delimitR/delimR_run_delimR_priors.R $dir_path $output_dir $SLURM_CPUS_PER_TASK $library
#  Once this has been created results plotted using this script
Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/delimitR/delimitR_plot_results.R $dir_path $output_dir $SLURM_CPUS_PER_TASK $library

## Plot each demographic scenario using demographic plots
mkdir -p demo_plots/

## Plot the demographic scenario whose SFS was most closely aligned to the PCA of the obs data
best_par=$(cat best_par_file.txt)

Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/delimitR/ParFileViewer.r $best_par
mv $best_par.pdf demo_plots/Best_${library}_par_file.pdf

# PLot all demographic scerenrios
num_scenarios=$(find $dir_path/$output_dir -type f -name "*.tpl" | wc -l)
num_scenarios=$(find . -type f -name "*.tpl" | wc -l)

for i in $(seq 1 ${num_scenarios})
 do
    echo "${library}_${i}/${library}_${i}_1.par"
    Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/delimitR/ParFileViewer.r ${library}_${i}/${library}_${i}_1.par 
    mv ${library}_${i}/${library}_${i}_1.par.pdf demo_plots/${library}_${i}_1.par.pdf
 done


cd ~