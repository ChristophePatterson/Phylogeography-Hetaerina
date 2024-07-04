#!/bin/bash

#SBATCH -c 64 
#SBATCH --mem=10G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:10G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 70:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

#Modules
module load r/4.2.1
module load bioinformatics
#Name of output SNP library
## Input directory
input_dir=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/VCF_chrom_r10000/)
model_name=(WGS_titia_chr1-12)

cd $input_dir

# Make input files for ruby script
Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/WGS/WGS_SNAPP.R $input_dir $model_name

# Run ruby script for creating xml SNAPP config file
# Uses https://github.com/mmatschiner/tutorials/blob/master/divergence_time_estimation_with_snp_data/README.md

cd results/SNAPP/

ruby /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/SNAPP/snapp_prep.rb -p $model_name.phy -t $model_name.txt \
-c $model_name.con.txt -l 1000000 -x $model_name.xml -o $model_name

/nobackup/tmjj24/apps/beast/bin/beast -threads $SLURM_CPUS_PER_TASK -overwrite $model_name.xml > $model_name.screen.log

/nobackup/tmjj24/apps/beast/bin/treeannotator -burnin 10 $model_name.trees $model_name.trees.Anon

cd ~