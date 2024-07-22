#!/bin/bash

#SBATCH -c 24 
#SBATCH --mem=10G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:10G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 48:00:00         # time limit in format dd-hh:mm:ss

# Specify the tasks to run:
#SBATCH --output=slurm-%x.%j.out

module load r/4.1.2

# Demultiplexed file locations
# Number of samples to use per population
select_N=("2")


input_directory=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/ipyrad)
cd $input_directory/ipyrad_N$select_N

# Get vcf file

Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/denovo/ipyrad/ipyrad_vcf_filter.R $select_N $input_directory/ipyrad_N$select_N

cd denovo_N2_outfiles
mkdir -p SNAPP
model_name=(denovo_N${select_N}_filter_1SNPperRAD)

ruby /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/SNAPP/snapp_prep.rb -p $model_name.phy -t $model_name.txt \
-c $model_name.con.txt -l 1000000 -x SNAPP/$model_name.xml -o $model_name

cd SNAPP

/nobackup/tmjj24/apps/beast/bin/beast -threads $SLURM_CPUS_PER_TASK -overwrite $model_name.xml > $model_name.screen.log

/nobackup/tmjj24/apps/beast/bin/treeannotator -burnin 10 $model_name.trees $model_name.trees.Anon

