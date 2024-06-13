#!/bin/bash

#SBATCH -c 24 
#SBATCH --mem=50G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:50G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 48:00:00         # time limit in format dd-hh:mm:ss

# Specify the tasks to run:
#SBATCH --array=1-2  # Create 6 tasks, numbers 1 to 6
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here

#Modules
module load r/4.2.1
module load bioinformatics
module load raxml/8.2.12"
module show raxml/8.2.12"

line_num=$(expr $SLURM_ARRAY_TASK_ID)
# line_num=$(expr 2)
echo "$line_num"

# Get library and genome names

Library_name=$(sed -n "${line_num}p" /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/library_combinations/library_name)
genome=$(sed -n "${line_num}p" /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/library_combinations/genome)

library_version=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/)

echo "Processing database $Library_name using $genome"

# Remove previous run files
rm -f $library_version/$Library_name/RAD_loci/*.phylip.*

Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/PHRAPL/PHRAPL_RAxML_run.R $Library_name $library_version

# Create PHRAPL directory
mkdir -p $library_version/$Library_name/PHRAPL
## Merge all trees into single file
cat $library_version/$Library_name/RAD_loci/RAxML_bestTree.*.tre > $library_version/$Library_name/PHRAPL/Hetaerina_all_ddRAD_titia_dg_phrapl.trees


