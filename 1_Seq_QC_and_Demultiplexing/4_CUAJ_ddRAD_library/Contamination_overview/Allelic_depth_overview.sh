#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=5G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:5G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 48:00:00         # time limit in format dd-hh:mm:ss

# Specify the tasks to run:
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here

################################
#### Varibles and set up #######
################################

# Get library and genome names
Library_name=(CUAJ_Masters_library)
genome=/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/draft_genomes/HetTit1.0.p.genome.fasta

echo "Processing database $Library_name using $genome"

#Output directory
output_dir=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/$Library_name)
plot_dir=(/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/1_Seq_QC_and_Demultiplexing/4_CUAJ_ddRAD_library/Contamination_overview/plots)

mkdir -p $plot_dir

Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/1_Seq_QC_and_Demultiplexing/4_CUAJ_ddRAD_library/Contamination_overview/Allelic_depth_overview.R \
 $Library_name $output_dir $plot_dir