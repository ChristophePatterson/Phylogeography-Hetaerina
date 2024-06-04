#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=5G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:5G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 01:00:00         # time limit in format dd-hh:mm:ss
# Specify the tasks to run:
#SBATCH --output=/nobackup/tmjj24/ddRAD/CUAJ_sequence_processing/excess_slurm_out/slurm-%x.%j.out

module load bioinformatics
module load samtools
module load r/4.1.2

## Merge all bamstats files
cat /nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/bwa_stats_array_SDC/bwa_HetTit1.0_dg/*.bamstats > /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/bamstats/Allsamples_HetTit1.0_dg.bamstats
cat /nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/bwa_stats_array_SDC/bwa_HetAmer1.0_dg/*.bamstats > /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/bamstats/Allsamples_HetAmer1.0_dg.bamstats

Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/5_2_bamFile_stats_and_list_creation.R