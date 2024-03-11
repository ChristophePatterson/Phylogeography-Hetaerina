#!/bin/bash

#SBATCH -c 12 
#SBATCH --mem=10G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:10G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 06:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

module load bioinformatics
module load trimmomatic/0.39

dir_input=(/nobackup/tmjj24/ddRAD/RAW_Durham_data/All_fastq_files)
dir_output=(/nobackup/tmjj24/ddRAD/Durham_sequence_processing/trimmed)

mkdir -p $dir_output

input_file_1=(i5_A6_i7_A6_1_EKDL230000697-1A_HNHC3DSX5_L3_1.fq.gz)
input_file_2=(i5_A6_i7_A6_1_EKDL230000697-1A_HNHC3DSX5_L3_2.fq.gz)

# PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...


trimmomatic PE -phred33 -threads 12 \
	$dir_input/$input_file_1 $dir_input/$input_file_2 \
	$dir_output/trimmed.paired.$input_file_1 $dir_output/trimmed.paired.$input_file_2 \
	$dir_output/trimmed.unpaired.$input_file_1 $dir_output/trimmed.unpaired.$input_file_2 \
	ILLUMINACLIP:/home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/ddRAD_Durham/adapters/adapters.fasta:2:30:10