#!/bin/bash

#SBATCH -c 12 
#SBATCH --mem=25G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:25G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 12:00:00         # time limit in format dd-hh:mm:ss

# Specify the tasks to run:
#SBATCH --array=0   # Create 32 tasks, numbers 1 to 32
#SBATCH --output=slurm-%x.%j.out

module load bioinformatics
module load trimmomatic/0.39

dir_input=(/nobackup/tmjj24/ddRAD/RAW_Pool_X_data/All_fastq_files)
dir_output=(/nobackup/tmjj24/ddRAD/Pool_X_sequence_processing/trimmed_adapt_sw4-20_minlen15)

mkdir -p $dir_output

cd $dir_input

FILES_1=(*1.fq.gz)
FILES_2=(*2.fq.gz)

input_file_1=${FILES_1[$SLURM_ARRAY_TASK_ID]}
input_file_2=${FILES_2[$SLURM_ARRAY_TASK_ID]}

echo "This is job $SLURM_ARRAY_TASK_ID"

echo "This job should use $input_file_1 and $input_file_2"

mkdir -p $dir_output

# PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...

trimmomatic PE -phred33 -threads 12 \
	$dir_input/$input_file_1 $dir_input/$input_file_2 \
	$dir_output/trimmed.paired.$input_file_1 $dir_output/trimmed.unpaired.$input_file_1 \
	$dir_output/trimmed.paired.$input_file_2 $dir_output/trimmed.unpaired.$input_file_2 \
	ILLUMINACLIP:/home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/ddRAD_Durham/adapters/adapters.fasta:2:30:10 \
	SLIDINGWINDOW:4:20 MINLEN:15




