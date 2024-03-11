#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=30G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:30G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 12:00:00         # time limit in format dd-hh:mm:ss

# Specify the tasks to run:
#SBATCH --array=0-1   # Create 32 tasks, numbers 1 to 32
#SBATCH --output=/nobackup/tmjj24/ddRAD/CUAJ_sequence_processing/slurm-%x.%j.out

# Each separate task can be identified based on the SLURM_ARRAY_TASK_ID
# environment variable:

module load bioinformatics
module load stacks/2.60

# Output and input directory

input_directory=(/nobackup/tmjj24/ddRAD/CUAJ_sequence_processing/filtered_clones_trimmed_adapt_sw4-20_minlen15/)
dir_output=(/nobackup/tmjj24/ddRAD/CUAJ_sequence_processing/demultiplexed/demultiplexed_filtered_clones_trimmed_adapt_sw4-20_minlen15_$SLURM_ARRAY_TASK_ID/)

mkdir -p $dir_output

# Creates vector of adapter barcodes

BARCODES=(/home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/ddRAD_CUAJ/barcodes/*)

# List of files
FILES_1=($input_directory/*1.1.fq.gz)
FILES_2=($input_directory/*2.2.fq.gz)

# Extracts specific paired files depended on which slurm job it is
FILE_1=${FILES_1[$SLURM_ARRAY_TASK_ID]}
FILE_2=${FILES_2[$SLURM_ARRAY_TASK_ID]}
BARCODE_1=${BARCODES[$SLURM_ARRAY_TASK_ID]}

echo "This is job $SLURM_ARRAY_TASK_ID"

echo "This job should use $FILE_1 and $FILE_2 and barcodes $BARCODE_1"

# Load modules
module load bioinformatics
module load stacks/2.60

# Runs process_radtags
process_radtags \
-1 $FILE_1 \
-2 $FILE_2 \
--inline-inline \
--renz_1 pstI --renz_2 ecoRI \
-b $BARCODE_1 \
-o $dir_output \
-c -q -r 

cd ~









