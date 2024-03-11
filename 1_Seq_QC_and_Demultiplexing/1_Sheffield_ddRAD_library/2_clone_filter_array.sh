#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=150G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:150G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 24:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

# Specify the tasks to run:
#SBATCH --array=0-1   # Create 32 tasks, numbers 1 to 32

# Each separate task can be identified based on the SLURM_ARRAY_TASK_ID
# environment variable:


dir_path=(/nobackup/tmjj24/ddRAD/Sheffield_sequence_processing/trimmed_adapt_sw4-20_minlen15)
OUT_DIR=(/nobackup/tmjj24/ddRAD/Sheffield_sequence_processing/clone_filtered)

mkdir -p $OUT_DIR

cd $dir_path

FILES_1=(*.paired.*R1_001.fastq.gz)
FILES_2=(*.paired.*R2_001.fastq.gz)

FILE_1=${FILES_1[$SLURM_ARRAY_TASK_ID]}
FILE_2=${FILES_2[$SLURM_ARRAY_TASK_ID]}

echo "This is job $SLURM_ARRAY_TASK_ID"

echo "This job should use $FILE_1 and $FILE_2"

# Load modules
module load bioinformatics
module load stacks/2.60

clone_filter \
-1 $FILE_1 \
-2 $FILE_2 \
-i gzfastq \
-o $OUT_DIR \
--null_inline \
--oligo_len_2 7 

cd ~

