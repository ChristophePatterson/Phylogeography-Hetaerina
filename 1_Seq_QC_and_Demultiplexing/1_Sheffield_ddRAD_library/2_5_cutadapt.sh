#!/bin/bash

#SBATCH -c 15 
#SBATCH --mem=10G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:10G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 24:00:00         # time limit in format dd-hh:mm:ss

# Specify the tasks to run:
#SBATCH --array=0-1   # Create 32 tasks, numbers 1 to 32
#SBATCH --output=slurm-%x.%j.out

# Each separate task can be identified based on the SLURM_ARRAY_TASK_ID
# environment variable:

module load bioinformatics
module load cutadapt/4.1


# Output and input directory

input_directory=(/nobackup/tmjj24/ddRAD/Sheffield_sequence_processing/clone_filtered)
dir_output=(/nobackup/tmjj24/ddRAD/Sheffield_sequence_processing/clone_filtered_cutadapt_cores)


mkdir -p $dir_output

# List of files
FILES_1=($input_directory/*.1.fq.gz)
FILES_2=($input_directory/*.2.fq.gz)

# Extracts specific paired files depended on which slurm job it is
FILE_1=${FILES_1[$SLURM_ARRAY_TASK_ID]}
FILE_2=${FILES_2[$SLURM_ARRAY_TASK_ID]}

FILE_1_simple=(${FILE_1##*/})
FILE_2_simple=(${FILE_2##*/})

echo $FILE_1_simple
echo $FILE_2_simple

cutadapt -a GAATTC -G GAATTC -o $dir_output/cutadapt.$FILE_1_simple -p $dir_output/cutadapt.$FILE_2_simple --action=retain -j 0 $FILE_1 $FILE_2

cd ~







