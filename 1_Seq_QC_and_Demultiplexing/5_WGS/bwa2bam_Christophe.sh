#!/bin/bash

#SBATCH -c 64
#SBATCH --mem=199G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:199G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 69:00:00         # time limit in format dd-hh:mm:ss

# Specify the tasks to run:
#SBATCH --array=0-5   # Create 4 tasks
#SBATCH --output=/nobackup/tmjj24/ddRAD/WGS_titia/slurm-%x.%j.out

module purge
module load bioinformatics
module load bwa
module load samtools

########################
#### BEFORE RUNNING ####
########################

########################

# Output and input directory

## Merging sample RLPEa01 which contains sequence data from multiple lanes
#cat /nobackup/tmjj24/ddRAD/WGS_titia/data/Christophe/X204SC24013540-Z01-F001/01.RawData/RLPEb01_R/RLPEb01_R_EKDN240013190-1A_223LLLLT4_L4_1.fq.gz /nobackup/tmjj24/ddRAD/WGS_titia/data/Christophe/X204SC24013540-Z01-F001/01.RawData/RLPEb01_R/RLPEb01_R_EKDN240013190-1A_223FHVLT4_L3_1.fq.gz > /nobackup/tmjj24/ddRAD/WGS_titia/data/Christophe/all-data/RLPEb01_R_EKDN240013190-1A_223LLLLT4_L4_1.fq.gz
#cat /nobackup/tmjj24/ddRAD/WGS_titia/data/Christophe/X204SC24013540-Z01-F001/01.RawData/RLPEb01_R/RLPEb01_R_EKDN240013190-1A_223LLLLT4_L4_2.fq.gz /nobackup/tmjj24/ddRAD/WGS_titia/data/Christophe/X204SC24013540-Z01-F001/01.RawData/RLPEb01_R/RLPEb01_R_EKDN240013190-1A_223FHVLT4_L3_2.fq.gz > /nobackup/tmjj24/ddRAD/WGS_titia/data/Christophe/all-data/RLPEb01_R_EKDN240013190-1A_223LLLLT4_L4_2.fq.gz


input_directory=(/nobackup/tmjj24/ddRAD/WGS_titia/data/Christophe/all-data)
dir_output=(/nobackup/tmjj24/ddRAD/WGS_titia/bam/Christophe)

## Find all fq files and move to specific directory (done prior to running)
## find /nobackup/tmjj24/ddRAD/WGS_titia/data/Christophe -name "*.fq.gz" -exec cp {} /nobackup/tmjj24/ddRAD/WGS_titia/data/Christophe/all-data \;

# Draft genome to use
genome=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/draft_genomes/HetTit1.0.p.genome.fasta)

# index genome only needs doing once
# samtools faidx $genome

mkdir -p $dir_output

# Creates vector of adapter barcodes

# List of files
FILES_1=($input_directory/*1.fq.gz)
FILES_2=($input_directory/*2.fq.gz)

# Extracts specific paired files depended on which slurm job it is
FILE_1=${FILES_1[$SLURM_ARRAY_TASK_ID]}
FILE_2=${FILES_2[$SLURM_ARRAY_TASK_ID]}

FILE_1_simple=(${FILE_1##*/})
FILE_1_simple=$(echo "$FILE_1_simple" | sed 's/_EKDN.*//')

echo "This is job $SLURM_ARRAY_TASK_ID"

echo "This job should use $FILE_1 and $FILE_2"

echo "Samples will be saved as $dir_output/$FILE_1_simple.bam"

echo "Output will be temporary saved to $TMPDIR"

# Runs bwa
bwa mem -t 64 \
  $genome \
  $FILE_1 $FILE_2 > $TMPDIR/$FILE_1_simple.sam
  
samtools sort --threads 64 -o $dir_output/$FILE_1_simple.bam $TMPDIR/$FILE_1_simple.sam

samtools index $dir_output/$FILE_1_simple.bam

cd ~
