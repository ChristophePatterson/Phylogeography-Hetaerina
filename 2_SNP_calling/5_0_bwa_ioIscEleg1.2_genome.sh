#!/bin/bash

#SBATCH -c 5 
#SBATCH --mem=10G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:10G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 12:00:00         # time limit in format dd-hh:mm:ss

# Specify the tasks to run:
#SBATCH --array=0-387   # Create 382 tasks, numbers 0 to 381
#SBATCH --output=slurm_out/slurm-%x.%j.out

module purge
module load bioinformatics
module load bwa
module load samtools

########################
#### BEFORE RUNNING ####
########################

# Transfer all fastq files into a single directory and rm or "rem" files

# Find all *fq.gz files and transfer to single directory*
# find /nobackup/tmjj24/ddRAD/Durham_sequence_processing/demultiplexed_filtered_clones_trimmed_adapt_sw4-20_minlen15_all -name "*.fq.gz" -exec cp {} /nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/demultiplex_all \;
# find /nobackup/tmjj24/ddRAD/CUAJ_sequence_processing/demultiplexed -name "*.fq.gz" -exec cp {} /nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/demultiplex_all \;


# Remove all *.rem.* files and dig2con files
# rm /nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/demultiplex_all/*.rem.*.fq.gz
# rm /nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/demultiplex_all/DIG2CON*

########################

# Output and input directory

input_directory=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/demultiplex_all)
dir_output=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/bwa_ioIscEleg1.2_dg)

# Draft genome to use
genome=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/draft_genomes/GCA_921293095.2_ioIscEleg1.2_genomic.fna)

# index genome only needs doing once
# samtools faidx $genome
# bwa index $genome

mkdir -p $dir_output

# Creates vector of adapter barcodes

# List of files
FILES_1=($input_directory/*.1.fq.gz)
FILES_2=($input_directory/*.2.fq.gz)

# Extracts specific paired files depended on which slurm job it is
FILE_1=${FILES_1[$SLURM_ARRAY_TASK_ID]}
FILE_2=${FILES_2[$SLURM_ARRAY_TASK_ID]}

FILE_1_simple=(${FILE_1##*/})
FILE_1_simple=(${FILE_1_simple%%.1.fq.gz*})

FILE_2_simple=(${FILE_2##*/})
FILE_2_simple=(${FILE_1_simple%%.1.fq.gz*})

echo "This is job $SLURM_ARRAY_TASK_ID"

echo "This job should use $FILE_1 and $FILE_2"

# Runs bwa
bwa mem -t $SLURM_CPUS_PER_TASK \
  $genome \
  $FILE_1 $FILE_2 | \
  samtools sort --threads $SLURM_CPUS_PER_TASK -o $dir_output/$FILE_1_simple.paired.bam

samtools index $dir_output/$FILE_1_simple.paired.bam

cd ~








