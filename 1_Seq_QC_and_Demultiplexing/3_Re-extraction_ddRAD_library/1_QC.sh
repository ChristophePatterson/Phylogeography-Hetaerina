#!/bin/bash

#SBATCH -c 32 
#SBATCH --mem=50G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:50G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 06:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

module load bioinformatics
module load fastqc
module load multiqc/1.14 

dir_output=(/nobackup/tmjj24/ddRAD/Pool_X_sequence_processing/QC)
dir_output_raw=($dir_output/raw_multiqc_Pool_X)
dir_output_trim=($dir_output/trim_multiqc_Pool_X)
dir_output_clone=($dir_output/clone_multiqc_Pool_X)
dir_output_demulti=($dir_output/demulti_multiqc_Pool_X)

dir_raw=(/nobackup/tmjj24/ddRAD/RAW_Pool_X_data/All_fastq_files)
dir_trim=(/nobackup/tmjj24/ddRAD/Pool_X_sequence_processing/trimmed_adapt_sw4-20_minlen15/)
dir_clone=(/nobackup/tmjj24/ddRAD/Pool_X_sequence_processing/filtered_clones_trimmed_adapt_sw4-20_minlen15/)
dir_demulti=(/nobackup/tmjj24/ddRAD/Pool_X_sequence_processing/demultiplexed_filtered_clones_trimmed_adapt_sw4-20_minlen15_0_all/)

mkdir -p $dir_output
mkdir -p $dir_output_raw
mkdir -p $dir_output_trim
mkdir -p $dir_output_clone
mkdir -p $dir_output_demulti

# Run fastqc on all raw fastq files

fastqc -o $dir_output_raw -t 32 $dir_raw/*.fq.gz
fastqc -o $dir_output_trim -t 32 $dir_trim/*.paired.*.fq.gz
fastqc -o $dir_output_clone -t 32 $dir_clone/*.fq.gz
fastqc -o $dir_output_demulti -t 32 $dir_demulti/*.fq.gz

Combine reports into sinlge multiqc report
multiqc -o $dir_output_raw $dir_output_raw/*fastqc.zip --interactive
multiqc -o $dir_output_trim $dir_output_trim/*fastqc.zip --interactive
multiqc -o $dir_output_clone $dir_output_clone/*fastqc.zip --interactive
multiqc -o $dir_output_demulti $dir_output_demulti/*fastqc.zip --interactive
