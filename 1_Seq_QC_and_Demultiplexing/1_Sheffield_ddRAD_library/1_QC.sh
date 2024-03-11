#!/bin/bash

#SBATCH -c 32 
#SBATCH --mem=50G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:50G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 06:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

module load bioinformatics
module load fastqc
module load multiqc/1.14 

dir_output=(/nobackup/tmjj24/ddRAD/Sheffield_sequence_processing/QC)
dir_output_raw=($dir_output/raw_multiqc_sheffield)
dir_output_trim=($dir_output/trim_multiqc_sheffield)
dir_output_clone=($dir_output/clone_multiqc_sheffield)
dir_output_cutadp=($dir_output/cutadp_multiqc_sheffield)
dir_output_demulti_shef=($dir_output/demultiplexed_shef_all)

dir_raw=(/nobackup/tmjj24/ddRAD/RAW_sheffield_data/All_Sheffield_seq)
dir_trim=(/nobackup/tmjj24/ddRAD/Sheffield_sequence_processing/trimmed_adapt_sw4-20_minlen15)
dir_clone=(/nobackup/tmjj24/ddRAD/Sheffield_sequence_processing/clone_filtered)
dir_cutadp=(/nobackup/tmjj24/ddRAD/Sheffield_sequence_processing/clone_filtered_cutadapt_cores)
dir_demulti_shef=(/nobackup/tmjj24/ddRAD/Sheffield_sequence_processing/demultiplexed_shef_all)

mkdir -p $dir_output
mkdir -p $dir_output_raw
mkdir -p $dir_output_trim
mkdir -p $dir_output_clone
mkdir -p $dir_output_cutadp
mkdir -p $dir_output_demulti_shef

# Run fastqc on all raw fastq files

# fastqc -o $dir_output_raw -t 32 $dir_raw/*.fastq.gz
# fastqc -o $dir_output_trim -t 32 $dir_trim/*.paired.*.fastq.gz
# fastqc -o $dir_output_clone -t 32 $dir_clone/*.fq.gz
# fastqc -o $dir_output_cutadp -t 32 $dir_cutadp/*.fq.gz
fastqc -o $dir_output_demulti_shef -t 32 $dir_demulti_shef/*.fq.gz

Combine reports into sinlge multiqc report
# multiqc -o $dir_output_raw $dir_output_raw/*fastqc.zip --interactive
# multiqc -o $dir_output_trim $dir_output_trim/*fastqc.zip --interactive
# multiqc -o $dir_output_clone $dir_output_clone/*fastqc.zip --interactive
# multiqc -o $dir_output_cutadp $dir_output_cutadp/cutadapt.*fastqc.zip --interactive
multiqc -o $dir_output_demulti_shef $dir_output_demulti_shef/*fastqc.zip --ignore "*rem*" --interactive

