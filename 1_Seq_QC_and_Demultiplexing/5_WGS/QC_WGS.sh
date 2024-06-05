#!/bin/bash

#SBATCH -c 64 
#SBATCH --mem=50G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:50G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 06:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

module load bioinformatics
module load fastqc/0.11.9
module load multiqc/1.15 

dir_output_raw=(/nobackup/tmjj24/ddRAD/WGS_titia/data/Christophe/QC)
dir_raw=(/nobackup/tmjj24/ddRAD/WGS_titia/data/Christophe/all-data)

mkdir -p $dir_output_raw

# Run fastqc on all raw fastq files

fastqc -o $dir_output_raw -t 64 $dir_raw/*.fq.gz

# Combine reports into sinlge multiqc report
multiqc -o $dir_output_raw $dir_output_raw/*fastqc.zip --interactive

