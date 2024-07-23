#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=10G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:10G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 12:00:00         # time limit in format dd-hh:mm:ss

# Specify the tasks to run:
#SBATCH --output=slurm-%x.%j.out

module load r/4.1.2

# Demultiplexed file locations
input_directory=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/demultiplex_all)
dir_output=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/ipyrad)

# Number of samples to use per population
select_N=("3")

mkdir -p $dir_output
mkdir -p $dir_output/ipyrad_N$select_N
mkdir -p $dir_output/ipyrad_N$select_N/demultiplexed

Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/denovo/ipyrad/ipyrad_sample_N.R $input_directory $dir_output $select_N

# Copy over raw sequences files and rename into ipyrad paired end convension
while read line; do
    filename=$line
    echo $filename
    cp $input_directory/${filename}.1.fq.gz $dir_output/ipyrad_N$select_N/demultiplexed/${filename}_R1_.fq.gz
    cp $input_directory/${filename}.2.fq.gz $dir_output/ipyrad_N$select_N/demultiplexed/${filename}_R2_.fq.gz
done < $dir_output/samples_$select_N.txt

cd $dir_output/ipyrad_N$select_N

# Load ipyrad conda enviroment
echo "Activating conda enviroment"
source /home/${USER}/.bashrc
conda activate ipyrad

# Create new ipyrad dataset
## NEEDS TO BE THEN CUSTOM SET
ipyrad -n denovo_N$select_N

conda deactivate