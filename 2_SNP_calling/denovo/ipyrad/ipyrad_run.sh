#!/bin/bash

#SBATCH -c 24 
#SBATCH --mem=200G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:200G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 48:00:00         # time limit in format dd-hh:mm:ss

# Specify the tasks to run:
#SBATCH --output=slurm-%x.%j.out

# Demultiplexed file locations
input_directory=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/demultiplex_all)
dir_output=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/ipyrad)

# Number of samples to use per population
select_N=("20")
cd $dir_output/ipyrad_N$select_N

# Load ipyrad conda enviroment
echo "Activating conda enviroment"
source /home/${USER}/.bashrc
conda activate ipyrad

# MAKE SURE PARAM FILE HAS BEEN FULLING CONFIG

ipyrad -p params-denovo_N$select_N.txt -s 123 -c 20 -d -f

ipyrad -p params-denovo_N$select_N.txt -s 4567 -c 20 -d

conda deactivate


