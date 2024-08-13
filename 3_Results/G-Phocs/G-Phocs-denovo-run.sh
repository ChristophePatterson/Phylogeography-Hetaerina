#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=50G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:50G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 24:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

# Load modules
module load r/4.2.1

# Which species to run the programme on. Can either be "titia" or "americana"
spp=("titia")

sample_N=("3")
# File inputs and output lociation
top_dir=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/ipyrad)
input_dir=($top_dir/ipyrad_N${sample_N}/denovo_N${sample_N}_${spp}_outfiles)
output_dir=($top_dir/ipyrad_N${sample_N}/ipyrad_${spp}_N${sample_N}/G-Phocs)
mkdir -p $top_dir/ipyrad_N${sample_N}/ipyrad_${spp}_N${sample_N}/
mkdir -p $output_dir

# Name of loci file from ipyrad
loci_file=(denovo_N3_${spp})

# Change to input directory
cd $input_dir

## Get loci that have more than 8 samples with sequence data
grep 'locus.* 9 .*' $loci_file.gphocs > $output_dir/high_cov_loci.txt
num_loci=(`wc -l $output_dir/high_cov_loci.txt`)
# Create new subset loci file
echo $num_loci > $output_dir/ipyrad_${spp}_N${sample_N}_N${sample_N}.gphocs
# Then following lines are are the loci that match the randomly selected loci
# Yes the output has two N3 its to keep combatlatly with with config genreation code
grep -f $output_dir/high_cov_loci.txt --no-group-separator -A 10 $loci_file.gphocs >> $output_dir/ipyrad_${spp}_N${sample_N}_N${sample_N}.gphocs
# Replace "-" with "N"
sed -i 's/-/N/g' $output_dir/ipyrad_${spp}_N${sample_N}_N${sample_N}.gphocs

## Run G-Phocs
## Copy over pop assign file
cp $top_dir/samples_${sample_N}_${spp}_assign.txt $output_dir/GPhocs_samples_N${sample_N}.txt

cd $output_dir
Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/G-Phocs/g-phocs-config-generation.R ipyrad_${spp}_N${sample_N} $top_dir/ipyrad_N${sample_N}/ $sample_N $output_dir/model_runs/


