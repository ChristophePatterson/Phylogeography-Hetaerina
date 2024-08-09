#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=50G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:50G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 24:00:00         # time limit in format dd-hh:mm:ss

#SBATCH --output=slurm-%x.%j.out


# File inputs and output lociation
input_dir=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/ipyrad/ipyrad_N3/denovo_N3_outfiles)
output_dir=($input_dir/gphocs)
mkdir -p $output_dir

# Name of loci file from ipyrad
loci_file=(denovo_N3)

# Change to input directory
cd $input_dir

## Get loci that have more than 8 samples with sequence data
grep 'locus.* 9 .*' $loci_file.gphocs > $output_dir/high_cov_loci.txt
num_loci=(`wc -l $output_dir/high_cov_loci.txt`)
# Create new subset loci file
echo $num_loci > $output_dir/$loci_file.lociHC.gphocs
# Then following lines are are the loci that match the randomly selected loci
grep -f $output_dir/high_cov_loci.txt --no-group-separator -A 10 $loci_file.gphocs >> $output_dir/$loci_file.lociHC.gphocs

## Run G-Phocs

