#!/bin/bash

#SBATCH -c 24 
#SBATCH --mem=10G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:10G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 24:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here
module purge
module load gcc/native
module load gdal
module load geos/3.10.1
module load r/4.2.1
module load bioinformatics
module load paup/4.0a

#Output directory
## Input directory
input_dir=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/VCF_chrom_r10000)
phy_file=(WGS_titia_chr1-12)
cd $input_dir
output_dir=($input_dir/results/SVDQuartets)

mkdir -p $output_dir

Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/WGS/vcf_2_phy.R $input_dir $output_dir $phy_file

cd $output_dir

# Create paup input file
echo "log file=${output_dir}/${phy_file}_T20.log;" > $output_dir/Paup_excute_file.nex
echo "cd ${output_dir};" >> $output_dir/Paup_excute_file.nex
echo "set hashcomment=y;" >> $output_dir/Paup_excute_file.nex
echo "execute ${output_dir}/${phy_file}.nex;" >> $output_dir/Paup_excute_file.nex
echo "" >> $output_dir/Paup_excute_file.nex

echo "svdQuartets bootstrap=standard nthreads=${SLURM_CPUS_PER_TASK};" >> $output_dir/Paup_excute_file.nex
echo "saveTrees file = ${output_dir}/${phy_file}_T${SLURM_CPUS_PER_TASK}.tre supportValues = nodeLabels;" >> $output_dir/Paup_excute_file.nex
echo "quit;" >> $output_dir/Paup_excute_file.nex

# Run SVDquartets
# Run SVDquartet
paup $output_dir/Paup_excute_file.nex -n

cd ~








