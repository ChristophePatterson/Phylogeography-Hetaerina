#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=2G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:2G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 02:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here
module purge
module load gcc/native
module load gdal
module load geos/3.10.1
module load r/4.2.1
module load bioinformatics
module load raxml/8.2.12 

#Output directory
## Input directory
input_dir=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/VCF_chrom_r10000)
phy_file=(WGS_titia_chr1-12)
cd $input_dir
output_dir=($input_dir/results/dsuite)

mkdir -p $output_dir

# Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/WGS/vcf_2_phy.R $input_dir $output_dir $phy_file

cd $input_dir

# run dsuite across all snps
/nobackup/tmjj24/apps/Dsuite/Build/Dsuite Dquartets ${phy_file}_filtered.vcf.gz \
${phy_file}_sites.txt \
-o ${output_dir}/${phy_file}_filtered_allChr.txt 

# run dsuite for each chromosome
chr=(1)
cat ${phy_file}_chrom_start_end.txt | while read line; do
   /nobackup/tmjj24/apps/Dsuite/Build/Dsuite Dquartets ${phy_file}_filtered.vcf.gz \
   ${phy_file}_sites.txt \
   -r $line \
   -o ${output_dir}/${phy_file}_filtered_chr${chr}.txt 
   ((chr++))
done

cd ~








