#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=50G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:50G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 48:00:00         # time limit in format dd-hh:mm:ss

# Specify the tasks to run:
#SBATCH --array=1-6   # Create 6 tasks, numbers 1 to 6
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here

module purge
module load bioinformatics
module load samtools
module load bcftools
module load plink
module load python/3.9.9
module load r/4.2.1


################################
#### Varibles and set up #######
################################

#Name of output SNP library

line_num=$(expr $SLURM_ARRAY_TASK_ID)

echo "$line_num"

# Get library and genome names

Library_name=$(sed -n "${line_num}p" /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/library_combinations/library_name)
genome=$(sed -n "${line_num}p" /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/library_combinations/genome)

echo "Processing database $Library_name using $genome"

#Output directory
output_stats_dir=(/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/library_combinations)
vcf_dir=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_v2/$Library_name)

#Make output directories
mkdir -p $vcf_dir


#############################
#### SNP Filtering Stats ####
#############################

#Name of output SNP library
BCF_FILE=(${Library_name}.all)

echo "$BCF_FILE processing"

cd $vcf_dir

# All Varients
echo '0. All varients'
bcftools view -H $BCF_FILE.bcf | grep -v -c '^#' > $output_stats_dir/$Library_name.vcf.stats.txt

echo "Number of samples"
bcftools query -l $BCF_FILE.bcf | grep -v -c '^#' >> $output_stats_dir/$Library_name.vcf.stats.txt

# All snps
echo '1. All snps'
bcftools view -H $BCF_FILE.snps.bcf | grep -v -c '^#' >> $output_stats_dir/$Library_name.vcf.stats.txt

# SNPs genotyped in more than 

# Depth greater than 10, average depth for each sample greater than 10 and less than 200, quality score greater than 60

echo '2. Depth greater than 10, average depth for each sample greater than 10 and less than 200, quality score greater than 60'

bcftools view -H $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.bcf | grep -v -c '^#' >> $output_stats_dir/$Library_name.vcf.stats.txt

# Removes SNPs that are not present in more than 80% samples
echo "4. Removing SNPs that arn't genotyped in more than 80% samples"
bcftools view -H $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.bcf | grep -v -c '^#' >> $output_stats_dir/$Library_name.vcf.stats.txt

# Removing SNPs with a minor allele freq less than 0.05
echo "5. Removing alleles that have a minor allele count of less that 2"
bcftools view -H $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.bcf | grep -v -c '^#' >> $output_stats_dir/$Library_name.vcf.stats.txt

##### Use output from above to determine appropreate filtering step

echo '6. SNPS randomly thinned to one per 1000 bases'
bcftools view -H $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.bcf | grep -v -c '^#' >> $output_stats_dir/$Library_name.vcf.stats.txt

echo "$BCF_FILE Complete"

cd ~





