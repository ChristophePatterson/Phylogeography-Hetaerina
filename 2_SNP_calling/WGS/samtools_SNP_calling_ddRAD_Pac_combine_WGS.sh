#!/bin/bash

#SBATCH -c 4 
#SBATCH --mem=20G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:20G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 24:00:00         # time limit in format dd-hh:mm:ss

# Specify the tasks to run:
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here

module purge
module load bioinformatics
module load samtools
module load bcftools
module load plink
module load python/3.9.9
module load r/4.2.1

# Script to run SNP calling ddRAD pacific populations and then combine that with the WGS data

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 1: SNP call pacific populations for ddRAD data
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# List of Pacific samples bam files
bamFiles=(/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/WGS/ddRAD_PAC_samples.txt)

# If genomes have been downloaded from NCBI they need to reformatted into default fasta format
# awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' in.fna > out.fna

Library_name="WGS_ddRAD_Pacific_titia"

genome=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/draft_genomes/HetTit1.0.p.genome.fasta)

echo "Processing database $Library_name using $genome"

#Output directory
output_dir=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/$Library_name)

#Make output directories
mkdir -p $output_dir
cd $output_dir

#####################
#SNP Calling ####
#####################

cd $output_dir
mkdir -p $output_dir/ddRAD_temp
cd $output_dir/ddRAD_temp

# Variant calling

# add -v to bcftools call inorder to only retain varient site

bcftools mpileup -Ou \
--max-depth 10000 -q 20 -Q 20 \
 -P ILLUMINA --annotate FORMAT/DP,FORMAT/AD \
-f $genome \
-b $bamFiles | \
bcftools call -v -m -P 1e-6 -f GQ -G - \
-O b -o $Library_name.all.bcf


#######################
#######SNP Filtering ####
#######################

#Name of output SNP library
BCF_FILE=(${Library_name}.all)

echo "$BCF_FILE processing"

cd $output_dir/ddRAD_temp

# Indexs file, Im not sure if this is nessacary tbh
bcftools index $BCF_FILE.bcf
bcftools index -t $BCF_FILE.bcf

# All Varients
echo '0. All varients'
bcftools view -H $BCF_FILE.bcf | grep -v -c '^#'

echo "Number of samples"
bcftools query -l $BCF_FILE.bcf | grep -v -c '^#'

# Remove indels
bcftools view -V indels $BCF_FILE.bcf -O b > $BCF_FILE.snps.bcf

# All snps
echo '1. All single nucleotides'
bcftools view -H $BCF_FILE.snps.bcf | grep -v -c '^#'

# SNPs genotyped in more than 

# Depth greater than 10, average depth for each sample greater than 10 and less than 200, quality score greater than 60

echo '2. Depth greater than 10, average depth for each sample greater than 10 and less than 200, quality score greater than 60'

bcftools filter -S . -e 'FMT/DP<10' $BCF_FILE.snps.bcf | \
bcftools view -e 'AVG(FMT/DP)<10 || AVG(FMT/DP)>200 || QUAL<60' -O b > \
$BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.bcf

bcftools view -H $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.bcf | grep -v -c '^#'

bcftools view -O z $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.bcf > $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.vcf.gz

# Removes SNPs that are not present in more than 80% samples
echo "4. Removing SNPs that arn't genotyped in more than 80% samples"
bcftools view -e 'AN/2<N_SAMPLES*0.8' -O b $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.bcf > $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.bcf
bcftools view -H $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.bcf | grep -v -c '^#'

# Removing SNPs with a minor allele freq less than 0.05
echo "5. Removing alleles that have a minor allele count of less that 2"
bcftools view --min-ac 2 -O b $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.bcf > $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.bcf
bcftools view -H $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.bcf | grep -v -c '^#'

bcftools view -O z $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.bcf > $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.vcf.gz


##########################
##### LD calculation #####
##########################
# Requires SNP library to have been created

# Go to new directory
cd $output_dir/ddRAD_temp

mkdir ld_decay
# move into it
cd ld_decay

# bcftools view -O z $output_dir/$BCF_FILE.snps.bcf > $output_dir/$BCF_FILE.snps.vcf.gz

# calc ld with plink
plink --vcf $output_dir/$BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.01 --geno 0.1 --mind 0.5 \
-r2 gz --ld-window 1000 --ld-window-kb 1000 \
--ld-window-r2 0 \
--make-bed --out $BCF_FILE.LD

# Calculate the cororlation over distance 
/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/6_LD_dist_calc_python3.py -i $BCF_FILE.LD.ld.gz -o $BCF_FILE.LD

# Plot results using R

Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/6_LD_dist_calc_plot.R

cd $output_dir/ddRAD_temp

#convert to Plink

# Extracts and edits sample names
bcftools query -l $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.vcf.gz > sample_names.txt
# Removes last 12 characters from samples
cat sample_names.txt | rev | cut -c12- | rev > sample_names_short.txt

while read line; do
    filename=$(basename "$line")
    echo "$filename"
done < sample_names_short.txt > sample_names_short_base.txt

sed 's/_X//' sample_names_short_base.txt > sample_names_short_base_minus_X.txt

# Reheaders vcf
bcftools reheader -s sample_names_short_base_minus_X.txt $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.vcf.gz > $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.reheader.vcf.gz
mkdir -p plink
plink -vcf $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.reheader.vcf.gz --allow-extra-chr -recode --out plink/$BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2
plink --file plink/$BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2 --chr-set 75 --allow-extra-chr --make-bed --recode A --out plink/$BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2


echo "$BCF_FILE Complete"

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 2: Extract position of SNPs
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
cd $output_dir/ddRAD_temp

SNP_POS=(${BCF_FILE}_SNP_position)
zcat $BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.vcf.gz | bioawk -c vcf '{print $chrom "\t" $pos}' > $SNP_POS.txt
awk '{gsub("_", "\.");print}' $SNP_POS.txt > $output_dir/$SNP_POS.dot.txt
tr -s '\n'  '\t'< $SNP_POS.txt > $output_dir/$SNP_POS.tab.txt

# Move back into main directory
cd $output_dir/
mkdir -p WGS_temp
cd $output_dir/WGS_temp

# Pacific samples with WGS data
bamFiles_Pac=(/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/WGS/WGS_titia_PAC_bam_list.txt)
# Conduct SNP calling on WGS files
bcftools mpileup -Ou \
--max-depth 10000 -q 20 -Q 20 \
 -P ILLUMINA -a FORMAT/DP,FORMAT/AD \
 -R $output_dir/$SNP_POS.dot.txt \
-f $genome \
-b $bamFiles_Pac | \
bcftools call -m -P 1e-6 -f GQ -G - \
-O b -o $Library_name.WGS.bcf

# Index
bcftools index $output_dir/WGS_temp/$Library_name.WGS.bcf
bcftools index $output_dir/ddRAD_temp/$BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.bcf
# Convert to vcf.gz
bcftools view -O z $Library_name.WGS.bcf > $Library_name.WGS.vcf.gz

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 3: Combine into single vcf
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
cd $output_dir

bcftools merge WGS_temp/$Library_name.WGS.bcf ddRAD_temp/$BCF_FILE.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.bcf -O z > $output_dir/$Library_name.all.vcf.gz

# Remove non SNP variets
bcftools view -V indels $output_dir/$Library_name.all.vcf.gz -O z > $output_dir/$Library_name.vcf.gz

# Limit to one every 10000 RAD base pairs
bcftools +prune -n 1 -N rand -w 10000bp $output_dir/$Library_name.vcf.gz -O z -o $output_dir/$Library_name.rand10000.vcf.gz

## Filter and convert input files#
mkdir -p $output_dir/results/SVDQuartets

Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/WGS/WGS_ddRAD_filter.R $Library_name $output_dir $SLURM_CPUS_PER_TASK

mkdir -p $output_dir/results/RAxML
cd $output_dir/results/RAxML
raxmlHPC-PTHREADS -s $Library_name.phy -f a -m ASC_GTRGAMMA --asc-corr=lewis -n $phy_file.raxmlASC_GTRGAM_HPCaT15 -p 12345 -x 12345 -N 100 -T $SLURM_CPUS_PER_TASK
