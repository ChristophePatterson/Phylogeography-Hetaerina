#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=30G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:30G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 70:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

module load r/4.2.1
module load bioinformatics
module load bcftools/1.15

# Activate  conda enviroment
echo "Activating conda enviroment"
source /home/${USER}/.bashrc
conda activate easySFS

# Input VCF file
# Get library name
library=(VCF_chrom_r10000)
dir_path=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/$library/)
species="titia"
cd $dir_path
# Make output for SFS
mkdir -p $dir_path/SFS

## Use the VCF that contains only biallic SNPs and excludes samples from CUAJ created within the custom R script 7_SNP_filtering_vcfR.R
VCF=(WGS_titia_chr1-12_filtered)
# Create popfile
echo -e  "CUAJa03_R\nRLPEb01_R\nZANAa07_R\nZANAa05\n" > $dir_path/SFS/samples_names.txt
echo -e  "CUAJ01\nRLPE02\nZANA01\nZANA01\n" > $dir_path/SFS/samples_sites.txt
paste $dir_path/SFS/samples_names.txt $dir_path/SFS/samples_sites.txt > $dir_path/SFS/popfile_Pacific_${species}.txt

pop_file=($dir_path/SFS/popfile_Pacific_${species}.txt)

cp $VCF.vcf.gz $dir_path/SFS/$VCF.vcf.gz
## Subset vcf to only contain samples that have LEA assignment
## Created three distinct population files based on LEA snmf ancestory estimate
## Estimate proporstion of missing data
#### easySFS.py -i $dir_path/SFS/$VCF.vcf.gz -p $pop_file -a -f --preview > $dir_path/SFS/SFS_projection_preview.txt

# Assign values to sp1, sp2, and sp3 based on the value of library

sp_0=2
sp_1=2
sp_2=4

# Output the values of sp1, sp2, and sp3
echo "sp1: $sp_0"
echo "sp2: $sp_1"
echo "sp3: $sp_2"

# Create output directory
mkdir -p $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}

## Run SFS
easySFS.py -i $dir_path/SFS/$VCF.vcf.gz -p $pop_file \
 -a -f --proj $sp_0,$sp_1,$sp_2 -o $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}
#
conda deactivate
#
### Calculate number of SNPs within the vcf file
bcftools view -H $dir_path/SFS/$VCF.vcf.gz | wc -l > $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}/snp_num.txt
#
## Create delimR input files
Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/delimitR/delimR_inputfile.R $sp_0 $sp_1 $sp_2 $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}
#
## Create directory to run delimitR
#
mkdir -p delimitR/
output_dir=(delimitR/${sp_0}_${sp_1}_${sp_2}_nreps10000_MIG1_trees_3_sbatch)
mkdir -p $output_dir
#
#mkdir -p $dir_path/$output_dir
cp $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}/fastsimcoal2/*.obs $dir_path/$output_dir/
cp $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}/snp_num.txt $dir_path/$output_dir/snp_num.txt
cp $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}/all_traits.txt $dir_path/$output_dir/all_traits.txt
#
## Run delimR - creating a range of fastsimcoal2 config files and then executes them
Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/WGS/WGS_delimitR/delimit_R_run_WGS.R $sp_0 $sp_1 $sp_2 $dir_path $VCF $output_dir $SLURM_CPUS_PER_TASK $library
#
cd ~