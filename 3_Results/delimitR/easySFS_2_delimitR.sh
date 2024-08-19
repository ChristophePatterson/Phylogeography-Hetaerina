#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=30G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:30G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 70:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --array=3,4   # Create 32 tasks, numbers 1 to 32
#SBATCH --output=slurm-%x.%j.out

## 3,4,5,6

module load r/4.2.1
module load bioinformatics
module load bcftools/1.15

# Activate  conda enviroment
echo "Activating conda enviroment"
source /home/${USER}/.bashrc
conda activate easySFS

# Input VCF file
# Get library name
line_num=$(expr $SLURM_ARRAY_TASK_ID)
echo "$line_num"
library=$(sed -n "${line_num}p" /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/library_combinations/library_name)
dir_path=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/$library/)
species=$(sed -n "${line_num}p" /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/library_combinations/species)
cd $dir_path
# Make output for SFS
mkdir -p $dir_path/SFS

## Use the VCF that contains only biallic SNPs and excludes samples from CUAJ created within the custom R script 7_SNP_filtering_vcfR.R
VCF=($library.all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.biSNP0_20.noX.noCUAJ)
pop_file=(/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/LEA/LEA_pop_assign/popfile_${species}_noCUAJ.txt)

cp $VCF.vcf.gz $dir_path/SFS/$VCF.vcf.gz
## Subset vcf to only contain samples that have LEA assignment
## Created three distinct population files based on LEA snmf ancestory estimate
## Estimate proporstion of missing data
easySFS.py -i $dir_path/SFS/$VCF.vcf.gz -p $pop_file -a -f --preview > $dir_path/SFS/SFS_projection_preview.txt

# Assign values to sp1, sp2, and sp3 based on the value of library
case "$library" in
    "Hetaerina_americana_ddRAD_titia_dg")
        sp_0=10
        sp_1=11
        sp_2=18
        ;;
    "Hetaerina_titia_ddRAD_titia_dg")
        sp_0=26
        sp_1=50
        sp_2=20
        ;;
    "Hetaerina_americana_ddRAD_americana_dg")
        sp_0=10
        sp_1=12
        sp_2=18
        ;;
    "Hetaerina_titia_ddRAD_americana_dg")
        sp_0=20
        sp_1=25
        sp_2=20
        ;;
    *)  # Default case if library doesn't match any expected value
        echo "Unknown library"
        exit 1
        ;;
esac

# Output the values of sp1, sp2, and sp3
echo "sp1: $sp_0"
echo "sp2: $sp_1"
echo "sp3: $sp_2"

mkdir -p $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}

# Reduced testing size for quicker analysis
## sp_0=12
## sp_1=24
## sp_2=12
#
#
## 86,174,82 allows for maxium number of segregating sites within each population Pop order is SAtl, NAtl, Pac
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
output_dir=(delimitR/${sp_0}_${sp_1}_${sp_2}_nreps10000_MIG3_divWgene_sbatch)
mkdir -p $output_dir

echo "Saving output to"
echo $output_dir

ls $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}/fastsimcoal2/
#
#mkdir -p $dir_path/$output_dir
cp $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}/fastsimcoal2/*.obs $dir_path/$output_dir/
cp $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}/snp_num.txt $dir_path/$output_dir/snp_num.txt
cp $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}/all_traits.txt $dir_path/$output_dir/all_traits.txt
#
## Run delimR - creating a range of fastsimcoal2 config files and then executes them
Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/delimitR/delimR_run.R $sp_0 $sp_1 $sp_2 $dir_path $VCF $output_dir $SLURM_CPUS_PER_TASK $library
#
cd ~