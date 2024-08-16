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

library="WGS_titia"
link_filt="10000"

dir_path=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/${library}_r$link_filt/)
cd $dir_path
# Make output for SFS
mkdir -p $dir_path/SFS

## Use the VCF that contains only biallic SNPs and excludes samples from CUAJ created within the custom R script 7_SNP_filtering_vcfR.R
VCF=(${library}_chr1-12_filtered)
# Create popfile
echo -e  "CUAJa03\nRLPEb01\nZANAa07\nZANAa05\n" > $dir_path/SFS/${library}_samples_names.txt
echo -e  "CUAJ01\nRLPE02\nZANA01\nZANA01\n" > $dir_path/SFS/${library}_samples_sites.txt
# Assign values to sp1, sp2, and sp3 based on the value of library
case "$library" in
    "WGS_titia_PAC")
        echo -e  "CUAJa03\nRLPEb01\nZANAa07\nZANAa05\n" > $dir_path/SFS/${library}_samples_names.txt
        echo -e  "CUAJ01\nRLPE02\nZANA01\nZANA01\n" > $dir_path/SFS/${library}_samples_sites.txt
        sp_0=2
        sp_1=2
        sp_2=4
        ;;
    "WGS_titia")
        echo -e  "MIXTa01\nMIXTa02\nCUAJb01\nCUAJa01\nCUAJa03\nRLPEb01\nZANAa07\nZANAa05\n" > $dir_path/SFS/${library}_samples_names.txt
        echo -e  "MIXT01\nMIXT01\nCUAJ01_ATL\nCUAJ01_ATL\nCUAJ01_PAC\nRLPE02\nZANA01\nZANA01\n" > $dir_path/SFS/${library}_samples_sites.txt
        sp_0=4
        sp_1=4
        sp_2=2
        sp_3=2
        sp_4=4
        ;;
    *)  # Default case if library doesn't match any expected value
        echo "Unknown library"
        exit 1
        ;;
esac

# Convert to popfile
pop_file=($dir_path/SFS/${library}_popfile.txt)
paste $dir_path/SFS/${library}_samples_names.txt $dir_path/SFS/${library}_samples_sites.txt > $pop_file
# Copy over vcf for ease
cp $VCF.vcf.gz $dir_path/SFS/$VCF.vcf.gz
## Subset vcf to only contain samples that are wanted
bcftools view --samples-file $dir_path/SFS/${library}_samples_names.txt -O z $dir_path/SFS/$VCF.vcf.gz > $dir_path/SFS/${VCF}_temp.vcf.gz
# Remove SNPs that are no longer polymorphic
bcftools view --min-ac 2 -O z $dir_path/SFS/${VCF}_temp.vcf.gz > $dir_path/SFS/${VCF}_final.vcf.gz

# Remove temp files $dir_path/SFS/${VCF}_temp.vcf.gz

## Created three distinct population files based on LEA snmf ancestory estimate
## Estimate proporstion of missing data
### easySFS.py -i $dir_path/SFS/${VCF}_final.vcf.gz -p $pop_file -a -f --preview > $dir_path/SFS/SFS_projection_preview.txt

### # Assign values to sp1, sp2, and sp3 based on the value of library
mkdir -p $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}

case "$library" in
    "WGS_titia_PAC")
        easySFS.py -i $dir_path/SFS/${VCF}_final.vcf.gz -p $pop_file \
        -a -f --proj $sp_0,$sp_1,$sp_2 -o $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}
        echo -e "$sp_0\n$sp_1\n$sp_2" > $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}/projections.txt
        ;;
    "WGS_titia")
        easySFS.py -i $dir_path/SFS/${VCF}_final.vcf.gz -p $pop_file \
        -a -f --proj $sp_0,$sp_1,$sp_2,$sp_3,$sp_4 -o $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}
        echo -e "$sp_0\n$sp_1\n$sp_2\n$sp_3\n$sp_4" > $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}/projections.txt
        ;;
    *)  # Default case if library doesn't match any expected value
        echo "Unknown library"
        exit 1
        ;;
esac


### 
### #
### conda deactivate
### #
### ### Calculate number of SNPs within the vcf file
bcftools view -H $dir_path/SFS/${VCF}_final.vcf.gz | wc -l > $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}/snp_num.txt

### #
### ## Create delimR input files
Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/delimitR/delimR_inputfile.R $sp_0 $sp_1 $sp_2 $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}
#
## Create directory to run delimitR
#
mkdir -p delimitR/
output_dir=(delimitR/${sp_0}_${sp_1}_${sp_2}_nreps1000_MIG3_trees_3_sbatch)
mkdir -p $output_dir
#
#mkdir -p $dir_path/$output_dir
cp $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}/fastsimcoal2/*.obs $dir_path/$output_dir/
cp $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}/snp_num.txt $dir_path/$output_dir/snp_num.txt
cp $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}/all_traits.txt $dir_path/$output_dir/all_traits.txt
cp $dir_path/SFS/SFS_${sp_0}_${sp_1}_${sp_2}/projections.txt $dir_path/$output_dir/projections.txt
#
## Run delimR - creating a range of fastsimcoal2 config files and then executes them
Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/WGS/WGS_delimitR/delimit_R_run_WGS.R $sp_0 $sp_1 $sp_2 $dir_path $VCF $output_dir $SLURM_CPUS_PER_TASK $library
#
cd ~