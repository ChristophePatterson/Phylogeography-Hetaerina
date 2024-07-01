#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=50G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:50G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 48:00:00         # time limit in format dd-hh:mm:ss

# Specify the tasks to run:
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here

#Modules
module load r/4.2.1
module load bioinformatics
module load bcftools

#Name of output SNP library

# Get library and genome names
Library_name=(VCF_chrom_r10000)
genome=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/draft_genomes/HetTit1.0.p.genome.fasta)

library_version=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/)

echo "Processing database $Library_name using $genome"

## SNP library to use
SNP_file=(WGS_titia_chr1-12)

bcftools query -l $library_version/$Library_name/$SNP_file > $library_version/$Library_name/samples_temp.txt

# Get list of position for each SNP file
Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/WGS/8_WGS_loci_calling.R $Library_name $library_version $SNP_file

# index bcf
bcftools view -O b -o $library_version/$Library_name/$SNP_file.bcf $library_version/$Library_name/$SNP_file.vcf
bcftools index $library_version/$Library_name/$SNP_file.bcf

rm -r $library_version/$Library_name/fasta_files
mkdir $library_version/$Library_name/fasta_files

# Create temporary bcf with all calls
bcftools mpileup -Ou \
    --max-depth 10000 -q 20 -Q 20 \
    -P ILLUMINA --annotate FORMAT/DP,FORMAT/AD \
    -R $library_version/$Library_name/${Library_name}_RAD_posistions_delim.txt \
    -f $genome \
    -b $library_version/$Library_name/samples_temp.txt | \
    bcftools call -m -P 1e-6 -f GQ -G - \
    -O b -o $TMPDIR/$SNP_file.all.bcf



bcftools view -V indels $TMPDIR/$SNP_file.all.bcf -O b > $TMPDIR/$SNP_file.snps.bcf
bcftools index $TMPDIR/$SNP_file.snps.bcf


cat $library_version/$Library_name/samples_temp.txt | while read line; do
    #echo $library_version/$Library_name/fasta_files/$line.fa
    # Get sample name
    filename=$(basename "$line")
    # Remove extension
    filename=${filename:0:-4}
    echo $filename
    samtools faidx $genome -r $library_version/$Library_name/${Library_name}_RAD_posistions.txt | \
    bcftools consensus -s $line -e 'FORMAT/AD<5 & FORMAT/AD>200 & QUAL<50' -a "?" --iupac-codes -o $library_version/$Library_name/fasta_files/${filename}.fasta $TMPDIR/$SNP_file.snps.bcf 
done 

rm -r $library_version/$Library_name/RAD_loci
mkdir -p $library_version/$Library_name/RAD_loci

# Convert to phylip
# Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/8_RAD_loci_calling_phylip.R $Library_name $library_version





