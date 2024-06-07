
#SBATCH -c 1 
#SBATCH --mem=20G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:20G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 24:00:00         # time limit in format dd-hh:mm:ss

# Specify the tasks to run:
#SBATCH --array=1-6  # Create 6 tasks, numbers 1 to 6
#SBATCH --output=slurm-%x.%j.out

# Commands to execute start here

#Modules
module load r/4.2.1
module load bioinformatics
module load bcftools

#Name of output SNP library

line_num=$(expr $SLURM_ARRAY_TASK_ID)

echo "$line_num"

# Get library and genome names

Library_name=$(sed -n "${line_num}p" /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/library_combinations/library_name)
genome=$(sed -n "${line_num}p" /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/library_combinations/genome)

library_version=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/)

echo "Processing database $Library_name using $genome"

## SNP library to use
SNP_file=($Library_name.all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.biSNP0_20.noX.vcf.gz)

# Get list of position for each SNP file
Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/8_RAD_loci_calling.R $Library_name $library_version

# index  bcf
bcftools index $library_version/$Library_name/$Library_name.all.snps.bcf
#Create 
bcftools query -l $library_version/$Library_name/$Library_name.all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.bcf > $library_version/$Library_name/samples_temp.txt

mkdir $library_version/$Library_name/fasta_files

cat $library_version/$Library_name/samples_temp.txt | while read line; do
    #echo $library_version/$Library_name/fasta_files/$line.fa
    # Get sample name
    filename=$(basename "$line")
    # Remove extension
    filename=${filename:0:-11}
    echo $filename
    samtools faidx $genome -r $library_version/$Library_name/${Library_name}_RAD_posistions.txt | \
    bcftools consensus -s $line -e 'FORMAT/AD<5 & FORMAT/AD>200 & QUAL<50' -a "?" -o $library_version/$Library_name/fasta_files/${filename}.fasta $library_version/$Library_name/$Library_name.all.snps.bcf 
done 

Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/8_RAD_loci_calling_phylip.R $Library_name $library_version





