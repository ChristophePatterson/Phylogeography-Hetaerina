#!/bin/bash

#SBATCH -c 1 
#SBATCH --mem=5G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:5G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 01:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --array=1-407   # Create 32 tasks, numbers 1 to 32
# Specify the tasks to run:
#SBATCH --output=/nobackup/tmjj24/ddRAD/CUAJ_sequence_processing/excess_slurm_out/slurm-%x.%j.out

module load bioinformatics
module load samtools
module load r/4.1.2

## Copy over new CUAJ files to correct folder
# cp /nobackup/tmjj24/ddRAD/CUAJ_sequence_processing/demultiplexed/bwa_HetTit1.0_dg/* /nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/bwa_HetTit1.0_dg
# cp /nobackup/tmjj24/ddRAD/CUAJ_sequence_processing/demultiplexed/bwa_HetAmer1.0_dg/* /nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/bwa_HetAmer1.0_dg

# Varibles

bwa_dir=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing)
bwa_stats_dir=($bwa_dir/bwa_stats_array_SDC)
americana_dir=(bwa_HetAmer1.0_dg)
titia_dir=(bwa_HetTit1.0_dg)

line_num=$(expr $SLURM_ARRAY_TASK_ID)

FILE=$(sed -n "${line_num}p" /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/ddRAD_Durham_and_Sheffield_CUAJ/library_combinations/samples)

# Making directories to output stats
mkdir -p $bwa_stats_dir/$titia_dir
mkdir -p $bwa_stats_dir/$americana_dir

samtools depth $bwa_dir/$americana_dir/$FILE.paired.bam | awk '{sum+=$3; sumsq+=$3*$3} END { printf (sum/NR"\t"sqrt(sumsq/NR - (sum/NR)**2)"\t"NR"\n")}' > $bwa_stats_dir/$americana_dir/$FILE.sample_cov.txt
samtools view -c -F 260 $bwa_dir/$americana_dir/$FILE.paired.bam > $bwa_stats_dir/$americana_dir/$FILE.num.mapped.txt
samtools flagstat $bwa_dir/$americana_dir/$FILE.paired.bam | awk -F "[(|%]" 'NR == 5 {print $2}' > $bwa_stats_dir/$americana_dir/$FILE.mapping_proportion.txt

samtools depth $bwa_dir/$titia_dir/$FILE.paired.bam | awk '{sum+=$3; sumsq+=$3*$3} END { printf (sum/NR"\t"sqrt(sumsq/NR - (sum/NR)**2)"\t"NR"\n")}' > $bwa_stats_dir/$titia_dir/$FILE.sample_cov.txt
samtools view -c -F 260 $bwa_dir/$titia_dir/$FILE.paired.bam > $bwa_stats_dir/$titia_dir/$FILE.num.mapped.txt
samtools flagstat $bwa_dir/$titia_dir/$FILE.paired.bam | awk -F "[(|%]" 'NR == 5 {print $2}' > $bwa_stats_dir/$titia_dir/$FILE.mapping_proportion.txt

echo $FILE > $bwa_stats_dir/$titia_dir/$FILE.samplename.txt
echo $FILE > $bwa_stats_dir/$americana_dir/$FILE.samplename.txt

paste $bwa_stats_dir/$titia_dir/$FILE.samplename.txt $bwa_stats_dir/$titia_dir/$FILE.num.mapped.txt $bwa_stats_dir/$titia_dir/$FILE.sample_cov.txt $bwa_stats_dir/$titia_dir/$FILE.mapping_proportion.txt > $bwa_stats_dir/$titia_dir/$FILE.bamstats
paste $bwa_stats_dir/$americana_dir/$FILE.samplename.txt $bwa_stats_dir/$americana_dir/$FILE.num.mapped.txt $bwa_stats_dir/$americana_dir/$FILE.sample_cov.txt $bwa_stats_dir/$americana_dir/$FILE.mapping_proportion.txt > $bwa_stats_dir/$americana_dir/$FILE.bamstats

rm $bwa_stats_dir/$titia_dir/$FILE.*.txt
rm $bwa_stats_dir/$americana_dir/$FILE.*.txt

# THEN RUN THIS IN THE CONSOLE

# cat /nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/bwa_stats_array_SDC/bwa_HetTit1.0_dg/*.bamstats > /nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/bwa_stats_array_SDC/Allsamples_HetTit1.0_dg.bamstats
# cat /nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/bwa_stats_array_SDC/bwa_HetAmer1.0_dg/*.bamstats > /nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/bwa_stats_array_SDC/Allsamples_HetAmer1.0_dg.bamstats