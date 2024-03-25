#!/bin/bash

#SBATCH -c 64 
#SBATCH --mem=10G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:10G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 64:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --array=1,2   # Create 32 tasks, numbers 1 to 32
#SBATCH --output=slurm-%x.%j.out

#Modules
module load r/4.2.1
module load bioinformatics
#Name of output SNP library

# Get library name
line_num=$(expr $SLURM_ARRAY_TASK_ID)
echo "$line_num"
SNP_library=$(sed -n "${line_num}p" /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/2_SNP_calling/library_combinations/library_name)
#Output directory
library_version=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_v3/)
# module load beast
# May need to run the code
# /home/tmjj24/apps/beast/bin/packagemanager -add SNAPP
#Name of output SNP library
echo $SNP_library
out_dir=($library_version/$SNP_library/SNAPP)
mkdir -p $out_dir

# select number of samples to use - Can be no more than 4
select_N=(4)

model_name=(${SNP_library}_ind_${select_N})

cd $library_version/$SNP_library

# Make input files for ruby script
Rscript /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/3_Results/SNAPP/SNAPP_input.R $SNP_library $library_version $select_N

cd SNAPP
# Run ruby script for creating xml SNAPP config file
# Uses https://github.com/mmatschiner/tutorials/blob/master/divergence_time_estimation_with_snp_data/README.md

ruby /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/ddRAD_Durham_and_Sheffield/SNAPP/snapp_prep.rb -p $model_name.phy -t $model_name.txt \
-c $model_name.con.txt -l 1000000 -x $model_name.xml -o $model_name

/nobackup/tmjj24/apps/beast/bin/beast -threads $SLURM_CPUS_PER_TASK -overwrite $model_name.xml > $model_name.screen.log

/nobackup/tmjj24/apps/beast/bin/treeannotator -burnin 10 $model_name.trees $model_name.trees.Anon

cd ~