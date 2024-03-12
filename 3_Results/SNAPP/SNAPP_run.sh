#!/bin/bash

#SBATCH -c 24 
#SBATCH --mem=10G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:10G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 72:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=slurm-%x.%j.out

#Modules
module load r/4.2.1
module load bioinformatics
#Name of output SNP library

# Get library and genome names
SNP_library=("Hetaerina_all_ddRAD_titia_dg")
library_version=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_v2/)
# module load beast
# May need to run the code
# /home/tmjj24/apps/beast/bin/packagemanager -add SNAPP
#Name of output SNP library
out_dir=($library_version/$SNP_library/SNAPP)
mkdir -p $out_dir

select_K=(8)
select_N=(2)

model_name=(${SNP_library}_ind_${select_N}_K${select_K}_run$line_num)

cd $library_version/$SNP_library

# Make input file
Rscript /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/3_Results/SNAPP/SNAPP_input.R $SNP_library $library_version $select_K $select_N


# Run ruby script for creating xml SNAPP config file
# Uses https://github.com/mmatschiner/tutorials/blob/master/divergence_time_estimation_with_snp_data/README.md

ruby /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/ddRAD_Durham_and_Sheffield/SNAPP/snapp_prep.rb -p $model_name.phy -t $model_name.txt \
-c $model_name.txt -l 10000000 -x $model_name.xml -o $model_name

/home/tmjj24/apps/beast/bin/beast -threads 24 -overwrite $model_name.xml > $model_name.screen.log

/home/tmjj24/apps/beast/bin/treeannotator -burnin 10 $model_name.trees $model_name.trees.Anon

cd ~