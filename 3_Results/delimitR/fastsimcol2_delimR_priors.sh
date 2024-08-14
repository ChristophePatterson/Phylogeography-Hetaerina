#!/bin/bash

#SBATCH -c 12 
#SBATCH --mem=200G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:200G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 48:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --array=4,6   # Create 32 tasks, numbers 1 to 32
#SBATCH --output=slurm-%x.%j.out

module load r/4.2.1

# Activate  conda enviroment
echo "Activating conda enviroment"
source /home/${USER}/.bashrc
conda activate delimR

# Path to working directory
# Get library name
line_num=$(expr $SLURM_ARRAY_TASK_ID)
echo "$line_num"
library=$(sed -n "${line_num}p" /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/library_combinations/library_name)
dir_path=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/$library/)
species=$(sed -n "${line_num}p" /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/library_combinations/species)
cd $dir_path

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
    *)  # Default case if library doesn't match any expected value
        echo "Unknown library"
        exit 1
        ;;
esac

output_dir=(delimitR/${sp_0}_${sp_1}_${sp_2}_nreps10000_MIG3_divWgene_sbatch)

# Location of popgen file
# pop_file=(/home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/ddRAD_Durham_and_Sheffield_CUAJ/library_combinations/LEA_populations/popfile_titia_noCUAJ.txt)
pop_file=(/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/LEA/LEA_pop_assign/popfile_${species}_noCUAJ.txt)

# Select directory that delimitR has inputed files into
cd $output_dir

# Remove previous prior runs
rm -r Prior/

echo "Running Rscript"
# Rscript that produces reduced priors
Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/delimitR/delimR_run_delimR_priors.R $dir_path $output_dir $SLURM_CPUS_PER_TASK $library
#  Once this has been created results plotted using this script
Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/delimitR/delimitR_plot_results.R $dir_path $output_dir $SLURM_CPUS_PER_TASK $library

## Plot each demographic scenario using demographic plots
mkdir -p demo_plots/

num_scenarios=$(find $dir_path/$output_dir -type f -name "*.tpl" | wc -l)
num_scenarios=$(find . -type f -name "*.tpl" | wc -l)

for i in $(seq 1 ${num_scenarios})
do
    echo "${library}_${i}/${library}_${i}_1.par"
    Rscript /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/delimitR/ParFileViewer.r ${library}_${i}/${library}_${i}_1.par 
    mv ${library}_${i}/${library}_${i}_1.par.pdf demo_plots/${library}_${i}_1.par.pdf
done