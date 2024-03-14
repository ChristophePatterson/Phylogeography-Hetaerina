
# Create filtered SNP library
sbatch /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/2_SNP_calling/7_SNP_filtering_vcfR.sh

# After generating filtered SNP libraries these results files can be run in parrell
sbatch /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/3_Results/SVDQuartets/SVDQuartets_arrary.sh
sbatch /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/3_Results/RAxML/9_RAxML_array.sh
sbatch /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/3_Results/LEA/LEA_PCA.sh


