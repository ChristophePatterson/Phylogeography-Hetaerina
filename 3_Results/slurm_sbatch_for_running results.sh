
# Create filtered SNP library
sbatch /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/7_SNP_filtering_vcfR.sh

# After generating filtered SNP libraries these results files can be run in parrell
sbatch /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/SVDQuartets/SVDQuartets_arrary.sh
sbatch /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/RAxML/9_RAxML_array.sh
sbatch /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/LEA/LEA_PCA.sh

## After running LEA
sbatch /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/G-Phocs/g-phocs-run.sh
sbatch /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/SNAPP/SNAPP_run.sh

## Runnning delimitR
# Run 1
sbatch /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/delimitR/easySFS_2_delimitR.sh
# Then run 2 (after checking downprojecting is correct)
sbatch /home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/delimitR/easySFS_2_delimitR.sh
