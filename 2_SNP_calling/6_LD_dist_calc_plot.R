.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
# Output file location
# Get task ID number
taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
taskid <- as.numeric(taskid)

libraries <- read.table("/home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/ddRAD_Durham_and_Sheffield_CUAJ/library_combinations/library_name", header = F)$V1

SNP.library.name <- libraries[taskid]

library(tidyverse)

# set path
my_bins <- paste0(SNP.library.name, ".LD.ld_decay_bins")
print(my_bins)

# read in data
ld_bins <- read_tsv(my_bins)

# plot LD decay
p <- ggplot(ld_bins, aes(distance, avg_R2)) + geom_point(alpha = 0.3) + geom_smooth() +
  xlab("Distance (bp)") + ylab(expression(italic(r)^2))

# Save file
ggsave(paste0(SNP.library.name, ".LD.ld_decay_bins.pdf"), p, height = 10, width = 20)
# ggsave(paste0("/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_Pool_X_HetTit1.0_HetAmer1.0/", SNP.library.name, ".LD.ld_decay_bins.pdf"), p, height = 10, width = 20)
