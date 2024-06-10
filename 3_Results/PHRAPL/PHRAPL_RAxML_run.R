# Load R packages
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
library(ggplot2)
library(ape)
library(poppr)
library(adegenet)
library(phrapl)

# Get system arguments
args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[1]
# SNP.library.name <- "Hetaerina_all_ddRAD_titia_dg"
SNP.library.location <- args[2]
# SNP.library.location <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/"
# Directory
dir.path <- paste0(SNP.library.location,SNP.library.name,"/")

## Load RAxML (get location by running module show raxml and then opening the file to see where RAxML is)
RAXML_filepath <- "/nobackup/dbl0hpc/apps/miniconda3/envs/raxml/bin/"
# ?RunRaxml

source("/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/PHRAPL/RunRaxml_custom_function.R")

RunRaxml_v2(raxmlPath=RAXML_filepath, raxmlVersion="raxmlHPC-PTHREADS",inputPath=paste0(dir.path, "RAD_loci/"), mutationModel="GTRGAMMA", 
         iterations=2, seed=sample(1:10000000,1), outputSeeds=FALSE,Threads=24, discard=TRUE)


## Merge all trees into one
Mer