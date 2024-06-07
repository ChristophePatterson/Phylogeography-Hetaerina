
# Load R packages
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
library(ggplot2)
library(ape)
library(vcfR)
library(poppr)
library(adegenet)

# Get system arguments
args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[1]
SNP.library.location <- args[2]
# Directory
dir.path <- paste0(SNP.library.location,SNP.library.name,"/")
filter_para <- ".all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.biSNP0_20.noX"
# Read in vcf
vcf <- read.vcfR(paste0(dir.path,SNP.library.name,filter_para,".vcf.gz"))
#vcf <- read.vcfR("/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/Hetaerina_titia_ddRAD_titia_dg/Hetaerina_titia_ddRAD_titia_dg.all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.biSNP0_20.noX.vcf.gz")
#vcf <- read.vcfR("/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/Hetaerina_all_ddRAD_titia_dg/Hetaerina_all_ddRAD_titia_dg.all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.biSNP0_20.noX.vcf.gz")
## Get SNP posistions
RAD.loci.length <- 1000
RAD.loci <- paste0(vcf@fix[,"CHROM"], ":", as.numeric(vcf@fix[,"POS"])-(RAD.loci.length/2), "-" , as.numeric(vcf@fix[,"POS"])+(RAD.loci.length/2))

RAD.loci <- gsub(RAD.loci, pattern = "_", replacement = "\\.")

writeLines(con = paste0(dir.path,SNP.library.name,"_RAD_posistions.txt"), text = RAD.loci)


