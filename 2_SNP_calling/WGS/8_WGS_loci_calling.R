
# Load R packages
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2.1/")
library(vcfR)

# Get system arguments
args <- commandArgs(trailingOnly = TRUE)
SNP.library.name <- args[1]
SNP.library.location <- args[2]
# Directory
dir.path <- paste0(SNP.library.location,SNP.library.name,"/")
# Read in vcf
vcf.SNPs <- read.vcfR(paste0(dir.path, args[3],".vcf"), verbose = F)
# vcf.SNPs <- read.vcfR("/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/VCF_chrom_r10000/WGS_titia_chr1-12.vcf")
#vcf <- read.vcfR("/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/Hetaerina_all_ddRAD_titia_dg/Hetaerina_all_ddRAD_titia_dg.all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.biSNP0_20.noX.vcf.gz")
## Get SNP posistions
RAD.loci.length <- 1000
RAD.loci <- paste0(vcf.SNPs@fix[,"CHROM"], ":", as.numeric(vcf.SNPs@fix[,"POS"])-(RAD.loci.length/2), "-" , as.numeric(vcf.SNPs@fix[,"POS"])+(RAD.loci.length/2))

N_loci <- 2000
subset.loci <- sort(sample(1:length(RAD.loci), N_loci))
# Name of RAD loci
RAD.loci <- gsub(RAD.loci, pattern = "_", replacement = "\\.")
# Sample 2000 of them
RAD.loci <- RAD.loci[subset.loci]
# Write out SNP posistion with a buffer of qualt to RAD loci length
writeLines(con = paste0(dir.path,SNP.library.name,"_RAD_posistions.txt"), text = RAD.loci)

# bcftool mpileup regions format file
RAD.loci <- cbind(vcf.SNPs@fix[,"CHROM"], as.numeric(vcf.SNPs@fix[,"POS"])-(RAD.loci.length/2), as.numeric(vcf.SNPs@fix[,"POS"])+(RAD.loci.length/2))
RAD.loci <- RAD.loci[subset.loci,]

write.table(file = paste0(dir.path,SNP.library.name,"_RAD_posistions_delim.txt"), RAD.loci, sep = "\t", quote = FALSE,row.names = FALSE,col.names = FALSE)



