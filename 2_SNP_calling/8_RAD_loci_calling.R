
# Load R packages
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
library(ape)
library(vcfR)
library(poppr)
library(adegenet)
library(tidyverse)

# Get system arguments
args <- commandArgs(trailingOnly = TRUE)
SNP.library.name <- args[1]
# SNP.library.name <- "Hetaerina_all_ddRAD_titia_dg"
SNP.library.location <- args[2]
# SNP.library.location <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/"
select_N <- as.numeric(args[3])
# Directory
dir.path <- paste0(SNP.library.location,SNP.library.name,"/")
filter_para <- ".all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.biSNP0_20.noX"
# Read in vcf
vcf.SNPs <- read.vcfR(paste0(dir.path,SNP.library.name,filter_para,".vcf.gz"))
#vcf <- read.vcfR("/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/Hetaerina_titia_ddRAD_titia_dg/Hetaerina_titia_ddRAD_titia_dg.all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.biSNP0_20.noX.vcf.gz")
#vcf <- read.vcfR("/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/Hetaerina_all_ddRAD_titia_dg/Hetaerina_all_ddRAD_titia_dg.all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.biSNP0_20.noX.vcf.gz")
## Get SNP posistions
RAD.loci.length <- 1000
RAD.loci <- paste0(vcf.SNPs@fix[,"CHROM"], ":", as.numeric(vcf.SNPs@fix[,"POS"])-(RAD.loci.length/2), "-" , as.numeric(vcf.SNPs@fix[,"POS"])+(RAD.loci.length/2))

# Name of RAD loci
RAD.loci <- gsub(RAD.loci, pattern = "_", replacement = "\\.")
# Write out SNP posistion with a buffer of qualt to RAD loci length
writeLines(con = paste0(dir.path,SNP.library.name,"_RAD_posistions.txt"), text = RAD.loci)


# Reorder samples so they are in alphabetical order
vcf.SNPs <- vcf.SNPs[samples = sort(colnames(vcf.SNPs@gt)[-1])]

## Read in populations assignment
pop_assign <- read.table("/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/LEA/LEA_pop_assign/popfile_all_noCUAJ_tabdelim.txt")
colnames(pop_assign) <- c("sample", "pop")

# Remove samples from vcf that are not in pop_assign
vcf.SNPs <- vcf.SNPs[samples = colnames(vcf.SNPs@gt[,-1])[colnames(vcf.SNPs@gt[,-1])%in%pop_assign$sample]]
# Remove samples in pop_assign that are not in vcf
pop_assign <- pop_assign[pop_assign$sample%in%colnames(vcf.SNPs@gt[,-1]),]

any(!pop_assign$sample%in%colnames(vcf.SNPs@gt[,-1]))
any(!colnames(vcf.SNPs@gt[,-1])%in%pop_assign$sample)

# Reorder pop_asssign to be in same order as vcf
pop_assign <- pop_assign[match(colnames(vcf.SNPs@gt[,-1]), pop_assign$sample),]
# Check order
pop_assign$sample==colnames(vcf.SNPs@gt[,-1])

# Removing dots from SNP names because genind doesn't like this
vcf.SNPs@fix[,1] <- gsub(vcf.SNPs@fix[,1], pattern = "\\.", replacement = "_")
# Get Genind
my_genind_ti_SNPs <- vcfR2genind(vcf.SNPs, sep = "/", return.alleles = TRUE)
# Proporation of samples with missing data
sample.miss <- propTyped(my_genind_ti_SNPs , by = "ind")

# Match the vcf data to the pop_assignment
names(sample.miss)==pop_assign$sample

## Get highest coverage N individuals for each population
pop_assign$cov <- sample.miss

# Number of samples that are assigned to each group
print("Number of samples that are assigned to each group")
table(pop_assign$pop)
min(table(pop_assign$pop))

top.cov.samples <- lapply(unique(pop_assign$pop), function(x) {
    pop_assign_temp <- pop_assign[pop_assign$pop==x,]
    if(dim(pop_assign_temp)[1]>select_N){
        return((pop_assign_temp[order(pop_assign_temp$cov, decreasing = T),])[1:select_N,])
    }
    if(dim(pop_assign_temp)[1]<=select_N){
        return((pop_assign_temp[order(pop_assign_temp$cov, decreasing = T),])[1:dim(pop_assign_temp)[1],])
    }
})

pop_assign_N <- do.call("rbind", top.cov.samples)

# Because older bcf file names have full file path need to match sample names with full file path
full_file_samples <- readLines(paste0(dir.path,"/samples_temp.txt"))
grep_file_match <- paste0(pop_assign_N$sample, "|", collapse = "")
grep_file_match <- substr(grep_file_match, 1, nchar(grep_file_match)-1)
full_file_samples[grep(grep_file_match, full_file_samples)]

# Write out selected  samples
writeLines(con = paste0(dir.path,SNP.library.name,"_PHYLIP_subsample.txt"), text = full_file_samples[grep(grep_file_match, full_file_samples)])