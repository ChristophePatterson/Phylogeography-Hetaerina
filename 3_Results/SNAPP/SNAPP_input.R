# Get system arguments
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")

library(vcfR)
library(LEA)
library(poppr)

# Load in arguments
args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[1]
SNP.library.location <- args[2]
select_K <- args[3]
select_N <- args[4]
select_N <- 2

# Read in vcf
dir.path <- paste0(SNP.library.location,SNP.library.name,"/")
filter_para <- ".all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.biSNP0_20.noX.noCUAJ"

vcf.SNPs <- read.vcfR(paste0(dir.path,SNP.library.name,filter_para,".vcf.gz"),verbose = F)
# Reorder samples so they are in alphabetical order
vcf.SNPs <- vcf.SNPs[samples = sort(colnames(vcf.SNPs@gt)[-1])]

## Read in populations assignment
pop_assign <- read.table("/home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/3_Results/LEA/LEA_pop_assign/popfile_all_noCUAJ_tabdelim.txt")
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

# Get Genind
my_genind_ti_SNPs <- vcfR2genind(vcf.SNPs, sep = "/", return.alleles = TRUE)
# Proporation of samples with missing data
sample.miss <- propTyped(my_genind_ti_SNPs , by = "ind")

# Match the vcf data to the pop_assignment
names(sample.miss)==pop_assign$sample

## Get highest coverage N individuals for each population
pop_assign$cov <- sample.miss
pops <- list(unique(pop_assign$pop))

lapply(pops, function(x) pop_assign$cov[pop_assign$pop==x])





