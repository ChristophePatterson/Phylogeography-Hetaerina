#Hamilton LEA and PCR plots
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2.1/")
# Combined PCA and admixture plot

# CalcuLation of PCR on titia popuLations
library(ggplot2)
library(ape)
library(vcfR)
library(tidyverse)

# Get library names

args <- commandArgs(trailingOnly = TRUE)
SNP.library.name <- args[1]
SNP.library.location <- args[2]

# Get task ID number
taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
taskid <- as.numeric(taskid)

print(paste("########## This is taskID", taskid, "and will use the SNP library:", SNP.library.name,"##########"))

species <- unlist(strsplit("Hetaerina_titia_ddRAD_titia_dg", split = "_"))[1:2]

species <- paste(species[1], species[2], sep = "_")

# Directory
dir.path <- paste0(SNP.library.location,SNP.library.name,"/")
filter_para <- ".all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.biSNP0_20.noX"

dir.create(paste0(dir.path, "dsuite/"))

vcf.SNPs <- read.vcfR(paste0(dir.path,SNP.library.name,filter_para,".vcf.gz"),verbose = T)
# Reorder samples so they are in alphabetical order
vcf.SNPs <- vcf.SNPs[samples = sort(colnames(vcf.SNPs@gt)[-1])] 

vcf.SNPs@fix[,1] <- gsub(vcf.SNPs@fix[,1], pattern = "\\.", replacement = "_")
# Get X chromosome
X_chrom <- "HetTit1_0_p_scaff-12-96647824"

# read in sample information
analysis.name <- 
sample_map<- read.csv("/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/All samples held in Durham_v17.csv", check.names=F)
names(sample_map)

#Extracting sample information and linking it to tree labels
samples <- colnames(vcf.SNPs@gt)[-1]
sites <- data.frame(samples)
sites$samples
sites <- merge(sites, sample_map[,c("Unique.ID","sex","species","Year","Site.ID","Lat","Long","Country","Province","Ocean.drainage")], by.x = "samples", by.y = "Unique.ID", all.x = T)
sites$site.sub <- paste(gsub('[[:digit:]]+', '', sites$Site.ID))
sites$site.sub[sites$site.sub=="BZ"] <- sites$Site.ID[sites$site.sub=="BZ"]
sites$site.sub[sites$site.sub=="NA"] <- sites$Site.ID[sites$site.sub=="NA"]

## Relevent sample sites
Atl.sites <- c("CT", "TXRS", "CUAJ", "MIXT")
Pac.sites <- c("PUMA", "RLPE", "ZANA", "TULI")

#Subset sites
hybrid.sites <- sites[sites$site.sub%in%c(Atl.sites, Pac.sites),]
hybrid.sites$Ocean.drainage <- as.factor(hybrid.sites$Ocean.drainage)
hybrid.sites$Lat <- as.numeric(hybrid.sites$Lat)
hybrid.sites$Long <- as.numeric(hybrid.sites$Long)

colnames(vcf.SNPs@gt)
vcf.SNPs <- vcf.SNPs[,colnames(vcf.SNPs@gt)%in%c("FORMAT", hybrid.sites$samples)]

#Remove non-bialletlic alleles
vcf.SNPs <- vcf.SNPs[is.biallelic(vcf.SNPs),]
vcf.SNPs <- vcf.SNPs[is.polymorphic(vcf.SNPs,na.omit=T),]

vcf.SNPs
## Sample site information for dsuite
write.vcf(vcf.SNPs, file = paste0(dir.path, "dsuite/hybrid_samples.vcf.gz"))
hybrid.sites$site.dsuite <- hybrid.sites$Site.ID
hybrid.sites$site.dsuite[hybrid.sites$Site.ID=="CUAJ01"] <- hybrid.sites$samples[hybrid.sites$Site.ID=="CUAJ01"]
cbind(hybrid.sites$samples,hybrid.sites$site.dsuite)
write.table(cbind(hybrid.sites$samples,hybrid.sites$site.dsuite), file = paste0(dir.path, "dsuite/pop_file.txt"),
            col.names = F, row.names = F, quote = F, sep = "\t")