# Get system arguments
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")

library(vcfR)
library(LEA)
library(poppr)
library(tidyverse)
library(ape)

# Load in arguments
args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[2]
# SNP.library.name <- "Hetaerina_all_ddRAD_americana_dg"
SNP.library.location <- args[1]
# SNP.library.location  <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_v2/"

dir.path.SNAPP <- paste0(SNP.library.location, "/results/SNAPP/")
print(dir.path.SNAPP)
dir.create(dir.path.SNAPP)
# Read in vcf

vcf <- read.vcfR(paste0(SNP.library.location, SNP.library.name,".vcf"),verbose = F)
# Reorder samples so they are in alphabetical order
# Reformat names
colnames(vcf@gt) <- gsub("^.*/","",colnames(vcf@gt))
colnames(vcf@gt) <- substr(colnames(vcf@gt), 1, nchar(colnames(vcf@gt))-4)
colnames(vcf@gt)[1] <- "FORMAT"

# Erandi sample code to lab code conversion
Erandi_sample_names <- cbind(c("ZANAa05", "CUAJa03","CUAJb01","MIXTa01"), c("E008","E018","E019","E020"))
colnames(vcf@gt)[(colnames(vcf@gt)%in%Erandi_sample_names[,2])] <- Erandi_sample_names[,1][match(colnames(vcf@gt),Erandi_sample_names[,2])[(colnames(vcf@gt)%in%Erandi_sample_names[,2])]]

# Remove CUAJa03 Erandi sample which is duplicated
vcf <- vcf[,!colnames(vcf@gt)=="CUAJa03"]
colnames(vcf@gt) <- gsub(colnames(vcf@gt), pattern = "_R", replacement = "")
vcf <- vcf[,!colnames(vcf@gt)=="CUAJa02"]

sample_map <- colnames(vcf@gt)[-1]

vcf@gt[1:10,]
vcf.bi <- vcf[is.biallelic(vcf)]
vcf.ply <- vcf.bi[is.polymorphic(vcf.bi, na.omit = T)]
# Convert to geno
X_chrom <- "HetTit1.0.p.scaff-12-96647824"

table(vcf.ply@fix[,1]!=X_chrom)
# remove X chrom
geno <- extract.gt(vcf.ply[vcf.ply@fix[,1]!=X_chrom])

DP.mat <- extract.gt(vcf.ply, element = "DP")
Samp.DP <- apply(DP.mat, MARGIN = 2, function(x) mean(as.numeric(x), na.rm = T))
SNP.DP <- apply(DP.mat, MARGIN = 1, function(x) mean(as.numeric(x), na.rm = T))

c(mean(Samp.DP,na.rm =T), median(Samp.DP,na.rm =T), range(Samp.DP,na.rm = T))
c(mean(SNP.DP,na.rm =T), median(SNP.DP,na.rm =T), range(SNP.DP,na.rm = T))

# Extract genotype
gt <- extract.gt(vcf.ply[vcf.ply@fix[,1]!=X_chrom,], return.alleles = TRUE, convertNA = TRUE)

# Overwrite data with Ambiguity codes
gt[gt=="A/A"] <- "A"
gt[gt=="T/T"] <- "T"
gt[gt=="G/G"] <- "G"
gt[gt=="C/C"] <- "C"
gt[gt=="A/G"] <- "R"
gt[gt=="G/A"] <- "R"
gt[gt=="C/T"] <- "Y"
gt[gt=="T/C"] <- "Y"
gt[gt=="A/C"] <- "M"
gt[gt=="C/A"] <- "M"
gt[gt=="G/T"] <- "K"
gt[gt=="T/G"] <- "K"
gt[gt=="C/G"] <- "S"
gt[gt=="G/C"] <- "S"
gt[gt=="A/T"] <- "W"
gt[gt=="T/A"] <- "W"
gt[gt=="."] <- "?"

## Set up files following https://github.com/mmatschiner/tutorials/blob/master/divergence_time_estimation_with_snp_data/README.md
## (1) Phylip file

ape::write.dna(t(gt), file = paste0(dir.path.SNAPP, SNP.library.name,".phy"),
               format ="interleaved", nbcol = -1, colsep = "")
## (2) Species table
species.df <- data.frame(species = paste0("titia_", colnames(gt)), sample = colnames(gt))
write.table(species.df, file = paste0(dir.path.SNAPP, SNP.library.name,".txt"),
            row.names = F,quote = F, sep = "\t")
# (3) Theta priors from Stranding et al 2022
# All Hetaerina samples
# normal(offset,mean,sigma)
contrant.df <- data.frame(x = "normal(0,3,0.1)", y = "crown", 
                          z = paste0(paste0("titia_", colnames(gt)), sep = ",",collapse = ""))

#remove trailing comma
contrant.df$z <- substr(contrant.df$z, start = 1 , stop = nchar(contrant.df$z)-1)

write.table(contrant.df, file = paste0(dir.path.SNAPP, SNP.library.name,".con.txt"),
            row.names = F, col.names = F, quote = F, sep = "\t")
