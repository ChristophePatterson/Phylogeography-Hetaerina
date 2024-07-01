
# Load R packages
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
library(ggplot2)
library(ape)
library(vcfR)
library(poppr)
library(adegenet)


#######################################################################################################
#######################################################################################################
###### vcfR SNP filtering ######
#######################################################################################################
#######################################################################################################

# Get system arguments
args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[1]
SNP.library.location <- args[2]
# Directory
dir.path <- paste0(SNP.library.location,SNP.library.name,"/")
filter_para <- ".all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000"

# Change working directory
setwd(dir.path)

## Running analysis name
analysis.name <- "H_titia_complete_snp_noX"

vcf <- read.vcfR(paste0(dir.path,SNP.library.name,filter_para,".vcf.gz"))

# Remames columns
# gsub("^.*/","",colnames(vcf@gt))
colnames(vcf@gt) <- gsub("^.*/","",colnames(vcf@gt))
colnames(vcf@gt) <- substr(colnames(vcf@gt), 1, nchar(colnames(vcf@gt))-11)
colnames(vcf@gt)[1] <- "FORMAT"
# Rename samples that end in _X
colnames(vcf@gt) <- gsub(pattern = "_X", x = colnames(vcf@gt), replacement = "")

# Removing dots from SNP names because genind doesn't like this
vcf@fix[,1] <- gsub(vcf@fix[,1], pattern = "\\.", replacement = "_")
# Removes SNPs that do not vary between samples (e.g heterozgous across all samples)
vcf.bi <- vcf[is.polymorphic(vcf, na.omit = T),]
# Removes SNPs that aren't biallelic
vcf.bi <- vcf.bi[is.biallelic(vcf.bi),]

# Removes potentially containminated samples
dodgy <- c("HtiTi12", "NA0101","CA0101","CUAJa02.Dur","CUAJa02.Shef","HXRCaAM03", "CUAJa03")
vcf.bi <- vcf.bi[,!colnames(vcf.bi@gt)%in%dodgy]

my_genind_ti <- vcfR2genind(vcf.bi, sep = "/", return.alleles = TRUE)

#Histogram of sample and SNP coverage
SNP.miss <- propTyped(my_genind_ti, by = "loc")
sample.miss <- propTyped(my_genind_ti, by = "ind")
hist(sample.miss)
hist(SNP.miss)

min_cov <- min(sample.miss)
max_cov <- max(sample.miss)

# Set cut of value for sample removal
snp_sub <- 0.20
snp_sub_text <-"0_20"

#Remove SNPs and samples with lower than cut off values
my_genind_ti_SNPs <- missingno(my_genind_ti, type = "geno", cutoff = snp_sub)

#Histogram of sample and SNP coverage
SNP.miss <- propTyped(my_genind_ti_SNPs, by = "loc")
sample.miss <- propTyped(my_genind_ti_SNPs, by = "ind")
hist(sample.miss)
hist(SNP.miss)

# Checks if CUAJa02 is included
loc_CUAJa02 <- which(rownames(my_genind_ti_SNPs@tab)=="CUAJa02")
sample.miss[loc_CUAJa02]

#Removes non polymorphic SNPs
poly_loci <- isPoly(my_genind_ti_SNPs, by = "locus")
print("Check Point 1")
my_genind_ti_SNPs <- my_genind_ti_SNPs[loc = names(which(poly_loci==T))]

#Number of SNPs and samples and proportion of missing sites
my_genind_ti_SNPs
dim(my_genind_ti_SNPs@tab)

# Mean percentage of missing data per sample
Sample.miss <- propTyped(my_genind_ti_SNPs, by = "ind")
range(Sample.miss*100)
mean(Sample.miss)*100

# Mean percentage missing data per SNP
SNP.miss <- propTyped(my_genind_ti_SNPs, by = "loc")
range(SNP.miss)
mean(SNP.miss)*100

print("Check Point 2")
# Get limit vcf file to those identified by filtering
vcf.SNPs <- vcf.bi[samples = rownames(my_genind_ti_SNPs@tab)]
poly.snps <- names(my_genind_ti_SNPs@all.names)
all_snps <- paste(vcf.SNPs@fix[,1],vcf.SNPs@fix[,2], sep = "_")
vcf.SNPs <- vcf.SNPs[all_snps%in%poly.snps,]

#### Write and read back in filtered vcf
# Write all snps
write.vcf(vcf.SNPs, file = paste0(dir.path,SNP.library.name,filter_para,".biSNP",snp_sub_text,".vcf.gz"),mask = F,APPEND = F)
# Write all snps without X
# Removes SNPs on the X chromosome
# For H. tita library
if(grepl(SNP.library.name, pattern="titia_dg")){
  X_chrom <- "HetTit1_0_p_scaff-12-96647824"
}
# For H. americana lots of small X chromosomes
if(grepl(SNP.library.name, pattern="americana_dg")){
  X_chrom <- c("JAKTNV010000016_1","JAKTNV010000031_1","JAKTNV010000013_1","JAKTNV010000060_1","JAKTNV010000023_1","JAKTNV010000051_1","JAKTNV010000012_1")
}
vcf.SNPs@fix[1:5,]

# For H. americana lots of small X chromosomes
if(grepl(SNP.library.name, pattern="elegans_dg")){
  X_chrom <- "OV121106.1"
}

vcf.SNPs <- vcf.SNPs[!vcf.SNPs@fix[,1]%in%X_chrom,]
table(vcf.SNPs@fix[,1]%in%X_chrom)

write.vcf(vcf.SNPs, file = paste0(dir.path,SNP.library.name,filter_para,".biSNP",snp_sub_text,".noX.vcf.gz"),mask = F,APPEND = F)

# Write all snps without X and CUAJ samples
# Get CUAJ samples if there are any
if(any(grepl(colnames(vcf.SNPs@gt), pattern = "CUAJ"))){
  nonCUAJ <- colnames(vcf.SNPs@gt)[-grep(colnames(vcf.SNPs@gt), pattern = "CUAJ")]
}else{nonCUAJ <- colnames(vcf.SNPs@gt)}
# removes first column with "format" from list of samples
print("Check Point 3")
nonCUAJ <- nonCUAJ[-1]
vcf.SNPs.nCUAJ <- vcf.SNPs[samples = nonCUAJ]
write.vcf(vcf.SNPs.nCUAJ, file = paste0(dir.path,SNP.library.name,filter_para,".biSNP",snp_sub_text,".noX.noCUAJ.vcf.gz"),mask = F,APPEND = F)

#######################################################################################################
### Create in put for RAxML
#######################################################################################################
dir.path.RAxML <- paste0(dir.path, "RAxML/")
dir.create(dir.path.RAxML)

# Extract genotypes
gt <- extract.gt(vcf.SNPs.nCUAJ, return.alleles = TRUE, convertNA = TRUE)

#Replacing  homozygous and heterozygous calls with IUPAC ambiguity codes, remove heterozgous sites
gt[gt=="A/A"] <- "A"
gt[gt=="T/T"] <- "T"
gt[gt=="G/G"] <- "G"
gt[gt=="C/C"] <- "C"
gt[gt=="A/G"] <- NA
gt[gt=="G/A"] <- NA
gt[gt=="C/T"] <- NA
gt[gt=="T/C"] <- NA
gt[gt=="A/C"] <- NA
gt[gt=="C/A"] <- NA
gt[gt=="G/T"] <- NA
gt[gt=="T/G"] <- NA
gt[gt=="C/G"] <- NA
gt[gt=="G/C"] <- NA
gt[gt=="A/T"] <- NA
gt[gt=="T/A"] <- NA
gt[gt=="."] <- NA

# Removes sites that are no longer variant
SNP.allele.num <- vector()
for(i in 1:dim(gt)[1]){
  temp.gt <- unique(gt[i,])
  temp.gt <- temp.gt[!is.na(temp.gt)]
  SNP.allele.num[i] <- length(temp.gt)
}

print("Check point 2")
gt <- gt[SNP.allele.num!=1,]

#Convert NA into "?"
gt[is.na(gt)] <- "?"

# Write .phy file for RAxML
ape::write.dna(t(gt), file = paste0(dir.path.RAxML,SNP.library.name,".noX.noCUAJ.phy"), format ="interleaved")

#######################################################################################################
###### CREATE NEXUS FOR SVDQUARTETs 
#######################################################################################################

dir.path.SVDQuartet <- paste0(dir.path, "SVDQuartets/")
dir.create(dir.path.SVDQuartet)
# Extract genotypes
gt <- extract.gt(vcf.SNPs.nCUAJ, return.alleles = TRUE, convertNA = TRUE)
#Replacing  homozygous and heterozygous calls with IUPAC ambiguity codes,
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
gt[gt=="."] <- NA

#Convert NA into "?"
gt[is.na(gt)] <- "?"

# Creating Paup excute file
nexus.list <- paste0("log  file=",dir.path.SVDQuartet,SNP.library.name, "_T20.log;")
nexus.list[2] <- paste0("cd ",dir.path.SVDQuartet,";")
nexus.list[3] <- "set hashcomment=y;"
nexus.list[4] <- paste0("execute ",dir.path.SVDQuartet,SNP.library.name, ".noX.noCUAJ.nex;")
nexus.list[5] <- ""

nexus.list[6] <- paste("svdQuartets bootstrap=standard nthreads=20;")
nexus.list[7] <- paste0("saveTrees file = ",dir.path.SVDQuartet,SNP.library.name, "_T20.tre supportValues = nodeLabels;")
nexus.list[8] <- "quit;"

writeLines(nexus.list, con = paste0(dir.path.SVDQuartet,"Paup_excute_file.nex"))

# Write .nexus file for SVD_quartets
write.nexus.data(t(gt), file = paste0(dir.path.SVDQuartet,SNP.library.name,".noX.noCUAJ.nex"), 
                 format = "DNA", missing = "?", interleaved = FALSE)






