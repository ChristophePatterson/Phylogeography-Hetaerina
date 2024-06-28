# Convert vcf file into format okay for RAxML
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
#Libraries needed
library(vcfR)
library(adegenet)
library(ape)

# Get system arguments
args <- commandArgs(trailingOnly = TRUE)

# Create output location
libary.dir <- args[2]

dir.path.vcf <- args[1]
vcf.name <- args[3]
dir.path.RAxML <- libary.dir <- args[2]
dir.create(dir.path.RAxML)

# REad in vcf file
vcf <- read.vcfR(paste0(dir.path.vcf, "/", vcf.name, ".vcf"))

# Reformat names
gsub("^.*/","",colnames(vcf@gt))
colnames(vcf@gt) <- gsub("^.*/","",colnames(vcf@gt))
colnames(vcf@gt) <- substr(colnames(vcf@gt), 1, nchar(colnames(vcf@gt))-4)
colnames(vcf@gt)[1] <- "FORMAT"
colnames(vcf@gt)

# Erandi sample code to lab code conversion
Erandi_sample_names <- cbind(c("ZANAa05", "CUAJa03","CUAJb01","MIXTa01"), c("E008","E018","E019","E020"))
colnames(vcf@gt)[(colnames(vcf@gt)%in%Erandi_sample_names[,2])] <- Erandi_sample_names[,1][match(colnames(vcf@gt),Erandi_sample_names[,2])[(colnames(vcf@gt)%in%Erandi_sample_names[,2])]]

# Remove CUAJa03 Erandi sample which is duplicated
vcf <- vcf[,!colnames(vcf@gt)=="CUAJa03"]
colnames(vcf@gt) <- gsub(colnames(vcf@gt), pattern = "_R", replacement = "")

#Remove non-bialletlic alleles
vcf.bi <- vcf[is.biallelic(vcf),]
dim(vcf.bi)

X_chrom <- "HetTit1.0.p.scaff-12-96647824"
print(table(vcf.bi@fix[,1]!=X_chrom))
# Extract genotypes
gt <- extract.gt(vcf.bi[vcf.bi@fix[,1]!=X_chrom], return.alleles = TRUE, convertNA = TRUE)

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

gt <- t(gt)

# Checks if each sample is in sample database
samples <- rownames(gt)
sites <- data.frame(samples)

#Remove low covage samples
snp_sub <- 0.2
snp_sub_text <-"0_2"

# Subsetting low coverage samples except from occisa samples
sites$coverage <- rowSums(is.na(as.matrix(gt)))/dim(as.matrix(gt))[2]
sites$retained <- sites$coverage<=snp_sub
gt <- gt[sites$retained,]
sites <- sites[sites$retained==T,]

gt <- t(gt)

print("Check point 1")

# Remove sites that are completly NA
gt <- gt[rowSums(is.na(gt)),]

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

print("Check point 3")

# Write .phy file for RAxML
ape::write.dna(t(gt), file = paste0(dir.path.RAxML,"/", vcf.name,".phy"), format ="interleaved")
