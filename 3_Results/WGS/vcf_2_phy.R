# Convert vcf file into format okay for RAxML
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2.1/")
#Libraries needed
library(vcfR)
library(ape)

# Get system arguments
args <- commandArgs(trailingOnly = TRUE)

# Create output location
libary.dir <- args[2]

dir.path.vcf <- args[1]
# dir.path.vcf <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/VCF_chrom_r10000"
# library.dir <- paste0(dir.path.vcf, "/results/RAxML")
vcf.name <- args[3]
#vcf.name <- "WGS_titia_chr1-12"
dir.path.RAxML <- library.dir
dir.create(dir.path.RAxML)

# REad in vcf file
vcf <- read.vcfR(paste0(dir.path.vcf, "/", vcf.name, ".vcf"), verbose = F)

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
#print(colnames(vcf@gt))
vcf <- vcf[,!colnames(vcf@gt) %in% c("CUAJa03", "CUAJa02_R")]
#colnames(vcf@gt) <- gsub(colnames(vcf@gt), pattern = "_R", replacement = "")

#Remove non-bialletlic alleles
vcf.bi <- vcf[is.biallelic(vcf),]
vcf.bi <- vcf.bi[is.polymorphic(vcf.bi,na.omit=T),]
dim(vcf.bi)

X_chrom <- "HetTit1.0.p.scaff-12-96647824"
print(table(vcf.bi@fix[,1]!=X_chrom))
# Extract genotypes
gt <- extract.gt(vcf.bi[vcf.bi@fix[,1]!=X_chrom,], return.alleles = TRUE, convertNA = TRUE)
print("Printing gt")
gt[1:10,1:8]

#Replacing  homozygous and heterozygous calls with IUPAC ambiguity codes, remove heterozgous sites
gt[gt=="A/A"] <- "A"
gt[gt=="T/T"] <- "T"
gt[gt=="G/G"] <- "G"
gt[gt=="C/C"] <- "C"
gt[gt=="A/G"] <- "?"
gt[gt=="G/A"] <- "?"
gt[gt=="C/T"] <- "?"
gt[gt=="T/C"] <- "?"
gt[gt=="A/C"] <- "?"
gt[gt=="C/A"] <- "?"
gt[gt=="G/T"] <- "?"
gt[gt=="T/G"] <- "?"
gt[gt=="C/G"] <- "?"
gt[gt=="G/C"] <- "?"
gt[gt=="A/T"] <- "?"
gt[gt=="T/A"] <- "?"
gt[gt=="."] <- "?"

# Remove sites that are completly NA
gt <- gt[rowSums(gt=="?")!=dim(gt)[1],]

# Removes sites that are no longer variant
SNP.allele.num <- vector()
for(i in 1:dim(gt)[1]){
  temp.gt <- unique(gt[i,])
  temp.gt <- temp.gt[temp.gt!="?"]
  SNP.allele.num[i] <- length(temp.gt)
}
gt <- gt[SNP.allele.num!=1,]

# Write .phy file for RAxML
ape::write.dna(t(gt), file = paste0(dir.path.RAxML,"/", vcf.name,".phy"), format ="interleaved")
