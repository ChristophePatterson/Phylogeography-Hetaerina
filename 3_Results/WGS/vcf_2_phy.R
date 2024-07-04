# Convert vcf file into format okay for RAxML
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2.1/")
#Libraries needed
library(vcfR)
library(ape)

# Get system arguments
args <- commandArgs(trailingOnly = TRUE)

# Create output location
libary.dir <- args[2]
vcf.name <- args[3]
dir.path.vcf <- args[1]
# dir.path.vcf <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/VCF_chrom_r10000"
library.dir.RAxML <- paste0(dir.path.vcf, "/results/RAxML/")
library.dir.SVD <- paste0(dir.path.vcf, "/results/SVDQuartets/")

#vcf.name <- "WGS_titia_chr1-12"
dir.create(library.dir.RAxML)
dir.create(library.dir.SVD)

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
vcf.bi <- vcf[,!colnames(vcf@gt) == "CUAJa03"]
#colnames(vcf@gt) <- gsub(colnames(vcf@gt), pattern = "_R", replacement = "")

#Remove non-bialletlic alleles
vcf.bi <- vcf.bi[is.biallelic(vcf.bi),]
vcf.bi <- vcf.bi[is.polymorphic(vcf.bi,na.omit=T),]
dim(vcf.bi)

X_chrom <- "HetTit1.0.p.scaff-12-96647824"
print(table(vcf.bi@fix[,1]!=X_chrom))
vcf.bi <- vcf.bi[vcf.bi@fix[,1]!=X_chrom,]
#Write out filtered VCF
write.vcf(vcf.bi, file = paste0(dir.path.vcf, "/", vcf.name, "_filtered.vcf.gz"))
# Write table of samples
sites <- gsub("_R", "",colnames(vcf.bi@gt)[-1])
sites <- substr(sites,1, nchar(sites)-2)
write.table(cbind(colnames(vcf.bi@gt)[-1], colnames(vcf.bi@gt)[-1]),
            file = paste0(dir.path.vcf, "/", vcf.name, "_sites.txt"),
            col.names = F, row.names = F, quote = F, sep = "\t")

# Chromosome position
chr.pos <- list()
for(i in 1:length(unique(vcf.bi@fix[,1]))){
  chr <- unique(vcf.bi@fix[,1])[i]
  chr.pos[[i]] <- which(vcf.bi@fix[,1]==chr)[c(1, length(which(vcf.bi@fix[,1]==chr)))]
}

write.table(do.call("rbind",chr.pos),
            file = paste0(dir.path.vcf, "/", vcf.name, "_chrom_start_end.txt"),
            col.names = F, row.names = F, quote = F, sep = ",")

# Remove CUAJa02 for phylogeny analysis
vcf.bi <- vcf[,!colnames(vcf@gt) %in% c("CUAJa03", "CUAJa02_R")]
#colnames(vcf@gt) <- gsub(colnames(vcf@gt), pattern = "_R", replacement = "")

#Remove non-bialletlic alleles
vcf.bi <- vcf[is.biallelic(vcf),]
vcf.bi <- vcf.bi[is.polymorphic(vcf.bi,na.omit=T),]
dim(vcf.bi)

X_chrom <- "HetTit1.0.p.scaff-12-96647824"
print(table(vcf.bi@fix[,1]!=X_chrom))

# Extract genotypes
gt <- extract.gt(vcf.bi[vcf.bi@fix[,1]!=X_chrom,], return.alleles = TRUE, convertNA = TRUE)

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

# Write .nexus file for SVD_quartets
write.nexus.data(t(gt), file = paste0(library.dir.SVD, vcf.name,".nex"), 
                 format = "DNA", missing = "?", interleaved = FALSE)

# Run without heterozgousity
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
ape::write.dna(t(gt), file = paste0(library.dir.RAxML, vcf.name,".phy"), format ="interleaved")
