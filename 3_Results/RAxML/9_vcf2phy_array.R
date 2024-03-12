# Convert vcf file into format okay for RAxML

#Libraries needed
library(vcfR)
library(adegenet)
library(ape)

# Get task ID number
taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
taskid <- as.numeric(taskid)

libraries <- read.table("/home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/ddRAD_Durham_and_Sheffield/library_combinations/library_name", header = F)$V1

SNP.library.name <- libraries[taskid]

print(paste("########## This is taskID", taskid, "and will use the SNP library:", SNP.library.name,"##########"))
SNP.library.filter <- ".snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.vcf.gz"

# write(x = paste0("taskID is ", taskid), file = paste0("taskID is ", taskid, " and will use library", SNP.library.name))

# Output file location

# Create output location
libary.dir <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_Pool_X_HetTit1.0_HetAmer1.0/"

dir.path.vcf <- paste0(libary.dir, SNP.library.name,"/")
dir.path.RAxML <- paste0(dir.path.vcf,"/RAxML/")
dir.create(dir.path.RAxML)

# REad in vcf file
vcf <- read.vcfR(paste0(dir.path.vcf,SNP.library.name,".snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.vcf.gz"))

#Remove non-bialletlic alleles
vcf.bi <- vcf[is.biallelic(vcf),]
dim(vcf.bi)

# Extract genotypes
gt <- extract.gt(vcf.bi, return.alleles = TRUE, convertNA = TRUE)

#Shortening sample names
colnames(gt) <- gsub("^.*/","",colnames(gt))
colnames(gt) <- substr(colnames(gt), 1, nchar(colnames(gt))-11)

#Removes Pool X notation
colnames(gt) <- gsub(pattern = "_X", x = colnames(gt), replacement = "")

# Replacing correct occisa label
colnames(gt)[colnames(gt)=="TILIaOc01"] <- "TULIaOc01"

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

gt <- t(gt)

# Checks if each sample is in sample database
sample_map<- read.csv("/home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/ddRAD_Durham_and_Sheffield/library_combinations/All samples held in Durham_v15_titia_correct.csv", check.names=F)
samples <- rownames(gt)
sites <- data.frame(samples)
#Convert lat and long to number
sample_map$Long <- as.numeric(sample_map$Long)
sample_map$Lat <- as.numeric(sample_map$Lat)

for(i in 1:length(sites$samples)){
  if(sites$samples[i]%in%sample_map$Unique.ID) {
    sites$species[i] <- sample_map$species[which(sample_map$Unique.ID==sites$samples[i])][1]
    sites$site.name[i] <- sample_map$Site.ID[which(sample_map$Unique.ID==sites$samples[i])][1]
    sites$site.county[i] <- sample_map$Country[which(sample_map$Unique.ID==sites$samples[i])][1]
    sites$Province[i] <- sample_map$Province[which(sample_map$Unique.ID==sites$samples[i])][1]
    sites$long[i] <- round(sample_map$Long[which(sample_map$Unique.ID==sites$samples[i])][1], 2)
    sites$lat[i] <- round(sample_map$Lat[which(sample_map$Unique.ID==sites$samples[i])][1], 2)
    sites$ocean.drainage[i] <- sample_map$Ocean.drainage[which(sample_map$Unique.ID==sites$samples[i])][1]
  } else sites[i,] <- c(sites$samples[i], rep(NA,dim(sites)[2]-1))
}


#Remove low covage samples
snp_sub <- 0.2
snp_sub_text <-"0_2"

# Subsetting low coverage samples except from occisa samples
sites$coverage <- rowSums(is.na(as.matrix(gt)))/dim(as.matrix(gt))[2]
sites$retained <- sites$coverage<=snp_sub|sites$species=="occisa"
table(sites$retained, sites$species)
gt <- gt[sites$retained,]
sites <- sites[sites$retained==T,]

gt <- t(gt)

print("Check point 1")

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
ape::write.dna(t(gt), file = paste0(dir.path.RAxML,SNP.library.name,".phy"), format ="interleaved")
