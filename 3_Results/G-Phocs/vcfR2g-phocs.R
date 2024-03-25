# Convert vcf file into format okay for G-phocs
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
#Libraries needed
library(vcfR)
library(adegenet)
library(ape)
library(ggplot2)
library(tidyverse)

# Get system arguments
args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[1]
SNP.library.location <- args[2]
output_dir <- as.numeric(args[4])
# Number of samples to use
select_N <- as.numeric(args[3])

## SNP.library.name <- "Hetaerina_titia_ddRAD_titia_dg"
## SNP.library.location  <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_v2/"
## select_N <- 3

dir.path <- paste0(SNP.library.location,SNP.library.name,"/")
dir.path.GPhocs <- paste0(dir.path, "G-Phocs/")
dir.create(dir.path.GPhocs)
dir.create(paste0(dir.path.GPhocs, "gphocs-loci"))

# SNP library
SNP.library.filter <- ".all.snps.NOGTDP10.MEANGTDP10_200.Q60"

vcf.SNPs <- read.vcfR(paste0(dir.path,SNP.library.name,SNP.library.filter,".vcf.gz"), verbose = F)
# Reorder samples so they are in alphabetical order
vcf.SNPs <- vcf.SNPs[samples = sort(colnames(vcf.SNPs@gt)[-1])]
# Rename samples
colnames(vcf.SNPs@gt) <- gsub("^.*/","",colnames(vcf.SNPs@gt))
colnames(vcf.SNPs@gt) <- substr(colnames(vcf.SNPs@gt), 1, nchar(colnames(vcf.SNPs@gt))-11)
colnames(vcf.SNPs@gt)[1] <- "FORMAT"
# Rename samples that end in _X
colnames(vcf.SNPs@gt) <- gsub(pattern = "_X", x = colnames(vcf.SNPs@gt), replacement = "")

## Read in populations assignment
pop_assign <- read.table("/home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/3_Results/LEA/LEA_pop_assign/popfile_all_noCUAJ_tabdelim.txt")
colnames(pop_assign) <- c("sample", "pop")

# Remove samples from vcf that are not in pop_assign
vcf.SNPs <- vcf.SNPs[samples = (colnames(vcf.SNPs@gt)[-1])[(colnames(vcf.SNPs@gt)[-1])%in%pop_assign$sample]]
# Remove samples in pop_assign that are not in vcf
pop_assign <- pop_assign[pop_assign$sample%in%(colnames(vcf.SNPs@gt)[-1]),]

any(!pop_assign$sample%in%(colnames(vcf.SNPs@gt)[-1]))
any(!(colnames(vcf.SNPs@gt)[-1])%in%pop_assign$sample)

# Reorder pop_asssign to be in same order as vcf
pop_assign <- pop_assign[match((colnames(vcf.SNPs@gt)[-1]), pop_assign$sample),]
# Check order
pop_assign$sample==(colnames(vcf.SNPs@gt)[-1])

# Calculate genotyep
# extract.gt(vcf.SNPs[1:10])
sample.miss <- apply(extract.gt(vcf.SNPs), MARGIN=2, function(x) sum(is.na(x)))

# Match the vcf data to the pop_assignment
names(sample.miss)==pop_assign$sample

## Get highest coverage N individuals for each population
pop_assign$cov <- sample.miss

top.cov.samples <- group_by(pop_assign, by = pop) %>%
  summarise(fst.sample = sample[order(cov)[1]],
            Snd.sample = sample[order(cov)[2]],
            thrd.sample = sample[order(cov)[3]],
            Frthsample = sample[order(cov)[4]])
top.cov.samples 
top.N.samples <- unlist(c(top.cov.samples[2:(select_N+1)]))

vcf.SNPs.N <- vcf.SNPs[samples = top.N.samples]

SNP.miss <- apply(extract.gt(vcf.SNPs.N), MARGIN=1, function(x) sum(is.na(x)))
# Remove sites that are not seen in all samples
vcf.SNPs.N <- vcf.SNPs.N[which(!SNP.miss>=1),]

if(grepl(SNP.library.name, pattern="titia_dg")){
  X_chrom <- "HetTit1.0.p.scaff-12-96647824"
}
# For H. americana lots of small X chromosomes
if(grepl(SNP.library.name, pattern="americana_dg")){
  X_chrom <- c("JAKTNV010000016.1","JAKTNV010000031.1","JAKTNV010000013.1","JAKTNV010000060.1","JAKTNV010000023.1","JAKTNV010000051.1","JAKTNV010000012.1")
}
# Remove X chrom
vcf.SNPs.N <- vcf.SNPs.N[!vcf.SNPs.N@fix[,"CHROM"]%in%X_chrom,]

############################
## G-phocs file creation
############################
loci.distance <- 1000
min.loci.length <- 100
max.loci <- 2000
min.cov <- 0.5
vcf.list <- list()

# Identifying individual loci
contig.dimensions <- dim(vcf.SNPs.N)
# Calcuates difference between SNP locations on contig and asked if greater than the set loci distance
contig.diff <- diff(as.numeric(vcf.SNPs.N@fix[,"POS"]))>loci.distance|diff(as.numeric(vcf.SNPs.N@fix[,"POS"]))<0
# Gets position of first loci
current.contig <- vcf.SNPs.N@fix[1,"CHROM"]
loci.position <- vcf.SNPs.N@fix[1,"POS"]
current.loci <- paste(vcf.SNPs.N@fix[1,1], loci.position, sep = ":")
# Creates storage area for loci start position
loci.list <- vector()
# Add the first loci to the first position
loci.list[[1]] <- paste(current.loci, 1, sep = "-")
# Prints name of first loci
print("Starting for loop that identifies each individual loci")

i<-3
# Loops through every SNP within the contig
for(i in 2:contig.dimensions[1]){
  # Checks if difference between this and the previous contig is less than the set loci distance
  # If yes adds the current name of the current loci
  if(!contig.diff[i-1]){loci.list[[i]] <- paste(current.loci, 
                                                as.numeric(vcf.SNPs.N@fix[i,"POS"])-as.numeric(loci.position)+1,
                                                sep = "-")}
  if(contig.diff[i-1]){
    #IF not updates the current loci to a new loci name
    loci.position <- vcf.SNPs.N@fix[i,"POS"]
    current.loci <- paste(vcf.SNPs.N@fix[i,1], loci.position, sep = ":")
    # print(paste("New.contig:", current.loci))
    loci.list[[i]] <- paste(current.loci, 
                            as.numeric(vcf.SNPs.N@fix[i,"POS"])-as.numeric(loci.position)+1,
                            sep = "-")
  }
}

print("Complete for loop that identifies each individual loci")
#Check loci poisiton is same length as SNP
length(vcf.SNPs.N@fix[,"POS"])==length(loci.list)

# Investage how new loci have been called in relation to the position of each bp
# Difference requires NA in last postion because there is now "next" bp to calculate distance from
loci.df <- data.frame(postion = vcf.SNPs.N@fix[,"POS"], bp.name = loci.list, difference = c(diff(as.numeric(vcf.SNPs.N@fix[,"POS"])), NA))

# Extract loci name and position within loci for each SNP for H. titia scaf names
if(grepl(SNP.library.name, pattern="titia_dg")){
  temp.loci.name <- do.call(rbind, strsplit(loci.df$bp.name, "-"))
  # Add loci name
  loci.df$loci.name <- do.call("c", lapply(1:dim(temp.loci.name)[1], function(x) paste(temp.loci.name[x,1:3] ,collapse="-")))
  loci.df$loci.position <- do.call("c", lapply(1:dim(temp.loci.name)[1], function(x) paste(temp.loci.name[x,4] ,collapse="-")))
}
# Extract loci name and position within loci for each SNP for H. americana scaf names
if(grepl(SNP.library.name, pattern="americana_dg")){
  temp.loci.name <- do.call(rbind, strsplit(loci.df$bp.name, "-"))
  loci.df$loci.name <- do.call("c", lapply(1:dim(temp.loci.name)[1], function(x) paste(temp.loci.name[x,1] ,collapse="-")))
  loci.df$loci.name <- gsub(loci.df$loci.name, pattern = ":", replacement = "-")
  loci.df$loci.position <- do.call("c", lapply(1:dim(temp.loci.name)[1], function(x) paste(temp.loci.name[x,2] ,collapse="-")))
}

# Number of loci identified
length(unique(loci.df$loci.name))

loci.stats <- data.frame(loci.name = unique(loci.df$loci.name))
# Length of loci
loci.stats$loci.length <- as.numeric(sapply(loci.stats$loci.name, FUN = function(x) max(as.numeric(loci.df$loci.position[loci.df$loci.name==x]))))
# Proportion of loci with called data
loci.stats$loci.called.bp <- sapply(loci.stats$loci.name, FUN = function(x) length(as.numeric(loci.df$loci.position[loci.df$loci.name==x])))
#Number of gaps in loci
loci.stats$loci.gaps <- sapply(loci.stats$loci.name, FUN = function(x) table(as.numeric(loci.df$difference[loci.df$loci.name==x])>1&loci.df$difference[loci.df$loci.name==x]<loci.distance)[2])
loci.stats$loci.gaps[is.na(loci.stats$loci.gaps)] <- 0
# Largest gap within loci (except final distance with is distance between this and the next loci)
# This should identify, if presence the a gap between the forward and the reverse read
# Will return -INF when loci is 1 base long
loci.stats$loci.largest.gap <- sapply(loci.stats$loci.name, FUN = function(x) max(as.numeric(loci.df$difference[loci.df$loci.name==x])[-length(loci.df$difference[loci.df$loci.name==x])]))
loci.stats$loci.largest.gap[is.infinite(loci.stats$loci.largest.gap)] <- 0
#Distance between loci and next
loci.stats$dist.next.loci <- sapply(loci.stats$loci.name, FUN = function(x) as.numeric(loci.df$difference[loci.df$loci.name==x])[length(loci.df$difference[loci.df$loci.name==x])])
loci.stats$loci.site.cov <- loci.stats$loci.called.bp/loci.stats$loci.length
#summary of loci stats
summary(loci.stats)
# Number of loci above min length
print("Loci over min length")
table(loci.stats$loci.length>min.loci.length)
# Number of loci with above X percentage bases called
print("Loci with greater than 50% missing data")
table((loci.stats$loci.site.cov)>min.cov)
# Number of loci with both criterea
table((loci.stats$loci.site.cov>min.cov)&(loci.stats$loci.length>min.loci.length))

print("Printed loci stats")
# Update ID column
# Extract loci name and position within loci for each SNP for H. titia scaf names
if(grepl(SNP.library.name, pattern="titia_dg")){
  vcf.SNPs.N@fix[,"ID"] <- loci.list
  loci.table <- do.call(rbind, strsplit(vcf.SNPs.N@fix[,"ID"], "-"))
  loci.table <- cbind(do.call("c", lapply(1:dim(loci.table)[1], function(x) paste(loci.table[x,1:3],collapse="-"))), loci.table[,4])
}
# Extract loci name and position within loci for each SNP for H. americana scaf names
if(grepl(SNP.library.name, pattern="americana_dg")){
  vcf.SNPs.N@fix[,"ID"] <- loci.list
  loci.table <- do.call(rbind, strsplit(vcf.SNPs.N@fix[,"ID"], "-"))
  loci.table[,1] <- gsub(loci.table[,1], pattern = ":", replacement="-")
}

# Removes loci that have less than min.loci.length called and missing more than 50% data
valid.loci <- loci.stats$loci.name[which((loci.stats$loci.site.cov>min.cov)&(loci.stats$loci.length>min.loci.length))]
vcf.loci <- vcf.SNPs.N[loci.table[,1]%in%valid.loci,]

# Extract loci name and position within loci for each SNP for H. titia scaf names
if(grepl(SNP.library.name, pattern="titia_dg")){
  loci.table.valid <- do.call(rbind, strsplit(vcf.loci@fix[,"ID"], "-"))
  loci.table.valid <- cbind(do.call("c", lapply(1:dim(loci.table.valid)[1], function(x) paste(loci.table.valid[x,1:3],collapse="-"))), loci.table.valid[,4])
}
# Extract loci name and position within loci for each SNP for H. americana scaf names
if(grepl(SNP.library.name, pattern="americana_dg")){
  loci.table.valid <- do.call(rbind, strsplit(vcf.loci@fix[,"ID"], "-"))
  loci.table.valid[,1] <- gsub(loci.table.valid[,1], pattern = ":", replacement="-")
}


p <- ggplot(data.frame(position = vcf.loci@fix[,"POS"], loci = loci.table.valid[,1], chrom=vcf.loci@fix[,"CHROM"])) +
  geom_point(aes(as.numeric(position), as.numeric(position), col = loci), show.legend = F) +
  facet_wrap(~chrom)

ggsave(p, filename=paste0(dir.path.GPhocs,"g-phocs-loci-check.png"))
# Get names of all loci
unique.loci <- unique(loci.table.valid[,1])

if(length(unique.loci) > max.loci){
  print("More than max number of loci: randomly selecting loci to reduce number")
  unique.loci <- sample(unique.loci, max.loci)
}
loci <- unique.loci[1]
for(loci in unique.loci){
  loci.temp <- vcf.loci[loci.table.valid[,1]==loci]
  diff(as.numeric(loci.temp@fix[,"POS"]))
  # Check how much missing 
  loci.text <- gsub(loci, pattern = ":", replacement = "-")
  gt <- extract.gt(loci.temp, return.alleles = TRUE, convertNA = TRUE)
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
  gt[is.na(gt)] <- "N"
  gt[gt=="."] <- "N"
  
  # Combine each loci into Gphocs format
  g_phocs <- c(paste(loci.text, dim(gt)[2], dim(gt)[1], sep = " "),
               do.call("c",lapply(1:dim(gt)[2], function(x) paste(rownames(t(gt))[x], paste(t(gt)[x,], collapse = ""), collapse = "\t"))),
               "")
  #Save as individual file
  writeLines(g_phocs, paste0(dir.path.GPhocs, "gphocs-loci/",loci.text,".gphocs"))
  
}
print("Completed for loop to create individual loci files")
# Write top line
writeLines(c(as.character(length(unique.loci)), ""), paste0(dir.path.GPhocs,"gphocs-loci/Gphocs_header.txt"))

## Create df of samples for use in generating contig files
pop_assign_N <- pop_assign[pop_assign$sample%in%top.N.samples,]

write.table(pop_assign_N[,1:2], file = paste0(dir.path.GPhocs,"GPhocs_samples_N",select_N,".txt"), col.names=F, row.names = F, quote=F)


