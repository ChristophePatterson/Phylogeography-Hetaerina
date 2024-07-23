# Calculation of PCR on titia populations
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")

library(ggplot2)
library(vcfR)
#install.packages("rlang")
# install.packages("rgdal")
# install.packages("BiocManager")
# BiocManager::install("LEA")
library(LEA)

args <- commandArgs(trailingOnly = TRUE)

select_N <- as.numeric(args[1])
input_directory<-args[2]
# input_directory<-"/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/ipyrad/ipyrad_N2"
# select_N <- as.numeric("2")


# Read in vcf
vcf_file <- paste0(input_directory,"/denovo_N",select_N,"_outfiles/denovo_N",select_N)
vcf <- read.vcfR(paste0(vcf_file,".vcf"))

#number of unique RAD loci
# Stats
myMiss <- apply(extract.gt(vcf), MARGIN = 1, function(x){ sum(is.na(x)) })
myMiss <- data.frame(sample = names(myMiss), per.gt = myMiss)

vcf <- vcf[which(myMiss$per.gt<4),]

# filter to remove polybilleic alleles
vcf <- vcf[is.biallelic(vcf),]
vcf <- vcf[is.polymorphic(vcf, na.omit = T),]

# subset snp to one SNP per loci
SNP_per_RAD <- do.call("c",lapply(unique(vcf@fix[,1]), FUN = function(x){ 
  sample(which(vcf@fix[,1]==x), size = 1)
}))

# Filer to top coverage samples
vcf.SNPs <- vcf[SNP_per_RAD,]

# Creaete genind object
my_genind_ti_SNPs <- vcfR2genind(vcf.SNPs, sep = "/", return.alleles = TRUE)

# Create geno object
geno <- vcfR2loci(vcf.SNPs, return.alleles = F)
geno.mat <- as.matrix(geno)
table(geno.mat)


# Convert into geno format 2 for alternative homozygous, 1 for he
# 2 for alternative homozygous
geno.mat[geno.mat=="1/1"] <- 2
# 1 for heterozygous
geno.mat[geno.mat=="0/1"] <- 1
geno.mat[geno.mat=="1/0"] <- 1
# 0 for homozygous genome assembly allele
geno.mat[geno.mat=="0/0"] <- 0
# Checks if any data points are entirely heterozygous
is.only.het <- apply(geno.mat, MARGIN = 2, function(x) gsub(paste0(unique(x), collapse = ""), pattern = "NA",replacement = "")=="1")
if(any(is.only.het)){geno.mat <- geno.mat[,-which(is.only.het)]}
# Missing data coded as 9
geno.mat[is.na(geno.mat)] <- 9
table(geno.mat)

# Convert to data frame to save
geno.df <- data.frame(t(geno.mat))

# Write out geno file
write.table(x = geno.df, file = paste0(vcf_file,".geno"),
            col.names = F, row.names = F, quote = F, sep = "")

# Write out vcf 
write.vcf(vcf.SNPs, file = paste0(vcf_file,"_filter_1SNPperRAD.vcf.gz"),mask = F,APPEND = F)

## Write out SNAPP

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

# Extract genotype
gt <- extract.gt(vcf.SNPs, return.alleles = TRUE, convertNA = TRUE)

colnames(gt)==pop_assign$sample
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
ape::write.dna(t(gt), file = paste0(vcf_file,"_filter_1SNPperRAD.phy"),
               format ="interleaved", nbcol = -1, colsep = "")
## (2) Species table
species.df <- data.frame(species = pop_assign$pop, sample = pop_assign$sample)
write.table(species.df, file = paste0(vcf_file,"_filter_1SNPperRAD.txt"),
            row.names = F,quote = F, sep = "\t")
# (3) Theta priors from Stranding et al 2022
# All Hetaerina samples
# normal(offset,mean,sigma)
contrant.df <- data.frame(x = "normal(0,33.08,5.53)", y = "crown", 
                          z = paste0(unique(pop_assign$pop), sep = ",",collapse = ""))
# Americana/calverti
contrant.df[2,] <- data.frame(x = "normal(0,3.76,1.87)", y = "crown", 
                              z = paste0(unique(pop_assign$pop[grepl(pop_assign$pop, pattern = "americana|calverti")]), sep = ",",collapse = ""))
#remove trailing comma
contrant.df$z <- substr(contrant.df$z, start = 1 , stop = nchar(contrant.df$z)-1)

write.table(contrant.df, file = paste0(vcf_file,"_filter_1SNPperRAD.con.txt"),
            row.names = F, col.names = F, quote = F, sep = "\t")

