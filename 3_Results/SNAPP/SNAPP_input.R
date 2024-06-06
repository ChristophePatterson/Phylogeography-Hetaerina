# Get system arguments
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")

library(vcfR)
library(LEA)
library(poppr)
library(tidyverse)
library(ape)

# Load in arguments
args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[1]
# SNP.library.name <- "Hetaerina_all_ddRAD_americana_dg"
SNP.library.location <- args[2]
# SNP.library.location  <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_v2/"
select_N <- as.numeric(args[3])
# select_N <- 2

# Read in vcf
dir.path <- paste0(SNP.library.location,"/",SNP.library.name,"/")
dir.path.SNAPP <- paste0(dir.path, "SNAPP/")
dir.create(dir.path.SNAPP)
filter_para <- ".all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.biSNP0_20.noX.noCUAJ"

vcf.SNPs <- read.vcfR(paste0(dir.path,SNP.library.name,filter_para,".vcf.gz"),verbose = F)
# Reorder samples so they are in alphabetical order
vcf.SNPs <- vcf.SNPs[samples = sort(colnames(vcf.SNPs@gt)[-1])]

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

# Removing dots from SNP names because genind doesn't like this
vcf.SNPs@fix[,1] <- gsub(vcf.SNPs@fix[,1], pattern = "\\.", replacement = "_")

# Get Genind
my_genind_ti_SNPs <- vcfR2genind(vcf.SNPs, sep = "/", return.alleles = TRUE)
# Proporation of samples with missing data
sample.miss <- propTyped(my_genind_ti_SNPs , by = "ind")

# Match the vcf data to the pop_assignment
names(sample.miss)==pop_assign$sample

## Get highest coverage N individuals for each population
pop_assign$cov <- sample.miss

# Number of samples that are assigned to each group
print("Number of samples that are assigned to each group")
table(pop_assign$pop)
min(table(pop_assign$pop))

top.cov.samples <- group_by(pop_assign, by = pop) %>%
  summarise(fst.sample = sample[order(cov)[1]],
            Snd.sample = sample[order(cov)[2]],
            thrd.sample = sample[order(cov)[3]],
            Frthsample = sample[order(cov)[4]],
            Fithsample = sample[order(cov)[5]],
            sixsample = sample[order(cov)[6]])
top.cov.samples

# Get samples to N for each pop
top_N_samples <- unlist(c(top.cov.samples[,2:(select_N+1)]))
# pop_assign
pop_assign_N <- pop_assign[pop_assign$sample%in%top_N_samples,]
vcf.SNPs.N <- vcf.SNPs[samples = top_N_samples]
# Ensure order is correct
pop_assign_N <- pop_assign_N[match(colnames(vcf.SNPs.N@gt[,-1]), pop_assign_N$sample),]
pop_assign_N$sample%in%colnames(vcf.SNPs.N@gt[,-1])

# Check if SNPs are polymorphic
vcf.SNPs.N <- vcf.SNPs.N[(which(is.polymorphic(vcf.SNPs.N, na.omit = T))),]

# Extract genotype
gt <- extract.gt(vcf.SNPs.N, return.alleles = TRUE, convertNA = TRUE)

colnames(gt)==pop_assign_N$sample
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

# Rename file with population name
colnames(gt) <- paste(pop_assign_N$pop, pop_assign_N$sample, sep = "_")
#Write nexus
write.nexus.data(t(gt), file = paste0(dir.path.SNAPP, SNP.library.name,"_ind_",select_N,".nex"), 
                 format = "DNA", missing = "?", interleaved = F)

## Set up files following https://github.com/mmatschiner/tutorials/blob/master/divergence_time_estimation_with_snp_data/README.md
## (1) Phylip file
colnames(gt) <- pop_assign_N$sample

ape::write.dna(t(gt), file = paste0(dir.path.SNAPP, SNP.library.name,"_ind_",select_N,".phy"),
               format ="interleaved", nbcol = -1, colsep = "")
## (2) Species table
species.df <- data.frame(species = pop_assign_N$pop, sample = pop_assign_N$sample)
write.table(species.df, file = paste0(dir.path.SNAPP, SNP.library.name,"_ind_",select_N,".txt"),
            row.names = F,quote = F, sep = "\t")
# (3) Theta priors from Stranding et al 2022
# All Hetaerina samples
# normal(offset,mean,sigma)
contrant.df <- data.frame(x = "normal(0,33.08,5.53)", y = "crown", 
                          z = paste0(unique(pop_assign_N$pop), sep = ",",collapse = ""))
# Americana/calverti
contrant.df[2,] <- data.frame(x = "normal(0,3.76,1.87)", y = "crown", 
                              z = paste0(unique(pop_assign_N$pop[grepl(pop_assign_N$pop, pattern = "americana|calverti")]), sep = ",",collapse = ""))
#remove trailing comma
contrant.df$z <- substr(contrant.df$z, start = 1 , stop = nchar(contrant.df$z)-1)

write.table(contrant.df, file = paste0(dir.path.SNAPP, SNP.library.name,"_ind_",select_N,".con.txt"),
            row.names = F, col.names = F, quote = F, sep = "\t")
