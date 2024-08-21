## Whole genome sequence data combined with ddRAD data
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
library(vcfR)
library(LEA)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(introgress)
library(patchwork)
library(adegenet)
library(hierfstat)
library(ape)

args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[1]
# SNP.library.name <- "WGS_titia"
# input_dir <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/WGS_titia_r10000/"
input_dir <- args[2]
ncores <- args[3]
output_dir <- paste0(input_dir,"/results/")
plot.dir <- paste0("/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/4_Manuscript/plots/WGS_titia/",SNP.library.name,"/")


# Output file location
SNP.library.name <- "WGS_ddRAD_Pacific_titia"
species <- "titia"
dir.path <- output_dir
filter_para <- ".rand10000"
snp_sub_text <- "0_20"
snp_sub <- 0.20
# Plot output file location
dir.create(plot.dir)

## Running analysis name
analysis.name <- SNP.library.name

vcf <- read.vcfR(paste0(input_dir, "/", SNP.library.name, filter_para ,".vcf.gz"))

colnames(vcf@gt)
# Reformat names
gsub("^.*/","",colnames(vcf@gt))
colnames(vcf@gt) <- gsub("^.*/","",colnames(vcf@gt))
colnames(vcf@gt) <- substr(colnames(vcf@gt), 1, nchar(colnames(vcf@gt))-4)
colnames(vcf@gt) <- gsub("\\.paired","",colnames(vcf@gt))
colnames(vcf@gt)[1] <- "FORMAT"
colnames(vcf@gt)

Erandi_sample_names <- cbind(c("ZANAa05_R", "CUAJa03_R","CUAJb01_R","MIXTa01_R"), c("E008","E018","E019","E020"))
colnames(vcf@gt)[(colnames(vcf@gt)%in%Erandi_sample_names[,2])] <- Erandi_sample_names[,1][match(colnames(vcf@gt),Erandi_sample_names[,2])[(colnames(vcf@gt)%in%Erandi_sample_names[,2])]]

## Which samples are duplicated
dups_ddRAD <- colnames(vcf@gt)[!colnames(vcf@gt)%in%gsub("_R", "", colnames(vcf@gt)[grep("_R", colnames(vcf@gt))])]
vcf <- vcf[samples = dups_ddRAD]
# Rename WGS samples
colnames(vcf@gt) <- gsub("_R","",colnames(vcf@gt))

# Removes . from varient names as genind does like that
vcf@fix[,1] <- gsub(vcf@fix[,1], pattern = "\\.", replacement = "_")

# Removes SNPs on the X chromosome
X_chrom <- "HetTit1_0_p_scaff-12-96647824"
table(vcf@fix[,1]!=X_chrom)
vcf.SNPs <- vcf[vcf@fix[,1]!=X_chrom,]

# Stats
myMiss <- apply(extract.gt(vcf.SNPs), MARGIN = 2, function(x){ sum(is.na(x)) })
myMiss <- data.frame(sample = names(myMiss), per.gt = myMiss/nrow(vcf)*100)

summary(myMiss)
ggplot(myMiss) +
  geom_violin(aes(x = "X", y = per.gt), outlier.colour = NA) +
  geom_jitter(aes(x = "X", y = per.gt, col = grepl("_R", sample)), height = 0, width = 0.2)

# Remove low coverage samples
myMiss[order(myMiss$per.gt),]
vcf.SNPs <- vcf.SNPs[samples = c("FORMAT",myMiss$sample[myMiss$per.gt<=(snp_sub*100)])]

# Remove nonpolymorphic SNPsvcf@gt[1:10,]
vcf.SNPs <- vcf.SNPs[is.biallelic(vcf.SNPs)]
vcf.SNPs <- vcf.SNPs[is.polymorphic(vcf.SNPs, na.omit = T)]

# Get coverage statistics 
DP.mat <- extract.gt(vcf.SNPs, element = "DP")
Samp.DP <- apply(DP.mat, MARGIN = 2, function(x) mean(as.numeric(x), na.rm = T))
SNP.DP <- apply(DP.mat, MARGIN = 1, function(x) mean(as.numeric(x), na.rm = T))

c(mean(Samp.DP,na.rm =T), median(Samp.DP,na.rm =T), range(Samp.DP,na.rm = T))
c(mean(SNP.DP,na.rm =T), median(SNP.DP,na.rm =T), range(SNP.DP,na.rm = T))

## Write out vcf
write.vcf(vcf.SNPs, file = paste0(input_dir, "/", SNP.library.name, filter_para,"_filtered.noX.vcf.gz"),mask = F,APPEND = F)

# Creaete genind object
my_genind_ti_SNPs <- vcfR2genind(vcf.SNPs, sep = "/", return.alleles = TRUE)
sample_map <- colnames(vcf.SNPs@gt)[-1]
# Create geno object
geno <- vcfR2loci(vcf.SNPs, return.alleles = F)
geno.mat <- as.matrix(geno)

# Convert into geno format 2 for alternative homozygous, 1 for he
# 2 for alternative homozygous
geno.mat[geno.mat=="1/1"] <- 2
# 1 for heterozygous
geno.mat[geno.mat=="0/1"] <- 1
# 0 for homozygous genome assembly allele
geno.mat[geno.mat=="0/0"] <- 0
# Checks if any data points are entirely heterozygous
is.only.het <- apply(geno.mat, MARGIN = 2, function(x) gsub(paste0(unique(x), collapse = ""), pattern = "NA",replacement = "")=="1")
if(any(is.only.het)){geno.mat <- geno.mat[,-which(is.only.het)]}
# Missing data coded as 9
geno.mat[is.na(geno.mat)] <- 9

# Convert to data frame to save
geno.df <- data.frame(t(geno.mat))

paste0(output_dir, SNP.library.name, ".geno")
write.table(x = geno.df, file = paste0(output_dir, SNP.library.name, ".geno"),
            col.names = F, row.names = F, quote = F, sep = "")

#Read geno object
geno <- read.geno(paste0(output_dir, SNP.library.name, ".geno"))
dim(geno)

# Get sample information
sample_map<- read.csv("/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/All samples held in Durham_v17.csv", check.names=F)
names(sample_map)

#Extracting sample information and linking it to tree labels
samples <- rownames(my_genind_ti_SNPs@tab)
sites <- data.frame(samples)
sites$samples

sites <- merge(sites, sample_map[,c("Unique.ID","sex","species","Year","Site.ID","Lat","Long","Country","Province","Ocean.drainage")], by.x = "samples", by.y = "Unique.ID", all.x = T)
any(is.na(sites))
sites <- sites[match(rownames(my_genind_ti_SNPs@tab), sites$samples),]

#Conduct PCA
geno2lfmm(paste0(output_dir, SNP.library.name, ".geno"),
          paste0(output_dir, SNP.library.name, ".geno.lfmm"))

#PCA
pc <- pca(paste0(output_dir, SNP.library.name, ".geno.lfmm"), scale = TRUE)

pc.sum <- summary(pc)

pca.data <- data.frame(pc$projections)
colnames(pca.data) <- paste0("pca", 1:dim(pc$projections)[2])
pca.labs <- paste("pca", 1:dim(pc$projections)[2], " (",round(as.numeric(pc.sum[2,1:dim(pc$projections)[2]])*100, 1), "%)", sep = "")
pca.data$samples <- rownames(my_genind_ti_SNPs@tab)

## Conduct snmf
max.K <- 8
obj.at <- snmf(paste0(output_dir, SNP.library.name, ".geno"), K = 1:max.K, ploidy = 2, entropy = T,
               CPU = as.numeric(ncores), project = "new", repetitions = 20, alpha = 100)
titia.snmf <- load.snmfProject(file = paste0(output_dir, SNP.library.name, ".snmfProject"))
titia.snmf.sum <- summary(titia.snmf)

# Cross entropy plot
ce <- cbind(1:max.K, t(titia.snmf.sum$crossEntropy))
colnames(ce) <- c("K", "min","mean","max")
ce <- data.frame(ce)

summary(titia.snmf)

which.min(ce$mean)
ce.plot <- ggplot(ce) +
  geom_point(aes(K, mean), size = 2, shape = 19) +
  geom_errorbar(aes(x = K, ymin=min, ymax=max), width=.2,
                position=position_dodge(.9)) +
  ggtitle(paste("Minimum mean cross-entropy value K =", which.min(ce$mean))) +
  ylab("Cross-entropy") +
  theme_bw()

ggsave(paste0(plot.dir,SNP.library.name,"_cross_entropy.pdf"), plot = ce.plot)

# Choose best K
K = 3
best <- which.min(cross.entropy(titia.snmf, K = K))
qmatrix = Q(titia.snmf, K = K, run = best)

qtable <- cbind(rep(sites$samples,K), rep(1:K, each = length(sites$samples)), c(qmatrix[,1:K]))
qtable <-  data.frame(qtable)
colnames(qtable) <- c("sample","Qid", "Q")
qtable

## Pac.clust <- which.max(qtable$Q[qtable$sample=="ZANAa05"])
## SAtl.clust <- which.min(qtable$Q[qtable$sample=="ZANAa05"])
## NAtl.clust <- 3
# Cols
het.cols <- c("#AF0F09","#E5D9BA","#3E3C3A")
## het.cols[Pac.clust] <- "#AF0F09"
## het.cols[SAtl.clust] <- "#E5D9BA"
## het.cols[NAtl.clust] <- "#3E3C3A"

# Reorder bar plot
qtable$sample <-  factor(qtable$sample, levels = sites$sample[order(sites$Long)])

p.bar <- ggplot(qtable)+
  geom_bar(stat="identity", aes(sample, as.numeric(Q), fill = Qid,), position = "stack", width = 1, col = "black", show.legend = F) +
  scale_fill_manual(values = het.cols) +
  ylab(paste0("Ancestry assignment (K = ",K,")")) +
  xlab("Sample") +
  theme_bw() +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") 
# theme(axis.text.x = element_text(size=6))
ggsave(paste0(plot.dir, SNP.library.name, "snmf_K", K, "_bar__plot.pdf"), p.bar)

pca.data$assign <- apply(pca.data,MARGIN=1, function(x) which.max(qtable$Q[which(qtable$sample==x["samples"])]))
pca.data
## PCA plot
p <- ggplot(pca.data) +
  geom_point(aes(as.numeric(pca1), as.numeric(pca2), fill=as.factor(assign)), size = 3, shape = 21, show.legend=F) +
  geom_label_repel(data = pca.data[pca.data$sample=="CUAJa03",], aes(as.numeric(pca1), as.numeric(pca2), label = samples),
                     size = 4, nudge_x = 20, nudge_y = -5, min.segment.length = 0) +
  scale_fill_manual(values = het.cols) +
  xlab(pca.labs[1]) +
  ylab(pca.labs[2]) +
  theme_bw() +
  theme(text = element_text(size = 20))

ggsave(paste0(plot.dir, "PCA_",SNP.library.name,"_plot.pdf"), p)

ggsave(paste0(plot.dir, "PCA_sNMF_",SNP.library.name,"_plot.pdf"), plot = (p / p.bar), width = 10, height = 10)

### Convert to RAxML and SVDQuartets
# Extract genotypes
gt <- extract.gt(vcf.SNPs, return.alleles = TRUE, convertNA = TRUE)

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
write.nexus.data(t(gt), file = paste0(output_dir, "SVDQuartets/", SNP.library.name,".nex"), 
                 format = "DNA", missing = "?", interleaved = FALSE)

# Run without heterozgousity
gt <- extract.gt(vcf.SNPs, return.alleles = TRUE, convertNA = TRUE)
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

gt[435:445,]
table(gt)

# Write .phy file for RAxML
ape::write.dna(t(gt), file = paste0(output_dir, "RAxML/", SNP.library.name,".phy"), format ="interleaved")


# Assign populations
# Two different versions being considered
rownames(my_genind_ti_SNPs@tab)==sites$samples
sites$site.sub <- paste(gsub('[[:digit:]]+', '', sites$Site.ID))
sites$site.sub
my_genind_ti_SNPs@pop <- as.factor(sites$site.sub)

# Calcultes Fst values between populations
print("Starting Fst calculation")
pop.diff <- genet.dist(my_genind_ti_SNPs, diploid = T, method = "WC84")
pop.diff

saveRDS(pop.diff, file = paste0(output_dir, "Fst_table_", SNP.library.name, ".rds"))
write.table(as.matrix(pop.diff), paste0(output_dir, "Fst_table_", SNP.library.name, ".txt"))

my_stats <- basic.stats(my_genind_ti_SNPs)
my_stats$overall

my_stats

pdf(paste0(output_dir, SNP.library.name, "basic_stats_plot.pdf"))
boxplot(my_stats$Ho)
boxplot(my_stats$Hs)
boxplot(my_stats$perloc)
dev.off()
