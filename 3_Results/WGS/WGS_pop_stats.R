.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
library(vcfR)
library(LEA)
library(ggplot2)
library(adegenet)
library(poppr)
library(hierfstat)

# Get system arguments
args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[1]
# SNP.library.name <- "WGS_titia"
# input_dir <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/WGS_titia_r10000/"
input_dir <- args[2]
output_dir <- paste0(input_dir,"/results/")
plot.dir <- paste0("/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/4_Manuscript/plots/WGS_titia/",SNP.library.name,"/")

ncores <- args[3]

dir.create(output_dir)
dir.create(plot.dir)
vcf <- read.vcfR(paste0(input_dir ,SNP.library.name,"_chr1-12_filtered.vcf.gz"))

colnames(vcf@gt)
sample_map <- colnames(vcf@gt)[-1]

vcf@gt[1:10,]
vcf.bi <- vcf[is.biallelic(vcf)]
vcf.ply <- vcf.bi[is.polymorphic(vcf.bi, na.omit = T)]
## vcf.ply@gt[1:10,1:6]
# Convert to geno
X_chrom <- "HetTit1.0.p.scaff-12-96647824"

# Removes . from varient names as genind does like that
vcf.ply@fix[,1] <- gsub(vcf.ply@fix[,1], pattern = "\\.", replacement = "_")

# Creaete genind object
my_genind_ti_SNPs <- vcfR2genind(vcf.ply, sep = "/", return.alleles = TRUE)

# Assign populations
# Two different versions being considered
rownames(my_genind_ti_SNPs@tab)
my_genind_ti_SNPs@pop <- as.factor(substr(rownames(my_genind_ti_SNPs@tab), 1, 4))

my_genind_ti_SNPs@pop
#my_genind_ti_SNPs@pop <- as.factor(paste(pca.q.df$assign, pca.q.df$Country_Ocean, sep = "_"))

# Calcultes Fst values between populations
pop.diff <- genet.dist(my_genind_ti_SNPs, diploid = T, method = "WC84")
pop.diff

my_stats <- basic.stats(my_genind_ti_SNPs)
my_stats$overall

my_stats

write.table(pop.diff, paste0("Fst_table_", SNP.library.name, ".txt"))

pdf(paste0(SNP.library.name, "basic_stats_plot.pdf"))
boxplot(my_stats$Ho)
boxplot(my_stats$Hs)
boxplot(my_stats$perloc)
dev.off()