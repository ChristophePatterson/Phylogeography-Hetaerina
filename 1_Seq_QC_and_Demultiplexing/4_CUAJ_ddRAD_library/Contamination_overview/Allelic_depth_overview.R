.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
library(vcfR)
library(ggplot2)
library(patchwork)

# Get system arguments
args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[1]
# SNP.library.name <- "CUAJ_Masters_library"
SNP.library.location <- args[2]
# SNP.library.location <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/CUAJ_Masters_library"
plot.dir <- args[3]
# plot.dir <- "/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/1_Seq_QC_and_Demultiplexing/4_CUAJ_ddRAD_library/Contamination_overview/plots"

filter_para <- ".all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000"

## Running analysis name
analysis.name <- "CUAJ_contamination_check"

vcf <- read.vcfR(paste0(SNP.library.location,"/",SNP.library.name,filter_para,".vcf.gz"))# vcf <- read.vcfR("data/bcftools/H_titia_190samp_Allscaf_snps.NOGTDP3.MEANGTDP3.Q20.SAMP160.MAF005.vcf.gz")

# Remames columns
gsub("^.*/","",colnames(vcf@gt))
colnames(vcf@gt) <- gsub("^.*/","",colnames(vcf@gt))
colnames(vcf@gt) <- gsub(".paired.bam","",colnames(vcf@gt))
colnames(vcf@gt) <- gsub(".bam","",colnames(vcf@gt))
colnames(vcf@gt)[1] <- "FORMAT"
#Limit SNPs to those that are biallelic
colnames(vcf@gt) <- gsub(pattern = "_X", x = colnames(vcf@gt), replacement = "")

# Extract allele depth
AD <- extract.gt(vcf, element = "AD")

# Extract the AD depth for the reference 
AD.0 <- apply(AD, MARGIN = 2, function(x) as.numeric(stringr::str_split_i(x, pattern = ",", i =1)))
# Extract the AD depth for the alternative 
AD.1 <- apply(AD, MARGIN = 2, function(x) as.numeric(stringr::str_split_i(x, pattern = ",", i =2)))

# What is the ratio of alternative to reference
AD.diff <- data.frame(AD.1/(AD.0+AD.1))

# Melt into single dataframe
AD.diff.melt <- reshape2::melt(AD.diff, value.name = "AD.ratio")
AD.diff.melt$sample.simple <- as.character(AD.diff.melt$variable)
AD.diff.melt$sample.simple[!grepl(AD.diff.melt$variable, pattern = "CUAJa02|CUAJa03")] <- "Other"


p <- ggplot(AD.diff.melt, aes(x = AD.ratio, group = variable, col = sample.simple,
                         )) +
  geom_freqpoly(show.legend = T, bins = 20, linewidth = 1) +
  xlim(c(0,1))

# Extract the AD depth for the reference 
AD.df <- lapply(X = colnames(vcf@gt)[-1], function(x){
  ID.sample <- x
  AD.0.CUAJa03 <- as.numeric(stringr::str_split_i(AD[,ID.sample], pattern = ",", i =1))
  # Extract the AD depth for the alternative 
  AD.1.CUAJa03 <- as.numeric(stringr::str_split_i(AD[,ID.sample], pattern = ",", i =2))
  
  AD.CUAJ <- data.frame(sample = ID.sample, SNP = rownames(AD), gt = extract.gt(vcf[sample = ID.sample])[,1], 
                        AD.0 = AD.0.CUAJa03, AD.1 = AD.1.CUAJa03,
                        AD.0.ratio = AD.0.CUAJa03/(AD.0.CUAJa03+AD.1.CUAJa03),
                        AD.1.ratio = AD.1.CUAJa03/(AD.0.CUAJa03+AD.1.CUAJa03))
  
  AD.CUAJ$het <- is.het(matrix(AD.CUAJ$gt),na_is_false = FALSE)
  return(AD.CUAJ)
  
}

)

AD.df <- do.call("rbind", AD.df)

# Label samples dependent of whether they came from the containminated pool or not
dodgy_lib <- c("CUAJb19", "CUAJb18", "CUAJb01", "CUAJb06", "CUAJb21", "CUAJb02", "CUAJb07", "CUAJb13", 
            "CUAJb15", "CUAJa03", "CUAJb14", "CUAJb12", "CUAJb16", "CUAJb10", "CUAJb08", "CUAJb03", 
            "CUAJb20", "CUAJb09", "CUAJb17", "CUAJb11")
AD.df$dodgy <- AD.df$sample%in%dodgy_lib
    

p <- ggplot(AD.df) +
  geom_violin(aes(AD.0/(AD.0+AD.1), y = sample, col = het, fill = het),
              bw=0.01, position="identity", alpha = 0.5, scale = "width")
        
ggsave(paste0(plot.dir,"/",SNP.library.name,"_AD_violin.pdf"), p, 
       height = 40, width = 5)

q1 <- ggplot(AD.df[AD.df$dodgy,]) +
  geom_point(aes((AD.0),(AD.1), col = het)) +
  facet_wrap(~sample, scales = "free") +
  xlim(c(0,100)) +
  ylim(c(0,100)) +
  ggtitle("Samples from contaminated ddRAD pool") +
  theme(strip.background = element_rect(fill="#fe5757"))

q2 <- ggplot(AD.df[!AD.df$dodgy,]) +
  geom_point(aes((AD.0),(AD.1), col = het)) +
  facet_wrap(~sample, scales = "free") +
  xlim(c(0,100)) +
  ylim(c(0,100)) +
  ggtitle("Selection of samples from other ddRAD pool") +
  theme(strip.background = element_rect(fill="lightgreen"))

ggsave(paste0(plot.dir,"/",SNP.library.name,"_AD_scatter.pdf"), (q1 + q2), 
    height = 20, width = 25)



