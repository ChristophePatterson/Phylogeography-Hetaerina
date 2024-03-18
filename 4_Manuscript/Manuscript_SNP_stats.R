## Quick SNP stats for all libraries

# Names of all libries
lib.names <- c("Hetaerina_all_ddRAD_americana_dg","Hetaerina_all_ddRAD_titia_dg",
               "Hetaerina_americana_ddRAD_americana_dg","Hetaerina_americana_ddRAD_titia_dg",
               "Hetaerina_titia_ddRAD_americana_dg", "Hetaerina_titia_ddRAD_titia_dg")

# Extract info for SNP filtering
SNP.stats <- list()
## Read vcf stats
for(i in 1:length(lib.names)){
  SNP.library.name <- lib.names[i]
  SNP.stats.temp <- readLines(paste0("2_SNP_calling/library_combinations/",SNP.library.name, ".vcf.stats.txt"))
  SNP.stats[[i]] <- c(SNP.library.name, SNP.stats.temp)
}
# Combine into single data frame
SNP.stats <- as.data.frame(do.call("rbind", SNP.stats))
# Names
names(SNP.stats) <- c("library", "all_base_calls", "N.samples", "single_site_calls", "D10_100_Q20", "Samp0_8", "MAF2", "r1000")

# Extract SNP  and sample data from the final SNP datasets
dir.path <- paste0("4_Manuscript/data/SNP_libraries/",SNP.library.name,"/")
filter_para <- ".all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.biSNP0_20"

list.files(results.dir)
for(i in 1:length(SNP.stats$library)){
  SNP.library.name <- SNP.stats$library[i]
  dir.path <- paste0("4_Manuscript/data/SNP_libraries/",SNP.library.name,"/")
  vcf.temp <- read.vcfR(paste0(dir.path,SNP.library.name,filter_para,".vcf.gz"))
  gt <- extract.gt(vcf.temp)
  SNP.stats$N.samples.end[i] <- dim(gt)[2]
  SNP.stats$N.SNPs[i] <- dim(gt)[1]
  # Stats
  myMiss <- apply(gt, MARGIN = 2, function(x){ sum(is.na(x)) })
  SNP.stats$Sample.mean.miss[i] <- mean(myMiss/nrow(gt))
  SNP.stats$Sample.med.miss[i] <- median(myMiss/nrow(gt))
  SNP.stats$Sample.max.miss[i] <- max(myMiss/nrow(gt))
  SNP.stats$Sample.min.miss[i] <- min(myMiss/nrow(gt))
  
}

write.table(SNP.stats, "4_Manuscript/plots/SNP_stats_all.txt", row.names = F, quote = F)

