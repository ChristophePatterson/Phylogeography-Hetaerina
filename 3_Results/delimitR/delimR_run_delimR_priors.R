.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
library(delimitR)
library(vcfR)

# Change python path
Sys.setenv(PATH = paste("/nobackup/tmjj24/apps/miniconda3/envs/delimR/bin", Sys.getenv("PATH"), sep=":"))

# Get sample sizes
args <- commandArgs(trailingOnly = TRUE)

# Set output and working directory

# dir_output <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_Pool_X_HetTit1.0_HetAmer1.0/Hetaerina_titia_ddRAD_titia_dg/delimitR"
upper_dir <- args[1]
lower_dir <- args[2]
cores <- as.numeric(args[3])
# dir_output <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_Pool_X_HetTit1.0_HetAmer1.0/Hetaerina_titia_ddRAD_titia_dg/delimitR"
dir_output <- paste0(upper_dir, "/delimitR/", lower_dir)

#Change working directory
print(dir_output)
setwd(dir_output)

# SFS observations - name of file without extension
obsprefix <- args[4]
obsspecies<- 3
nclasses <- 3
traitsfile <- "all_traits.txt"

print("Full Prior")
FullPrior <- makeprior(prefix=obsprefix,
                       nspec=obsspecies,
                       nclasses=nclasses,
                       getwd(),
                       traitsfile = traitsfile,
                       threshold=100, 
                       thefolder = 'Prior',
                       ncores = cores)


# We want to remove rows that have zero variance. 
# For example, if no SNPs are ever observed in a certain bin of the SFS, across all models and all simulations, 
# this bin adds nothing to our analysis, and we must remove it before building the RF classifier. 
print("Run Prior_reduced")

ReducedPrior <- Prior_reduced(FullPrior)
saveRDS(ReducedPrior , "ReducedPrior.rds")

myRF <- RF_build_abcrf(ReducedPrior,FullPrior,1000)
saveRDS(myRF, "RF_build_abcrf.rds")
myRF

plot(myRF, training = ReducedPrior)

# Now look at the actual data
observedSFS <- paste0(args[4],"_MSFS")

myobserved <- prepobserved(
  observedSFS,
  FullPrior,
  ReducedPrior,
  nclasses,
  obsspecies,
  traitsfile=traitsfile,
  threshold = 100)

saveRDS(myobserved, "myobserved.rds")

prediction <- RF_predict_abcrf(myRF, myobserved, ReducedPrior, FullPrior, 1000)
prediction

#write results to file
write.csv(as.matrix(prediction), file="prediction.csv")
# Save out of bag error rate
write.table(myRF$model.rf$confusion.matrix, file = "out_of_bag_error.txt")


pca <- prcomp(ReducedPrior, scale. = T)
pca.data <- data.frame(cbind(pca.1 = pca$x[,1], pca.2 = pca$x[,2], model = rep(1:13, each = 10000)))

p <- ggplot(pca.data) +
  geom_point(aes(pca.1, pca.2, col = as.factor(model)), shape = 3) +
  geom_point(data = data.frame(predict(pca, myobserved)), aes(PC1, PC2, col = "observed"), size = 5, col = "black") +
  labs(col = "Model", x = "pca1", y = "pca2")

ggsave("PCA_plot.png", p , width = 10, height = 10)

