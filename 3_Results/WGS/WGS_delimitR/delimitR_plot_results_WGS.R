.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
library(delimitR)
library(vcfR)
library(ggplot2)
library(patchwork)

# Change python path
Sys.setenv(PATH = paste("/nobackup/tmjj24/apps/miniconda3/envs/delimR/bin", Sys.getenv("PATH"), sep=":"))

# Get sample sizes
args <- commandArgs(trailingOnly = TRUE)

upper_dir <- args[1]
lower_dir <- args[2]
ncores <- args[3]
SNP.library <- args[4]
# SNP.library <-"WGS_titia_PAC"

str(ncores)
print(paste("Using",ncores, "cores"))
ncores <- as.numeric(ncores)
# dir_output <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/WGS_titia_PAC_r10000/delimitR/2_2_4_nreps1000_MIG3_trees_3_sbatch"
dir_output <- paste0(upper_dir, lower_dir)

#Change working directory
print(dir_output)
setwd(dir_output)

# SFS observations - name of file without extension
# obsprefix <- "Hetaerina_americana_ddRAD_titia_dg"
obsprefix <- args[4]
obsspecies<- 3
nclasses <- 3
traitsfile <- "all_traits.txt"

# Read in observed data
myobserved <- readRDS("myobserved.rds")
## Read in reduced prior
ReducedPrior <- readRDS("ReducedPrior.rds")
# Read in Random forest
myRF <- readRDS("RF_build_abcrf.rds")

# Write out of the box error rate
write.table(myRF$model.rf$confusion.matrix, file = paste0(SNP.library,"_out_of_bag_error.txt"))

# Get number of scereios   
num.scenarios <- length(grep(".est",list.files()))
#how many simulations where conducted
num.sims <- dim(ReducedPrior)[1]/num.scenarios

# Get PCA of fastsimcoal simulations
pca <- prcomp(ReducedPrior, scale. = T)
pca.data <- data.frame(cbind(pca.1 = pca$x[,1], pca.2 = pca$x[,2], pca.3 = pca$x[,3], pca.4 = pca$x[,4], model = rep(1:num.scenarios, each = num.sims)))
#  Convert observed data into PCA
myobserved.pca <- data.frame(predict(pca, myobserved))

# Plot PCA data
p <- ggplot(pca.data) +
  geom_point(aes(pca.1, pca.2, col = as.factor(model)), shape = 3) +
  geom_point(data = myobserved.pca, aes(PC1, PC2, col = "observed"), size = 5, col = "black") +
  labs(col = "Model", x = "pca1", y = "pca2")

q <- ggplot(pca.data) +
  geom_point(aes(pca.3, pca.4, col = as.factor(model)), shape = 3) +
  geom_point(data = myobserved.pca, aes(PC3, PC4, col = "observed"), size = 5, col = "black") +
  labs(col = "Model", x = "pca3", y = "pca4")

ggsave(paste0("DelimitR_PCA12_plot_",SNP.library,".png"), p  , width = 14, height = 10)
### ggsave(paste0("DelimitR_PCA34_plot_",SNP.library,".png"), q  , width = 14, height = 10)

# Calculate which par file produced the most similar SFS to that of the obs data
prediction <- read.csv(paste0(obsprefix,"_prediction.csv"))

# Which prediction is the highest 
best.model <- prediction$selected.model

# Get pca for best model
pca.data.best.model <- pca.data #[pca.data$model==best.model,]

pca.data.best.model$dist.observered <- apply(pca.data.best.model, MARGIN=1, function(x) {
  sqrt(sum((myobserved.pca[1,1:2] - x[1:2])^2))
})

dist.p <- ggplot(pca.data.best.model) +
              geom_point(aes(pca.1, pca.2, col = dist.observered ))  +
              geom_point(data = myobserved.pca, aes(PC1, PC2, col = "observed"), size = 5, col = "black")

# For the best model which is closest in PCA space
best.par <- which.min(pca.data.best.model$dist.observered[pca.data.best.model$model==best.model])
writeLines(text = paste0(obsprefix, "_", best.model, "/", obsprefix, "_", best.model, "_", best.par,".par"), con = "best_par_file.txt")


