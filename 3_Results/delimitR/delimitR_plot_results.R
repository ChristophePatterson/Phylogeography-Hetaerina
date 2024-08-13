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
# SNP.library <-"Hetaerina_titia_ddRAD_titia_dg"

str(ncores)
print(paste("Using",ncores, "cores"))
ncores <- as.numeric(ncores)
# dir_output <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_v2/Hetaerina_titia_ddRAD_titia_dg/delimitR/82_174_76_nreps10000_MIG3_divWgene_sbatch"
dir_output <- paste0(upper_dir, lower_dir)

#Change working directory
print(dir_output)
setwd(dir_output)

# SFS observations - name of file without extension
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


pca <- prcomp(ReducedPrior, scale. = T)
pca.data <- data.frame(cbind(pca.1 = pca$x[,1], pca.2 = pca$x[,2], pca.3 = pca$x[,3], pca.4 = pca$x[,4], model = rep(1:num.scenarios, each = num.sims)))
myobserved.pca <- data.frame(predict(pca, myobserved))
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

