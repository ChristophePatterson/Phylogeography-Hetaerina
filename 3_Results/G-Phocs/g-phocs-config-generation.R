# Code for the creation of config files for G-Phocs
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
library(stringr)
# Get system arguments
args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[1]
# SNP.library.name <- "Hetaerina_titia_ddRAD_titia_dg"
SNP.library.location <- args[2]
# SNP.library.location  <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_v2/"
dir.path <- paste0(SNP.library.location,SNP.library.name,"/")
dir.path.GPhocs <- paste0(dir.path, "G-Phocs/")

# Number of samples to use for each population
sample.num <- as.numeric(args[3])
# sample.num <- 2

## List of each species (must be how that is defined within the samples data set species column)
species <- str_split_1(SNP.library.name, pattern = "_")[2]

# Different models
# Iso = No migration
# migA = All migration bands (Except from Pac to N/S)
# migR = Mig bands between most recent diverged samples
model.types <- c("Iso", "migA","migR")

alphas <- c(1)
betas <- c(20,200)

# Create data frame of all different combination of models
model.df <- expand.grid("model"=model.types,"species"=species, "alpha" = alphas, "beta"= betas)

# Samples used in gphocs
pop_assign <- read.table(paste0(dir.path.GPhocs,"GPhocs_samples_N", sample.num,".txt"))
colnames(pop_assign) <- c("sample", "pop")


# Creat name for each model
model.names <- paste0("G-Phocs-a", model.df$alpha, "-b",model.df$beta, "-", model.df$species, "-", model.df$model, "-run0")
# Write file that has the names of each unique model name
write.table(x = model.names, file = paste0(dir.path.GPhocs,"/","G-Phocs-model-names.txt"),
              row.names = F, col.names = F, quote = F, sep = "\t")


set.seed(28)
i <- 1
dir.create(paste0(dir.path.GPhocs,"model_runs/"))

for(i in 1:dim(model.df)[1]){
    #print(model.df[i,])
    # read in correct default config file
    model.name <- model.names[i]
    print(model.name)
    # Read in correct default config file
    config.file <- readLines(paste0("/home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/3_Results/G-Phocs/G-phocs-configs/config-defaults/",
                        "G-phocs-a1-b200-",model.df$species[i],"_max_",model.df$model[i],".config"))
    # Get data set of samples to be use

    # Start replacing lines in config file
    # Input sequence
    config.file[3]  <- paste0("\t\tseq-file\t\t",dir.path.GPhocs,SNP.library.name,".gphocs")
    # Output trace file
    config.file[4]  <- paste0("\t\ttrace-file\t\t",dir.path.GPhocs,"model_runs/",model.name,".trace")
    # Replace alpha and beta
    config.file[13]  <- paste0("\t\ttau-theta-alpha\t\t",model.df$alpha[i])
    config.file[14]  <- paste0("\t\ttau-theta-beta\t\t",model.df$beta[i])
    # If species = americana over write samples using americana/calv populations
    if(model.df$species[i]=="americana"){
    config.file[34] <- paste0(c("\t\t\t\tsamples\t\t", paste0(pop_assign$sample[pop_assign$pop=="calverti-Mex"], sep = " d ", collapse = "")), collapse = "")
    config.file[39] <- paste0(c("\t\t\t\tsamples\t\t", paste0(pop_assign$sample[pop_assign$pop=="americana-USA"], sep = " d ", collapse = "")), collapse = "")
    config.file[44] <- paste0(c("\t\t\t\tsamples\t\t", paste0(pop_assign$sample[pop_assign$pop=="americana-Mex"], sep = " d ", collapse = "")), collapse = "")
    }
    # If species = titia over write samples using titia populations
    if(model.df$species[i]=="titia"){
    config.file[34] <- paste0(c("\t\t\t\tsamples\t\t", paste0(pop_assign$sample[pop_assign$pop=="titia-Pac"], sep = " d ", collapse = "")), collapse = "")
    config.file[39] <- paste0(c("\t\t\t\tsamples\t\t", paste0(pop_assign$sample[pop_assign$pop=="titia-NAtl"], sep = " d ", collapse = "")), collapse = "")
    config.file[44] <- paste0(c("\t\t\t\tsamples\t\t", paste0(pop_assign$sample[pop_assign$pop=="titia-SAtl"], sep = " d ", collapse = "")), collapse = "")
    }
    
    # Write code output
    writeLines(config.file, con = paste0(dir.path.GPhocs,"model_runs/",model.name,".config"),sep = "\n")

    # Create slurm script
    
}

# Function that creates and excutes slurm scripts for each model
slurm_submit_gphocs <- function(prefix,models,pathtogphocs,ncores){
  for(model in models){
    # Creates a sbatch file for each of the individu
    sbatch_file <- c("#!/bin/bash",
                      paste0("#SBATCH -c ", ncores),
                      "#SBATCH --mem=10G",
                      "#SBATCH --gres=tmp:10G",
                      "#SBATCH -t 72:00:00",
                      paste0("#SBATCH --output=slurm-%x.%j.out"),
                      "",
                      print(paste(pathtogphocs, " ",prefix,model,".config -n ",ncores, " > ",prefix,model,".log", sep = "")))
  writeLines(sbatch_file, paste(prefix, model,"_run.sh", sep = ""))
  system(paste("sbatch ", paste(prefix, model,"_run.sh", sep = "")), wait = F)
  }
}

slurm_submit_gphocs(prefix = paste0(dir.path.GPhocs, "model_runs/"), models = model.names, pathtogphocs = "/nobackup/tmjj24/apps/G-PhoCS/bin/G-PhoCS", ncores = 24)


