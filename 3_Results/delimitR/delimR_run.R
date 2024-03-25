.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")

# Get sample sizes
args <- commandArgs(trailingOnly = TRUE)

# Set output and working directory

# dir_output <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_Pool_X_HetTit1.0_HetAmer1.0/Hetaerina_titia_ddRAD_titia_dg/delimitR"
sp_0 <- args[1]
sp_1 <- args[2]
sp_2 <- args[3]
# dir_output <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_Pool_X_HetTit1.0_HetAmer1.0/Hetaerina_titia_ddRAD_titia_dg/delimitR"
dir_output <- paste0(args[4], "/", args[6])
vcf <- args[5]
threads <- args[7]

print(dir_output)
setwd(dir_output)

library(delimitR)
library(vcfR)

# SFS observations - name of file without extension
observedSFS <- args[8]

#location of our traits file (2 column file which maps alleles to species, must be in wd)
traitsfile <- "all_traits.txt"

#guide tree
observedtree <- '((0,1),2);'

#migration matrix (must be symmetrical)
migmatrix <- matrix(c(FALSE, TRUE, TRUE,
                    TRUE, FALSE, TRUE,
                    TRUE, TRUE, FALSE),
                    nrow = 3, ncol = 3, byrow = TRUE)
#test divergence with gene flow?
# divwgeneflow <- FALSE
divwgeneflow <- TRUE

#test secondary contact?
seccontact <- TRUE

#what is the maximum number of migration events to consider on your guide tree?
# maxedges <- 2
maxedges <- 3

#how many species are in your guide tree?
obsspecies<- 3

#the number of "alleles" retained after downsampling SFS
# obssamplesize <- as.numeric(c(86,174, 82))
obssamplesize <- as.numeric(c(sp_0,sp_1, sp_2))

#The user must specify the number of linkage blocks to simulate
#For unlinked SNPs, this is equal to the number of SNPs used to build your observed SFS
# vcf <- "Hetaerina_titia_ddRAD_titia_dg.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.bi.0_20.vcf.gz"
obssnps <- readLines("snp_num.txt")

obsprefix <- args[8]

#The user must specify priors on population sizes
#The first vector is for population 0, the second for population 1, and the third for population 2
#Note that these are in terms of the number of haploid individuals (as specified in the fsc2 documentation)
obspopsizeprior <- list(c(100000,2000000),c(100000,2000000),c(100000,2000000))

#priors for divergence times given in terms of the number of generations and must be supplied as a list of vectors
#Divergence time priors should be provided in order of coalescent interval
obsdivtimeprior <- list(c(500000,8000000),c(500000,12000000))

# Create rules about how divergence dates should be ordered
myrules <- c('Tdiv2$>Tdiv1$')

#prior on migration rates, program only allows one prior for all migration rates in the default model sets
obsmigrateprior <- list(c(0.000005,0.0005))

#set up your prior models for fastsimcoal2
setup_fsc2(tree=observedtree,
           nspec=obsspecies,
           samplesizes=obssamplesize,
           nsnps=obssnps,
           prefix=obsprefix,
           migmatrix=migmatrix,
           popsizeprior=obspopsizeprior,
           divtimeprior=obsdivtimeprior,
           myrules=myrules,
           migrateprior=obsmigrateprior,
           secondarycontact= seccontact,
           divwgeneflow= divwgeneflow,
           maxmigrations = maxedges)

## Function that writes out slurm scripts for each of the fastsimcol simulations

# Updated fastsimcoalsims function that includes multithreading by submitting an individual sbatch for each run
# This code instead on submitting each simualtion to fastsimcoal2 one line at a time, creates a slurm script or each simulation and excutes them all.
fastsimcoalsim_sbatch <- function(prefix,pathtofsc,nreps){
  listoftpl <- list()
  listofest <- list()
  tpllist <- system(paste("ls ", prefix,"*.tpl",sep=""),intern=T)
  estlist <- system(paste("ls ", prefix,"*.est",sep=""),intern=T)
  listoftpl <- c(listoftpl, tpllist)
  listofest <- c(listofest, estlist)
  count <- 1
  for(j in 1:length(listoftpl)){
    # Creates a sbatch file for each of the individu
    sbatch_file <- c("#!/bin/bash",
                      "#SBATCH -c 1" ,
                      "#SBATCH --mem=10G",
                      "#SBATCH --gres=tmp:10G",
                      "#SBATCH -t 72:00:00",
                      paste0("#SBATCH --output=", getwd(), "slurm-%x.%j.out"),
                      "",
                      print(paste(pathtofsc, " -t ", prefix, "_", count, ".tpl",  " -e ", prefix, "_", count, ".est -n 1 --msfs -q --multiSFS -x -E" ,nreps, sep = "")))
  writeLines(sbatch_file, paste(prefix, "_", count , "_run.sh", sep = ""))
  system(paste(paste("sbatch ", prefix, "_", count , "_run.sh", sep = "")), wait = F)
  Sys.sleep(5)
count <- count + 1
  }
}

# Run fastsimcoal as slurm jobs
fastsimcoalsim_sbatch(prefix=obsprefix,
                pathtofsc='/nobackup/tmjj24/apps/fsc26_linux64/fsc26',
                nreps=10000)

