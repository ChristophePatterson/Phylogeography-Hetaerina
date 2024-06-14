.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
library(ggplot2)
library(ape)
library(poppr)
library(adegenet)
library(phangorn)
library(phrapl)

# Get system arguments
args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[1]
# SNP.library.name <- "Hetaerina_all_ddRAD_titia_dg"
SNP.library.location <- args[2]
# SNP.library.location <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/"
# Directory
dir.path <- paste0(SNP.library.location,SNP.library.name,"/")

#trees<-read.tree(paste(path.package("phrapl"),"/extdata/trees.tre",sep=""))
trees <- read.tree(paste0(dir.path,"PHRAPL/", SNP.library.name, "_phrapl.trees"))

## # plot trees
## png(paste0("/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/PHRAPL/", SNP.library.name, "_phrapl_trees.png"), units = "px", width = 500, height= 500)
## par(mfrow = c(3, 3))
## for(i in 1:9){
##     plot(trees[[i]])
## }
## dev.off()
## par(mfrow = c(1, 1))

#assignFile<-read.table(paste(path.package("phrapl"),"/extdata/cladeAssignments.txt",sep=""),
#     header=TRUE,stringsAsFactors=FALSE)
#assignFile

assignFile <- read.table("/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/LEA/LEA_pop_assign/popfile_all_noCUAJ_tabdelim.txt")
colnames(assignFile) <- c("indiv", "popLabel")
# Get all sample names from trees adn remove duplicates
samples <- do.call("c", lapply(trees, function(x) x$tip.label))
samples <- samples[!duplicated(samples)]

# Remove samples from assignFile that are not in trees
# Any tree tips that are not in the assignfile
any(!(samples%in%assignFile$indiv))
assignFile <- assignFile[assignFile$indiv %in% samples,]

# Convert assign pop to alpha
assignFile$popLabel_alpha <- LETTERS[match(as.factor(assignFile$popLabel), levels(as.factor(assignFile$popLabel)))]
assignFile <- assignFile[order(assignFile$popLabel_alpha),]
# Create pop vector (recommend subsampling in a way that keeps the sum of popAssignments at 16 or less)
subsample_N <- 2
popN <- length(unique(assignFile$popLabel))
popVector <- rep(subsample_N, length(unique(assignFile$popLabel)))
sum(popVector)

lapply(trees, function(x) {
    tree <- as.phylo(x)
    assigntips <- assignFile$popLabel_alpha[na.omit(match(assignFile$indiv, tree$tip.label))]
    return(table(assigntips))
}
)

assignFile$indiv[is.na(match(assignFile$indiv, trees[[38]]$tip.label))]
trees[[605]]$tip.label

# Check each tree has a least N samples per populations
missing_data_tree <- do.call("c",lapply(trees, function(x) {
    tree <- as.phylo(x)
    tip.assign.per.tree <- table(assignFile$popLabel_alpha[assignFile$indiv %in% tree$tip.label])
    return(any(as.numeric(tip.assign.per.tree) < subsample_N) | length(tip.assign.per.tree)<popN)
}
))
# Subset trees
trees <- trees[which(!missing_data_tree)]

assignmentsGlobal<-assignFile[,c(1,3)]
colnames(assignmentsGlobal) <- c("indiv","popLabel")
observedTrees <- trees  
popAssignments<-list(popVector) 
subsamplesPerGene<- 2
outgroup=FALSE
outgroupPrune=FALSE  

observedTrees <- PrepSubsampling(assignmentsGlobal=assignmentsGlobal,observedTrees=observedTrees,
    popAssignments=popAssignments,subsamplesPerGene=subsamplesPerGene,outgroup=outgroup,
    outgroupPrune=outgroupPrune)
