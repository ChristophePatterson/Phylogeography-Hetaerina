.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
#Libraries needed
library(ape)

# Get system arguments
args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[1]
# SNP.library.name <- "Hetaerina_titia_ddRAD_titia_dg"
SNP.library.location <- args[2]
N_loci <- as.numeric(args[3])
# SNP.library.location <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/"
# Directory
dir.path <- paste0(SNP.library.location,SNP.library.name,"/")

# Get sequences names
sequences <- list.files(paste0(dir.path,"/RAD_loci/"),pattern = "*.phylip$")
# Remove sequences that are not in the autosomes
if(grepl(SNP.library.name, pattern="titia_dg")){
  X_chrom <- "HetTit1.0.p.scaff-12-96647824"
  # Remove X chromosome
  sequences <- sequences[!sequences%in%X_chrom]
  #Remove non autosomes 
  sequences <- sequences[sapply(strsplit(sequences,"-") , "[", 2 )%in%1:11]
}
# For H. americana lots of small X chromosomes
if(grepl(SNP.library.name, pattern="americana_dg")){
  X_chrom <- c("JAKTNV010000016.1","JAKTNV010000031.1","JAKTNV010000013.1","JAKTNV010000060.1","JAKTNV010000023.1","JAKTNV010000051.1","JAKTNV010000012.1")
}

# Limit to N number of sequences
# If more than N loci sample N
if(length(sequences)>N_loci){
    loci.select <- sort(sample(1:length(sequences), N_loci))
    sequences <- sequences[loci.select]
}

sequences <- gsub(x = sequences, "\\.phylip", "")
# Read in DNA sequences
RAD.loci <- lapply(paste0(dir.path,"RAD_loci/", sequences, ".phylip"), read.dna, format = "interleaved")

# Loop through all loci and format in to BPP format
BPP <- lapply(1:length(RAD.loci), function(x) {
    c(paste(dim(RAD.loci[[x]])[1], dim(RAD.loci[[x]])[2], collapse = ""),
    "",
    paste0("^", rownames(RAD.loci[[x]]), " ", apply(as.character(RAD.loci[[x]]), MARGIN = 1, function(y) {
        paste0(toupper(y), collapse = "")})),
    "")
})
# Write out lines
writeLines(do.call("c", BPP), con = paste0(dir.path, "BPP/", SNP.library.name, "_N", N_loci, ".txt"))

# Write out names of configs information
assignFile <- read.table("/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/LEA/LEA_pop_assign/popfile_all_noCUAJ_tabdelim.txt")
colnames(assignFile) <- c("indiv", "popLabel")
# Get all sample names from trees adn remove duplicates
samples <- do.call("c", lapply(RAD.loci, function(x) rownames(x)))
samples <- samples[!duplicated(samples)]
# Remove samples from assignFile that are not in trees
assignFile <- assignFile[assignFile$indiv %in% samples,]

# Write out pop file
write.table(assignFile, file = paste0(dir.path,"/BPP/BPP_samples_N",N_loci,".txt"), col.names=F, row.names = F, quote=F)
