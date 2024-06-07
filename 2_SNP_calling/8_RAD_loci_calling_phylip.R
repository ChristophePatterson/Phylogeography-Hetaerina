
# Load R packages
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
library(ggplot2)
library(ape)
library(vcfR)
library(poppr)
library(adegenet)


# Get system arguments
args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[1]
SNP.library.location <- args[2]
# Directory
dir.path <- paste0(SNP.library.location,SNP.library.name,"/")

# Get sequences names
sequences <- list.files(paste0(dir.path,/"/fasta_files/"))
# Read in DNA sequences
alignments <- lapply(paste0(dir.path,"fasta_files/", sequences), read.dna, format = "fasta")

# Get loci names
RAD.loci.names <- rownames(alignments[[1]])

DNA.dims <- dim(alignments[[1]])[1]
# Loop through each loci and read out
RAD.loci <- lapply(1:DNA.dims[1], function(x) {
    # Get each loci from each sample
    align.tmp <- do.call(rbind.DNAbin, lapply(alignments, function(y) y[x,]))
    #Get sample names
    rownames(align.tmp) <- substr(sequences, 1, nchar(sequences)-6)
    # Which bases are missing from all samples
    missing.data <- which(colSums(as.character(align.tmp) == "?")!=DNA.dims[1])
    # Which samples are missing data from every base
    missing.seq <- which(rowSums(as.character(align.tmp) == "?")!=DNA.dims[2])
    # Remove missing samples and remove all sequences before and after regions where no samples are called
    align.tmp <- align.tmp[missing.seq,missing.data[1]:missing.data[length(missing.data)]]
    # Write out DNA
    write.dna(align.tmp, file = paste("/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/Hetaerina_titia_ddRAD_titia_dg/RAD_loci/",RAD.loci.names[x],".phy"),
    format="interleaved")
})





