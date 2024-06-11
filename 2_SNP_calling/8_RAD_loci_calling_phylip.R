
# Load R packages
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
library(ggplot2)
library(ape)
library(poppr)
library(adegenet)

# Get system arguments
args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[1]
# SNP.library.name <- "Hetaerina_all_ddRAD_titia_dg"
SNP.library.location <- args[2]
# SNP.library.location <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/"
# Directory
dir.path <- paste0(SNP.library.location,SNP.library.name,"/")

# Get sequences names
sequences <- list.files(paste0(dir.path,"/fasta_files/"))
# Read in DNA sequences
alignments <- lapply(paste0(dir.path,"fasta_files/", sequences), read.dna, format = "fasta")


# Get loci names
RAD.loci.names <- rownames(alignments[[1]])

Alignment.dims <- dim(alignments[[1]])
# Loop through each loci and read out
# min number of samples base has to be in
min.bases <- 0.2
Alignment.dims[1]
RAD.loci <- lapply(1:Alignment.dims[1], function(x) {
    # Get each loci from each sample
    align.tmp <- do.call(rbind.DNAbin, lapply(alignments, function(y) y[x,]))
    DNA.dims <- dim(align.tmp)
    #Get sample names
    rownames(align.tmp) <- substr(sequences, 1, nchar(sequences)-6)
    # Which bases are missing from all samples
    missing.bases <- which(colSums(as.character(align.tmp) != "?")>=(DNA.dims[1]*min.bases))
    align.tmp <- align.tmp[,missing.bases[1]:missing.bases[length(missing.bases)]]
    # Recalculate DNA dim
    DNA.dims <- dim(align.tmp)
    # Which samples are missing data from every base
    missing.samples <- which(rowSums(as.character(align.tmp) == "?")!=DNA.dims[2])
    # Remove missing samples and remove all sequences before and after regions where no samples are called
    align.tmp <- align.tmp[missing.samples,]
    # Write out DNA
    write.dna(align.tmp, file = paste0(dir.path, "/RAD_loci/",RAD.loci.names[x],".phylip"),
    format="interleaved")
    return(align.tmp)
})

# Convert to G-Phocs

