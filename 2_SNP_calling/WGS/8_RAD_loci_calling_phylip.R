
# Load R packages
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2.1/")
library(ape)


# Get system arguments
args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[1]
# SNP.library.name <- "VCF_chrom_r10000"
SNP.library.location <- args[2]
# SNP.library.location <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/"
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

# Get sequences names
sequences <- list.files(paste0(dir.path,"/RAD_loci/"),pattern = "*.phylip$")

sequences <- gsub(x = sequences, "\\.phylip", "")
# Read in DNA sequences
RAD.loci <- lapply(paste0(dir.path,"RAD_loci/", sequences, ".phylip"), read.dna, format = "interleaved")

# Rename samples
Erandi_sample_names <- cbind(c("ZANAa05", "CUAJa03","CUAJb01","MIXTa01"), c("E008","E018","E019","E020"))

for(i in 1:length(RAD.loci)){
    rownames(RAD.loci[[i]])[rownames(RAD.loci[[i]])%in%Erandi_sample_names[,2]] <- Erandi_sample_names[,1][match(rownames(RAD.loci[[i]]),Erandi_sample_names[,2])[(rownames(RAD.loci[[i]])%in%Erandi_sample_names[,2])]]
    # Remove CUAJa03 and CUAJa02
    RAD.loci[[i]] <- RAD.loci[[i]][!rownames(RAD.loci[[i]])%in%c("CUAJa03", "CUAJa02_R"),]
    rownames(RAD.loci[[i]]) <- gsub("_R", "", rownames(RAD.loci[[i]]))
}

# Loop through all loci and format in to gphocs format
g_phocs <- lapply(1:length(RAD.loci), function(x) {
    c(paste(sequences[x], dim(RAD.loci[[x]])[1], dim(RAD.loci[[x]])[2], collapse = ""),
    paste(rownames(RAD.loci[[x]]), apply(as.character(RAD.loci[[x]]), MARGIN = 1, function(y) {
        gsub(x = paste0(toupper(y), collapse = ""), pattern= "\\?", replacement = "N")})),
    "")
})

dir.create(paste0(dir.path, "G-Phocs/"))
# Write out lines
print("Saving gphocs file to")
print(paste0(dir.path, "G-Phocs/", SNP.library.name,".gphocs"))
writeLines(c(length(RAD.loci), "", do.call("c", g_phocs)), con = paste0(dir.path, "G-Phocs/", SNP.library.name,".gphocs"))

## Create g-phocs file just for pacific individuals
all.samples <- c("CUAJa01", "CUAJa03", "ZANAa05", "CUAJb01", "MIXTa01", "MIXTa02", "RLPEb01", "ZANAa07")
pac.samples <- c("CUAJa03", "ZANAa05", "RLPEb01", "ZANAa07")
atl.samples <- c("CUAJa01", "CUAJb01", "MIXTa01", "MIXTa02")

RAD.loci.pac <- lapply(RAD.loci, function(x) x[row.names(x)%in%pac.samples,])
RAD.loci.atl <- lapply(RAD.loci, function(x) x[row.names(x)%in%pac.samples,])

# Loop through all loci and format in to gphocs format
g_phocs <- lapply(1:length(RAD.loci.pac), function(x) {
    c(paste(sequences[x], dim(RAD.loci.pac[[x]])[1], dim(RAD.loci.pac[[x]])[2], collapse = ""),
    paste(rownames(RAD.loci.pac[[x]]), apply(as.character(RAD.loci.pac[[x]]), MARGIN = 1, function(y) {
        gsub(x = paste0(toupper(y), collapse = ""), pattern= "\\?", replacement = "N")})),
    "")
})

writeLines(c(length(RAD.loci), "", do.call("c", g_phocs)), con = paste0(dir.path, "G-Phocs/", SNP.library.name,"_pac.gphocs"))

# Write out names of configs information
assignFile <- cbind(pac.samples, substr(pac.samples, 1, nchar(pac.samples)-3))

# Write out pop file
write.table(assignFile, file = paste0(dir.path,"/G-Phocs/GPhocs_samples_pac.txt"), col.names=F, row.names = F, quote=F)


