.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")


args <- commandArgs(trailingOnly = TRUE)

# demultiplex file location
# input_folder <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/demultiplex_all"
input_folder <- args[1]

# Output folder
# output_folder <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/stacks"
output_folder <- args[2]
# Number of samples to use per population
# select_N <- 2
select_N <- as.numeric(args[3])

# Get list of demuliplexed files
sequences <- list.files(input_folder)
# Only retain forwards reads
sequences <- sequences[grep("\\.1\\.", sequences)]

seq_info <- file.info(paste0(input_folder, "/", list.files(input_folder)), extra_cols = F)
seq_info$sequences <- row.names(seq_info)

seq_info <-  seq_info[grep("\\.1\\.", seq_info$sequences),]

seq_info$sample <- substr(seq_info$sequences, nchar(input_folder)+2, nchar(seq_info$sequences)-8)

head(seq_info)

# Remove contaiminated samples
dodgy <- c("HtiTi12", "NA0101","CA0101","CUAJa02.Dur","CUAJa02.Shef","HXRCaAM03", 
            "CUAJb19", "CUAJb18", "CUAJb01", "CUAJb06", "CUAJb21", "CUAJb02", "CUAJb07", "CUAJb13", 
            "CUAJb15", "CUAJa03", "CUAJb14", "CUAJb12", "CUAJb16", "CUAJb10", "CUAJb08", "CUAJb03", 
            "CUAJb20", "CUAJb09", "CUAJb17", "CUAJb11")

seq_info <- seq_info[!seq_info$sample%in%dodgy,]

# Read in sample assignment
hetaerina_assign <- read.table("/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/LEA/LEA_pop_assign/popfile_all_noCUAJ.txt")
colnames(hetaerina_assign) <- c("sample", "assign")

print(paste("Any samples duplicated in assign LEA table? - ", any(duplicated(hetaerina_assign$sample))))

# Merge with raw sequence reads
seq_info <- merge(hetaerina_assign, seq_info, by.x = "sample", by.y="sample")

# Get samples with highest number of reads from each region

top.cov.samples <- lapply(unique(seq_info$assign), function(x) {
    pop_assign_temp <- seq_info[seq_info$assign==x,]
    if(dim(pop_assign_temp)[1]>select_N){
        print(paste0(x, "greater than N"))
        print(order(pop_assign_temp$size, decreasing = T))
        return((pop_assign_temp[order(pop_assign_temp$size, decreasing = T),])[1:select_N,])
    }
    if(dim(pop_assign_temp)[1]<=select_N){
        print(paste0(x, "less than and equal to N"))
        return((pop_assign_temp[order(pop_assign_temp$size, decreasing = T),])[1:dim(pop_assign_temp)[1],])
    }
})
# MErge into one
#top.cov.samples
pop_assign_N <- do.call("rbind", top.cov.samples)
#pop_assign_N

print(paste("Any samples duplicated in selected top coverage samples? - ", any(duplicated(pop_assign_N$sample))))


# Save popmap
write.table(pop_assign_N[,c("sample")], paste0(output_folder,"/", "samples_", select_N,".txt"), 
            row.names = F, col.names = F, quote = F, sep = "\t")

## Write popmap for each species
titia_assign <- c("titia-Pac", "titia-NAtl", "titia-SAtl")
ameri_assign <- c("calverti-Mex", "americana-Mex", "americana-USA")

# Write out titia
write.table(pop_assign_N[pop_assign_N$assign%in%titia_assign,c("sample")], paste0(output_folder,"/", "samples_", select_N,"_titia.txt"), 
            row.names = F, col.names = F, quote = F, sep = "\t")

write.table(pop_assign_N[pop_assign_N$assign%in%titia_assign,c("sample", "assign")], paste0(output_folder,"/", "samples_", select_N,"_titia_assign.txt"), 
            row.names = F, col.names = F, quote = F, sep = "\t")

write.table(pop_assign_N[pop_assign_N$assign%in%ameri_assign,c("sample")], paste0(output_folder,"/", "samples_", select_N,"_americana.txt"), 
            row.names = F, col.names = F, quote = F, sep = "\t")

write.table(pop_assign_N[pop_assign_N$assign%in%ameri_assign,c("sample", "assign")], paste0(output_folder,"/", "samples_", select_N,"_americana_assign.txt"), 
            row.names = F, col.names = F, quote = F, sep = "\t")