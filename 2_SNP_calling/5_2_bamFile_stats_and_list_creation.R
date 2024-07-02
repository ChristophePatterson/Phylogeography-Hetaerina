# Given x number of files within a directory that contain sequences names and nothing else
# create x number of new files that contain the full directory path for all bam files
# The directory cannot contain any other files
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2.1/")
library(tidyverse)
#install.packages("tidyverse")
#####################################################
# For Durham and Sheffield and CUAJ data November 2023###########
#####################################################

# Location of bamfile basic names
# directory <- "C:/Users/tmjj24/OneDrive - Durham University/Christophe/Work/Sequence analysis/Demultiplex_seq_processing_SDC/"
# directory <- "C:/Users/chris/OneDrive - Durham University/Christophe/Work/Sequence analysis/Demultiplex_seq_processing_SDC/"
directory <- "/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/bamstats/"

#OUtput directory
dir_out <- "/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/2_SNP_calling/library_combinations/bamfiles/"
dir.create(dir_out)

americana_dg <- read.table(paste0(directory,"Allsamples_HetAmer1.0_dg.bamstats"))
titia_dg <- read.table(paste0(directory,"Allsamples_HetTit1.0_dg.bamstats"))
elegans_dg <- read.table(paste0(directory,"Allsamples_ioIscEleg1.2_dg.bamstats"))
head(titia_dg)
colnames(americana_dg) <- c("sample", "num.reads", "coverage", "coverage.SD", "unknown", "prop.mapped")
colnames(titia_dg) <- c("sample", "num.reads", "coverage", "coverage.SD", "unknown", "prop.mapped")
colnames(elegans_dg) <- c("sample", "num.reads", "coverage", "coverage.SD", "unknown", "prop.mapped")

#Calculating 2sd of map reads
sd.2.americana <- mean(americana_dg$num.reads)-sd(americana_dg$num.reads)*2
americana.min.cov <- 5
americana_dg$cov10_NoR2sd <- americana_dg$coverage>=americana.min.cov&americana_dg$num.reads>=sd.2.americana
sd.2.titia <- mean(titia_dg$num.reads)-sd(titia_dg$num.reads)*2
titia.min.cov <- 5
titia_dg$cov10_NoR2sd <- titia_dg$coverage>=titia.min.cov&titia_dg$num.reads>=sd.2.titia
sd.2.elegans <- mean(elegans_dg$num.reads)-sd(elegans_dg$num.reads)
elegans.min.cov <- 15
elegans_dg$cov10_NoR2sd <- elegans_dg$coverage>=elegans.min.cov&elegans_dg$num.reads>=sd.2.elegans

duplicated(titia_dg$sample)
dim(titia_dg)

americana_dg$library <- "americana_dg"
titia_dg$library <- "titia_dg"
elegans_dg$library <- "elegans_dg"

df <- rbind(americana_dg, titia_dg, elegans_dg)

dodgy <- c("HtiTi12", "NA0101","CA0101","CUAJa02.Dur","CUAJa02.Shef","HXRCaAM03", 
            "CUAJb19", "CUAJb18", "CUAJb01", "CUAJb06", "CUAJb21", "CUAJb02", "CUAJb07", "CUAJb13", 
            "CUAJb15", "CUAJa03", "CUAJb14", "CUAJb12", "CUAJb16", "CUAJb10", "CUAJb08", "CUAJb03", 
            "CUAJb20", "CUAJb09", "CUAJb17", "CUAJb11")

df <- df[!df$sample%in%dodgy,]

table(df$library, df$cov10_NoR2sd)

## Summary stats of bamfiles

bamstats <- group_by(df, by = library) %>%
  summarise(Total.mapped.reads = sum(num.reads), 
            mean.mapped.reads = mean(num.reads),
            max.mapped.reads = max(num.reads),
            min.mapped.reads = min(num.reads),
            sd.2 = mean(num.reads)-sd(num.reads)*2,
            mean.perc.mapped = mean(prop.mapped))

sum(df$num.reads[df$library=="titia_dg"])

colnames(bamstats)[1] <- "library"

bamstats
write.table(x = bamstats, file = paste0(directory,"/bam_summary_stats_SDC.txt"), sep = "\t")

p <- ggplot(df) +
  geom_point(aes(coverage, num.reads, col = cov10_NoR2sd), show.legend = F, size = 3) +
  #geom_vline(xintercept = 10) +
  #geom_hline(yintercept = sd.2.inter) +
  #geom_label(aes(60, sd.2.inter, label = round(sd.2.inter))) +
  xlab("Mean coverage of mapped reads") +
  ylab("Total number of mapped reads") +
  geom_hline(data = bamstats, aes(yintercept = sd.2), col = "#F8766D") +
  geom_vline(xintercept = 5, col = "#F8766D") +
  geom_label(data = bamstats, aes(60, sd.2, label = round(sd.2)), col = "#F8766D") +
  geom_label(data = bamstats, aes(5, max(df$num.reads), label = "5x"), col = "#F8766D") +
  facet_wrap(~library) +
  ylim(0, max(df$num.reads)) +
  theme_bw()
p

df$CUAJ <- grepl(df$sample, pattern = "CUAJ")

ggsave(p, filename = paste0(directory, "/Number of Mapped reads and mean percentage coverage HetTit1.0 and HetAmer1.0 and elegans July 2024.jpeg"), height = 6, width = 12)

# Location of bamfiles
americana_directory <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/bwa_HetAmer1.0_dg"
titia_directory <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/bwa_HetTit1.0_dg"
elegans_directory <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/bwa_ioIscEleg1.2_dg"

# Hetaerina sample data
hetaerina_data <- read.csv(paste0("/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/All samples held in Durham_v17.csv"), check.names = F)
hetaerina_data <- hetaerina_data[,c("Unique.ID", "species")]
#Adding in new names for samples that were sequenced in pool X
hetaerina_data$Unique.ID[hetaerina_data$Unique.ID%in%c(dodgy, "CUAJa02")] <- paste0(hetaerina_data$Unique.ID[hetaerina_data$Unique.ID%in%c(dodgy, "CUAJa02")],"_X")


df <- merge(df,hetaerina_data, by.x = "sample", by.y = "Unique.ID")
dim(df)/2
df$sample[grep(df$sample, pattern = "_X")]

p <- ggplot(df) +
  geom_point(aes(coverage, num.reads, col = species, shape = cov10_NoR2sd), size = 3) +
  #geom_label(aes(coverage, num.reads, col = species, label = sample), size = 3)
  #geom_vline(xintercept = 10) +
  #geom_hline(yintercept = sd.2.inter) +
  #geom_label(aes(60, sd.2.inter, label = round(sd.2.inter))) +
  xlab("Mean coverage of mapped reads") +
  ylab("Total number of mapped reads") +
  geom_hline(data = bamstats, aes(yintercept = sd.2), col = "#F8766D") +
  geom_vline(xintercept = 5, col = "#F8766D") +
  geom_label(data = bamstats, aes(60, sd.2, label = round(sd.2)), col = "#F8766D") +
  geom_label(data = bamstats, aes(5, max(df$num.reads), label = "5x"), col = "#F8766D") +
  facet_wrap(~library) +
  ylim(0, max(df$num.reads)) +
  theme_bw()
p

ggsave(p, filename = paste0(directory, "/Number of Mapped reads and mean percentage coverage HetTit1.0_and_HetAmer1.0_and_elegans_2024 with spp.jpeg"), height = 6, width = 12)
ggsave(p, filename = paste0(directory, "/Number of Mapped reads and mean percentage coverage HetTit1.0_and_HetAmer1.0_and_elegans_2024 with spp.pdf"), height = 6, width = 12)


# Plot with pdf with sample names to identify outliers
p <- ggplot(df) +
  geom_label(aes(coverage, num.reads, label = sample), size = 3) +
  #geom_label(aes(coverage, num.reads, col = species, label = sample), size = 3)
  #geom_vline(xintercept = 10) +
  #geom_hline(yintercept = sd.2.inter) +
  #geom_label(aes(60, sd.2.inter, label = round(sd.2.inter))) +
  xlab("Mean coverage of mapped reads") +
  ylab("Total number of mapped reads") +
  geom_hline(data = bamstats, aes(yintercept = sd.2), col = "#F8766D") +
  geom_vline(xintercept = 5, col = "#F8766D") +
  geom_label(data = bamstats, aes(60, sd.2, label = round(sd.2)), col = "#F8766D") +
  geom_label(data = bamstats, aes(5, max(df$num.reads), label = "5x"), col = "#F8766D") +
  facet_wrap(~library) +
  ylim(0, max(df$num.reads)) +
  theme_bw()
p

ggsave(p, filename = paste0(directory, "/Number of Mapped reads and mean percentage coverage HetTit1.0_and_HetAmer1.0_and_elegans July 2024 with sample names.pdf"), height = 6, width = 12)

#Creates sumset of data for each draft genome bamfiles
hetaerina_data.americana_dg <- df[df$cov10_NoR2sd&df$library=="americana_dg",]
dim(hetaerina_data.americana_dg)
any(duplicated(hetaerina_data.americana_dg$sample))
hetaerina_data.titia_dg <- df[df$cov10_NoR2sd&df$library=="titia_dg",]
dim(hetaerina_data.titia_dg)
any(duplicated(hetaerina_data.titia_dg$sample))
hetaerina_data.elegans_dg <- df[df$cov10_NoR2sd&df$library=="elegans_dg",]
dim(hetaerina_data.elegans_dg)
any(duplicated(hetaerina_data.elegans_dg$sample))

unique(df$species)
is.na(df$species)
### DATA SET 1: Hetaerina_all_ddRAD_titia_dg
write.table(file = paste0(dir_out,"Hetaerina_all_ddRAD_titia_dg",".direct_path.txt"),
            paste0(titia_directory,"/",hetaerina_data.titia_dg$sample,".paired.bam"), 
            #removes all formatting from outputed file
            sep = "\t", quote = FALSE,row.names = FALSE,col.names = FALSE)

### DATA SET 2: Hetaerina_all_ddRAD_americana_dg
write.table(file = paste0(dir_out,"Hetaerina_all_ddRAD_americana_dg",".direct_path.txt"),
            paste0(americana_directory,"/",hetaerina_data.americana_dg$sample,".paired.bam"), 
            #removes all formatting from outputed file
            sep = "\t", quote = FALSE,row.names = FALSE,col.names = FALSE)


### DATA SET 3: Hetaerina_titia_ddRAD_titia_dg
bamfiles_titia.titia_dg <- hetaerina_data.titia_dg$sample[hetaerina_data.titia_dg$species=="titia"]

write.table(file = paste0(dir_out,"Hetaerina_titia_ddRAD_titia_dg",".direct_path.txt"),
            paste0(titia_directory,"/",bamfiles_titia.titia_dg,".paired.bam"), 
            #removes all formatting from outputed file
            sep = "\t", quote = FALSE,row.names = FALSE,col.names = FALSE)

### DATA SET 4: Hetaerina_titia_ddRAD_americana_dg
bamfiles_titia.americana_dg <-  hetaerina_data.americana_dg$sample[hetaerina_data.americana_dg$species=="titia"]

write.table(file = paste0(dir_out,"Hetaerina_titia_ddRAD_americana_dg",".direct_path.txt"),
            paste0(americana_directory,"/",bamfiles_titia.americana_dg,".paired.bam"), 
            #removes all formatting from outputed file
            sep = "\t", quote = FALSE,row.names = FALSE,col.names = FALSE)

### DATA SET 5: Hetaerina_americana_ddRAD_titia_dg
bamfiles_americana.titia_dg <-  hetaerina_data.titia_dg$sample[hetaerina_data.titia_dg$species=="americana"]

write.table(file = paste0(dir_out,"Hetaerina_americana_ddRAD_titia_dg",".direct_path.txt"),
            paste0(titia_directory,"/",bamfiles_americana.titia_dg,".paired.bam"), 
            #removes all formatting from outputed file
            sep = "\t", quote = FALSE,row.names = FALSE,col.names = FALSE)

### DATA SET 6: Hetaerina_americana_ddRAD_americana_dg
bamfiles_americana.americana_dg <- hetaerina_data.americana_dg$sample[hetaerina_data.americana_dg$species=="americana"]

write.table(file = paste0(dir_out,"Hetaerina_americana_ddRAD_americana_dg",".direct_path.txt"),
            paste0(americana_directory,"/",bamfiles_americana.americana_dg,".paired.bam"), 
            #removes all formatting from outputed file
            sep = "\t", quote = FALSE,row.names = FALSE,col.names = FALSE)

### DATA SET 7: Hetaerina_all_ddRAD_elegans_dg
write.table(file = paste0(dir_out,"Hetaerina_all_ddRAD_elegans_dg",".direct_path.txt"),
            paste0(elegans_directory,"/",hetaerina_data.elegans_dg$sample,".paired.bam"), 
            #removes all formatting from outputed file
            sep = "\t", quote = FALSE,row.names = FALSE,col.names = FALSE)

