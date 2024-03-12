#Hamilton LEA and PCR plots
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
# Combined PCA and admixture plot

# CalcuLation of PCR on titia popuLations
library(patchwork)
library(ggplot2)
library(ape)
library(vcfR)
library(tidyverse)
#install.packages("rlang")
# install.packages("rgdal")
#BiocManager::install("LEA")
library(poppr)
library(adegenet)
library(LEA)
library(patchwork)
library(ggrepel)
library(scatterpie)
library(ggnewscale)


# Get library names

args <- commandArgs(trailingOnly = TRUE)
SNP.library.name <- args[1]
SNP.library.location <- args[2]

# Get task ID number
taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
taskid <- as.numeric(taskid)

print(paste("########## This is taskID", taskid, "and will use the SNP library:", SNP.library.name,"##########"))

species <- unlist(strsplit("Hetaerina_titia_ddRAD_titia_dg", split = "_"))[1:2]

species <- paste(species[1], species[2], sep = "_")

# Directory
dir.path <- paste0(SNP.library.location,SNP.library.name,"/")
filter_para <- ".all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.biSNP0_20.noX"
# Plot output file location
plot.dir <- paste0("/home/tmjj24/plots/Chapter_3/", SNP.library.name, "/")
dir.create(plot.dir)
dir.create(paste0(plot.dir, "LEA_PCA/"))

vcf.SNPs <- read.vcfR(paste0(dir.path,SNP.library.name,filter_para,".vcf.gz"),verbose = T)
# Reorder samples so they are in alphabetical order
vcf.SNPs <- vcf.SNPs[samples = sort(colnames(vcf.SNPs@gt)[-1])] 

# Get Genind
my_genind_ti_SNPs <- vcfR2genind(vcf.SNPs, sep = "/", return.alleles = TRUE)

# Set cut of value for sample removal
snp_sub <- 0.20
snp_sub_text <-"0_20"

geno <- vcfR2loci(vcf.SNPs, return.alleles = F)
geno.mat <- as.matrix(geno)
dim(geno.mat)
geno.mat[geno.mat=="1/1"] <- 2
geno.mat[geno.mat=="0/1"] <- 1
geno.mat[geno.mat=="0/0"] <- 0

# Check none of the SNPs are entirely heterozgous and remove them if they are
is.only.het <- apply(geno.mat, MARGIN = 2, function(x) gsub(paste0(unique(x), collapse = ""), pattern = "NA",replacement = "")=="1")
if(any(is.only.het)){geno.mat <- geno.mat[,-which(is.only.het)]}
# Make missing SNPs equal to "9"
geno.mat[is.na(geno.mat)] <- 9

geno.df <- data.frame(t(geno.mat))
dim(geno.df)

write.table(x = geno.df, file = paste0(dir.path,SNP.library.name,"_snp",snp_sub_text,".geno"),
            col.names = F, row.names = F, quote = F, sep = "")

#Read geno object
geno <- read.geno(paste0(dir.path,SNP.library.name,"_snp",snp_sub_text,".geno"))
dim(geno)

sample_map<- read.csv("/home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/3_Results/All samples held in Durham_v17.csv", check.names=F)
names(sample_map)

#Extracting sample information and linking it to tree labels
samples <- rownames(my_genind_ti_SNPs@tab)
sites <- data.frame(samples)
sites$samples

sites <- merge(sites, sample_map[,c("Unique.ID","sex","species","Year","Site.ID","Lat","Long","Country","Province","Ocean.drainage")], by.x = "samples", by.y = "Unique.ID", all.x = T)
any(is.na(sites))

names(sites)

# Creates data set for Country and Ocean drainage
sites$Country_Ocean <- paste(sites$Ocean.drainage, sites$Country, sep = "_")
sites$Country_Ocean <- factor(sites$Country_Ocean, 
                              levels = c("Gulf_United States", "Pacific_Mexico", "Pacific_Costa Rica", 
                                         "Carribbean_Costa Rica", "Carribbean_Belize", "Gulf_Mexico", 
                                         "Carribbean_Mexico", "Atlantic_United States", "occisa","americana"))
any(is.na(sites$Country_Ocean))
# creates a table of each site where the affinty to each grouping to the average of all samples within that site
#list of unique sites
sites$Site.ID
sites$site.sub <- paste(gsub('[[:digit:]]+', '', sites$Site.ID))
sites$site.sub[sites$site.sub=="BZ"] <- sites$Site.ID[sites$site.sub=="BZ"]
sites$site.sub[sites$site.sub=="NA"] <- sites$Site.ID[sites$site.sub=="NA"]
sites$site.sub.county <- paste(sites$site.sub, sites$Country, sep = "_")
sites$site.sub.county.drain <- paste(sites$site.sub, sites$Country, sites$Ocean.drainage, sep = "_")
sites$species.country.drainage <- paste(paste("H.",sites$species), sites$Country, sites$Ocean.drainage, sep = " - ")

cbPalette <- c("#009E73", "#0072B2", "#E69F00",  "#D55E00","#56B4E9", 
               "#CC79A7","#999999", "#F0E442", 
               "deeppink","darkred", "purple","black", "grey50")

write.table(sites, file = paste0(dir.path, "coorrd_Titia_complete.txt"))

sites <- read.table(paste0(dir.path, "coorrd_Titia_complete.txt"))

#Calculates structure for samples from K=1 to k=15

max.K <- 10
# MAY NEED TO PAUSE ONEDRIVE
# File names are becoming too Long

obj.at <- snmf(paste0(dir.path,SNP.library.name,"_snp",snp_sub_text,".geno"), K = 1:max.K, ploidy = 2, entropy = T,
               CPU = 3, project = "new", repetitions = 100, alpha = 100)
hetaerina.snmf <- load.snmfProject(file = paste0(dir.path,SNP.library.name,"_snp",snp_sub_text,".snmfProject"))
hetaerina.snmf.sum <- summary(hetaerina.snmf)

plot(hetaerina.snmf, col = "blue4", cex = 1.4, pch = 19)

cross.entropy(hetaerina.snmf, K = 6)

ce <- cbind(1:max.K, t(hetaerina.snmf.sum$crossEntropy))
colnames(ce) <- c("K", "min","mean","max")
ce <- data.frame(ce)

summary(hetaerina.snmf)

which.min(ce$mean)
ce.plot <- ggplot(ce) +
  geom_point(aes(K, mean), size = 2, shape = 19) +
  geom_errorbar(aes(x = K, ymin=min, ymax=max), width=.2,
                position=position_dodge(.9)) +
  ggtitle(paste("Minimum mean cross-entropy value K =", which.min(ce$mean))) +
  ylab("Cross-entropy") +
  theme_bw()

ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_K1-",max.K,"_snp",snp_sub_text,"cross_entropy.pdf"), plot = ce.plot)
ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_K1-",max.K,"_snp",snp_sub_text,"cross_entropy.jpg"), plot = ce.plot)
#Choose K
K <- which.min(ce$mean)
K <- 3
best <- which.min(cross.entropy(hetaerina.snmf, K = K))
qmatrix = Q(hetaerina.snmf, K = K, run = best)

pop <- unique(sites$site.sub.county)

pop
#Number of unique sites
Npop = length(pop)
Npop
# Creating qmatrix
qpop = matrix(NA, ncol = K, nrow = Npop)
qpop
coord.pop = matrix(NA, ncol = 2, nrow = Npop)
for (i in 1:length(unique(pop))){
  tmp.pop <- unique(pop)[i]
  print(sites$site.sub.county[which(sites$site.sub.county==tmp.pop)])
  print(paste("There are ", length(which(sites$site.sub.county==tmp.pop)), "samples"))
  if(length(which(sites$site.sub.county==tmp.pop)) == 1) {
    qpop[i,] <- qmatrix[sites$site.sub.county == tmp.pop,]
    coord.pop[i,] <- apply(sites[sites$site.sub.county == tmp.pop,][,c("Lat","Long")], 2, mean)
  } else {
    qpop[i,] = apply(qmatrix[sites$site.sub.county == tmp.pop,], 2, mean)
    coord.pop[i,] = apply(sites[sites$site.sub.county == tmp.pop,][,c("Lat","Long")], 2, mean)
  }
}

print("Check point 1")
q.coord.pop <- data.frame(pop, coord.pop, qpop)

q.coord.pop

print("Check point 2")

colnames(q.coord.pop) <- c("site", "Lat", "Lon", LETTERS[1:K])

print("Check point 3")
q.coord.pop$Lat <- as.numeric(q.coord.pop$Lat)
print("Check point 4")
q.coord.pop$Lon <- as.numeric(q.coord.pop$Lon)
print("Check point 5")


print("##### RUNNINNG PCA ######")
#Conduct PCA
geno2lfmm(paste0(dir.path,SNP.library.name,"_snp",snp_sub_text,".geno"), 
          paste0(dir.path,SNP.library.name,"_snp",snp_sub_text,".lfmm"), force = TRUE)
#PCA
pc <- pca(paste0(dir.path,SNP.library.name,"_snp",snp_sub_text,".lfmm"), scale = TRUE)

pc.sum <- summary(pc)
# Links PCA data to 
pca.comp <- data.frame(pc$projections[,1:6])
colnames(pca.comp) <- paste("pca", 1:6, sep = "")
pca.labs <- paste("pca", 1:4, " (",round(pc.sum[2,1:4]*100, 1), "%)", sep = "")
pca.comp$sample <- rownames(my_genind_ti_SNPs@tab)

pca.data <- cbind(sites, pca.comp)
#Random order for ploting
pca.data <- pca.data[sample(1:length(pca.data$samples)),]
#plot(pca.comp)


########
#PLOTS##
########

world.e <- data.frame(Long = c(-120,-60), Lat = c(8,38))
mex.e <- data.frame( Long = c(-105,-88), Lat = c(23,15))
mex.e.zoom <- data.frame( Long = c(-96,-94), Lat = c(18,16))
CR.e <- data.frame(Long = c(-85,-82), Lat = c(8,11))


source("/home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/3_Results/theme_black.R")
source("/home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/Chapter_3/3_Results/World/hydrobasins_extract_code.R")

# Colour blind pallette
# cbPalette <- c("#F0E442","#D55E00","#0072B2","#999999", "#E69F00" , "#56B4E9", "#009E73", "#CC79A7", "black")

p <- ggplot() +
  geom_polygon(data = worldmap, aes(long, lat, group = group), fill = "#FFE6D4", col = "black") +
  #geom_point(data = sites, aes(x = Long, y = Lat)) +
  #geom_point(data = q.coord.pop, aes(x = Long, y = Lat)) 
  geom_scatterpie(data = q.coord.pop, aes(x = Lon, y = Lat, group = site, r = 1), cols = LETTERS[1:K]) +
  #theme_bw()  +
  #xlim(c(-130,-60)) +
  #ylim(c(-0,60)) +
  #coord_sf(expand = FALSE) +
  scale_fill_manual(values = cbPalette) +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic")) +
  labs(x = "Longitude", y = "Latitude", col = "Species") +
  theme(panel.background = element_rect(fill = NA),
        panel.ontop = TRUE) +
  theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
        panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
  theme(legend.position = "none") +
  coord_map("bonne", Lat0 = 50, xlim = world.e$Long, ylim = world.e$Lat) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))  +
  ggtitle(paste("(a) K = ", K))

q <- ggplot() +
  geom_raster(data = hill.df.Mex, aes(lon, lat, fill = hill), alpha = 1) +
  #geom_polygon(data = worldmap, aes(Long, Lat, group = group), col = "black", fill = rgb(0,0,0,alpha = 0.3)) +
  scale_fill_gradientn(colours=c("#5C2700","#FFE6D4")) +
  new_scale_fill() +
  geom_sf(data = hydrobasins_geo, col = "black", fill = NA) +
  geom_sf(data = hydrorivers_geo, col = "#5C2700", lineend = "round") +
  geom_scatterpie(data = q.coord.pop, aes(x = Lon, y = Lat, r = 0.2 , group = site), cols = LETTERS[1:K]) +
  theme(legend.position="none") +
  #theme_bw()  +
  #xlim(c(-130,-60)) +
  #ylim(c(-0,60)) +
  #coord_sf(expand = FALSE) +
  scale_fill_manual(values = cbPalette) +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic")) +
  labs(x = "Longitude", y = "Latitude", col = "Species") +
  theme(panel.background = element_rect(fill = NA),
        panel.ontop = TRUE) +
  theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
        panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
  theme(legend.position = "none") +
  coord_sf(xlim = mex.e$Long, ylim = mex.e$Lat) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))+
  ggtitle("(b)")

q.zoom <- ggplot() +
  geom_raster(data = hill.df.Mex.z, aes(lon, lat, fill = hill), alpha = 1) +
  #geom_polygon(data = worldmap, aes(Long, Lat, group = group), col = "black", fill = rgb(0,0,0,alpha = 0.3)) +
  scale_fill_gradientn(colours=c("#5C2700","#FFE6D4")) +
  new_scale_fill() +
  geom_sf(data = hydrobasins_geo, col = "black", fill = NA) +
  geom_sf(data = hydrorivers_geo, col = "#5C2700", lineend = "round") +
  geom_scatterpie(data = q.coord.pop, aes(x = Lon, y = Lat, r = 0.08 , group = site), cols = LETTERS[1:K]) +
  theme(legend.position="none") +
  #theme_bw()  +
  #xlim(c(-130,-60)) +
  #ylim(c(-0,60)) +
  #coord_sf(expand = FALSE) +
  scale_fill_manual(values = cbPalette) +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic")) +
  labs(x = "Longitude", y = "Latitude", col = "Species") +
  theme(panel.background = element_rect(fill = NA),
        panel.ontop = TRUE) +
  theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
        panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
  theme(legend.position = "none") +
  coord_sf(xlim = mex.e.zoom$Long, ylim = mex.e.zoom$Lat) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))+
  ggtitle("(b)")

r <- ggplot() +
  geom_raster(data = hill.df.CR, aes(lon, lat, fill = hill), alpha = 1) +
  scale_fill_gradientn(colours=c("#5C2700","#FFE6D4")) +
  theme(legend.position = "none") +
  new_scale_fill() +
  #geom_polygon(data = worldmap, aes(Long, Lat, group = group), col = "black", fill = rgb(0,0,0,alpha = 0.3)) +
  #scale_fill_gradientn(colours=c("#d95f0e","#FFE6D4")) +
  new_scale_color() +
  geom_sf(data = hydrobasins_geo, col = "black", fill = NA) +
  geom_sf(data = hydrorivers_geo, col = "#5C2700", lineend = "round") +
  #geom_point(data = sites, aes(x = Long, y = Lat)) +
  #geom_point(data = q.coord.pop, aes(x = Long, y = Lat)) 
  geom_scatterpie(data = q.coord.pop, aes(x = Lon, y = Lat, r = 0.2 , group = site), cols = LETTERS[1:K]) +
  theme(legend.position="none") +
  #theme_bw()  +
  #xlim(c(-130,-60)) +
  #ylim(c(-0,60)) +
  #coord_sf(expand = FALSE) +
  scale_fill_manual(values = cbPalette) +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic")) +
  labs(x = "Longitude", y = "Latitude", col = "Species") +
  theme(panel.background = element_rect(fill = NA),
        panel.ontop = TRUE) +
  theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
        panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
  theme(legend.position = "none") +
  coord_sf(xlim = CR.e$Long, ylim = CR.e$Lat) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
  ggtitle("(c)")

best <- which.min(cross.entropy(hetaerina.snmf, K = K))
qmatrix = Q(hetaerina.snmf, K = K, run = best)
qtable <- cbind(rep(sites$samples,K), rep(sites$Lat,K), rep(1:K, each = length(sites$samples)), c(qmatrix[,1:K]))
qtable <-  data.frame(qtable)
colnames(qtable) <- c("sample","Lat","Qid", "Q")

qtable$sample <-  factor(qtable$sample, levels = sites$sample[order(paste(sites$Country_Ocean, sites$Lat))])
qtable$Q <- as.numeric(qtable$Q)
qtable$sample_Lat <- paste(qtable$Lat, qtable$sample, sep = "_")

#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "black")
#cbPalette <- c("#F0E442","#D55E00","#0072B2","#999999", "#E69F00" , "#56B4E9", "#009E73", "#CC79A7", "black")

v <- ggplot(qtable)+
  geom_bar(stat="identity", aes(sample, Q, fill = Qid,), position = "stack", width = 1, col = "black") +
  scale_fill_manual(values = cbPalette) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
  ylab(label = paste("K =", K))+
  ggtitle("(d)")

#Creates multiple K barcharts

#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "black")

s <- list()

max.K <- min(c(max.K, 6))

for(i in 2:max.K){
  best <- which.min(cross.entropy(hetaerina.snmf, K = i))
  qmatrix = Q(hetaerina.snmf, K = i, run = best)
  qtable <- cbind(rep(sites$samples,i), rep(sites$Lat,i), rep(1:i, each = length(sites$samples)), c(qmatrix[,1:i]))
  qtable <-  data.frame(qtable)
  colnames(qtable) <- c("sample","Lat","Qid", "Q")
  
  qtable$sample <-  factor(qtable$sample, levels = sites$sample[order(paste(sites$Country_Ocean, sites$Lat))])
  qtable$Q <- as.numeric(qtable$Q)
  qtable$sample_Lat <- paste(qtable$Lat, qtable$sample, "_")
  
  s[[i]] <- ggplot(qtable)+
    geom_bar(stat="identity", aes(sample, Q, fill = Qid,), position = "stack", width = 1, col = "black") +
    scale_fill_manual(values = cbPalette) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none") +
    theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) +
    ylab(label = paste("K =",i))
  
}


s[[max.K]] <- ggplot(qtable)+
  geom_bar(stat="identity", aes(sample, Q, fill = Qid,), position = "stack", width = 1, col = "black") +
  scale_fill_manual(values = cbPalette) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
  ylab(label = paste("K =",max.K))


plot2 <- s[[2]]
for(i in 3:max.K){
  plot2 <- plot2 / s[[i]]
}
s <- plot2

## PCA plots

# cbPalette <- c("#999999", "#009E73", "#0072B2", "#E69F00", "#F0E442", "#56B4E9", "#D55E00", "#CC79A7", "black")
plot(1:length(cbPalette), col = cbPalette, pch = 19, cex =8)
x <- ggplot(pca.data) +
  #geom_point(aes(pca1, pca2, colour = Country_Ocean), size = 4) +
  #scale_colour_manual(values = cbPalette, name = "Region") +
  geom_text(aes(pca1, pca2, colour = species.country.drainage, label = samples), alpha = .8) +
  #geom_label_repel(aes(pca1, pca2, label = samples, colour = species), size = 0.5
  #                 ,box.padding = 0.1,label.padding = 0.1, min.segment.length = 10) +
  xlab(pca.labs[1]) +
  ylab(pca.labs[2]) + 
  theme_bw() 

ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_",snp_sub_text,"_pca1_2.pdf"), plot=x, height=10, width=10)
ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_",snp_sub_text,"_pca1_2.jpg"), plot=x, height=10, width=10)

#theme(legend.position = "none")
y <- ggplot(pca.data) +
  geom_point(aes(pca1, pca2, colour = species.country.drainage), size = 6, alpha = 0.5) +
  scale_colour_manual(values = cbPalette, name = "Region") +
  xlab(pca.labs[1]) +
  ylab(pca.labs[2]) +
  theme(legend.position = c(0.21, 0.8)) 
z <- ggplot(pca.data) +
  geom_point(aes(pca3, pca4, colour = species.country.drainage), size = 6, alpha = 0.5) + 
  scale_colour_manual(values = cbPalette, name = "Region") +
  xlab(pca.labs[3]) +
  ylab(pca.labs[4]) +
  theme(legend.position = "none") 
w <- ggplot(pca.data) +
  geom_point(aes(pca5, pca6, colour = species.country.drainage), size = 6, alpha = 0.5) + 
  scale_colour_manual(values = cbPalette, name = "Region") +
  #xlab(pca.labs[5]) +
  #ylab(pca.labs[6]) +
  theme(legend.position = "none")
ww <- ggplot(pca.data) +
  geom_label(aes(pca5, pca6, label = samples, colour = species.country.drainage), size = 6, alpha = 0.5) +
  scale_colour_manual(values = cbPalette, name = "Region") +
  xlab(pca.labs[5]) +
  ylab(pca.labs[6]) +
  theme(legend.position = "none")
a <- ggplot() +
  geom_polygon(data = worldmap, aes(long, lat, group = group), fill = "#FFE6D4", col = "black") +
  geom_point(data = pca.data, aes(Long, Lat, colour = species.country.drainage), size = 6) +
  scale_colour_manual(values = cbPalette, name = "Region") +
  scale_fill_manual(values = cbPalette) +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic")) +
  labs(x = "Longitude", y = "Latitude", col = "Species") +
  theme(panel.background = element_rect(fill = NA),
        panel.ontop = TRUE) +
  theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
        panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
  theme(legend.position = "none") +
  coord_map("bonne", Lat0 = 50, xlim = world.e$Long, ylim = world.e$Lat)


plot1 <- (p+q+r)/v + plot_layout(heights = c(4, 1))
plot2 <- y+z + w + a + plot_layout(nrow = 2)

ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_K",K,"snp",snp_sub_text,"_complete_DP10_basins_patchwork.pdf"), plot=plot1, height=7, width=15)
ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_K",K,"snp",snp_sub_text,"_complete_DP10_basins_patchwork.jpg"), plot=plot1, height=7, width=15)
ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_K",K,"snp",snp_sub_text,"_complete_DP10_basins_patchwork_Mex_z.pdf"), plot=q.zoom, height=10, width=10)
ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_K",K,"snp",snp_sub_text,"_complete_DP10_basins_patchwork_Mex_z.jpg"), plot=q.zoom, height=10, width=10)
ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_barplot_1-",max.K,"_snp",snp_sub_text,"_complete_DP10_basins_patchwork.pdf"), plot=s, height=20, width=15)
ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_LEA_barplot_1-",max.K,"_snp",snp_sub_text,"_complete_DP10_basins_patchwork.jpg"), plot=s, height=20, width=15)
ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_PCA_complete","snp",snp_sub_text,"_DP10_basins_patchwork.pdf"), height=9, width=10, plot=plot2)
ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_PCA_complete","snp",snp_sub_text,"_DP10_basins_patchwork.jpg"), height=9, width=10, plot=plot2)
ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_PCA_1-2-complete","snp",snp_sub_text,"_DP10_basins_patchwork.pdf"), height=8, width=5, plot=(y/a))
ggsave(paste0(plot.dir, "LEA_PCA/", SNP.library.name,"_PCA_1-2-complete","snp",snp_sub_text,"_DP10_basins_patchwork.jpg"), height=8, width=5, plot=(y/a))
