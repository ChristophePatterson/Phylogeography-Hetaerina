# Script for extracting and creating plots for final manuscript of Chapter 3

# Calculation of PCR on titia populations
library(ggplot2)
library(ape)
library(vcfR)
#install.packages("rlang")
# install.packages("rgdal")
# install.packages("BiocManager")
# BiocManager::install("LEA")
library(poppr)
library(adegenet)
library(LEA)
library(patchwork)
library(ggrepel)
library(tidyverse)
# library(PBSmapping)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

###### H titia phylogeographu plots ######

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Output file location
SNP.library.name <- "Hetaerina_titia_ddRAD_titia_dg"
species <- "titia"
dir.path <- paste0("4_Manuscript/data/SNP_libraries/",SNP.library.name,"/")
filter_para <- ".all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.biSNP0_20"
snp_sub_text <- "0_20"
# Plot output file location
plot.dir <- paste0("4_Manuscript/plots/",SNP.library.name,"/")
dir.create(plot.dir)


## Running analysis name
analysis.name <- "H_titia_complete_snp_noX"

vcf <- read.vcfR(paste0(dir.path,SNP.library.name,filter_para,".vcf.gz"))

colnames(vcf@gt)
# Removes . from varient names as genind does like that
vcf@fix[,1] <- gsub(vcf@fix[,1], pattern = "\\.", replacement = "_")

# Removes SNPs on the X chromosome
X_chrom <- "HetTit1_0_p_scaff-12-96647824"
table(vcf@fix[,1]!=X_chrom)
vcf.SNPs <- vcf[vcf@fix[,1]!=X_chrom,]

# Stats
myMiss <- apply(extract.gt(vcf.SNPs), MARGIN = 2, function(x){ sum(is.na(x)) })
myMiss <- data.frame(sample = names(myMiss), per.gt = myMiss/nrow(vcf)*100)

summary(myMiss)
ggplot(myMiss) +
  geom_violin(aes(x = "X", y = per.gt), outlier.colour = NA) +
  geom_jitter(aes(x = "X", y = per.gt), height = 0, width = 0.2)

# Get coverage statistics 
DP.mat <- extract.gt(vcf.SNPs, element = "DP")
Samp.DP <- apply(DP.mat, MARGIN = 2, function(x) mean(as.numeric(x), na.rm = T))
SNP.DP <- apply(DP.mat, MARGIN = 1, function(x) mean(as.numeric(x), na.rm = T))

c(mean(Samp.DP,na.rm =T), median(Samp.DP,na.rm =T), range(Samp.DP,na.rm = T))
c(mean(SNP.DP,na.rm =T), median(SNP.DP,na.rm =T), range(SNP.DP,na.rm = T))

# Creaete genind object
my_genind_ti_SNPs <- vcfR2genind(vcf.SNPs, sep = "/", return.alleles = TRUE)

# Create geno object
geno <- vcfR2loci(vcf.SNPs, return.alleles = F)
geno.mat <- as.matrix(geno)


# Convert into geno format 2 for alternative homozygous, 1 for he
# 2 for alternative homozygous
geno.mat[geno.mat=="1/1"] <- 2
# 1 for heterozygous
geno.mat[geno.mat=="0/1"] <- 1
# 0 for homozygous genome assembly allele
geno.mat[geno.mat=="0/0"] <- 0
# Checks if any data points are entirely heterozygous
is.only.het <- apply(geno.mat, MARGIN = 2, function(x) gsub(paste0(unique(x), collapse = ""), pattern = "NA",replacement = "")=="1")
if(any(is.only.het)){geno.mat <- geno.mat[,-which(is.only.het)]}
# Missing data coded as 9
geno.mat[is.na(geno.mat)] <- 9

# Convert to data frame to save
geno.df <- data.frame(t(geno.mat))

write.table(x = geno.df, file = paste0(dir.path,analysis.name,snp_sub_text,".geno"),
            col.names = F, row.names = F, quote = F, sep = "")

#Read geno object
geno <- read.geno(paste0(dir.path, analysis.name,snp_sub_text,".geno"))
dim(geno)


# Read in data on all samples
sample_map <- read.csv("4_Manuscript/data/Hetaerina_sample_data.csv")
#Extracting sample information and linking it to tree labels
samples <- rownames(my_genind_ti_SNPs@tab)
sites <- data.frame(samples)
sites <- merge(sites, sample_map, by.x = "samples", by.y = "Unique.ID")

# Creates data set for Country and Ocean drainage
sites$Country_Ocean <- paste(sites$Ocean.drainage, sites$Country, sep = "_")
sites$Country_Ocean <- factor(sites$Country_Ocean, 
                              levels = c("Gulf_United States", "Pacific_Mexico", "Pacific_Costa Rica", 
                                         "Carribbean_Costa Rica", "Carribbean_Belize", "Gulf_Mexico", 
                                         "Carribbean_Mexico", "Atlantic_United States", "occisa","americana"))
any(is.na(sites$Country_Ocean))
# creates a table of each site where the affinty to each grouping to the average of all samples within that site
#list of unique sites
sites$site.sub <- paste(gsub('[[:digit:]]+', '', sites$Site.ID))
# Because site is labelled NA (I DIDNT DO THIS!) need to override name back to NA01
sites$site.sub[sites$site.sub=="NA"] <- "NA01"
sites$site.sub[sites$site.sub=="BZ"] <- sites$Site.ID[sites$site.sub=="BZ"]
sites$site.sub.county <- paste(sites$site.sub, sites$Country, sep = "_")
sites$site.sub.county.drain <- paste(sites$site.sub, sites$Country, sites$Ocean.drainage, sep = "_")

# Colourblind colour pallette
cbPalette <- c("#999999", "#009E73", "#0072B2", "#E69F00", "#F0E442", "#56B4E9", "#D55E00", "#CC79A7")

# Write and read back in sample data
write.table(sites, file = paste0(dir.path, "coorrd_Titia_complete",analysis.name ,".txt"))
sites <- read.table(paste0(dir.path, "coorrd_Titia_complete",analysis.name ,".txt"))

sites[sites$Site.ID=="CUAJ01",]
#Calculates structure for samples from K=1 to k=10
max.K <- 10
# MAY NEED TO PAUSE ONEDRIVE
# obj.at <- snmf(paste0(dir.path, analysis.name,snp_sub_text,".geno"), K = 1:max.K, ploidy = 2, entropy = T,
#             CPU = 2, project = "new", repetitions = 20, alpha = 100)
# If swapping between home and work computer need to change the .smfproject file to \\homeblue02\tmjj24\My_Documents\Github\
#or to 
getwd()
titia.snmf <- load.snmfProject(file = paste0(dir.path, analysis.name,snp_sub_text,".snmfProject"))
titia.snmf.sum <- summary(titia.snmf)

plot(titia.snmf, col = "blue4", cex = 1.4, pch = 19)

# cross.entropy(titia.snmf, K = 6)

ce <- cbind(1:max.K, t(titia.snmf.sum$crossEntropy))
colnames(ce) <- c("K", "min","mean","max")
ce <- data.frame(ce)

summary(titia.snmf)

which.min(ce$mean)
ce.plot <- ggplot(ce) +
  geom_point(aes(K, mean), size = 2, shape = 19) +
  geom_errorbar(aes(x = K, ymin=min, ymax=max), width=.2,
                position=position_dodge(.9)) +
  ggtitle(paste("Minimum mean cross-entropy value K =", which.min(ce$mean))) +
  ylab("Cross-entropy") +
  theme_bw()
ce.plot

ggsave(paste0(plot.dir,"H_titia_LEA_complete_snp_K1-8",snp_sub_text,"cross_entropy.pdf"), plot = ce.plot)
ggsave(paste0(plot.dir,"H_titia_LEA_complete_snp_K1-8",snp_sub_text,"cross_entropy.jpg"), plot = ce.plot)
#Choose K
K <- which.min(ce$mean)
K <- 3
best <- which.min(cross.entropy(titia.snmf, K = K))
qmatrix = Q(titia.snmf, K = K, run = best)
# Get list of unique populations
pop <- unique(sites$site.sub.county.drain)
length(unique(sites$Lat))
length(unique(sites$site.sub.county.drain))

#Number of unique sites
Npop = length(pop)

qpop = matrix(NA, ncol = K, nrow = Npop)
coord.pop = matrix(NA, ncol = 2, nrow = Npop)
for (i in 1:length(unique(pop))){
  tmp.pop <- unique(pop)[i]
  print(sites$site.sub.county.drain[which(sites$site.sub.county.drain==tmp.pop)])
  if(length(which(sites$site.sub.county.drain==tmp.pop)) == 1) {
    qpop[i,] <- qmatrix[sites$site.sub.county.drain == tmp.pop,]
    coord.pop[i,] <- apply(sites[sites$site.sub.county.drain == tmp.pop,][,c("Lat","Long")], 2, mean)
  } else {
    qpop[i,] = apply(qmatrix[sites$site.sub.county.drain == tmp.pop,], 2, mean)
    coord.pop[i,] = apply(sites[sites$site.sub.county.drain == tmp.pop,][,c("Lat","Long")], 2, mean)
  }
}

# plot(coord.pop[,2:1])
q.coord.pop <- data.frame(pop, coord.pop, qpop)
colnames(q.coord.pop) <- c("site", "lat", "long", LETTERS[1:K])
q.coord.pop$lat <- as.numeric(q.coord.pop$lat)
q.coord.pop$long <- as.numeric(q.coord.pop$long)
q.coord.pop$num_samples <- apply(X = q.coord.pop, MARGIN = 1, function(x) length(which(sites$site.sub.county.drain==x["site"])))

# Figure out what number each region has been assigned to (As it varies between sNMF runs)
Pac.clust <- which.max(q.coord.pop[q.coord.pop$site=="ZANA_Mexico_Pacific",LETTERS[1:K]])
SAtl.clust <- which.max(q.coord.pop[q.coord.pop$site=="ESRB_Costa Rica_Carribbean",LETTERS[1:K]])
NAtl.clust <- which.max(q.coord.pop[q.coord.pop$site=="HCAR_United States_Gulf",LETTERS[1:K]])

#Conduct PCA
geno2lfmm(input.file = paste0(dir.path, analysis.name,snp_sub_text,".geno"), 
          output.file = paste0(dir.path, analysis.name,snp_sub_text,".lfmm"), force = TRUE)
#PCA
pc <- pca(paste0(dir.path, analysis.name,snp_sub_text,".lfmm"), scale = TRUE)

pc.sum <- summary(pc)
pc.sum
# Plot eigenvalues.
plot(pc, lwd=5, col="red",xlab=("PCs"),ylab="eigen")

tracy.widom(pc)
plot(tracy.widom(pc)$pvalues)
# Links PCA data to 
pca.comp <- data.frame(pc$projections[,1:6])
colnames(pca.comp) <- paste("pca", 1:6, sep = "")
pca.labs <- paste("pca", 1:6, " (",round(as.numeric(pc.sum[2,1:6])*100, 1), "%)", sep = "")
pca.comp$sample <- rownames(my_genind_ti_SNPs@tab)

pca.data <- cbind(sites, pca.comp)
#Random order for ploting
pca.data <- pca.data[sample(1:length(pca.data$samples)),]
# plot(pca.comp)
pca.data$Country_Ocean <- factor(pca.data$Country_Ocean, 
                                 levels = c("Gulf_United States", "Pacific_Mexico", "Pacific_Costa Rica", 
                                            "Carribbean_Costa Rica", "Carribbean_Belize", "Gulf_Mexico", 
                                            "Carribbean_Mexico", "Atlantic_United States", "occisa","americana"))
names(pca.data)

# Plot basic PCA data
ggplot(pca.data) +
  geom_point(aes(as.numeric(pca1), as.numeric(pca2), colour = Country_Ocean), size = 2)

## Number of SNP and samples

pc
dim(geno)
min(ce$mean)

# # # # # # # # # # #
###### PLOTS ######
# # # # # # # # # # #

# Outlines for each plot
world.e <- data.frame(Long = c(-120,-60), Lat = c(8,38))
mex.e <- data.frame( Long = c(-97.5,-86.3), Lat = c(19.5,14.3))
#mex.e <- data.frame( Long = c(-100,-91), Lat = c(22,14.4))
mex.e.zoom <- data.frame( Long = c(-95.5,-94.2), Lat = c(17.5,16.1))
CR.e <- data.frame(Long = c(-85.25,-82.5), Lat = c(8,11))


# World projection
crs_text <- "+proj=laea +x_0=0 +y_0=0 +lon_0=-74 +lat_0=40"
crs_text <- paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=",360+mean(sites$Long), " +lat_0=",mean(sites$Lat))

# Download terrrian and river data
source("4_Manuscript/hydrobasins_extract_code_sf.R")
library(scatterpie)
library(ggnewscale)

cbPalette <- c("#F0E442","#D55E00" , "#0072B2", "#999999","#E69F00" , "#56B4E9", "#009E73", "#CC79A7", "black")
# Colour scheme
# Order Het cols based on what regions where assign in sNMF
het.cols <- c("#AF0F09","#E5D9BA","#3E3C3A")
het.cols[Pac.clust] <- "#AF0F09"
het.cols[SAtl.clust] <- "#E5D9BA"
het.cols[NAtl.clust] <- "#3E3C3A"
#het.cols <- c("#3E3C3A","#AF0F09","#E5D9BA", "#694438")
#amer.cols <- c("#848A39","#726230","#920E02","#AA9599")

plot(1:4, col = het.cols, cex = 20, pch = 19)

## Creating bounding boxes for plots
#mexico
mex.e.e <- mex.e[1,]
mex.e.e[2,] <- c(mex.e[1,1],mex.e[2,2])
mex.e.e[3,] <- mex.e[2,]
mex.e.e[4,] <- c(mex.e[2,1],mex.e[1,2])
## Creating bounding boxes for plots
#mexico zoom
mex.e.zoom.e <- mex.e.zoom[1,]
mex.e.zoom.e[2,] <- c(mex.e.zoom[1,1],mex.e.zoom[2,2])
mex.e.zoom.e[3,] <- mex.e.zoom[2,]
mex.e.zoom.e[4,] <- c(mex.e.zoom[2,1],mex.e.zoom[1,2])
#costa rica
CR.e.e <- CR.e[1,]
CR.e.e[2,] <- c(CR.e[1,1],CR.e[2,2])
CR.e.e[3,] <- CR.e[2,]
CR.e.e[4,] <- c(CR.e[2,1],CR.e[1,2])

# Shifts all the Costa Rica points so that they do not overlap
# Min dist in decimal degrees that populations whould be assigned 
min.dist <- 0.35
# Subset population to only those in costa rica
set.seed(3)
q.coord.pop.CR <- q.coord.pop[grep(q.coord.pop$site, pattern = "Costa Rica"),]
# Assigned new location col
q.coord.pop.CR$long.new <- q.coord.pop.CR$long
q.coord.pop.CR$lat.new <- q.coord.pop.CR$lat

# Loop through each point
for(i in 1:length(q.coord.pop.CR$site)){
  #Sets distance
  dist.temp <- NA
  #Calcualted the distance between point i all other points
  temp.dist <- raster::pointDistance(c(q.coord.pop.CR$long.new[i], q.coord.pop.CR$lat.new[i]), cbind(q.coord.pop.CR$long.new, q.coord.pop.CR$lat.new), lonlat = T)/111/1000
  # Removes distance to itself which is always 0
  temp.dist[i] <- NA
  # Asks which is the closets point
  min.dist.i <- which(temp.dist==min(temp.dist, na.rm = T))
  #Prints check
  print(paste("For point", i, "the min distance is point ",min.dist.i,":",min(temp.dist, na.rm = T)))
  # Sets the distance for point to be shfited by to the min distance
  shift.dist <- min.dist
  # Sets break point counter
  break.point <- 1
  long.shift <- sample(seq(0,shift.dist,0.001),1)
  # Loops while the min distance is smaller than the min distance
  while(min(temp.dist, na.rm = T)<min.dist){
    # If point is on the pacific slope move southwest by random amount
    if(grepl(q.coord.pop.CR$site[i], pattern = "Pacific")){
      q.coord.pop.CR$long.new[i] <- q.coord.pop.CR$long.new[i]-long.shift
      q.coord.pop.CR$lat.new[i] <- q.coord.pop.CR$lat.new[i]-sample(seq(0,shift.dist,0.001),1)
    }
    # If point is on the Atlantic slope move Northeast by random amount
    if(grepl(q.coord.pop.CR$site[i], pattern = "Carribbean")){
      q.coord.pop.CR$lat.new[i] <- q.coord.pop.CR$lat.new[i]+sample(seq(0,shift.dist,0.001),1)
      q.coord.pop.CR$long.new[i] <- q.coord.pop.CR$long.new[i]+long.shift
    }
    # Recalculate distance to all points
    temp.dist <- raster::pointDistance(c(q.coord.pop.CR$long.new[i], q.coord.pop.CR$lat.new[i]), cbind(q.coord.pop.CR$long.new[-i], q.coord.pop.CR$lat.new[-i]), lonlat = T)/111/1000
    break.point <- break.point+1
    if(break.point>10000) break
  }
  
}

ggplot(q.coord.pop.CR) +
  geom_segment(aes(x = q.coord.pop.CR$long , y = q.coord.pop.CR$lat , xend = q.coord.pop.CR$long.new, yend = q.coord.pop.CR$lat.new)) +
  geom_point(aes(x = q.coord.pop.CR$long , y = q.coord.pop.CR$lat)) +
  geom_point(aes(x = q.coord.pop.CR$long.new , y = q.coord.pop.CR$lat.new), color = "red") +
  geom_polygon(data = CR.e.e, aes(x = Long, y = Lat), linewidth = 1,  fill = NA, colour = "black")

# Shifts all the Mexico points so that they do not overlap
# Min dist in decimal degrees that populations whould be assigned 
min.dist <- 0.41
# Subset population to only those in costa rica
q.coord.pop.Mx <- q.coord.pop[grepl(q.coord.pop$site, pattern = paste(c("Mexico", "Belize"), collapse = "|")),]
# Assigned new location col
q.coord.pop.Mx$long.new <- q.coord.pop.Mx$long
q.coord.pop.Mx$lat.new <- q.coord.pop.Mx$lat
set.seed(5)
# Loop through each point
for(i in 1:length(q.coord.pop.Mx$site)){
  #Sets distance
  dist.temp <- NA
  #Calcualted the distance between point i all other points
  temp.dist <- raster::pointDistance(c(q.coord.pop.Mx$long.new[i], q.coord.pop.Mx$lat.new[i]), cbind(q.coord.pop.Mx$long.new, q.coord.pop.Mx$lat.new), lonlat = T)/111/1000
  # Removes distance to itself which is always 0
  temp.dist[i] <- NA
  # Asks which is the closets point
  min.dist.i <- which(temp.dist==min(temp.dist, na.rm = T))
  #Prints check
  print(paste("For point", i, "the min distance is point ",min.dist.i,":",min(temp.dist, na.rm = T)))
  # Sets the distance for point to be shfited by to the min distance
  shift.dist <- min.dist
  # Sets break point counter
  break.point <- 1
  long.shift <- sample(seq(0,shift.dist,0.001),1)
  # Loops while the min distance is smaller than the min distance
  while(min(temp.dist, na.rm = T)<min.dist){
    # If point is on the pacific slope move southwest by random amount
    if(grepl(q.coord.pop.Mx$site[i], pattern = "Pacific")){
      q.coord.pop.Mx$long.new[i] <- q.coord.pop.Mx$long.new[i]-long.shift
      q.coord.pop.Mx$lat.new[i] <- q.coord.pop.Mx$lat.new[i]-sample(seq(0,shift.dist,0.001),1)
    }
    # If point is on the Atlantic slope move Northeast by random amount
    if(grepl(q.coord.pop.Mx$site[i], pattern = paste(c("Carribbean", "Gulf"), collapse = "|"))){
      q.coord.pop.Mx$lat.new[i] <- q.coord.pop.Mx$lat.new[i]+sample(seq(0,shift.dist,0.001),1)
      q.coord.pop.Mx$long.new[i] <- q.coord.pop.Mx$long.new[i]+long.shift
    }
    # Recalculate distance to all points
    temp.dist <- raster::pointDistance(c(q.coord.pop.Mx$long.new[i], q.coord.pop.Mx$lat.new[i]), cbind(q.coord.pop.Mx$long.new[-i], q.coord.pop.Mx$lat.new[-i]), lonlat = T)/111/1000
    break.point <- break.point+1
    if(break.point>10000) break
  }
  
}

ggplot(q.coord.pop.Mx) +
  geom_segment(aes(x = q.coord.pop.Mx$long , y = q.coord.pop.Mx$lat , xend = q.coord.pop.Mx$long.new, yend = q.coord.pop.Mx$lat.new)) +
  geom_point(aes(x = q.coord.pop.Mx$long , y = q.coord.pop.Mx$lat)) +
  geom_point(aes(x = q.coord.pop.Mx$long.new , y = q.coord.pop.Mx$lat.new), color = "red") + 
  geom_polygon(data = mex.e.e, aes(x = Long, y = Lat), linewidth = 1,  fill = NA, colour = "black")


# Custom plots

worldmap_proj <- st_transform(worldmap, crs = crs_text)

p <- ggplot() +
  geom_sf(data = worldmap_proj, fill = "#FFE6D4", col = "black", linewidth = 0.5) +
  geom_polygon(data = mex.e.e, aes(x = Long, y = Lat), size = 1.2, fill = NA, colour = "black") +
  geom_polygon(data = CR.e.e, aes(x = Long, y = Lat), size = 1.2,  fill = NA, colour = "black") +
  #geom_point(data = sites, aes(x = long, y = lat)) +
  #geom_point(data = q.coord.pop, aes(x = long, y = lat)) 
  geom_scatterpie(data = q.coord.pop, aes(x = long, y = lat, group = site, r = (sqrt(num_samples)/pi)*2), cols = LETTERS[1:K]) +
  geom_text(data = mex.e.e[4,],aes(Long + 3, Lat, label = "(b)"), size = 5) +
  geom_text(data = CR.e.e[3,],aes(Long + 2, Lat, label = "(c)"), size = 5) +
  #theme_bw()  +
  #xlim(c(-130,-60)) +
  #ylim(c(-0,60)) +
  #coord_sf(expand = FALSE) +
  scale_fill_manual(values = het.cols) +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic")) +
  labs(x = "Longitude", y = "Latitude", col = "Species") +
  theme(panel.background = element_rect(fill = NA),
        panel.ontop = TRUE) +
  theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
        panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
  theme(legend.position = "none") +
  coord_sf(xlim = range(q.coord.pop$long)+c(-5,+5), ylim =range(q.coord.pop$lat)+c(-2,+2),default_crs = st_crs(4326))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))  +
  ggtitle(paste("K = ", K))

p

q <- ggplot() +
  geom_raster(data = hill.df.Mex, aes(lon, lat, fill = hill), alpha = 1) +
  #geom_polygon(data = worldmap, aes(long, lat, group = group), col = "black", fill = rgb(0,0,0,alpha = 0.3)) +
  geom_polygon(data = mex.e.zoom.e, aes(x = Long, y = Lat), size = 1.2,  fill = NA, colour = "black") +
  scale_fill_gradientn(colours=c("#5C2700","#FFE6D4")) +
  new_scale_fill() +
  geom_sf(data = hydrobasins_geo, col = "black", fill = NA, linewidth = 0.5) +
  geom_sf(data = hydrorivers_geo, col = "#5C2700", lineend = "round") +
  geom_segment(aes(x = q.coord.pop.Mx$long , y = q.coord.pop.Mx$lat , xend = q.coord.pop.Mx$long.new, yend = q.coord.pop.Mx$lat.new),
               linewidth = 1.2, lineend =  "round") +
  geom_text(data = mex.e.zoom.e[4,],aes(Long+.8, Lat+.2, label = "(Fig.5b)"), size = 5) +
  geom_scatterpie(data = q.coord.pop.Mx, aes(x = long.new, y = lat.new, r = (sqrt(num_samples)/pi)/3 , group = site), cols = LETTERS[1:K]) +
  theme(legend.position="none") +
  #theme_bw()  +
  #xlim(c(-130,-60)) +
  #ylim(c(-0,60)) +
  #coord_sf(expand = FALSE) +
  scale_fill_manual(values = het.cols) +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic")) +
  labs(x = "Longitude", y = "Latitude", col = "Species") +
  theme(panel.background = element_rect(fill = NA),
        panel.ontop = TRUE) +
  theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
        panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
  theme(legend.position = "none") +
  coord_sf(xlim = mex.e$Long, ylim = mex.e$Lat) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))

q.zoom <- ggplot() +
  geom_raster(data = hill.df.Mex.z, aes(lon, lat, fill = hill), alpha = 1) +
  #geom_polygon(data = worldmap, aes(long, lat, group = group), col = "black", fill = rgb(0,0,0,alpha = 0.3)) +
  scale_fill_gradientn(colours=c("#5C2700","#FFE6D4")) +
  new_scale_fill() +
  geom_sf(data = hydrobasins_geo, col = "black", fill = NA, linewidth = 0.5) +
  geom_sf(data = hydrorivers_geo, col = "#5C2700", lineend = "round") +
  geom_scatterpie(data = q.coord.pop, aes(x = long, y = lat, r = (sqrt(num_samples)/pi)/6 , group = site), cols = LETTERS[1:K]) +
  theme(legend.position="none") +
  #theme_bw()  +
  #xlim(c(-130,-60)) +
  #ylim(c(-0,60)) +
  #coord_sf(expand = FALSE) +
  scale_fill_manual(values = het.cols) +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic")) +
  labs(x = "Longitude", y = "Latitude", col = "Species") +
  theme(panel.background = element_rect(fill = NA),
        panel.ontop = TRUE) +
  theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
        panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
  theme(legend.position = "none") +
  coord_sf(xlim = mex.e.zoom$Long, ylim = mex.e.zoom$Lat) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))


r <- ggplot() +
  geom_raster(data = hill.df.CR, aes(lon, lat, fill = hill), alpha = 1) +
  scale_fill_gradientn(colours=c("#5C2700","#FFE6D4")) +
  theme(legend.position = "none") +
  new_scale_fill() +
  #geom_polygon(data = worldmap, aes(long, lat, group = group), col = "black", fill = rgb(0,0,0,alpha = 0.3)) +
  #scale_fill_gradientn(colours=c("#d95f0e","#FFE6D4")) +
  new_scale_color() +
  geom_sf(data = hydrobasins_geo, col = "black", fill = NA, linewidth = 0.5) +
  geom_sf(data = hydrorivers_geo, col = "#5C2700", lineend = "round") +
  #geom_point(data = sites, aes(x = long, y = lat)) +
  #geom_point(data = q.coord.pop, aes(x = long, y = lat)) 
  geom_segment(aes(x = q.coord.pop.CR$long , y = q.coord.pop.CR$lat , xend = q.coord.pop.CR$long.new, yend = q.coord.pop.CR$lat.new)
               , linewidth = 1.2, lineend =  "round") +
  geom_scatterpie(data = q.coord.pop.CR, aes(x = long.new, y = lat.new, r = (sqrt(num_samples)/pi)/3 , group = site), cols = LETTERS[1:K]) +
  theme(legend.position="none") +
  #theme_bw()  +
  #xlim(c(-130,-60)) +
  #ylim(c(-0,60)) +
  #coord_sf(expand = FALSE) +
  scale_fill_manual(values = het.cols) +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic")) +
  labs(x = "Longitude", y = "Latitude", col = "Species") +
  theme(panel.background = element_rect(fill = NA),
        panel.ontop = TRUE) +
  theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
        panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
  theme(legend.position = "none") +
  coord_sf(xlim = CR.e$Long, ylim = CR.e$Lat) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))



best <- which.min(cross.entropy(titia.snmf, K = K))
qmatrix = Q(titia.snmf, K = K, run = best)
qtable <- cbind(rep(sites$samples,K), rep(sites$Lat,K), rep(1:K, each = length(sites$samples)), c(qmatrix[,1:K]))
qtable <-  data.frame(qtable)
colnames(qtable) <- c("sample","lat","Qid", "Q")


# New ocean drainage 
sites$Ocean.drainage.new <- sites$Ocean.drainage
sites$Ocean.drainage.new[sites$Ocean.drainage=="Atlantic"] <- "Gulf"
qtable$sample <-  factor(qtable$sample, levels = sites$sample[order(paste(sites$species,sites$Ocean.drainage.new, sites$Country, sites$Long, sites$site.name))])
qtable$Q <- as.numeric(qtable$Q)
qtable$sample_lat <- paste(qtable$lat, qtable$sample)

cbPalette <- c("#F0E442","#D55E00","#0072B2","#999999", "#E69F00" , "#56B4E9", "#009E73", "#CC79A7", "black")
plot(1:length(cbPalette), col = cbPalette, pch = 19, cex = 10)

CUAJ_sample_position <- range(grep("CUAJ", sites$sample[order(paste(sites$species,sites$Ocean.drainage.new, sites$Country, sites$Long, sites$site.name))]))+c(-0.5,0.5)

v <- ggplot(qtable)+
  geom_bar(stat="identity", aes(sample, Q, fill = Qid,), position = "stack", width = 1, col = "black") +
  geom_segment(aes(y = -0.03, yend = -0.03, 
                   x = CUAJ_sample_position[1], xend = CUAJ_sample_position[2]),
               linewidth =2) +
  #geom_text(aes(y = -0.04, x = mean(CUAJ_sample_position), label = "CUAJ")) +
  scale_fill_manual(values = het.cols) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),axis.line=element_blank()) +
  ylab(label = paste("K =", K))
# theme(axis.text.x = element_text(size=6))

v

# Are there any sample that have odd Q assigns
qtable[qtable$Q>0.2&qtable$Q<0.8,]
hist(qtable$Q)

# Assign PCA colour based off LEA anaylsis
qdf <- cbind(sites[,c("samples","Lat","Long")],Q(titia.snmf, K = K, run = best))
barchart(titia.snmf, K = 3, run = best)
for(i in 1:length(qdf$samples)){
  qdf$assign[i] <- which.max(qdf[i,c(4:(3+K))])
}

qdf$assign
pca.q.df <- merge(pca.data, qdf[,c("samples","assign")], by = "samples")
pca.q.df$site.sub[is.na(pca.q.df$site.sub)] <- "NA01"

pca.q.df$assign <- as.factor(pca.q.df$assign)

# Write out dataframe of sample and qassign
write.table(pca.q.df, file = paste0(dir.path, "LEA_qassign_pca_samples.txt"))
# Without CUAJ
write.table(pca.q.df[sites$site.sub!="CUAJ",], file = paste0(dir.path, "LEA_qassign_pca_samples_noCUAJ.txt"))

# Create simple LEA_pop_file
# Which q assign is Pacfic
LEA_popfile <- pca.q.df[,c("samples", "assign")]
# Convert to charater as factors do not work with assigning text
LEA_popfile$assign <- as.character(LEA_popfile$assign)
LEA_popfile$assign[LEA_popfile$assign==pca.q.df$assign[sites$site.sub=="ZANA"][1]] <- "titia-Pac"
LEA_popfile$assign[LEA_popfile$assign==pca.q.df$assign[sites$site.sub=="ESRB"][1]] <- "titia-SAtl"
LEA_popfile$assign[LEA_popfile$assign==pca.q.df$assign[sites$site.sub=="HCAR"][1]] <- "titia-NAtl"


write.table(LEA_popfile, paste0(dir.path, "popfile_", species,".txt"),
            quote = F,row.names = F, col.names = F)

write.table(LEA_popfile[sites$site.sub!="CUAJ",], paste0(dir.path, "popfile_", species,"_noCUAJ.txt"),
            quote = F,row.names = F, col.names = F)

# Names rather than number in PCA plot
pca.q.df$assign_name <- LEA_popfile$assign
# PCA plots
pca.q.df
pca.label <- c()
pca.label[NAtl.clust] <- "North Atlantic"
pca.label[Pac.clust] <- "Pacific"
pca.label[SAtl.clust] <- "South Atlantic"

y <- ggplot(pca.q.df) +
  geom_point(aes(as.numeric(pca1), as.numeric(pca2), fill = assign), size = 6, shape = 21, color = "black") +
  scale_fill_manual(values = het.cols, name = "Ancestory assignment", labels = pca.label) +
  xlab(pca.labs[1]) +
  ylab(pca.labs[2]) +
  theme(legend.position = c(0.3, 0.85), text = element_text(size = 15), axis.title = element_text(size = 25))

y
a <- ggplot() +
  geom_sf(data = worldmap, fill = "#FFE6D4", col = "black") +
  geom_point(data = pca.q.df, aes(Long, Lat, fill = assign), shape = 21, size = 6) +
  scale_colour_manual(values = het.cols, name = "Ancestory assignment") +
  scale_fill_manual(values = het.cols) +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic")) +
  labs(x = "Longitude", y = "Latitude", col = "Species") +
  theme(panel.background = element_rect(fill = NA),
        panel.ontop = TRUE) +
  theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
        panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
  theme(legend.position = "none") +
  coord_sf(xlim = range(q.coord.pop$long)+c(-5,+5), ylim =range(q.coord.pop$lat)+c(-2,+2),default_crs = st_crs(4326))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))
a

# Combine into single figure plots
library(patchwork)
plot1 <- (p+q+r)/v + plot_layout(heights = c(4, 2)) + theme(plot.background = element_blank()) + 
  plot_annotation(tag_levels = 'a',tag_prefix = "(", tag_suffix = ")") + 
  plot_layout(heights = c(1,0.8))
plot2 <- y #+ plot_annotation(tag_levels = 'a',tag_prefix = "(", tag_suffix = ")")

ggsave(paste0(plot.dir,"LEA_K",K,"snp",snp_sub_text,"_complete_DP10_basins_patchwork_BES_poster_v4.pdf"), plot=plot1, height=7, width=15)
ggsave(paste0(plot.dir,"LEA_K",K,"snp",snp_sub_text,"_complete_DP10_basins_patchwork_BES_poster_v4.png"), plot=plot1, height=7, width=15)
ggsave(paste0(plot.dir,"PCA1-2_snp",snp_sub_text,"_complete_DP10_basins_patchwork_BES_poster_v4.pdf"), plot=plot2, height=6, width=6)
ggsave(paste0(plot.dir,"PCA1-2_snp",snp_sub_text,"_complete_DP10_basins_patchwork_BES_poster_v4.png"), plot=plot2, height=6, width=6)


# # # # # # # # # # ## # # # # # # # # # #
####### Fst and population structure  ####
# # # # # # # # # # ## # # # # # # # # # #
# install.packages("hierfstat")
library(hierfstat)
pca.q.df <- read.table(paste0(dir.path, "LEA_qassign_pca_samples.txt"))
# ASsign populations
# Two different versions being considered
pca.q.df$samples==rownames(my_genind_ti_SNPs@tab)
my_genind_ti_SNPs@pop <- as.factor(pca.q.df$assign)
#my_genind_ti_SNPs@pop <- as.factor(paste(pca.q.df$assign, pca.q.df$Country_Ocean, sep = "_"))

# Calcultes Fst values between populations
pop.diff <- genet.dist(my_genind_ti_SNPs, diploid = T, method = "WC84")
pop.diff

# Same calculation for each sample site
## my_genind_ti_SNPs@pop <- as.factor(pca.q.df$site.sub.county.drain)
# Calcultes Fst values between populations
## pop.diff.site <- genet.dist(my_genind_ti_SNPs, diploid = T, method = "WC84")
## pop.diff.site

pop.diff
my_stats <- basic.stats(my_genind_ti_SNPs)
my_stats$overall

boxplot(my_stats$Ho)
boxplot(my_stats$Hs)
boxplot(my_stats$perloc)

# # # # # # # # # # # # # # # #
####### Kinship analysis ######
# # # # # # # # # # # # # # # #
#install.packages("popkin")
library(popkin)


# Load geno file and change missing values to "NA
geno.kin <- geno
geno.kin[geno.kin=="9"] <- NA

# Choose subpopulations
#Check dimensions are correct
dim(sites)[1]==dim(geno)[1]
subpops <- pca.q.df$Ocean.drainage
# Reorder samples
kin.plot.order <- order(paste(pca.q.df$Ocean.drainage,pca.q.df$Country, pca.q.df$site.sub))
geno.kin <- geno.kin[kin.plot.order,]
subpops <- subpops[kin.plot.order]
subpops.site.sub <- pca.q.df$site.sub[kin.plot.order]

kinship <- popkin(t(geno.kin), subpops = subpops)
plot_popkin(
  kinship,
  labs = subpops,
  # shared bottom and left margin value, to make space for labels
  mar = 1
)


# # # # # # # # # # # # # # 
####### Hybrid index ######
# # # # # # # # # # # # # # 
library(tidyverse)
library(shadowtext)
## install.packages("genetics")
## url <- "http://cran.nexr.com/src/contrib/introgress_1.2.3.tar.gz"
## pkgFile <- "introgress_1.2.3.tar.gz"
## download.file(url = url, destfile = pkgFile)
## install.packages(pkgs=pkgFile, type="source", repos=NULL)
library(introgress)

# Set cut of value for sample removal
snp_sub <- 0.20
snp_sub_text <-"0_20"

# Read in VCF and sample site database
vcf.SNPs <- read.vcfR(paste0(dir.path,SNP.library.name,filter_para,".vcf.gz"))
vcf.SNPs@fix[,1] <- gsub(vcf.SNPs@fix[,1], pattern = "\\.", replacement = "_")
# Get X chromosome
X_chrom <- "HetTit1_0_p_scaff-12-96647824"
vcf.X <- vcf.SNPs[vcf.SNPs@fix[,1]==X_chrom,] 
sites <- read.table(paste0(dir.path, "coorrd_Titia_complete",analysis.name ,".txt"))

vcf.X

## Relevent sample sites
Atl.sites <- c("CT", "TXRS", "CUAJ", "MIXT")
Pac.sites <- c("PUMA", "RLPE", "ZANA", "TULI")


#Subset sites
hybrid.sites <- sites[sites$site.sub%in%c(Atl.sites, Pac.sites),]
hybrid.sites$Ocean.drainage <- as.factor(hybrid.sites$Ocean.drainage)
hybrid.sites$Lat <- as.numeric(hybrid.sites$Lat)
hybrid.sites$Long <- as.numeric(hybrid.sites$Long)

#Check location of subset sites
# Check subset samples are from the correct region
ggplot() +
  geom_sf(data = worldmap, fill = "#FFE6D4", col = "black", linewidth = 0.5) +
  geom_label(data = hybrid.sites, aes(Long, Lat, label = site.sub, fill = Ocean.drainage)) +
  coord_sf(default_crs = st_crs(4326))

# Extract genotypes from vcf
mat <- extract.gt(vcf.SNPs)

#remove non hybrid sites
mat <-mat[,colnames(mat)%in%hybrid.sites$samples]
colnames(mat)==hybrid.sites$samples

conv.mat<-mat
conv.mat[conv.mat == "0/0"]<-0
conv.mat[conv.mat == "0/1"]<-1
conv.mat[conv.mat == "1/1"]<-2
conv.mat<-as.data.frame(conv.mat)
#convert to numeric
for (i in 1:ncol(conv.mat)){
  conv.mat[,i]<-as.numeric(as.character(conv.mat[,i]))
}

# Get sample positions that are not CUAJ
atl.samps <- which(hybrid.sites$site.sub%in%Atl.sites&hybrid.sites$site.sub!="CUAJ")
pac.samps <- which(hybrid.sites$site.sub%in%Pac.sites&hybrid.sites$site.sub!="CUAJ")

#calc AF for the samples you will use to call fixed differences
atl.af<-(rowSums(conv.mat[,atl.samps], na.rm=T)/(rowSums(is.na(conv.mat[,atl.samps]) == FALSE)))/2
pac.af<-(rowSums(conv.mat[,pac.samps], na.rm=T)/(rowSums(is.na(conv.mat[,pac.samps]) == FALSE)))/2
hist(atl.af-pac.af)
#find fixed SNPs
diff<-abs(pac.af - atl.af)
hist(diff)
# Set fixed allele threshold
fixed_thres <- 0.8
#how many SNPs are fixed
table(is.na(diff) == FALSE & diff > fixed_thres)

#subsample original matrix to only fixed diff SNPs
gen.mat<-mat[is.na(diff) == FALSE & diff > fixed_thres,]
dim(gen.mat)

#subsample matrix converted for AF calcs to only fixed SNPS
conv.mat<-conv.mat[is.na(diff) == FALSE & diff > fixed_thres,]
dim(conv.mat)

# Prevalence of 1 being the allele in the Atlantic samples (i.e. not the draft genome genotype)
hist((rowSums(conv.mat[,atl.samps], na.rm=T)/(rowSums(is.na(conv.mat[,atl.samps]) == FALSE)))/2)

#write a logical test to convert alleles so that a single number represents one parental ancestry
for (i in 1:nrow(gen.mat)){
  #if 1 is the atl allele (ie < .2 frequency in the pac samples used for identifying informative SNPs)
  if((sum(conv.mat[i,atl.samps], na.rm=T)/(sum(is.na(conv.mat[i,atl.samps]) == FALSE)))/2 < (1-fixed_thres)){
    #swap all '0/0' cells with '2/2'
    gen.mat[i,][gen.mat[i,] == "0/0"]<-"2/2"
    #swap all '1/1' cells with '0/0'
    gen.mat[i,][gen.mat[i,] == "1/1"]<-"0/0"
    #finally convert all '2/2' cells (originally 0/0) into '1/1'
    gen.mat[i,][gen.mat[i,] == "2/2"]<-"1/1"
    #no need to touch hets
  }
}

#convert R class NAs to the string "NA/NA"
gen.mat[is.na(gen.mat) == TRUE]<-"NA/NA"

#make locus info df
locus.info<-data.frame(locus=rownames(gen.mat),
                       type=rep("C", times=nrow(gen.mat)),
                       lg=vcf.SNPs@fix[,1][is.na(diff) == FALSE & diff > fixed_thres],
                       marker.pos=vcf.SNPs@fix[,2][is.na(diff) == FALSE & diff > fixed_thres])

locus.info$marker.pos<-as.numeric(as.character(locus.info$marker.pos))

# Relabel lg with chromosome number
locus.info$lg <- sapply(strsplit(locus.info$lg,"-") , "[", 2 )
#Label x chromosome
locus.info$lg[locus.info$lg==12] <- "X"
locus.info$lg <- fct_relevel(locus.info$lg, c(1:11,"X"))

# Creating column that is equvlent to distance over the whole genonome
nCHR <- length(unique(locus.info$lg))
locus.info$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(locus.info$lg)){
  nbp[i] <- max(locus.info[locus.info$lg == i,]$marker.pos)
  locus.info[locus.info$lg == i,"BPcum"] <- locus.info[locus.info$lg == i,"marker.pos"] + s
  s <- s + nbp[i]
}

head(locus.info)
table(locus.info$lg!="X")

#we now have a gt matrix in proper format for introgress
#convert genotype data into a matrix of allele counts
count.matrix<-prepare.data(admix.gen=gen.mat, loci.data=locus.info,
                           parental1="1",parental2="0", pop.id=F,
                           ind.id=F, fixed=T)

hi.index.sim<-est.h(introgress.data=count.matrix,loci.data=locus.info,
                    fixed=T, p1.allele="1", p2.allele="0")

#png(paste0(plot.dir, "hybrid_geno_count_plot.png"), width = 2000, height = 1000, units = "px")
#mk.image(introgress.data=count.matrix, loci.data=locus.info,
#        hi.index=hi.index.sim, ylab.image="Individuals",
#        xlab.h="population of Pacific ancestry", pdf=F,
#        col.image=c(rgb(1,0,0,alpha=.5),rgb(0,0,0,alpha=.8),rgb(0,0,1,alpha=.5)))
#dev.off()

# Combine data into tidy format for ploting in ggplot
i <- 1
locus.geno.type <- list()
for(i in 1:dim(gen.mat)[2]){
  locus.geno.type[[i]] <- cbind(sample = rep(colnames(gen.mat)[i], dim(gen.mat)[1]),
                                hybrid.index = hi.index.sim[i,2],
                                genotype = gen.mat[,i],
                                locus.name = row.names(gen.mat),
                                locus.info)
}

locus.geno.type <- do.call("rbind", locus.geno.type)
head(locus.geno.type)

locus.geno.type$sample <- as.factor(locus.geno.type$sample)
#Remove single site from chromo14 

# Reorder samples by hybrid index
hybrid.sites$new.lat <- hybrid.sites$Lat
#Annoying that TXRS is plotted below CUAJ when using lat
hybrid.sites$new.lat[hybrid.sites$site.sub=="TXRS"] <- hybrid.sites$new.lat[hybrid.sites$site.sub=="TXRS"]+1
CUAJ.lat <- hybrid.sites$new.lat[hybrid.sites$samples=="CUAJa02"]
locus.geno.type$sample <- fct_relevel(locus.geno.type$sample, colnames(gen.mat)[order(as.numeric(hi.index.sim[,2]), decreasing = T)])
locus.geno.type$sample <- fct_relevel(locus.geno.type$sample, colnames(gen.mat)[order((hybrid.sites$new.lat-CUAJ.lat)+(-hi.index.sim[,2])*0.001, decreasing = F)])
hybrid.sites$samples <- fct_relevel(hybrid.sites$sample, hybrid.sites$samples[order((hybrid.sites$new.lat-CUAJ.lat)+(-hi.index.sim[,2])*0.001, decreasing = F)])


mean(hybrid.sites$Long[hybrid.sites$site.sub=="CUAJ"])-hybrid.sites$Long

# Get locations of vertical lines where chromosomes are plotted
vline_chrom <- group_by(locus.geno.type[locus.geno.type$sample==unique(locus.geno.type$sample)[1],], by = lg) %>%
  summarise(max.length = which.max(BPcum),
            med.bp = length(BPcum))
#Reorder line database
vline_chrom <- vline_chrom[order(as.numeric(vline_chrom$by)),]
vline_chrom$cumsum.chrom <- cumsum(vline_chrom$max.length)
vline_chrom$text_chrom_pos <- vline_chrom$cumsum.chrom-(vline_chrom$max.length/2)
#Remove 14 line
vline_chrom <- vline_chrom[as.numeric(vline_chrom$by)<=12|vline_chrom$by=="X",]
#horizontal
p.h <- ggplot(locus.geno.type[]) +
  geom_raster(aes(x = as.factor(BPcum), y = sample, fill = genotype)) +
  geom_vline(xintercept = vline_chrom$cumsum.chrom+1) +
  geom_text(data = vline_chrom, aes(x  = text_chrom_pos, y = -1, label = by), size = 4) +
  scale_fill_manual(values = c("#AF0F09","darkorange","#E5D9BA","white")) +
  coord_cartesian(ylim = c(1, length(unique(locus.geno.type$sample))), # This focuses the x-axis on the range of interest
                  clip = 'off') +
  theme(plot.margin = unit(c(3,1,3,1), "lines"),legend.position = c(0.5,1.05),legend.direction = "horizontal",
        text = element_text(size = 15),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6)) +
  theme(legend.key=element_rect(colour="black")) +
  labs(fill = "Genotype")
p.h
#Vertical
p.v <- ggplot(locus.geno.type) +
  geom_raster(aes(y = as.factor(BPcum), x = sample, fill = genotype)) +
  geom_hline(yintercept = vline_chrom$cumsum.chrom+1) +
  geom_text(data = vline_chrom, aes(y  = text_chrom_pos, x = -1, label = by), ) +
  scale_fill_manual(values = c("#AF0F09","darkorange","#E5D9BA","white")) +
  coord_cartesian(xlim = c(1, length(unique(locus.geno.type$sample))), # This focuses the x-axis on the range of interest
                  clip = 'off') +
  theme(plot.margin = unit(c(3,1,3,3), "lines"),legend.position = c(0.5,1.05),legend.direction = "horizontal",
        axis.text.y = element_blank(), axis.title.y = element_blank()) +
  theme(legend.key=element_rect(colour="black"),
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #scale_y_discrete(limits=rev) +
  labs(fill = "Genotype")

# Just X chromosome
p.x <- ggplot(locus.geno.type[locus.geno.type$lg=="X",]) +
  geom_raster(aes(x = as.factor(BPcum), y = sample, fill = genotype)) +
  scale_fill_manual(values = c("#AF0F09","darkorange","#E5D9BA","white")) +
  coord_cartesian(ylim = c(1, length(unique(locus.geno.type$sample))), # This focuses the x-axis on the range of interest
                  clip = 'off') +
  theme(plot.margin = unit(c(3,1,3,1), "lines"),legend.position = c(0.5,1.05),legend.direction = "horizontal",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.key=element_rect(colour="black")) +
  labs(fill = "Genotype",x = "X") 

ggsave(file = paste0(plot.dir,"introgress_X_chrom.png"), p.x, width = 6.5, height = 6)

# Are any sites from CUAJa02 X chromosome heterozgous
locus.geno.type$locus.name[locus.geno.type$sample=="CUAJa02"&locus.geno.type$lg=="X"&locus.geno.type$genotype=="0/1"]


# Calculate heterozgousity of sample
# Recalculates heterozygousity but without the X chromosome which for males will only ever be homozgous
count.matrix.noX <- prepare.data(admix.gen=gen.mat[locus.info$lg!="X",], loci.data=locus.info[locus.info$lg!="X",],
                                 parental1="1",parental2="0", pop.id=F,
                                 ind.id=F, fixed=T)
# Calculate heterozgousity of samples
het<-calc.intersp.het(introgress.data=count.matrix.noX)

hi.index.sim.noX<-est.h(introgress.data=count.matrix.noX,loci.data=locus.info[locus.info$lg!="X",],
                    fixed=T, p1.allele="1", p2.allele="0")

# Inbuilt hybrid plot that used base R plotting
introgress::triangle.plot(hi.index=hi.index.sim.noX, int.het=het, pdf = F)

# attach het and hybrid to sample dataframe
hybrid.sites$het.fst <- het
hybrid.sites$hybrid.index <- hi.index.sim.noX[,2]

# Create triangle for hydrid index plot
tri <- cbind.data.frame(x = c(seq(0,0.5,length.out=10),seq(0.5,1,length.out=10)), y = c(seq(0,1,length.out=10), seq(1,0,length.out=10)))

q <- ggplot(hybrid.sites) +
  geom_line(data = tri, aes(x = x, y = y)) +
  geom_point(aes(x = hybrid.index, y = het.fst, fill = Ocean.drainage), shape = 21, size = 4) +
  scale_fill_manual(values = het.cols[c(1,3)]) +
  xlab("Hybrid Index") +
  ylab("Heterozygousity of >0.8 Fst SNPs") +
  labs(fill = "Slope") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 15))

ggsave(file = paste0(plot.dir,"introgress_tri_plot_v2.png"), q, width = 6.5, height = 6)

# Are there any samples that have a hybrid index suggested of introgression
hybrid.sites[hybrid.sites$hybrid.index>0.2&hybrid.sites$hybrid.index<0.8,]

#load in depth river shapefile
load("4_Manuscript/data/World_data/hydroRIVERS/HydroRIVERS_v10_na_sub_hybrid.rda")

hybrid.sites$nudge.dir <- NA

hybrid.sites$nudge.dir[hybrid.sites$Ocean.drainage=="Gulf"] <-0.3
hybrid.sites$nudge.dir[hybrid.sites$Ocean.drainage=="Pacific"] <- (-0.1)

# Remove samples outside of bounding box
mex.e$Long
r <- ggplot() +
  geom_raster(data = hill.df.Mex.z, aes(lon, lat, fill = hill), alpha = 1,show.legend = F) +
  scale_fill_gradientn(colours=c("#5C2700","#FFE6D4")) +
  new_scale_fill() +
  #
  geom_sf(data = hydrobasins_geo, col = "black", fill = NA, linewidth = 0.75) +
  geom_sf(data = hydrorivers.sub, col = "#5C2700", lineend = "round", size = 0.5) +
  geom_point(data = hybrid.sites[!duplicated(hybrid.sites$Lat)&hybrid.sites$Lat<17.4&hybrid.sites$Long<(-94.0),], 
             aes(Long, Lat), col = "black", size = 1.5) +
  #geom_label_repel(data = hybrid.sites[!duplicated(hybrid.sites$Lat),], 
  #                 aes(Long, Lat, label = Site.ID, fill = Ocean.drainage),
  #                 force_pull = 1,alpha = 1,force = 75, size = 3,
  #                 min.segment.length = 0, show.legend = F, seed = 4, nudge_x = -0.25, col = "black") +
  geom_label_repel(data = hybrid.sites[!duplicated(hybrid.sites$Lat)&hybrid.sites$Lat<17.4&hybrid.sites$Long<(-94.0),], 
                   aes(Long, Lat, label = Site.ID, fill = Ocean.drainage),
                   force_pull = 1,alpha = 1,force = 50, size = 4,
                   min.segment.length = 0, show.legend = F, seed = 12345, nudge_y = hybrid.sites$nudge.dir[!duplicated(hybrid.sites$Lat)&hybrid.sites$Lat<17.4&hybrid.sites$Long<(-94.0)],
                   nudge_x = -0.2, col =  "black") +
  scale_fill_manual(values = het.cols[c(1,3)]) +
  coord_sf(xlim = mex.e.zoom$Long, ylim = mex.e.zoom$Lat) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(text = element_text(size = 15), axis.text = element_text(size = 8))

# Inset map
in.map <- ggplot() +
  geom_sf(data = worldmap_proj, fill = "#FFE6D4", col = "black", linewidth = 0.5) +
  geom_point(data = hybrid.sites, aes(Long, Lat, col = Ocean.drainage), show.legend=F) +
  coord_sf(xlim = mex.e$Long, ylim=mex.e$Lat, default_crs = st_crs(4326)) +
  geom_polygon(data = mex.e.zoom.e, aes(x = Long, y = Lat), size = 1, fill = NA, colour = "black") +
  theme_void() +
  theme(plot.background = element_rect(fill = 'white', colour = 'black'))
  
r.map <- r + patchwork::inset_element(in.map, left = 0.5, right = 0.95, top = 1, bottom = 0.5)

# Stats
#Number of samples in analysis
length(hybrid.sites$samples)
# Number of higth fst SNPs
dim(gen.mat)
#CUAJa02 heterozgous
#with X
table(gen.mat[,"CUAJa02"])/(dim(gen.mat)[1]-table(gen.mat[,"CUAJa02"])["NA/NA"])*100
#Without X
table(hybrid.sites$Site.ID)
table(gen.mat[locus.info$lg!="X","CUAJa02"])/(dim(gen.mat[locus.info$lg!="X",])[1]-table(gen.mat[locus.info$lg!="X","CUAJa02"])["NA/NA"])*100

hybrid.sites$het.fst[hybrid.sites$samples=="CUAJa02"]
hybrid.sites$hybrid.index[hybrid.sites$samples=="CUAJa02"]

plot.hybrid <- p.h / (r+q) + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") 
ggsave(file = paste0(plot.dir,"introgress_tri_plot_with_genotype_horizontal_v2.png"), plot.hybrid, width = 12, height = 11)
ggsave(file = paste0(plot.dir,"introgress_tri_plot_with_genotype_horizontal_v2.pdf"), plot.hybrid, width = 12, height = 11)

plot.hybrid <- p.v + (q/r) + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
ggsave(file = paste0(plot.dir,"introgress_tri_plot_with_genotype_vertical_v2.png"), plot.hybrid, width = 10, height = 15)

## What is the coverage of the X chromosome of the hybrid individuals
X_chrom <- "HetTit1_0_p_scaff-12-96647824"

# Get depth of reads

# What is the coverage of samples across the genome

# Extrat raw read coverage
depth_df <- extract.gt(vcf.SNPs, element = "DP", as.numeric = T)
# Melth into tidy format
depth_df <- reshape2::melt(depth_df, id = c("SNP", "SAMPLE"), value.name = "depth_var")

# Extract chromosome position and chromosome bnumber for each SNP
depth_df$chrom <- str_split_i(depth_df$Var1, pattern = "-", i = 2)
depth_df$chrom_length <- str_split_i(str_split_i(depth_df$Var1, pattern = "-", i = 3), pattern = "_", i = 1)
depth_df$position <- str_split_i(str_split_i(depth_df$Var1, pattern = "-", i = 3), pattern = "_", i = 2)

# Extract depth for the hybrid individual
depth_CUAJ <- depth_df[depth_df$Var2=="CUAJa02",]
depth_CUAJ$chrom <- as.factor(depth_CUAJ$chrom)
depth_CUAJ$chrom <- factor(depth_CUAJ$chrom, levels = sort(as.numeric(unique(depth_CUAJ$chrom))))

# Extract SNP genotype for hybrid
vcf_CUAJ <- vcf.SNPs[sample = "CUAJa02"]
# Get heterozgousity of each SNP
depth_CUAJ$is_het <- is.het(extract.gt(vcf_CUAJ, element = "GT"), na_is_false = F)

# Plot Read depth over heterozgousity for each SNP for each Chromosome
dp.plot <- ggplot(depth_CUAJ[depth_CUAJ$chrom%in%1:12&!is.na(depth_CUAJ$is_het),]) +
  geom_boxplot(aes(x = chrom, y = depth_var), fill = NA, outlier.alpha = 0, size = 2) +
  geom_jitter(aes(x = chrom, y = depth_var, fill = is_het), shape = 21, height = 0, width = 0.2, size = 3) +
  scale_fill_manual(values = c("#E5D9BA","darkorange")) +
  labs(x = "Chromosome", y = "SNP raw read depth", fill = "Is heterozgous") +
  ylim(0,100) +
  theme_bw() +
  theme(legend.position = "bottom")
dp.plot

ggsave(file = paste0(plot.dir,"Raw_read_depth_and heterozgousity for CUAJa02.png"), dp.plot, width = 15, height = 10)
ggsave(file = paste0(plot.dir,"Raw_read_depth_and heterozgousity for CUAJa02.pdf"), dp.plot, width = 15, height = 10)





