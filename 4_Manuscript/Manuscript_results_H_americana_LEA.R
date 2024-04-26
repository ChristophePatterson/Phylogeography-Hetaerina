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


#######################################################################################################
#######################################################################################################

###### H titia phylogeographu plots ######

#######################################################################################################
#######################################################################################################


# Output file location
SNP.library.name <- "Hetaerina_americana_ddRAD_titia_dg"
species <- "americana"
dir.path <- paste0("4_Manuscript/data/SNP_libraries/",SNP.library.name,"/")
filter_para <- ".all.snps.NOGTDP10.MEANGTDP10_200.Q60.SAMP0.8.MAF2.rand1000.biSNP0_20"
snp_sub_text <- "0_20"
# Plot output file location
plot.dir <- paste0("4_Manuscript/plots/",SNP.library.name,"/")
dir.create(plot.dir)

## Running analysis name
analysis.name <- paste0("H_",species, "_complete_snp_noX")

vcf <- read.vcfR(paste0(dir.path,SNP.library.name,filter_para,".vcf.gz"))

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
sample_map<- read.csv("4_Manuscript/data/Hetaerina_sample_data.csv")
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
write.table(sites, file = paste0(dir.path, "coorrd_", species,"_complete_",analysis.name ,".txt"))
sites <- read.table(paste0(dir.path, "coorrd_", species,"_complete_",analysis.name ,".txt"))

#Calculates structure for samples from K=1 to k=10
max.K <- 10
# MAY NEED TO PAUSE ONEDRIVE
##obj.at <- snmf(paste0(dir.path, analysis.name,snp_sub_text,".geno"), K = 1:max.K, ploidy = 2, entropy = T,
##             CPU = 2, project = "new", repetitions = 20, alpha = 100)
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

ggsave(paste0(plot.dir,"H_", species,"_LEA_complete_snp_K1-8",snp_sub_text,"cross_entropy.pdf"), plot = ce.plot)
ggsave(paste0(plot.dir,"H_", species,"_LEA_complete_snp_K1-8",snp_sub_text,"cross_entropy.jpg"), plot = ce.plot)
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

# Figure out what number each region has been assigned to (As it varies between sNMF runs)
Namer.clust <- which.max(q.coord.pop[q.coord.pop$site=="TN_United States_Gulf",LETTERS[1:K]])
SAmer.clust <- which.max(q.coord.pop[q.coord.pop$site=="STDM_Mexico_Pacific",LETTERS[1:K]])
cal.clust <- which.max(q.coord.pop[q.coord.pop$site=="RG_Mexico_",LETTERS[1:K]])

#Conduct PCA
geno2lfmm(paste0(dir.path, analysis.name,snp_sub_text,".geno"), 
          paste0(dir.path, analysis.name,snp_sub_text,".lfmm"), force = TRUE)
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


####################
###### PLOTS ######
####################

# Outlines for each plot
world.e <- data.frame(Long = c(-120,-70), Lat = c(12,40))
mex.e <- data.frame( Long = c(-106,-92), Lat = c(22,14.5))
mex.e.zoom <- data.frame( Long = c(-96,-94), Lat = c(18,16))
CR.e <- data.frame(Long = c(-85,-82), Lat = c(8,11))


# World projection
crs_text <- "+proj=laea +x_0=0 +y_0=0 +lon_0=-74 +lat_0=40"
crs_text <- paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=",360+mean(sites$Long), " +lat_0=",mean(sites$Lat))

# Download terrrian and river data
source("4_Manuscript/hydrobasins_extract_code_sf.R")
library(scatterpie)
library(ggnewscale)

cbPalette <- c("#F0E442","#D55E00" , "#0072B2", "#999999","#E69F00" , "#56B4E9", "#009E73", "#CC79A7", "black")
# Colour scheme

het.cols <- c("#726230","#AA9599","#548A39")
het.cols[cal.clust] <- "#AA9599"
het.cols[Namer.clust] <- "#726230"
het.cols[SAmer.clust] <- "#548A39"

#het.cols <- c("#3E3C3A","#AF0F09","#E5D9BA", "#694438")
#amer.cols <- c("#948A39","#326230","#AA9599","#120E02")

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


# Shifts all the Mexico points so that they do not overlap
# Min dist in decimal degrees that populations whould be assigned 
min.dist <- 0.8
# Subset population to only those in costa rica
q.coord.pop.Mx <- q.coord.pop[grepl(q.coord.pop$site, pattern = paste(c("Mexico", "Belize"), collapse = "|")),]
# Assigned new location col
q.coord.pop.Mx$long.new <- q.coord.pop.Mx$long
q.coord.pop.Mx$lat.new <- q.coord.pop.Mx$lat
set.seed(3)
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
  geom_scatterpie(data = q.coord.pop, aes(x = long, y = lat, group = site, r = 1), cols = LETTERS[1:K]) +
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
  scale_fill_gradientn(colours=c("#5C2700","#FFE6D4")) +
  new_scale_fill() +
  geom_sf(data = hydrobasins_geo, col = "black", fill = NA, linewidth = 0.5) +
  geom_sf(data = hydrorivers_geo, col = "#5C2700", lineend = "round") +
  geom_segment(aes(x = q.coord.pop.Mx$long , y = q.coord.pop.Mx$lat , xend = q.coord.pop.Mx$long.new, yend = q.coord.pop.Mx$lat.new),
               linewidth = 1.2, lineend =  "round") +
  geom_scatterpie(data = q.coord.pop.Mx, aes(x = long.new, y = lat.new, r = 0.4 , group = site), cols = LETTERS[1:K]) +
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
v <- ggplot(qtable)+
  geom_bar(stat="identity", aes(sample, Q, fill = Qid,), position = "stack", width = 1, col = "black") +
  scale_fill_manual(values = het.cols) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
  ylab(label = paste("K =", K)) 
# theme(axis.text.x = element_text(size=6))

v
# hybrid triangle plot

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

write.table(LEA_popfile, paste0(dir.path, "popfile_", species,".txt"),
            quote = F,row.names = F, col.names = F)

write.table(LEA_popfile[sites$site.sub!="CUAJ",], paste0(dir.path, "popfile_", species,"_noCUAJ.txt"),
            quote = F,row.names = F, col.names = F)


# PCA plots
y <- ggplot(pca.q.df) +
  geom_point(aes(as.numeric(pca1), as.numeric(pca2), fill = assign), size = 6, shape = 21, color = "black") +
  scale_fill_manual(values = het.cols, name = "Ancestory assignment", 
                    c(expression(paste(italic("H. americana"), " north")), 
                      expression(paste(italic("H. americana"), " south")), 
                      expression(paste(italic("H. calverti"))))) +
  xlab(pca.labs[1]) +
  ylab(pca.labs[2]) +
  theme(legend.position = c(0.8, 0.2))
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
plot1 <- (p+q+y)/v + plot_layout(heights = c(4, 2)) + theme(plot.background = element_blank()) + 
  plot_annotation(tag_levels = 'a',tag_prefix = "(", tag_suffix = ")") + 
  plot_layout(heights = c(1,0.8))
plot2 <- y + plot_annotation(tag_levels = 'a',tag_prefix = "(", tag_suffix = ")")

ggsave(paste0(plot.dir,"LEA_K",K,"snp",snp_sub_text,"_complete_DP10_basins_patchwork.png"), plot=plot1, height=7, width=15)
ggsave(paste0(plot.dir,"LEA_K",K,"snp",snp_sub_text,"_complete_DP10_basins_patchwork.pdf"), plot=plot1, height=7, width=15)


############################################
####### Fst and population structure  ######
############################################
# install.packages("hierfstat")
library(hierfstat)
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
