# install.packages(c("ape","tidyverse","ggplot2","ggtree","patchwork","treeio","phytools","phangorn"))
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
# library(ape)
library(tidyverse)
library(ggtree)
library(patchwork)
library(treeio)
library(phytools)
library(phangorn)
library(rgdal)
library(mapproj)

# Output file location
# Output file location
# Get system arguments
args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[1]
SNP.library.location <- args[2]


# Create output location
libary.dir <- paste0(SNP.library.location,SNP.library.name) 
plot.dir <- paste0("/home/tmjj24/plots/Chapter_3/", SNP.library.name, "/")
dir.create(plot.dir)
dir.create(paste0(plot.dir, "SVDQuartets/"))

tre <- read.nexus(paste0(libary.dir,"/SVDQuartets/", SNP.library.name,"_T20.tre"))

print("Rooting tree at midpoint")
tre.midroot <- midpoint(tre, node.labels = "support")

#Extracting sample information and linking it to tree labels
sample_map<- read.csv("/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/3_Results/All samples held in Durham_v17.csv", check.names=F)

samples <- tre$tip.label
sites <- data.frame(samples)
sites$samples
#Convert lat and long to number
sample_map$Long <- as.numeric(sample_map$Long)
sample_map$Lat <- as.numeric(sample_map$Lat)

sites <- merge(sites, sample_map, by.x = "samples", by.y = "Unique.ID")
any(is.na(sites))

sites$samples[sites$species=="occisa"]
names(sites)
sites$species_drainage <- as.factor(paste(sites$species,sites$Ocean.drainage, sep = "_"))
tmp.levels <- c(levels(sites$species_drainage)[levels(sites$species_drainage)!="occisa_Pacific"], "occisa_Pacific")

sites$species_drainage <- factor(sites$species_drainage, levels = tmp.levels)
sites$species_drainage_country <- paste(sites$species,sites$Ocean.drainage, sites$Country,sep = "_")
sites$drainage_country <- paste(sites$Ocean.drainage, sites$Country,sep = "_")
sites$species.country.drainage <- paste(substr(sites$species, start = 1, stop = 2), substr(sites$Country, start = 1, stop = 3), substr(sites$Ocean.drainage, start = 1, stop = 3), sep = ".")

#Reorders factors to position occisa at the end of the level. This allows occisa to be
# removed from the map plot without reordering the colours used
sites$species.country.drainage[sites$species=="occisa"] <- "occisa"
sites$species.country.drainage <- as.factor(sites$species.country.drainage)
sites$species.country.drainage <- factor(sites$species.country.drainage, levels = c(levels(sites$species.country.drainage)[levels(sites$species.country.drainage)!="occisa"],
                                                                                    levels(sites$species.country.drainage)[levels(sites$species.country.drainage)=="occisa"]))

sites$species.fig <- NA
#sites$species.fig[sites$species=="occisa"] <- "H. occisa"
sites$species.fig[sites$species=="americana"] <- "H. americana/calverti"
sites$species.fig[sites$species=="titia"] <- "H. titia"

unique(sites$species_drainage_country)

#pdf(file = "plots/H_ti_Oc_DP10_ran1_1000_10MSVQ_region_colour.pdf", width = 5.5, height = 10)

cbPalette <- c("#009E73", "#0072B2", "#E69F00",  "#D55E00","#56B4E9", 
               "#CC79A7","#999999", "#F0E442", 
               "deeppink","darkred", "purple","black")

tre.midroot

tree.plot <- ggtree(tre.midroot,size = 0.8) %<+% sites
tree.plot$data$support <- as.numeric(tree.plot$data$label)
tree.plot$data$support.50 <- tree.plot$data$support
tree.plot$data$support.100 <- NA
tree.plot$data$support.100[tree.plot$data$support==100] <- "*"
tree.plot <- tree.plot +
  #geom_tree(aes(colour = support)) +
  geom_nodepoint(aes(colour = support.50, alpha = support.50), shape = 19, size = 4) +
  geom_tippoint(aes(fill = species.country.drainage), shape = 21, size = 3) +
  #geom_tiplab(data = tree.plot$data, aes(colour = species_drainage),offset = 0, size = 2,show.legend = FALSE) 
  #geom_tippoint(aes(colour = species_drainage)) +
  geom_text(aes(label = support.100), nudge_x = -0.002, nudge_y = 0, color = "red", size = 10) +
  scale_fill_manual(values = cbPalette) +
  scale_colour_gradient(low = "white", high = "black") +
  coord_cartesian(clip = 'off') +
  theme_tree(plot.margin = margin(2,2,2,2, "cm")) +
  guides(fill=guide_legend(title="Species and Drainage"), 
         colour=guide_legend(title="Bootstrap (Out of 100 reps)"),
         alpha="none") +
  ggtitle("SVDQuartet") +
  geom_treescale()


world.e <- data.frame(Long = c(-120,-60), Lat = c(8,38))
worldmap <- rgdal::readOGR("/home/tmjj24/World_map_data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
#Cropping the extent of worldmap to reduce plotting time
e <- raster::extent(-130, -50, 0, 50)
worldmap <- raster::crop(worldmap, e)

b <- ggplot() +
  geom_polygon(data = worldmap, aes(long, lat, group = group), fill = "#FFE6D4", col = "black") +
  geom_point(data = tree.plot$data[!is.na(tree.plot$data$species.fig),], aes(Long, Lat, fill = species.country.drainage), size = 3, shape = 21) +
  scale_fill_manual(values = cbPalette) +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic")) +
  labs(x = "Longitude", y = "Latitude", col = "Species") +
  theme(panel.background = element_rect(fill = NA),
        panel.ontop = TRUE) +
  theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
        panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
  theme(legend.position = "none") +
  ggtitle("Hetaerina sample distribution") +
  coord_map("bonne", lat0 = 50, xlim = world.e$Long, ylim = world.e$Lat) +
  facet_wrap(~species.fig, drop = T, nrow = 2) +
  theme(strip.text.x = element_text(size = 20),title = element_text(size = 25))

print("Saving image to:")
print(paste0(plot.dir, "SVDQuartets/","SVDQuartet_",SNP.library.name,".png"))
ggsave(plot = b/tree.plot,  filename = paste0(plot.dir, "SVDQuartets/","SVDQuartet_",SNP.library.name,".png"), height = 20, width = 15)
ggsave(plot = b/tree.plot,  filename = paste0(plot.dir, "SVDQuartets/", "SVDQuartet_",SNP.library.name,".pdf"), height = 20, width = 15)

       