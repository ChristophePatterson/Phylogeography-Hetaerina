# Manuscript_results_SVDQuartets and RAxML

library(ape)
library(tidyverse)
library(ggtree)
library(patchwork)
library(sf)
library(treeio)
library(phytools)
library(phangorn)
# Get colour pallettes
het.cols <- data.frame(assign = c(paste0("americana_",1:3),paste0("titia_",1:3)),
                       cols = c("#848A39","#AA9599","#920E02","#AF0F09","#E5D9BA","#3E3C3A"))
plot.dir <- "4_Manuscript/plots/Hetaerina_all_ddRAD_titia_dg/"
# Read in titia RAxML tree
SNP.library.name <- "Hetaerina_titia_ddRAD_titia_dg"
species <- "titia"
# Create output location
libary.dir <- paste0("4_Manuscript/data/SNP_libraries/",SNP.library.name) 
phy.amend <- ".raxmlASC_GTRGAM_HPCaT15_noCUAJ"

tre <- read.tree(paste0(libary.dir,"/RAxML/RAxML_bipartitions.",SNP.library.name, phy.amend))

# Rooting tree to the mid point
tre.midroot <- midpoint(tre, node.labels = "support")


#Extracting sample information and linking it to tree labels
## Read in LEA assignment
qtitia <- read.table("4_Manuscript/data/SNP_libraries/Hetaerina_titia_ddRAD_titia_dg/LEA_qassign_pca_samples.txt")
qtitia$assign_spp <- paste0(qtitia$species,"_", qtitia$assign)
qameri <- read.table("4_Manuscript/data/SNP_libraries/Hetaerina_americana_ddRAD_titia_dg/LEA_qassign_pca_samples.txt")
qameri$assign_spp <- paste0(qameri$species,"_", qameri$assign)
qtable <- rbind(qtitia, qameri)

# Merge with LEA assignment data
samples <- tre.midroot$tip.label
sites <- data.frame(samples)
#Convert lat and long to number

sites <- merge(sites, qtable, by.x = "samples", by.y = "samples")
table(sites$assign_spp)
any(is.na(sites$assign_spp))

# Merge tree data with sample data
tree.plot <- ggtree(tre.midroot,size = 1.5) %<+% sites
# Create support value
tree.plot$data$support <- as.numeric(tree.plot$data$label)
# Create binary 100 percent bootstrap tree
tree.plot$data$support.100 <- tree.plot$data$support
tree.plot$data$support.100[tree.plot$data$support<95] <- NA
tree.plot$data$support.100[tree.plot$data$support>=95] <- "*"
# tree.plot$data$support.100[tree.plot$data$support>=100] <- "*"


tree.plot$data[tree.plot$data$isTip&is.na(tree.plot$data$assign_spp),]

tree.plot$data$assign_spp <- as.factor(tree.plot$data$assign_spp)
tree.plot$data$assign_spp.col <- het.cols$cols[match(tree.plot$data$assign_spp, het.cols$assign)]

tree.plot.3 <- tree.plot +
  geom_segment(aes(x = x, xend = max(x)+0.01, y=y, yend=y), col = tree.plot$data$assign_spp.col) +
  geom_tippoint(aes(x = max(x)+0.01, col = assign_spp),col = tree.plot$data$assign_spp.col[tree.plot$data$isTip],
                shape =15, size = 1) +
  #geom_nodepoint(aes(colour = support, alpha = support), shape = 19, size = 4) +
  #geom_tiplab(aes(label=label), align = T, linetype = "dotted") +
  #geom_tiplab(data = tree.plot$data, aes(colour = species_drainage),offset = 0, size = 2,show.legend = FALSE) 
  #geom_tippoint(aes(colour = species_drainage)) +
  geom_text(aes(label = support.100), nudge_x = -0.002, nudge_y = 0, color = "deepskyblue", size = 10) +
  # scale_fill_manual(values = cbPalette) +
  scale_colour_gradient(low = "white", high = "black") +
  coord_cartesian(clip = 'off') +
  theme_tree(plot.margin = margin(2,2,2,2, "cm")) +
  guides(fill=guide_legend(title="Species and Drainage"), 
         colour=guide_legend(title="Bootstrap (Out of 100 reps)"),
         alpha="none") +
  geom_treescale() +
  theme(title = element_text(size = 25), legend.text = element_text(size=20), legend.title = element_text(size=25),
        plot.margin = unit(c(0,0,0,0), "cm"))

# Read in americana tree
SNP.library.name <- "Hetaerina_americana_ddRAD_titia_dg"
species <- "americana"
# Create output location
libary.dir <- paste0("4_Manuscript/data/SNP_libraries/",SNP.library.name) 
phy.amend <- ".raxmlASC_GTRGAM_HPCaT15_noCUAJ"

tre <- read.tree(paste0(libary.dir,"/RAxML/RAxML_bipartitions.",SNP.library.name, phy.amend))

# Rooting tree to the mid point
tre.midroot <- midpoint(tre, node.labels = "support")


#Extracting sample information and linking it to tree labels
## Read in LEA assignment
qtitia <- read.table("4_Manuscript/data/SNP_libraries/Hetaerina_titia_ddRAD_titia_dg/LEA_qassign_pca_samples.txt")
qtitia$assign_spp <- paste0(qtitia$species,"_", qtitia$assign)
qameri <- read.table("4_Manuscript/data/SNP_libraries/Hetaerina_americana_ddRAD_titia_dg/LEA_qassign_pca_samples.txt")
qameri$assign_spp <- paste0(qameri$species,"_", qameri$assign)
qtable <- rbind(qtitia, qameri)

# Merge with LEA assignment data
samples <- tre.midroot$tip.label
sites <- data.frame(samples)
#Convert lat and long to number

sites <- merge(sites, qtable, by.x = "samples", by.y = "samples")
table(sites$assign_spp)
any(is.na(sites$assign_spp))

# Merge tree data with sample data
tree.plot <- ggtree(tre.midroot,size = 1.5) %<+% sites
# Create support value
tree.plot$data$support <- as.numeric(tree.plot$data$label)
# Create binary 100 percent bootstrap tree
tree.plot$data$support.100 <- tree.plot$data$support
tree.plot$data$support.100[tree.plot$data$support<95] <- NA
tree.plot$data$support.100[tree.plot$data$support>=95] <- "*"
# tree.plot$data$support.100[tree.plot$data$support>=100] <- "*"


tree.plot$data[tree.plot$data$isTip&is.na(tree.plot$data$assign_spp),]

tree.plot$data$assign_spp <- as.factor(tree.plot$data$assign_spp)
tree.plot$data$assign_spp.col <- het.cols$cols[match(tree.plot$data$assign_spp, het.cols$assign)]

tree.plot.4 <- tree.plot +
  geom_segment(aes(x = x, xend = max(x)+0.01, y=y, yend=y), col = tree.plot$data$assign_spp.col) +
  geom_tippoint(aes(x = max(x)+0.01, col = assign_spp),col = tree.plot$data$assign_spp.col[tree.plot$data$isTip],
                shape =15, size = 1) +
  #geom_nodepoint(aes(colour = support, alpha = support), shape = 19, size = 4) +
  #geom_tiplab(aes(label=label), align = T, linetype = "dotted") +
  #geom_tiplab(data = tree.plot$data, aes(colour = species_drainage),offset = 0, size = 2,show.legend = FALSE) 
  #geom_tippoint(aes(colour = species_drainage)) +
  geom_text(aes(label = support.100), nudge_x = -0.002, nudge_y = 0, color = "deepskyblue", size = 10) +
  # scale_fill_manual(values = cbPalette) +
  scale_colour_gradient(low = "white", high = "black") +
  coord_cartesian(clip = 'off') +
  theme_tree(plot.margin = margin(2,2,2,2, "cm")) +
  guides(fill=guide_legend(title="Species and Drainage"), 
         colour=guide_legend(title="Bootstrap (Out of 100 reps)"),
         alpha="none") +
  geom_treescale() +
  theme(title = element_text(size = 25), legend.text = element_text(size=20), legend.title = element_text(size=25),
        plot.margin = unit(c(0,0,0,0), "cm"))

tree.plot.4

# MAP of samples
world.e <- data.frame(Long = c(-120,-60), Lat = c(8,38))
worldmap <- read_sf("4_Manuscript/data/World_data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
#Cropping the extent of worldmap to reduce plotting time
e <- c(xmin = -130, xmax = -60, ymin = 0, ymax = 45)
worldmap <- st_crop(worldmap, e)
crs_text <- "+proj=laea +x_0=0 +y_0=0 +lon_0=-74 +lat_0=40"
crs_text <- paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=",360+mean(sites$Long), " +lat_0=",mean(sites$Lat))
worldmap_proj <- st_transform(worldmap, crs = crs_text)

qtable$species.fig <- qtable$species
qtable$species.fig[qtable$species.fig=="titia"] <- "H. titia"
qtable$species.fig[qtable$species.fig=="americana"] <- "H. americana/cavlerti"

b <- ggplot() +
  geom_sf(data = worldmap_proj, fill = "#FFE6D4", col = "black") +
  geom_point(data = qtable, aes(Long, Lat, fill = assign_spp), size = 5, shape = 21) +
  scale_fill_manual(values = het.cols$cols) +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic")) +
  labs(x = "Longitude", y = "Latitude", col = "Species") +
  theme(panel.background = element_rect(fill = NA),
        panel.ontop = TRUE) +
  theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
        panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
  theme(legend.position = "none") +
  coord_sf(xlim = range(qtable$Long)+c(-5,+5), ylim =range(qtable$Lat)+c(-2,+2),default_crs = st_crs(4326))+
    facet_wrap(~species.fig, drop = T, nrow = 2) +
  theme(title = element_text(size = 25), legend.text = element_text(size=20), legend.title = element_text(size=25),
        strip.text.x = element_text(size = 20))
b
RAxML_plot <-b+(tree.plot.4/tree.plot.3) + plot_layout(width = c(1,2)) + plot_annotation(tag_levels = 'a',tag_prefix = "(", tag_suffix = ")")

# plot.tree.insert <- tree.plot.2 + inset_element(b, left = 0.05, bottom = 0.6, right = 0.6, top = 0.95) 
ggsave(plot = RAxML_plot, filename = paste0(plot.dir,"RAxML_tree_LEA_assignment.png"), height = 10, width = 15)
ggsave(plot = RAxML_plot, filename = paste0(plot.dir,"RAxML_tree_LEA_assignment.pdf"), height = 10, width = 15)
