# Manuscript_results_SVDQuartets and RAxML

library(ape)
library(tidyverse)
library(ggtree)
library(patchwork)
library(sf)
library(treeio)
library(phytools)
library(phangorn)
library(scatterpie)
# Get colour pallettes
het.cols <- data.frame(assign = c(paste0("americana_",1:3),paste0("titia_",1:3)),
                       cols = c("#726230","#AA9599","#548A39","#AF0F09","#E5D9BA","#3E3C3A"))
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
qameri <- read.table("4_Manuscript/data/SNP_libraries/Hetaerina_americana_ddRAD_americana_dg/LEA_qassign_pca_samples.txt")
qameri$assign_spp <- paste0(qameri$species,"_", qameri$assign)
names(qtitia)
names(qameri$sample)
qtable <- rbind(qtitia, qameri[,names(qtitia)])

# Merge with LEA assignment data
samples <- tre.midroot$tip.label
sites <- data.frame(samples)
#Convert lat and long to number

sites <- merge(sites, qtable, by.x = "samples", by.y = "samples")
table(sites$assign_spp)
any(is.na(sites$assign_spp))

# Order hetr cols correctly
het.cols$cols[het.cols$assign==paste0("titia_",sites$assign[sites$site.sub=="ZANA"][1])] <- "#AF0F09"
het.cols$cols[het.cols$assign==paste0("titia_",sites$assign[sites$site.sub=="ESRB"][1])] <- "#E5D9BA"
het.cols$cols[het.cols$assign==paste0("titia_",sites$assign[sites$site.sub=="HCAR"][1])] <- "#3E3C3A"
sites$assign_spp_name[sites$assign_spp==(sites$assign_spp[which(sites$site.sub=="ZANA")][1])] <- "Pacific"
sites$assign_spp_name[sites$assign_spp==(sites$assign_spp[which(sites$site.sub=="ESRB")][1])] <- "South Atlantic"
sites$assign_spp_name[sites$assign_spp==(sites$assign_spp[which(sites$site.sub=="HCAR")][1])] <- "North Atlantic"

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

titia_cluster_label <- group_by(tree.plot$data, assign_spp_name) %>%
  summarise(mn.y = median(y), mn.x = max(tree.plot$data$x)+(max(tree.plot$data$x)*0.07))
max(tree.plot$data$x)
min(tree.plot$data$x)  
tree.plot.3 <- tree.plot +
  geom_segment(aes(x = x, xend = max(x)+max(x*0.02), y=y, yend=y), col = tree.plot$data$assign_spp.col) +
  geom_tippoint(aes(x = max(x)+max(x*0.02), col = assign_spp),col = tree.plot$data$assign_spp.col[tree.plot$data$isTip],
                shape =15, size = 3) +
  #geom_nodepoint(aes(colour = support, alpha = support), shape = 19, size = 4) +
  #geom_tiplab(aes(label=label), align = T, linetype = "dotted") +
  #geom_tiplab(data = tree.plot$data, aes(colour = species_drainage),offset = 0, size = 2,show.legend = FALSE) 
  #geom_tippoint(aes(colour = species_drainage)) +
  geom_text(aes(label = support.100), nudge_x = -0.002, nudge_y = 0, color = "deepskyblue", size = 15) +
  geom_text(data = titia_cluster_label, aes(x = mn.x, y = mn.y, label = assign_spp_name),
            angle= -90, size = 8) +
  # scale_fill_manual(values = cbPalette) +
  scale_colour_gradient(low = "white", high = "black") +
  coord_cartesian(clip = 'off') +
  theme_tree(plot.margin = margin(2,2,2,2, "cm")) +
  guides(fill=guide_legend(title="Species and Drainage"), 
         colour=guide_legend(title="Bootstrap (Out of 100 reps)"),
         alpha="none") +
  geom_treescale(width = 0.05, fontsize = 6, linesize = 1) +
  theme(title = element_text(size = 25), legend.text = element_text(size=20), legend.title = element_text(size=25),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggtitle(label = expression(paste("(a) ", italic("H. titia"))))
  


tree.plot.3
# Read in americana tree
SNP.library.name <- "Hetaerina_americana_ddRAD_americana_dg"
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
qameri <- read.table("4_Manuscript/data/SNP_libraries/Hetaerina_americana_ddRAD_americana_dg/LEA_qassign_pca_samples.txt")
qameri$assign_spp <- paste0(qameri$species,"_", qameri$assign)
qtable <- rbind(qtitia, qameri[,names(qtitia)])

# Merge with LEA assignment data
samples <- tre.midroot$tip.label
sites <- data.frame(samples)
#Convert lat and long to number

sites <- merge(sites, qtable, by.x = "samples", by.y = "samples")
table(sites$assign_spp)
any(is.na(sites$assign_spp))

het.cols$cols[het.cols$assign==paste0("americana_",sites$assign[sites$site.sub=="RG"][1])] <- "#AA9599"
het.cols$cols[het.cols$assign==paste0("americana_",sites$assign[sites$site.sub=="STDM"][1])] <- "#548A39"
het.cols$cols[het.cols$assign==paste0("americana_",sites$assign[sites$site.sub=="TN"][1])] <- "#726230"
sites$assign_spp_name[sites$assign_spp==(sites$assign_spp[which(sites$site.sub=="RG")][1])] <- "italic(H.~calverti)"
sites$assign_spp_name[sites$assign_spp==(sites$assign_spp[which(sites$site.sub=="STDM")][1])] <- "~italic(H.~americana)~South"
sites$assign_spp_name[sites$assign_spp==(sites$assign_spp[which(sites$site.sub=="TN")][1])] <- "~italic(H.~americana)~North"



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


americana_cluster_label <- group_by(tree.plot$data, assign_spp_name) %>%
  summarise(mn.y = median(y), mn.x = max(tree.plot$data$x)+(max(tree.plot$data$x)*0.07))

tree.plot.4 <- tree.plot +
  geom_segment(aes(x = x, xend = max(x)+max(x)*0.02, y=y, yend=y), col = tree.plot$data$assign_spp.col) +
  geom_tippoint(aes(x = max(x)+max(x)*0.02, col = assign_spp),col = tree.plot$data$assign_spp.col[tree.plot$data$isTip],
                shape =15, size = 3) +
  #geom_nodepoint(aes(colour = support, alpha = support), shape = 19, size = 4) +
  #geom_tiplab(aes(label=label), align = T, linetype = "dotted") +
  #geom_tiplab(data = tree.plot$data, aes(colour = species_drainage),offset = 0, size = 2,show.legend = FALSE) 
  #geom_tippoint(aes(colour = species_drainage)) +
  geom_text(aes(label = support.100), nudge_x = -0.002, nudge_y = 0, color = "deepskyblue", size = 15) +
  geom_text(data = americana_cluster_label, aes(x = mn.x, y = mn.y, label = assign_spp_name),
            angle= -90, size = 8, parse = T) +
  # scale_fill_manual(values = cbPalette) +
  scale_colour_gradient(low = "white", high = "black") +
  coord_cartesian(clip = 'off') +
  theme_tree(plot.margin = margin(2,2,2,2, "cm")) +
  guides(fill=guide_legend(title="Species and Drainage"), 
         colour=guide_legend(title="Bootstrap (Out of 100 reps)"),
         alpha="none") +
  geom_treescale(width = 0.05, fontsize = 6, linesize = 1, x = 0.05) +
  theme(title = element_text(size = 25), legend.text = element_text(size=20), legend.title = element_text(size=25),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggtitle(label = expression(paste("(b) ", italic("H. americana"), " sensu lato")))

tree.plot.4


RAxML_plot <- tree.plot.3+tree.plot.4
ggsave(plot = RAxML_plot, filename = paste0(plot.dir,"RAxML_tree_LEA_assignment_no_map.png"), height = 20, width = 15)
ggsave(plot = RAxML_plot, filename = paste0(plot.dir,"RAxML_tree_LEA_assignment_no_map.pdf"), height = 20, width = 15)


# MAP of samples
world.e <- data.frame(Long = c(-120,-60), Lat = c(8,38))
worldmap <- read_sf("4_Manuscript/data/World_data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
#Cropping the extent of worldmap to reduce plotting time
e <- c(xmin = -130, xmax = -60, ymin = 0, ymax = 45)
worldmap <- st_crop(worldmap, e)
crs_text <- "+proj=laea +x_0=0 +y_0=0 +lon_0=-74 +lat_0=40"
crs_text <- paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=",360+mean(sites$Long), " +lat_0=",mean(sites$Lat))
worldmap_proj <- st_transform(worldmap, crs = crs_text)

# Get list of unique populations
pop <- unique(qtable$site.sub.county.drain)

#Number of unique sites
Npop = length(pop)

qtable$species.fig <- qtable$species
qtable$species.fig[qtable$species.fig=="titia"] <- "H. titia"
qtable$species.fig[qtable$species.fig=="americana"] <- "H. americana/cavlerti"

qmap <- as.data.frame(matrix(NA, nrow = Npop, ncol = length(unique(qtable$assign_spp)) + 1))
colnames(qmap) <- c("pop", unique(qtable$assign_spp))
qmap[,1] <- pop
spp <- unique(qtable$assign_spp)

# get counts of species for each population
tmp.count <- vector()
assign.count <- apply(qmap, MARGIN = 1, function(x) {
  for(i in 1:length(spp)){
    tmp.spp <- spp[i]
    qtmp <- qtable[qtable$site.sub.county.drain==x[1],]
    tmp.count[i] <- length(which(qtmp$assign_spp==tmp.spp))
    }
  return(tmp.count)
  })
# Combine
assign.count <- cbind(pop, t(assign.count))
colnames(assign.count) <- c("pop", spp)
assign.count <- as.data.frame(assign.count)
# Calculate mean lat and lon
assign.count$Lat <- apply(assign.count, MARGIN = 1, function(x) mean(qtable$Lat[qtable$site.sub.county.drain==x[1]]))
assign.count$Lon <- apply(assign.count, MARGIN = 1, function(x) mean(qtable$Lon[qtable$site.sub.county.drain==x[1]]))
# COnvert to numeric
for(colspp in which(colnames(assign.count)%in%spp)){
  assign.count[,colspp] <- as.numeric(assign.count[,colspp])
}

## Split into titia and americana
assign.count.titia <- assign.count[,c(1, 8, 9 , grep("titia", colnames(assign.count)))]
assign.count.titia$num_samples <- rowSums(assign.count.titia[,grep("titia", colnames(assign.count.titia))])
#
assign.count.ameri <- assign.count[,c(1, 8, 9 , grep("ameri", colnames(assign.count)))]
assign.count.ameri$num_samples <- rowSums(assign.count.ameri[,grep("ameri", colnames(assign.count.ameri))])

p.tit <- ggplot() +
  geom_sf(data = worldmap_proj, fill = "#FFE6D4", col = "black") +
  geom_scatterpie(data = assign.count.titia, aes(x = Lon, y = Lat, r = (sqrt(num_samples)/pi)*2), cols = paste0("titia_",1:3), show.legend = F) +
  scale_fill_manual(values = het.cols$cols[match(paste0("titia_",1:3), het.cols$assign)]) +
  theme(strip.text = element_text(face = "italic"),
      legend.text = element_text(face = "italic")) +
  labs(x = "Longitude", y = "Latitude", col = "Species") +
  theme(panel.background = element_rect(fill = NA, color = "black"),
        panel.ontop = TRUE) +
  theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
        panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
  theme(legend.position = "none") +
  coord_sf(xlim = range(qtable$Long)+c(-5,+5), ylim =range(qtable$Lat)+c(-2,+2),default_crs = st_crs(4326))+
  # facet_wrap(~species.fig, drop = T, nrow = 2) +
  theme(title = element_text(size = 25), legend.text = element_text(size=20), legend.title = element_text(size=25),
        strip.text.x = element_text(size = 20))
  # ggtitle(label = expression(paste(italic("H. titia"))))

p.am <- ggplot() +
  geom_sf(data = worldmap_proj, fill = "#FFE6D4", col = "black") +
  geom_scatterpie(data = assign.count.ameri, aes(x = Lon, y = Lat, r = (sqrt(num_samples)/pi)*2), cols = paste0("americana_",1:3), show.legend = F) +
  scale_fill_manual(values = het.cols$cols[match(paste0("americana_",1:3), het.cols$assign)]) +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic")) +
  labs(x = "Longitude", y = "Latitude", col = "Species") +
  theme(panel.background = element_rect(fill = NA, color = "black"),
        panel.ontop = TRUE) +
  theme(panel.grid.major = element_line(color = rgb(0,0,0,alpha = 0.1)),
        panel.grid.minor = element_line(color = rgb(0,0,0,alpha = 0.1))) +
  theme(legend.position = "none") +
  coord_sf(xlim = range(qtable$Long)+c(-5,+5), ylim =range(qtable$Lat)+c(-2,+2),default_crs = st_crs(4326))+
  # facet_wrap(~species.fig, drop = T, nrow = 2) +
  theme(title = element_text(size = 25), legend.text = element_text(size=20), legend.title = element_text(size=25),
        strip.text.x = element_text(size = 20))
  #ggtitle(label = expression(paste(italic("H. americana/calverti"))))
  
  
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
    facet_wrap(~species.fig, drop = T, nrow = 1) +
  theme(title = element_text(size = 25), legend.text = element_text(size=20), legend.title = element_text(size=25),
        strip.text.x = element_text(size = 20))
b

tree.plot.4.p <- tree.plot.4 + theme(title = element_blank())
tree.plot.3.p <- tree.plot.3 + theme(title = element_blank())

RAxML_plot <- (tree.plot.3 + tree.plot.4)/ plot_spacer() / (p.tit + p.am) + plot_layout(heights = c(9,0.25,3), width = c(1,1)) # + plot_annotation(tag_levels = 'a',tag_prefix = "(", tag_suffix = ")")
#theme(text=element_text(size = 10))

# plot.tree.insert <- tree.plot.2 + inset_element(b, left = 0.05, bottom = 0.6, right = 0.6, top = 0.95) 
ggsave(plot = RAxML_plot, filename = paste0(plot.dir,"RAxML_tree_LEA_assignment.png"), height = 22, width = 13)
ggsave(plot = RAxML_plot, filename = paste0(plot.dir,"RAxML_tree_LEA_assignment.pdf"), height = 22, width = 13)

