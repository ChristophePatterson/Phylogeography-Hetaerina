# Get system arguments
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")

library(vcfR)
library(ggtree)
library(ggplot2)
library(ape)
library(treeio)

# Load in arguments
args <- commandArgs(trailingOnly = TRUE)

SNP.library.name <- args[2]
# SNP.library.name <- "Hetaerina_all_ddRAD_americana_dg"
SNP.library.location <- args[1]
SNP_plot_name <- args[3]
# SNP.library.location  <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_v2/"

dir.path.SNAPP <- paste0(SNP.library.location, "/results/SNAPP/")

dir.plot <- paste0("/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/4_Manuscript/plots/",
                    SNP_plot_name,"/", SNP_plot_name,"/")

SNAPP.model <- list.files(dir.path.SNAPP)
SNAPP.model <- SNAPP.model[grep(pattern = ".trees.Anon", SNAPP.model)]

SNAPP.model <- gsub(SNAPP.model, pattern = ".trees.Anon",replacement = "")
SNAPP.tree <- read.mega(paste0(dir.path.SNAPP, SNAPP.model, ".trees.Anon"))

# get second deepest node
xmin <- sort(do.call("c", lapply(SNAPP.tree@data$height_0.95_HPD, max)), decreasing = T)[2]

SNAPP.tree@data$height_0.95_HPD

SNAPP.tree@data$height_0.95_HPD_1 <- sapply(SNAPP.tree@data$height_0.95_HPD, "[[", 1)
SNAPP.tree@data$height_0.95_HPD_2 <- sapply(SNAPP.tree@data$height_0.95_HPD, "[[", 2)
SNAPP.tree@data[13,]$height_0.95_HPD

xmin <- sort(SNAPP.tree@data$height_0.95_HPD_2, decreasing = T)[2]


p <- ggtree(SNAPP.tree) +
  geom_segment(aes(x = -height_0.95_HPD_1, xend = -height_0.95_HPD_2, y = y, yend = y), col = "deepskyblue", size = 4, alpha = .5) +
  geom_tiplab(aes(x = 0, label = gsub("titia_", "  ", label))) +
  geom_nodepoint(aes(col = posterior)) +
  scale_color_gradient(low = "white", high = "black") +
  # geom_nodelab(aes(x = branch, label = posterior)) +
  geom_nodelab(aes(x=-CAheight_mean, label=round(CAheight_mean,4)*1000000), size=5, vjust = -0.9, hjust = 1.1) +
  #scale_color_continuous(low="darkgreen", high="red") +
  theme(legend.position=c(.1, .8)) +
  #scale_x_range() +
  theme_tree2() +
  coord_cartesian(clip = 'off', xlim = c(-xmin, 0)) +
  scale_x_continuous(labels=function(x)x*1000000) +
  theme_tree2(plot.margin = margin(1,2.2,1,1, "cm")) +
  theme(panel.grid.major=element_line(colour = "grey70"),
        panel.grid.major.y = element_blank(), legend.position = "bottom") +
  xlab("Years") 
#ggtitle(gsub(SNP.library.name,replacement = " ", pattern = "_"))
# Reverts xaxis and relabel x axis ticks to be positive
p <- revts(p) 
# Note: Transparent ranges do not always show up in plot preview window

ggsave(paste(dir.plot, "WGS_SNAPP.png"), p)



# Read in posterior distribution
SNAPP.tree.nex <- read.nexus(paste0(dir.path.SNAPP, SNAPP.model, ".trees"))
# Remove burnin
burn.in <- 0.1
SNAPP.tree.nex[1:10]
SNAPP.tree.nex <- SNAPP.tree.nex[round(length(SNAPP.tree.nex)*burn.in):length(SNAPP.tree.nex)]
SNAPP.tree.nex[1]

# Get divergence times for Pac and atlantic titia
plot(SNAPP.tree.nex[1][[1]],use.edge.length = F)
tiplabels()
nodelabels()

SNAPP.tree.nex[1][[1]]$edge
SNAPP.tree.nex[1][[1]]$tip.label

CUAJa03_node <- vector()
CUAJa03_node_height <- vector()
for(i in 1:length(SNAPP.tree.nex)){
  tip.num <- which(SNAPP.tree.nex[1][[1]]$tip.label=="titia_CUAJa03")
  if(length(tip.num)>1){print("ERROR")}
  CUAJa03_node_temp <- SNAPP.tree.nex[i][[1]]$edge[,2]==tip.num
  CUAJa03_node[i] <- SNAPP.tree.nex[i][[1]]$edge[CUAJa03_node_temp,1]
  CUAJa03_node_height[i] <- phytools::nodeheight(SNAPP.tree.nex[i][[1]], node = CUAJa03_node[i])-phytools::nodeheight(SNAPP.tree.nex[i][[1]], node = 1)
}
hist(CUAJa03_node_height*1000000)
plot(CUAJa03_node_height*1000000, col = CUAJa03_node)
lines(CUAJa03_node_height*1000000)