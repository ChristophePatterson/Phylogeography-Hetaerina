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
# SNP.library.location  <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_v2/"

dir.path.SNAPP <- paste0(SNP.library.location, "/results/SNAPP/")

dir.plot <- "/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/4_Manuscript/plots/WGS_titia/"

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
