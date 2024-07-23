##################
####### SNAPP #######
##################
library(ggtree)
library(ggplot2)
library(ape)
library(treeio)
library(patchwork)

# read in anonotated tree
# Hetaerina_all_ddRAD_americana_dg OR Hetaerina_all_ddRAD_titia_dg
SNP.library.name <- "ipyrad"
# snap_Am_ti_default_div_est_N3 OR snap_Am_ti_Amdg_default_div_est_N3
# SNAPP.model <- "snap_Am_ti_default_div_est_N3"
dir.path <- paste0("4_Manuscript/data/SNP_libraries/",SNP.library.name,"/")
plot.dir <- paste0("4_Manuscript/plots/",SNP.library.name,"/")
num_samples <-3
dir.create(plot.dir)
SNAPP.model <- list.files(dir.path)[grep(pattern = "Anon",list.files(dir.path))]
SNAPP.model <- SNAPP.model[grep(pattern = paste0("N", num_samples),SNAPP.model)]

SNAPP.model <- gsub(SNAPP.model, pattern = ".trees.Anon",replacement = "")
SNAPP.tree <- read.mega(paste0(dir.path,SNAPP.model, ".trees.Anon"))

SNAPP.tree@data$height_0.95_HPD


SNAPP.tree@phylo$tip.label[SNAPP.tree@phylo$tip.label=="americana-Mex"] <- "americana-S"
SNAPP.tree@phylo$tip.label[SNAPP.tree@phylo$tip.label=="americana-USA"] <- "americana-N"
SNAPP.tree@phylo$tip.label[SNAPP.tree@phylo$tip.label=="calverti-Mex"] <- "calverti"
SNAPP.tree@data$posterior[SNAPP.tree@data$posterior<0.95] <- NA


p <- ggtree(SNAPP.tree) +
  geom_range(range = 'height_0.95_HPD', col = "deepskyblue", size = 4, alpha = .5) +
  geom_tiplab(aes(label=stringr::str_split_i(label, "_", 1))) +
  geom_nodepoint(aes(col=posterior>0.95), show.legend = F, size = 3) +
  scale_color_manual(values = c("black","white"), na.value = NA) +
  # geom_nodelab(aes(x = branch, label = posterior)) +
  #geom_nodelab(aes(x=-CAheight_mean, label=round(height,2)), size=5, vjust = -0.9, hjust = 1.1) +
  #scale_color_continuous(low="darkgreen", high="red") +
  theme(legend.position=c(.1, .8)) +
  #scale_x_range() +
  theme_tree2() +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin = margin(1,2.2,1,1, "cm")) +
  theme(panel.grid.major=element_line(colour = "grey70"),
        panel.grid.major.y = element_blank()) +
  xlab("Mya") 
#ggtitle(gsub(SNP.library.name,replacement = " ", pattern = "_"))
# Reverts xaxis and relabel x axis ticks to be positive
p <- revts(p) + scale_x_continuous(breaks = seq(-50, 0, 5),
                                   labels = abs(seq(50, 0, -5)))
# Note: Transparent ranges do not always show up in plot preview window
p

## STATS
divergence.stats <- na.omit(cbind(SNAPP.tree@data$CAheight_mean, 
                                  sapply(SNAPP.tree@data$height_0.95_HPD,"[[",1),
                                  sapply(SNAPP.tree@data$height_0.95_HPD,"[[",2)))
divergence.stats

# Read in posterior distribution
SNAPP.tree.nex <- read.nexus(paste0(dir.path, SNAPP.model, ".trees"))
# Remove burnin
burn.in <- 0.1
SNAPP.tree.nex[1:10]
SNAPP.tree.nex <- SNAPP.tree.nex[round(length(SNAPP.tree.nex)*burn.in):length(SNAPP.tree.nex)]
SNAPP.tree.nex[1]
# Go through each tree in turn and get the height
split.post.ti.am <- list()

max(phytools::nodeheight(SNAPP.tree.nex[[1]], node = 3))
# Identify nodes 
plot(SNAPP.tree.nex[2][[1]])
#tiplabels()
nodelabels(bg = NA, frame = "none")


# Get crown node
split.post.ti.am <- list()
for(i in 1:length(SNAPP.tree.nex)){
  # split.post.ti.am[[i]] <- max(phytools::nodeHeights(SNAPP.tree.nex[[i]])[,2])
  split.post.ti.am[[i]] <- phytools::nodeheight(SNAPP.tree.nex[i][[1]], node = 1)
}

#  Merge into one
split.post.ti.am <- do.call("c", split.post.ti.am)

#Plot as a histogram next to the proir
ti.am <- ggplot() +
  geom_density(aes(rnorm(1000000, mean = 33.08, sd = 5.53), fill = "prior"), alpha = 0.3) +
  geom_histogram(aes(split.post.ti.am, fill = "posterior", y = after_stat(density)), alpha = 0.6, ) +
  ggtitle(expression(italic("H. titia & H. americana/calverti"))) +
  scale_fill_manual(values = c("deepskyblue", "grey50")) +
  theme_bw() +
  xlab("mya") +
  ylab("freq") +
  labs(fill = "") +
  theme(legend.position = "bottom") +
  scale_x_reverse(limits = c(max(rnorm(1000000, mean = 33.08, sd = 5.53)),0))

#  Extract node height for split between H. am and H. cal
split.post.cal.am <- list()
for(i in 1:length(SNAPP.tree.nex)){
  split.post.cal.am[[i]] <- phytools::nodeheight(SNAPP.tree.nex[i][[1]], node = 8)
}
split.post.cal.am
#  Merge into one
split.post.cal.am <- do.call("c", split.post.cal.am)
plot(split.post.ti.am-split.post.cal.am)

# Plot 
am.cal <- ggplot() +
  geom_density(aes(rnorm(1000000, mean = 3.76, sd = 1.87), fill = "prior"), alpha = 0.3) +
  geom_histogram(aes(split.post.ti.am-split.post.cal.am, fill = "posterior", y = after_stat(density)), alpha = 0.6) +
  ggtitle(expression(italic("H. americana & H. calverti"))) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("deepskyblue", "grey50")) +
  xlab("mya") +
  ylab("freq") +
  labs(fill = "") +
  scale_x_reverse(limits = c(max(split.post.ti.am-split.post.cal.am),0))


split.post.tit.tit <- list()
# Get divergence times for Pac and atlantic titia
for(i in 1:length(SNAPP.tree.nex)){
  split.post.tit.tit[[i]] <- phytools::nodeheight(SNAPP.tree.nex[i][[1]], node = 10)
}

#  Merge into one
split.post.tit.tit <- do.call("c", split.post.tit.tit)
plot(split.post.ti.am-split.post.tit.tit)

lines(split.post.ti.am-split.post.tit.tit)

#Mean and SD of  split
median(split.post.ti.am-split.post.tit.tit)
sd(split.post.ti.am-split.post.tit.tit)

#plot
tit.tit <- ggplot() +
  #geom_histogram(aes(rnorm(length(split.post.ti.am), mean = 3.76, sd = 1.87), fill = "prior"), alpha = 0.5, ) +
  geom_histogram(aes(split.post.ti.am-split.post.tit.tit, fill = "posterior"), alpha = 0.5, show.legend = F) +
  ggtitle(expression(paste("Pacific & Atlantic ",italic("H. titia")))) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("deepskyblue", "grey50")) +
  xlab("mya") +
  ylab("freq") +
  labs(fill = "") +
  scale_x_reverse(limits = c(max(split.post.ti.am-split.post.cal.am),0))

split.post.ti.am


# Merge into one plot with tree
SNAPP.plot <- p + (tit.tit/am.cal/ti.am) + plot_layout(widths = c(2,1)) + 
  plot_annotation(tag_levels = 'a',tag_prefix = "(", tag_suffix = ")")

ggsave(plot = SNAPP.plot, filename = paste0(plot.dir, "SNAPP_", SNAPP.model,".png"), width = 11, height = 8)
ggsave(plot = SNAPP.plot, filename = paste0(plot.dir, "SNAPP_", SNAPP.model, ".pdf"), width = 11, height = 8)
  

# Identify nodes 
plot(SNAPP.tree.nex[2][[1]])
#tiplabels()
nodelabels(bg = NA, frame = "none")

## Final split between titia NA and titia SA
split.post.titN.titS <- list()
for(i in 1:length(SNAPP.tree.nex)){
  split.post.titN.titS[[i]] <- phytools::nodeheight(SNAPP.tree.nex[i][[1]], node = 88)
}

#  Merge into one
split.post.titN.titS <- do.call("c", split.post.titN.titS)
plot(split.post.ti.am-split.post.titN.titS)


### Calculate ESS
library(coda)
# Create trace file, need to minus tree height to get node age
trace <- mcmc(cbind(split.post.ti.am  = split.post.ti.am, 
                    split.post.cal.am = split.post.ti.am-split.post.cal.am, 
                    split.post.tit.tit= split.post.ti.am-split.post.tit.tit,
                    split.post.titN.titS= split.post.ti.am-split.post.titN.titS))

coda::effectiveSize(trace)
coda::traceplot(trace)
coda::crosscorr.plot(trace)

coda::HPDinterval(trace)
colMeans(trace)

trace.df <- as.data.frame(trace)


ggplot(trace.df) +
  geom_point(aes(split.post.tit.tit, split.post.cal.am, col = split.post.cal.am<split.post.tit.tit)) +
  geom_abline(slope = 1, intercept = 0)

ggplot(trace.df) +
  geom_histogram(aes(split.post.cal.am-split.post.tit.tit, fill = split.post.cal.am<split.post.tit.tit))

summary(trace.df$split.post.cal.am-trace.df$split.post.tit.tit)

sum(trace.df$split.post.cal.am>=trace.df$split.post.tit.tit)/length(trace.df$split.post.tit.tit)
HPDinterval(mcmc(trace.df$split.post.cal.am-trace.df$split.post.tit.tit))


