## Read in trace file
library(coda)
library(ggtree)
library(ggplot2)
library(patchwork)
library(tidyverse)

# # # # # # ## # # # # # ## # # # # # #
##### Heterospecific plots ####
# # # # # # ## # # # # # ## # # # # # #
### Produce plots that use hetero specific SNPs for g-phocs
input.dir.top <- "4_Manuscript/data/G-Phocs/"
N_select <- 2000
model.names <- paste0(c("americana_dg/G-Phocs-a1-b20-titia-Iso-N","titia_dg/G-Phocs-a1-b20-americana-Iso-N",
                 "americana_dg/G-Phocs-a1-b20-titia-migR-N","titia_dg/G-Phocs-a1-b20-americana-migR-N"),N_select) 
model.names <- c(model.names, c("denovo/G-Phocs-a1-b20-titia-Iso-N3","denovo/G-Phocs-a1-b20-americana-Iso-N3",
                                "denovo/G-Phocs-a1-b20-titia-migR-N3","denovo/G-Phocs-a1-b20-americana-migR-N3"))


u = 2.8*10^-9 
plot.dir <- "4_Manuscript/plots/G-Phocs/"
run <- "denovo/G-Phocs-a1-b20-titia-Iso-N3"
# Blank lists for filling in data
p <- list()
effect.size <- list()
run.df <- list()
p.map <- list()
HPD.trace.list <- list()
trace.df <- data.frame("tau_amer" = NA, "tau_amer.calv" = NA, "tau_titia.Atl" = NA,"tau_titia" = NA,
                       "theta_amer.N" = NA,"theta_amer.S"=NA, "theta_calv.1"= NA,
                       "theta_titia.Pac"=NA, "theta_titia.NAtl"=NA, "theta_titia.SAtl"=NA)
# Loop through every model
for(run in model.names){ 
  rm(trace.name,trace,trace.raw,trace.means.raw,trace.melt,trace.mig,HPD.trace,tree.text)
  trace.name <- run
  run <- stringr::str_split_i(run, pattern = "/", 2)
  print(trace.name)
  trace <- read.table(file = paste0(input.dir.top,
                                    trace.name, ".trace"), header = T)
  # Remove burn in
  burn.in <- round(length(trace$Sample)*0.1, 0)
  trace <- trace[burn.in:length(trace$Sample),]
  
  # Plot tracer
  # plot(trace[seq(1,length(trace$Sample),length.out = 5000),grep(names(trace), pattern = "tau")])
  
  effect.size[[trace.name]] <- effectiveSize(trace[,-1])
  effect.size[[trace.name]]
  
  # Allow visulasaion of trace file in ggplot
  trace.melt <- reshape2::melt(trace[,c(1,grep(names(trace), pattern = "tau"))], id.vars="Sample")
  table(trace.melt$variable)
  trace.melt <- reshape2::melt(trace[seq(1,length(trace$Sample),length.out = 5000),], id.vars="Sample")
  
  #ggplot(trace.melt) +
  #  geom_line(aes(y = value, x = Sample)) +
  #  facet_wrap(~variable, scales = "free", ncol = 2)
  
  # Raw means
  trace.raw <- trace
  trace.means.raw <- colMeans(trace.raw[,-1], na.rm = T)
  
  # Migration rates
  head(trace)
  trace.mig <- data.frame(mig = names(trace.raw[,grep(names(trace), pattern = "m_")]))
  if(dim(trace.mig)[1]>=1){
    trace.mig$mig.out <- gsub(sapply(strsplit(trace.mig$mig,"\\.\\."),"[[",1), pattern = "m_", replacement = "")
    trace.mig$mig.in <- sapply(strsplit(trace.mig$mig,"\\.\\."),"[[",2)
    if(dim(trace.mig)[1]>1){
      trace.mig$mig.rate <- colMeans(trace.raw[,grep(names(trace.raw), pattern = "m_")], na.rm = T)/0.001
    } else if(dim(trace.mig)[1]==1){trace.mig$mig.rate <- mean(trace.raw[,grep(names(trace), pattern = "m_")], na.rm = T)/0.001}
    
  }
  trace.mig$mig.in.theta <- do.call("c",lapply(trace.mig$mig.in, function(x) trace.means.raw[which(names(trace.means.raw)==paste0("theta_",x))]))
  trace.mig$mig.in.theta <- ((trace.mig$mig.in.theta/1000)/(4*u))
  trace.mig$mig.rate.prop <- trace.mig$mig.rate*u
  trace.mig$mig.rate.num.ind.per.gen <- trace.mig$mig.in.theta*trace.mig$mig.rate.prop
  trace.mig$mig.rate.num.ind.per.gen.pct <- trace.mig$mig.rate.num.ind.per.gen*100
  
  
  # Convert output to mya and number of individuals (millions)
  trace[,grep(names(trace), pattern = "tau")] <- ((trace[,grep(names(trace), pattern = "tau")]/1000)/u)*10^-6
  trace[,grep(names(trace), pattern = "m_")] <- ((trace[,grep(names(trace), pattern = "m_")]/0.001)*u)
  trace[,grep(names(trace), pattern = "theta")] <- ((trace[,grep(names(trace), pattern = "theta")]/1000)/(4*u))*10^-6
  
  
  for(i in colnames(trace.df)[colnames(trace.df)%in%names(trace)]){
    trace.df[run,i] <- mean(trace[,i])
  }
  
  #ggplot(trace.melt[]) +
  # geom_violin(aes(x = value, y = variable),orientation = "y") +
  # facet_wrap(~variable, scales = "free")
  
  # calculate HPD for all tau traces
  HPD.trace <- HPDinterval(as.mcmc(trace[,grep(names(trace), pattern = "tau|theta")]))
  HPD.trace.list[[trace.name]] <- HPDinterval(as.mcmc(trace[,grep(names(trace), pattern = "tau|theta")]))
  HPD.trace.list[[trace.name]] <- cbind(HPD.trace.list[[trace.name]], mean = colMeans(as.mcmc(trace[,grep(names(trace), pattern = "tau|theta")])))
  
  # Add to trace.df
  trace.df
  
  trace.means <- colMeans(trace[,-1], na.rm = T)
  names(trace.means)
  trace.means
  # Read in tree file
  if(grepl(run, pattern = "titia")){
    tree.text <- paste0("((titia.NAtl:",trace.means["tau_titia.Atl"],",titia.SAtl:",trace.means["tau_titia.Atl"],")titia.Atl:",trace.means["tau_titia"]-trace.means["tau_titia.Atl"]
                        ,",titia.Pac:",trace.means["tau_titia"],")titia:",trace.means["tau_titia"],";")
  }
  if(grepl(run, pattern = "americana")){
    tree.text <- paste0("((amer.N:",trace.means["tau_amer"],",amer.S:",trace.means["tau_amer"],")amer:",trace.means["tau_amer.calv"]-trace.means["tau_amer"]
                        ,",calv.1:",trace.means["tau_amer.calv"],")amer.calv:",trace.means["tau_amer.calv"],";")
  }
  
  tree <- ape::read.tree(text = tree.text)
  #plot(tree)
  # Merge HPD and tree
  tree.plot <- ggtree(tree) %<+% data.frame(node.labes = gsub(row.names(HPD.trace),pattern = "tau_", replacement= ""),HPD.trace)
  
  # Add in HPD in format that ggtree understands (list in a tibble)
  height_0.95_HPD_list <- list()
  for(i in 1:length(tree.plot$data$parent)){
    if(tree.plot$data$lower[i]>0&!is.na(tree.plot$data$lower[i])){
      height_0.95_HPD_list[[i]] <- c(tree.plot$data$lower[i],tree.plot$data$upper[i])
    }else height_0.95_HPD_list[[i]] <- c(tree.plot$data$x[i],tree.plot$data$x[i])
  }
  
  ## Add effective population size
  theta.means <- trace.means[grep(pattern = "theta",names(trace.means))]
  theta.means <- theta.means[match(paste0("theta_", tree.plot$data$label),names(theta.means))]
  
  cbind(paste0("theta_", tree.plot$data$label),theta.means)
  
  tree.plot$data$theta <- theta.means
  
  # Add to tree data
  tree.plot$data$height_0.95_HPD <- height_0.95_HPD_list
  tree.height <- max(ape::node.depth.edgelength(tree))
  # Get tree corrdinates for mig rates
  trace.mig$mig.out.x <- do.call("c",lapply(trace.mig$mig.out, function(x)  tree.plot$data$x[which(tree.plot$data$label==x)]-(tree.plot$data$branch.length[which(tree.plot$data$label==x)]/2)))
  trace.mig$mig.out.y <- do.call("c",lapply(trace.mig$mig.out, function(x)  tree.plot$data$y[which(tree.plot$data$label==x)]))
  trace.mig$mig.in.x <- do.call("c",lapply(trace.mig$mig.in, function(x)  tree.plot$data$x[which(tree.plot$data$label==x)]-(tree.plot$data$branch.length[which(tree.plot$data$label==x)]/2)))
  trace.mig$mig.in.y <- do.call("c",lapply(trace.mig$mig.in, function(x)  tree.plot$data$y[which(tree.plot$data$label==x)]))
  trace.mig$mig.out.x <-  tree.height-trace.mig$mig.out.x+rep(c(-0.1,0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.in.x <- tree.height-trace.mig$mig.in.x+rep(c(-0.1,0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.out.y <- trace.mig$mig.out.y+rep(c(0.1,-0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.in.y <- trace.mig$mig.in.y+rep(c(-0.1,0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.out.y.label <- trace.mig$mig.out.y-rep(c(-0.6,0.6),length.out = length(trace.mig$mig))
  trace.mig$mig.rate.perc <- (trace.mig$mig.rate*u)*100
  
  tree.plot$data$x.theta <- tree.height-tree.plot$data$x+(tree.plot$data$branch.length/2)
  tree.plot$data$y.theta <- tree.plot$data$y
  tree.plot$data$x.theta[tree.plot$data$label=="titia"|tree.plot$data$label=="amer.calv"] <- tree.height
  tree.plot$data$y.theta[tree.plot$data$label=="titia"|tree.plot$data$label=="amer.calv"] <- tree.plot$data$y[tree.plot$data$label=="titia"|tree.plot$data$label=="amer.calv"]-0.1
  
  trace.mig$mig.out.y.label
  # Plot tree
  
  p[[trace.name]] <- tree.plot +
    geom_tree(aes(linewidth = theta)) +
    #geom_rootedge(size = as.numeric(tree.plot$data$theta[4]*10),rootedge = 2) +
    #geom_segment(range = 'height_0.95_HPD', col = "deepskyblue", size = 4, alpha = .5) +
    geom_segment(aes(x = -lower, xend = -upper, y = , yend = y),col = "deepskyblue", linewidth = 10, alpha = .5,show.legend=F) +
    #geom_point(aes(x, y, size = theta)) +
    geom_tiplab(hjust = -0.2) +
    #coord_cartesian(clip = 'off') +
    #geom_nodelab(aes(x=x, label= paste(sprintf("%.1f", round(-x, digits = 1)))),size = 5, hjust= -0.4,vjust = 0.2,show.legend=F) +
    #geom_label(aes(x=x,y = y, label= paste0("Ne=",sprintf("%.1f", round(theta,1)))), size = 2.5, show.legend=F) +
    geom_label(aes(x=-x.theta,y = y.theta, label= paste0("Ne=",sprintf("%.1f", round(theta,1)))), size = 3, show.legend=F) +
    #geom_nodelab(aes(x=x, y=y, label=label, colour = "blue")) +
    #geom_treescale(fontsize=15, linesize=3, offset=2) +
    theme_tree2() +
    #coord_cartesian(clip = 'off') +
    theme_tree2(plot.margin = margin(1,1,1,1, "cm")) +
    theme(panel.grid.major=element_line(colour = "grey70"),
          panel.grid.major.y = element_blank()) +
    expand_limits(x = c(1,-5)) +
    xlab("Mya") +
    labs(size = "Effective population\nsize (millions)")
  
  # Option to add in migration rate to 
  if(dim(trace.mig)[1]>=1){
    p[[trace.name]] <- p[[trace.name]] + geom_segment(data = trace.mig, aes(x = -mig.out.x, y = mig.out.y, xend = -mig.in.x, yend = mig.in.y), arrow = arrow(length = unit(0.2,"cm"))) +
     geom_label(data = trace.mig, aes(x = -mig.out.x, y = mig.out.y.label, label = paste0(round(mig.rate.num.ind.per.gen, 3))), size = 2.5)
  }
  
  p[[trace.name]] <- revts(p[[trace.name]]) + scale_x_continuous(breaks = seq(-20, 0, 1),
                                                                 labels = abs(seq(20, 0, -1)))
  
  ggsave(plot = (p[[trace.name]]), filename = paste0(plot.dir,run,"-u-",
                                              gsub(u, replacement = "-", pattern ="\\."),".png"),
         width =12, height = 6)
  
}

names(p)
# Heterospecific Loci calling 
all_plot <- (p[[1]] + ggtitle(expression(paste(italic("H. titia"), " - No Migration")))) + 
  (p[[2]] + ggtitle(expression(paste(italic("H. americana/calverti"), " - No Migration")))) + 
  (p[[3]] + ggtitle(expression(paste(italic("H. titia"), " - Migration")))) +
  (p[[4]] + ggtitle(expression(paste(italic("H. americana/calverti"), " - Migration")))) +
  plot_annotation(title = "Heterospecific Loci calling", theme = theme(plot.title = element_text(size = 25)),
                  tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(ncol = 2)

ggsave(file = paste0(plot.dir, "Hetero_g-phocs-plots.png"), all_plot, height = 10, width = 12) 
ggsave(file = paste0(plot.dir, "Hetero_g-phocs-plots.pdf"), all_plot, height = 10, width = 12) 


# Merge all tau
trace.df <- trace.df[row.names(trace.df)!=1,]
trace.df$model <- row.names(trace.df)

trace.df$species <- sapply(strsplit(trace.df$model, "-"), "[", 5)
trace.df$mig.type <- sapply(strsplit(trace.df$model, "-"), "[", 6)
trace.df$beta <- sapply(strsplit(trace.df$model, "-"), "[", 4)
names(trace.df)

HPD.trace.list

trace.df
write.table(trace.df, paste0("4_Manuscript/data/G-Phocs/","G-Phocs_results_heterospecific_genome_", gsub(u, replacement = "-", pattern ="\\."),".txt"),
            row.names = F, quote=F)
write.csv(trace.df, paste0("4_Manuscript/data/G-Phocs/","G-Phocs_results_heterospecific_genome_", gsub(u, replacement = "-", pattern ="\\."),".csv"),
            row.names = F, quote=F)



# # # # # # ## # # # # # ## # # # # # #
#### ALL titia draft genome runs ####
# # # # # # # # # # # # # # ## # # # # # #

input.dir.top <- paste0("4_Manuscript/data/G-Phocs/")
draft_genome <- "titia_dg"
input.dir <- paste0(input.dir.top, draft_genome, "/")
model.names <- list.files(input.dir)
model.names <- model.names[grep(model.names, pattern = ".trace")]
model.names <- gsub(x = model.names, pattern = ".trace",replacement =  "")

trace.files <- paste0(model.names,".trace")
#trace.files <- trace.files[1]
#trace.name <- trace.files

# Broad use mutation rate used in other insect papers
u = 2.8*10^-9 
# Number of differences betweeen H. titia and I. elegans divide by divergence time 121 million years
# u = 2.2e-09
# Output folder for plots 
plot.dir <- "4_Manuscript/plots/G-Phocs/"

plot.dir <- paste0(plot.dir,draft_genome, "_",
                   gsub(u, replacement = "-", pattern ="\\."),"/")
dir.create(plot.dir)

# Blank lists for filling in data
p <- list()
effect.size <- list()
run.df <- list()
p.map <- list()
HPD.trace.list <- list()
trace.df <- data.frame("tau_amer" = NA, "tau_amer.calv" = NA, "tau_titia.Atl" = NA,"tau_titia" = NA)

# Loop through every model
for(run in model.names){ 
  rm(trace.name,trace,trace.raw,trace.means.raw,trace.melt,trace.mig,HPD.trace,tree.text)
  trace.name <- run
  print(trace.name)
  trace <- read.table(paste0(input.dir,
                             trace.name, ".trace"), header = T)
  # Remove burn in
  burn.in <- round(length(trace$Sample)*0.1, 0)
  trace <- trace[burn.in:length(trace$Sample),]
  
  # Plot tracer
  # plot(trace[seq(1,length(trace$Sample),length.out = 5000),grep(names(trace), pattern = "tau")])
  
  effect.size[[trace.name]] <- effectiveSize(trace[,-1])
  effect.size[[trace.name]]
  
  # Allow visulasaion of trace file in ggplot
  trace.melt <- reshape2::melt(trace[,c(1,grep(names(trace), pattern = "tau"))], id.vars="Sample")
  table(trace.melt$variable)
  trace.melt <- reshape2::melt(trace[seq(1,length(trace$Sample),length.out = 5000),], id.vars="Sample")
  
  #ggplot(trace.melt) +
  #  geom_line(aes(y = value, x = Sample)) +
  #  facet_wrap(~variable, scales = "free", ncol = 2)
  
  # Raw means
  trace.raw <- trace
  trace.means.raw <- colMeans(trace.raw[,-1], na.rm = T)
  
  # Migration rates
  head(trace)
  trace.mig <- data.frame(mig = names(trace.raw[,grep(names(trace), pattern = "m_")]))
  if(dim(trace.mig)[1]>=1){
    trace.mig$mig.out <- gsub(sapply(strsplit(trace.mig$mig,"\\.\\."),"[[",1), pattern = "m_", replacement = "")
    trace.mig$mig.in <- sapply(strsplit(trace.mig$mig,"\\.\\."),"[[",2)
    if(dim(trace.mig)[1]>1){
      trace.mig$mig.rate <- colMeans(trace.raw[,grep(names(trace.raw), pattern = "m_")], na.rm = T)/0.001
    } else if(dim(trace.mig)[1]==1){trace.mig$mig.rate <- mean(trace.raw[,grep(names(trace), pattern = "m_")], na.rm = T)/0.001}
    
  }
  trace.mig$mig.in.theta <- do.call("c",lapply(trace.mig$mig.in, function(x) trace.means.raw[which(names(trace.means.raw)==paste0("theta_",x))]))
  trace.mig$mig.in.theta <- ((trace.mig$mig.in.theta/1000)/(4*u))
  trace.mig$mig.rate.prop <- trace.mig$mig.rate*u
  trace.mig$mig.rate.num.ind.per.gen <- trace.mig$mig.in.theta*trace.mig$mig.rate.prop
  trace.mig$mig.rate.num.ind.per.gen.pct <- trace.mig$mig.rate.num.ind.per.gen*100
  
  
  # Convert output to mya and number of individuals (millions)
  trace[,grep(names(trace), pattern = "tau")] <- ((trace[,grep(names(trace), pattern = "tau")]/1000)/u)*10^-6
  trace[,grep(names(trace), pattern = "m_")] <- ((trace[,grep(names(trace), pattern = "m_")]/0.001)*u)
  trace[,grep(names(trace), pattern = "theta")] <- ((trace[,grep(names(trace), pattern = "theta")]/1000)/(4*u))*10^-6
  
  
  for(i in colnames(trace.df)[colnames(trace.df)%in%names(trace)]){
    trace.df[run,i] <- mean(trace[,i])
  }
  
  #ggplot(trace.melt[]) +
  # geom_violin(aes(x = value, y = variable),orientation = "y") +
  # facet_wrap(~variable, scales = "free")
  
  # calculate HPD for all tau traces
  HPD.trace <- HPDinterval(as.mcmc(trace[,grep(names(trace), pattern = "tau")]))
  HPD.trace.list[[trace.name]] <- matrix(c(row.names(HPD.trace)[1:2],
                                         c(paste(HPD.trace[1,1:2], collapse = " - "), paste(HPD.trace[2,1:2], collapse = " - "))),
                                         nrow=2, byrow = T)
  # Add to trace.df
  trace.df
  
  trace.means <- colMeans(trace[,-1], na.rm = T)
  names(trace.means)
  trace.means
  # Read in tree file
  if(grepl(run, pattern = "titia")){
    tree.text <- paste0("((titia.NAtl:",trace.means["tau_titia.Atl"],",titia.SAtl:",trace.means["tau_titia.Atl"],")titia.Atl:",trace.means["tau_titia"]-trace.means["tau_titia.Atl"]
                        ,",titia.Pac:",trace.means["tau_titia"],")titia:",trace.means["tau_titia"],";")
  }
  if(grepl(run, pattern = "americana")){
    tree.text <- paste0("((amer.N:",trace.means["tau_amer"],",amer.S:",trace.means["tau_amer"],")amer:",trace.means["tau_amer.calv"]-trace.means["tau_amer"]
                        ,",calv.1:",trace.means["tau_amer.calv"],")amer.calv:",trace.means["tau_amer.calv"],";")
  }
  
  tree <- ape::read.tree(text = tree.text)
  #plot(tree)
  # Merge HPD and tree
  tree.plot <- ggtree(tree) %<+% data.frame(node.labes = gsub(row.names(HPD.trace),pattern = "tau_", replacement= ""),HPD.trace)
  
  # Add in HPD in format that ggtree understands (list in a tibble)
  height_0.95_HPD_list <- list()
  for(i in 1:length(tree.plot$data$parent)){
    if(tree.plot$data$lower[i]>0&!is.na(tree.plot$data$lower[i])){
      height_0.95_HPD_list[[i]] <- c(tree.plot$data$lower[i],tree.plot$data$upper[i])
    }else height_0.95_HPD_list[[i]] <- c(tree.plot$data$x[i],tree.plot$data$x[i])
  }
  
  ## Add effective population size
  theta.means <- trace.means[grep(pattern = "theta",names(trace.means))]
  theta.means <- theta.means[match(paste0("theta_", tree.plot$data$label),names(theta.means))]
  
  cbind(paste0("theta_", tree.plot$data$label),theta.means)
  
  tree.plot$data$theta <- theta.means
  
  # Add to tree data
  tree.plot$data$height_0.95_HPD <- height_0.95_HPD_list
  tree.height <- max(ape::node.depth.edgelength(tree))
  # Get tree corrdinates for mig rates
  trace.mig$mig.out.x <- do.call("c",lapply(trace.mig$mig.out, function(x)  tree.plot$data$x[which(tree.plot$data$label==x)]-(tree.plot$data$branch.length[which(tree.plot$data$label==x)]/2)))
  trace.mig$mig.out.y <- do.call("c",lapply(trace.mig$mig.out, function(x)  tree.plot$data$y[which(tree.plot$data$label==x)]))
  trace.mig$mig.in.x <- do.call("c",lapply(trace.mig$mig.in, function(x)  tree.plot$data$x[which(tree.plot$data$label==x)]-(tree.plot$data$branch.length[which(tree.plot$data$label==x)]/2)))
  trace.mig$mig.in.y <- do.call("c",lapply(trace.mig$mig.in, function(x)  tree.plot$data$y[which(tree.plot$data$label==x)]))
  trace.mig$mig.out.x <-  tree.height-trace.mig$mig.out.x+rep(c(-0.1,0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.in.x <- tree.height-trace.mig$mig.in.x+rep(c(-0.1,0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.out.y <- trace.mig$mig.out.y+rep(c(0.1,-0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.in.y <- trace.mig$mig.in.y+rep(c(-0.1,0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.out.y.label <- trace.mig$mig.out.y-rep(c(-0.6,0.6),length.out = length(trace.mig$mig))
  trace.mig$mig.rate.perc <- (trace.mig$mig.rate*u)*100
  
  tree.plot$data$x.theta <- tree.height-tree.plot$data$x+(tree.plot$data$branch.length/2)
  tree.plot$data$y.theta <- tree.plot$data$y
  tree.plot$data$x.theta[tree.plot$data$label=="titia"|tree.plot$data$label=="amer.calv"] <- tree.height
  tree.plot$data$y.theta[tree.plot$data$label=="titia"|tree.plot$data$label=="amer.calv"] <- tree.plot$data$y[tree.plot$data$label=="titia"|tree.plot$data$label=="amer.calv"]-0.1
  
  trace.mig$mig.out.y.label
  # Plot tree
  
  p[[trace.name]] <- tree.plot +
    geom_tree(aes(linewidth = theta)) +
    #geom_rootedge(size = as.numeric(tree.plot$data$theta[4]*10),rootedge = 2) +
    #geom_segment(range = 'height_0.95_HPD', col = "deepskyblue", size = 4, alpha = .5) +
    geom_segment(aes(x = -lower, xend = -upper, y = , yend = y),col = "deepskyblue", linewidth = 10, alpha = .5,show.legend=F) +
    #geom_point(aes(x, y, size = theta)) +
    geom_tiplab(hjust = -0.2) +
    #coord_cartesian(clip = 'off') +
    # geom_nodelab(aes(x=x, label= paste(sprintf("%.1f", round(-x, digits = 1)))),size = 5, hjust= -0.4,vjust = 0.2,show.legend=F) +
    #geom_label(aes(x=x,y = y, label= paste0("Ne=",sprintf("%.1f", round(theta,1)))), size = 2.5, show.legend=F) +
    geom_label(aes(x=-x.theta,y = y.theta, label= paste0("Ne=",sprintf("%.1f", round(theta,1)))), size = 2.5, show.legend=F) +
    #geom_nodelab(aes(x=x, y=y, label=label, colour = "blue")) +
    #geom_treescale(fontsize=15, linesize=3, offset=2) +
    theme_tree2() +
    #coord_cartesian(clip = 'off') +
    theme_tree2(plot.margin = margin(1,1,1,1, "cm")) +
    theme(panel.grid.major=element_line(colour = "grey70"),
          panel.grid.major.y = element_blank()) +
    expand_limits(x = c(1,-10.2)) +
    xlab("Mya") +
    ggtitle(gsub(trace.name, pattern = "-run0", replacement = "")) +
    labs(size = "Effective population\nsize (millions)")
  
  if(dim(trace.mig)[1]>=1){
    p[[trace.name]] <- p[[trace.name]] + geom_segment(data = trace.mig, aes(x = -mig.out.x, y = mig.out.y, xend = -mig.in.x, yend = mig.in.y), arrow = arrow(length = unit(0.2,"cm"))) +
      geom_label(data = trace.mig, aes(x = -mig.out.x, y = mig.out.y.label, label = paste0(round(mig.rate.num.ind.per.gen, 3))), size = 2.5)
  }
  
  p[[trace.name]] <- revts(p[[trace.name]]) + scale_x_continuous(breaks = seq(-20, 0, 1),
                                                                 labels = abs(seq(20, 0, -1)))
  
  ggsave(plot = (p[[trace.name]]), filename = paste0(plot.dir,run,"-u-",
                                              gsub(u, replacement = "-", pattern ="\\."),".png"),
         width =12, height = 6)
  
}

# Merge all tau
trace.df <- trace.df[row.names(trace.df)!=1,]
trace.df$model <- row.names(trace.df)

trace.df$species <- sapply(strsplit(trace.df$model, "-"), "[", 5)
trace.df$mig.type <- sapply(strsplit(trace.df$model, "-"), "[", 6)
trace.df$beta <- sapply(strsplit(trace.df$model, "-"), "[", 4)
names(trace.df)

for(i in 1:length(trace.df$model)){
  HPD.temp <- HPD.trace.list[[trace.df$model[i]]]
  trace.df$HPD[i] <- paste0(paste(HPD.temp[1,],"=", HPD.temp[2,]), collapse = " : ")
}

trace.df
write.table(trace.df, paste0("4_Manuscript/data/G-Phocs/","G-Phocs_results_",draft_genome,"_", gsub(u, replacement = "-", pattern ="\\."),".txt"),
            row.names = F, quote=F)

effect.size

# Mean and range of divergence estiamtes across models
group_by(trace.df[trace.df$species=="titia",], mig.type) %>%
  summarise(mn.tau_titia = mean(tau_titia), min.tau_titia = min(tau_titia), max.tau_titia = max(tau_titia),
            mn.tau_titia.Atl = mean(tau_titia.Atl), min.tau_titia.Atl = min(tau_titia.Atl), max.tau_titia.Atl = max(tau_titia.Atl))

group_by(trace.df[trace.df$species=="americana",], mig.type) %>%
  summarise(mn.tau_amer.calv = mean(tau_amer.calv), min.tau_amer.calv = min(tau_amer.calv), max.tau_amer.calv = max(tau_amer.calv),
            mn.tau_amer = mean(tau_amer), min.tau_amer = min(tau_amer), max.tau_amer = max(tau_amer))

man.plot <- (p[[paste0("G-Phocs-a1-b20-titia-Iso-N", N_select)]]+ xlim(x = c(-8,1)) + theme(legend.position="none") + 
               p[[paste0("G-Phocs-a1-b20-americana-Iso-N", N_select)]]+ xlim(x = c(-8,1)) + theme(legend.position="none")) /
  (p[[paste0("G-Phocs-a1-b20-titia-migR-N", N_select)]]+ xlim(x = c(-8,1)) + theme(legend.position="none") + 
     p[[paste0("G-Phocs-a1-b20-americana-migR-N", N_select)]] + xlim(x = c(-8,1)) + theme(legend.position="none")) + 
  plot_annotation(title = expression(paste("Loci calling using ", italic("H. titia"), " draft genome")), theme = theme(plot.title = element_text(size = 25)),
tag_levels = 'a',tag_prefix = "(", tag_suffix = ")")
  

ggsave(plot = man.plot, filename = paste0(plot.dir,"Concat_plots_",draft_genome,"-u-",
                                          gsub(u, replacement = "-", pattern ="\\."),".png"),
       width =12, height = 10)


# # # # # # # # # # # # # # # # # # # # # #
#### ALL americana draft genome runs ####
# # # # # # # # # # # # # # # # # # # # # # # #

input.dir.top <- paste0("4_Manuscript/data/G-Phocs/")
draft_genome <- "americana_dg"
input.dir <- paste0(input.dir.top, draft_genome, "/")
model.names <- list.files(input.dir)
model.names <- model.names[grep(model.names, pattern = ".trace")]
model.names <- gsub(x = model.names, pattern = ".trace",replacement =  "")

trace.files <- paste0(model.names,".trace")
#trace.files <- trace.files[1]
#trace.name <- trace.files

# Broad use mutation rate used in other insect papers
u = 2.8*10^-9 
# Number of differences betweeen H. titia and I. elegans divide by divergence time 121 million years
# u = 2.2e-09
# Output folder for plots 
plot.dir <- "4_Manuscript/plots/G-Phocs/"

plot.dir <- paste0(plot.dir,draft_genome, "_",
                   gsub(u, replacement = "-", pattern ="\\."),"/")
dir.create(plot.dir)

# Blank lists for filling in data
p <- list()
effect.size <- list()
run.df <- list()
p.map <- list()
HPD.trace.list <- list()
trace.df <- data.frame("tau_amer" = NA, "tau_amer.calv" = NA, "tau_titia.Atl" = NA,"tau_titia" = NA)

# Loop through every model
for(run in model.names){ 
  rm(trace.name,trace,trace.raw,trace.means.raw,trace.melt,trace.mig,HPD.trace,tree.text)
  trace.name <- run
  print(trace.name)
  trace <- read.table(paste0(input.dir,
                             trace.name, ".trace"), header = T)
  # Remove burn in
  burn.in <- round(length(trace$Sample)*0.1, 0)
  trace <- trace[burn.in:length(trace$Sample),]
  
  # Plot tracer
  # plot(trace[seq(1,length(trace$Sample),length.out = 5000),grep(names(trace), pattern = "tau")])
  
  effect.size[[trace.name]] <- effectiveSize(trace[,-1])
  effect.size[[trace.name]]
  
  # Allow visulasaion of trace file in ggplot
  trace.melt <- reshape2::melt(trace[,c(1,grep(names(trace), pattern = "tau"))], id.vars="Sample")
  table(trace.melt$variable)
  trace.melt <- reshape2::melt(trace[seq(1,length(trace$Sample),length.out = 5000),], id.vars="Sample")
  
  #ggplot(trace.melt) +
  #  geom_line(aes(y = value, x = Sample)) +
  #  facet_wrap(~variable, scales = "free", ncol = 2)
  
  # Raw means
  trace.raw <- trace
  trace.means.raw <- colMeans(trace.raw[,-1], na.rm = T)
  
  # Migration rates
  head(trace)
  trace.mig <- data.frame(mig = names(trace.raw[,grep(names(trace), pattern = "m_")]))
  if(dim(trace.mig)[1]>=1){
    trace.mig$mig.out <- gsub(sapply(strsplit(trace.mig$mig,"\\.\\."),"[[",1), pattern = "m_", replacement = "")
    trace.mig$mig.in <- sapply(strsplit(trace.mig$mig,"\\.\\."),"[[",2)
    if(dim(trace.mig)[1]>1){
      trace.mig$mig.rate <- colMeans(trace.raw[,grep(names(trace.raw), pattern = "m_")], na.rm = T)/0.001
    } else if(dim(trace.mig)[1]==1){trace.mig$mig.rate <- mean(trace.raw[,grep(names(trace), pattern = "m_")], na.rm = T)/0.001}
    
  }
  trace.mig$mig.in.theta <- do.call("c",lapply(trace.mig$mig.in, function(x) trace.means.raw[which(names(trace.means.raw)==paste0("theta_",x))]))
  trace.mig$mig.in.theta <- ((trace.mig$mig.in.theta/1000)/(4*u))
  trace.mig$mig.rate.prop <- trace.mig$mig.rate*u
  trace.mig$mig.rate.num.ind.per.gen <- trace.mig$mig.in.theta*trace.mig$mig.rate.prop
  trace.mig$mig.rate.num.ind.per.gen.pct <- trace.mig$mig.rate.num.ind.per.gen*100
  
  
  # Convert output to mya and number of individuals (millions)
  trace[,grep(names(trace), pattern = "tau")] <- ((trace[,grep(names(trace), pattern = "tau")]/1000)/u)*10^-6
  trace[,grep(names(trace), pattern = "m_")] <- ((trace[,grep(names(trace), pattern = "m_")]/0.001)*u)
  trace[,grep(names(trace), pattern = "theta")] <- ((trace[,grep(names(trace), pattern = "theta")]/1000)/(4*u))*10^-6
  
  
  for(i in colnames(trace.df)[colnames(trace.df)%in%names(trace)]){
    trace.df[run,i] <- mean(trace[,i])
  }
  
  #ggplot(trace.melt[]) +
  # geom_violin(aes(x = value, y = variable),orientation = "y") +
  # facet_wrap(~variable, scales = "free")
  
  # calculate HPD for all tau traces
  HPD.trace <- HPDinterval(as.mcmc(trace[,grep(names(trace), pattern = "tau")]))
  HPD.trace.list[[trace.name]] <- matrix(c(row.names(HPD.trace)[1:2],
                                           c(paste(HPD.trace[1,1:2], collapse = " - "), paste(HPD.trace[2,1:2], collapse = " - "))),
                                         nrow=2, byrow = T)
  # Add to trace.df
  trace.df
  
  trace.means <- colMeans(trace[,-1], na.rm = T)
  names(trace.means)
  trace.means
  # Read in tree file
  if(grepl(run, pattern = "titia")){
    tree.text <- paste0("((titia.NAtl:",trace.means["tau_titia.Atl"],",titia.SAtl:",trace.means["tau_titia.Atl"],")titia.Atl:",trace.means["tau_titia"]-trace.means["tau_titia.Atl"]
                        ,",titia.Pac:",trace.means["tau_titia"],")titia:",trace.means["tau_titia"],";")
  }
  if(grepl(run, pattern = "americana")){
    tree.text <- paste0("((amer.N:",trace.means["tau_amer"],",amer.S:",trace.means["tau_amer"],")amer:",trace.means["tau_amer.calv"]-trace.means["tau_amer"]
                        ,",calv.1:",trace.means["tau_amer.calv"],")amer.calv:",trace.means["tau_amer.calv"],";")
  }
  
  tree <- ape::read.tree(text = tree.text)
  #plot(tree)
  # Merge HPD and tree
  tree.plot <- ggtree(tree) %<+% data.frame(node.labes = gsub(row.names(HPD.trace),pattern = "tau_", replacement= ""),HPD.trace)
  
  # Add in HPD in format that ggtree understands (list in a tibble)
  height_0.95_HPD_list <- list()
  for(i in 1:length(tree.plot$data$parent)){
    if(tree.plot$data$lower[i]>0&!is.na(tree.plot$data$lower[i])){
      height_0.95_HPD_list[[i]] <- c(tree.plot$data$lower[i],tree.plot$data$upper[i])
    }else height_0.95_HPD_list[[i]] <- c(tree.plot$data$x[i],tree.plot$data$x[i])
  }
  
  ## Add effective population size
  theta.means <- trace.means[grep(pattern = "theta",names(trace.means))]
  theta.means <- theta.means[match(paste0("theta_", tree.plot$data$label),names(theta.means))]
  
  cbind(paste0("theta_", tree.plot$data$label),theta.means)
  
  tree.plot$data$theta <- theta.means
  
  # Add to tree data
  tree.plot$data$height_0.95_HPD <- height_0.95_HPD_list
  tree.height <- max(ape::node.depth.edgelength(tree))
  # Get tree corrdinates for mig rates
  trace.mig$mig.out.x <- do.call("c",lapply(trace.mig$mig.out, function(x)  tree.plot$data$x[which(tree.plot$data$label==x)]-(tree.plot$data$branch.length[which(tree.plot$data$label==x)]/2)))
  trace.mig$mig.out.y <- do.call("c",lapply(trace.mig$mig.out, function(x)  tree.plot$data$y[which(tree.plot$data$label==x)]))
  trace.mig$mig.in.x <- do.call("c",lapply(trace.mig$mig.in, function(x)  tree.plot$data$x[which(tree.plot$data$label==x)]-(tree.plot$data$branch.length[which(tree.plot$data$label==x)]/2)))
  trace.mig$mig.in.y <- do.call("c",lapply(trace.mig$mig.in, function(x)  tree.plot$data$y[which(tree.plot$data$label==x)]))
  trace.mig$mig.out.x <-  tree.height-trace.mig$mig.out.x+rep(c(-0.1,0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.in.x <- tree.height-trace.mig$mig.in.x+rep(c(-0.1,0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.out.y <- trace.mig$mig.out.y+rep(c(0.1,-0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.in.y <- trace.mig$mig.in.y+rep(c(-0.1,0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.out.y.label <- trace.mig$mig.out.y-rep(c(-0.6,0.6),length.out = length(trace.mig$mig))
  trace.mig$mig.rate.perc <- (trace.mig$mig.rate*u)*100
  
  tree.plot$data$x.theta <- tree.height-tree.plot$data$x+(tree.plot$data$branch.length/2)
  tree.plot$data$y.theta <- tree.plot$data$y
  tree.plot$data$x.theta[tree.plot$data$label=="titia"|tree.plot$data$label=="amer.calv"] <- tree.height
  tree.plot$data$y.theta[tree.plot$data$label=="titia"|tree.plot$data$label=="amer.calv"] <- tree.plot$data$y[tree.plot$data$label=="titia"|tree.plot$data$label=="amer.calv"]-0.1
  
  trace.mig$mig.out.y.label
  # Plot tree
  
  p[[trace.name]] <- tree.plot +
    geom_tree(aes(linewidth = theta)) +
    #geom_rootedge(size = as.numeric(tree.plot$data$theta[4]*10),rootedge = 2) +
    #geom_segment(range = 'height_0.95_HPD', col = "deepskyblue", size = 4, alpha = .5) +
    geom_segment(aes(x = -lower, xend = -upper, y = , yend = y),col = "deepskyblue", linewidth = 10, alpha = .5,show.legend=F) +
    #geom_point(aes(x, y, size = theta)) +
    geom_tiplab(hjust = -0.2) +
    #coord_cartesian(clip = 'off') +
    # geom_nodelab(aes(x=x, label= paste(sprintf("%.1f", round(-x, digits = 1)))),size = 5, hjust= -0.4,vjust = 0.2,show.legend=F) +
    #geom_label(aes(x=x,y = y, label= paste0("Ne=",sprintf("%.1f", round(theta,1)))), size = 2.5, show.legend=F) +
    geom_label(aes(x=-x.theta,y = y.theta, label= paste0("Ne=",sprintf("%.1f", round(theta,1)))), size = 2.5, show.legend=F) +
    #geom_nodelab(aes(x=x, y=y, label=label, colour = "blue")) +
    #geom_treescale(fontsize=15, linesize=3, offset=2) +
    theme_tree2() +
    #coord_cartesian(clip = 'off') +
    theme_tree2(plot.margin = margin(1,1,1,1, "cm")) +
    theme(panel.grid.major=element_line(colour = "grey70"),
          panel.grid.major.y = element_blank()) +
    expand_limits(x = c(1,-10.2)) +
    xlab("Mya") +
    ggtitle(gsub(trace.name, pattern = "-run0", replacement = "")) +
    labs(size = "Effective population\nsize (millions)")
  
  if(dim(trace.mig)[1]>=1){
    p[[trace.name]] <- p[[trace.name]] + geom_segment(data = trace.mig, aes(x = -mig.out.x, y = mig.out.y, xend = -mig.in.x, yend = mig.in.y), arrow = arrow(length = unit(0.2,"cm"))) +
      geom_label(data = trace.mig, aes(x = -mig.out.x, y = mig.out.y.label, label = paste0(round(mig.rate.num.ind.per.gen, 3))), size = 2.5)
  }
  
  p[[trace.name]] <- revts(p[[trace.name]]) + scale_x_continuous(breaks = seq(-20, 0, 1),
                                                                 labels = abs(seq(20, 0, -1)))
  
  ggsave(plot = (p[[trace.name]]), filename = paste0(plot.dir,run,"-u-",
                                                     gsub(u, replacement = "-", pattern ="\\."),".png"),
         width =12, height = 6)
  
}

# Merge all tau
trace.df <- trace.df[row.names(trace.df)!=1,]
trace.df$model <- row.names(trace.df)

trace.df$species <- sapply(strsplit(trace.df$model, "-"), "[", 5)
trace.df$mig.type <- sapply(strsplit(trace.df$model, "-"), "[", 6)
trace.df$beta <- sapply(strsplit(trace.df$model, "-"), "[", 4)
names(trace.df)

for(i in 1:length(trace.df$model)){
  HPD.temp <- HPD.trace.list[[trace.df$model[i]]]
  trace.df$HPD[i] <- paste0(paste(HPD.temp[1,],"=", HPD.temp[2,]), collapse = " : ")
}

trace.df
write.table(trace.df, paste0("4_Manuscript/data/G-Phocs/","G-Phocs_results_",draft_genome,"_", gsub(u, replacement = "-", pattern ="\\."),".txt"),
            row.names = F, quote=F)

effect.size

# Mean and range of divergence estiamtes across models
group_by(trace.df[trace.df$species=="titia",], mig.type) %>%
  summarise(mn.tau_titia = mean(tau_titia), min.tau_titia = min(tau_titia), max.tau_titia = max(tau_titia),
            mn.tau_titia.Atl = mean(tau_titia.Atl), min.tau_titia.Atl = min(tau_titia.Atl), max.tau_titia.Atl = max(tau_titia.Atl))

group_by(trace.df[trace.df$species=="americana",], mig.type) %>%
  summarise(mn.tau_amer.calv = mean(tau_amer.calv), min.tau_amer.calv = min(tau_amer.calv), max.tau_amer.calv = max(tau_amer.calv),
            mn.tau_amer = mean(tau_amer), min.tau_amer = min(tau_amer), max.tau_amer = max(tau_amer))

man.plot <- (p[[paste0("G-Phocs-a1-b20-titia-Iso-N", N_select)]]+ xlim(x = c(-8,1)) + theme(legend.position="none") + 
               p[[paste0("G-Phocs-a1-b20-americana-Iso-N", N_select)]]+ xlim(x = c(-8,1)) + theme(legend.position="none")) /
  (p[[paste0("G-Phocs-a1-b20-titia-migR-N", N_select)]]+ xlim(x = c(-8,1)) + theme(legend.position="none") + 
     p[[paste0("G-Phocs-a1-b20-americana-migR-N", N_select)]] + xlim(x = c(-8,1)) + theme(legend.position="none")) + 
  plot_annotation(title = expression(paste("Loci calling using ", italic("H. americana"), " draft genome")), theme = theme(plot.title = element_text(size = 25))
    ,tag_levels = 'a',tag_prefix = "(", tag_suffix = ")")

ggsave(plot = man.plot, filename = paste0(plot.dir,"Concat_plots_",draft_genome,"-u-",
                                          gsub(u, replacement = "-", pattern ="\\."),".png"),
       width =12, height = 10)


# # # # # # # # # # # # # # # # # # # # # # # #
#### ALL denovo runs ####
# # # # # # # # # # # # # # # # # # # # # # # #

input.dir.top <- paste0("4_Manuscript/data/G-Phocs/")
draft_genome <- "denovo"
input.dir <- paste0(input.dir.top, draft_genome, "/")
model.names <- list.files(input.dir)
model.names <- model.names[grep(model.names, pattern = ".trace")]
model.names <- gsub(x = model.names, pattern = ".trace",replacement =  "")

trace.files <- paste0(model.names,".trace")
#trace.files <- trace.files[1]
#trace.name <- trace.files

# Broad use mutation rate used in other insect papers
u = 2.8*10^-9 
# Number of differences betweeen H. titia and I. elegans divide by divergence time 121 million years
# u = 2.2e-09
# Output folder for plots 
plot.dir <- "4_Manuscript/plots/G-Phocs/"

plot.dir <- paste0(plot.dir,draft_genome, "_",
                   gsub(u, replacement = "-", pattern ="\\."),"/")
dir.create(plot.dir)

# Blank lists for filling in data
p <- list()
effect.size <- list()
run.df <- list()
p.map <- list()
HPD.trace.list <- list()
trace.df <- data.frame("tau_amer" = NA, "tau_amer.calv" = NA, "tau_titia.Atl" = NA,"tau_titia" = NA)

# Loop through every model
for(run in model.names){ 
  rm(trace.name,trace,trace.raw,trace.means.raw,trace.melt,trace.mig,HPD.trace,tree.text)
  trace.name <- run
  print(trace.name)
  trace <- read.table(paste0(input.dir,
                             trace.name, ".trace"), header = T)
  # Remove burn in
  burn.in <- round(length(trace$Sample)*0.1, 0)
  trace <- trace[burn.in:length(trace$Sample),]
  
  # Plot tracer
  # plot(trace[seq(1,length(trace$Sample),length.out = 5000),grep(names(trace), pattern = "tau")])
  
  effect.size[[trace.name]] <- effectiveSize(trace[,-1])
  effect.size[[trace.name]]
  
  # Allow visulasaion of trace file in ggplot
  trace.melt <- reshape2::melt(trace[,c(1,grep(names(trace), pattern = "tau"))], id.vars="Sample")
  table(trace.melt$variable)
  trace.melt <- reshape2::melt(trace[seq(1,length(trace$Sample),length.out = 5000),], id.vars="Sample")
  
  #ggplot(trace.melt) +
  #  geom_line(aes(y = value, x = Sample)) +
  #  facet_wrap(~variable, scales = "free", ncol = 2)
  
  # Raw means
  trace.raw <- trace
  trace.means.raw <- colMeans(trace.raw[,-1], na.rm = T)
  
  # Migration rates
  head(trace)
  trace.mig <- data.frame(mig = names(trace.raw[,grep(names(trace), pattern = "m_")]))
  if(dim(trace.mig)[1]>=1){
    trace.mig$mig.out <- gsub(sapply(strsplit(trace.mig$mig,"\\.\\."),"[[",1), pattern = "m_", replacement = "")
    trace.mig$mig.in <- sapply(strsplit(trace.mig$mig,"\\.\\."),"[[",2)
    if(dim(trace.mig)[1]>1){
      trace.mig$mig.rate <- colMeans(trace.raw[,grep(names(trace.raw), pattern = "m_")], na.rm = T)/0.001
    } else if(dim(trace.mig)[1]==1){trace.mig$mig.rate <- mean(trace.raw[,grep(names(trace), pattern = "m_")], na.rm = T)/0.001}
    
  }
  trace.mig$mig.in.theta <- do.call("c",lapply(trace.mig$mig.in, function(x) trace.means.raw[which(names(trace.means.raw)==paste0("theta_",x))]))
  trace.mig$mig.in.theta <- ((trace.mig$mig.in.theta/1000)/(4*u))
  trace.mig$mig.rate.prop <- trace.mig$mig.rate*u
  trace.mig$mig.rate.num.ind.per.gen <- trace.mig$mig.in.theta*trace.mig$mig.rate.prop
  trace.mig$mig.rate.num.ind.per.gen.pct <- trace.mig$mig.rate.num.ind.per.gen*100
  
  
  # Convert output to mya and number of individuals (millions)
  trace[,grep(names(trace), pattern = "tau")] <- ((trace[,grep(names(trace), pattern = "tau")]/1000)/u)*10^-6
  trace[,grep(names(trace), pattern = "m_")] <- ((trace[,grep(names(trace), pattern = "m_")]/0.001)*u)
  trace[,grep(names(trace), pattern = "theta")] <- ((trace[,grep(names(trace), pattern = "theta")]/1000)/(4*u))*10^-6
  
  
  for(i in colnames(trace.df)[colnames(trace.df)%in%names(trace)]){
    trace.df[run,i] <- mean(trace[,i])
  }
  
  #ggplot(trace.melt[]) +
  # geom_violin(aes(x = value, y = variable),orientation = "y") +
  # facet_wrap(~variable, scales = "free")
  
  # calculate HPD for all tau traces
  HPD.trace <- HPDinterval(as.mcmc(trace[,grep(names(trace), pattern = "tau")]))
  HPD.trace.list[[trace.name]] <- matrix(c(row.names(HPD.trace)[1:2],
                                           c(paste(HPD.trace[1,1:2], collapse = " - "), paste(HPD.trace[2,1:2], collapse = " - "))),
                                         nrow=2, byrow = T)
  # Add to trace.df
  trace.df
  
  trace.means <- colMeans(trace[,-1], na.rm = T)
  names(trace.means)
  trace.means
  # Read in tree file
  if(grepl(run, pattern = "titia")){
    tree.text <- paste0("((titia.NAtl:",trace.means["tau_titia.Atl"],",titia.SAtl:",trace.means["tau_titia.Atl"],")titia.Atl:",trace.means["tau_titia"]-trace.means["tau_titia.Atl"]
                        ,",titia.Pac:",trace.means["tau_titia"],")titia:",trace.means["tau_titia"],";")
  }
  if(grepl(run, pattern = "americana")){
    tree.text <- paste0("((amer.N:",trace.means["tau_amer"],",amer.S:",trace.means["tau_amer"],")amer:",trace.means["tau_amer.calv"]-trace.means["tau_amer"]
                        ,",calv.1:",trace.means["tau_amer.calv"],")amer.calv:",trace.means["tau_amer.calv"],";")
  }
  
  tree <- ape::read.tree(text = tree.text)
  #plot(tree)
  # Merge HPD and tree
  tree.plot <- ggtree(tree) %<+% data.frame(node.labes = gsub(row.names(HPD.trace),pattern = "tau_", replacement= ""),HPD.trace)
  
  # Add in HPD in format that ggtree understands (list in a tibble)
  height_0.95_HPD_list <- list()
  for(i in 1:length(tree.plot$data$parent)){
    if(tree.plot$data$lower[i]>0&!is.na(tree.plot$data$lower[i])){
      height_0.95_HPD_list[[i]] <- c(tree.plot$data$lower[i],tree.plot$data$upper[i])
    }else height_0.95_HPD_list[[i]] <- c(tree.plot$data$x[i],tree.plot$data$x[i])
  }
  
  ## Add effective population size
  theta.means <- trace.means[grep(pattern = "theta",names(trace.means))]
  theta.means <- theta.means[match(paste0("theta_", tree.plot$data$label),names(theta.means))]
  
  cbind(paste0("theta_", tree.plot$data$label),theta.means)
  
  tree.plot$data$theta <- theta.means
  
  # Add to tree data
  tree.plot$data$height_0.95_HPD <- height_0.95_HPD_list
  tree.height <- max(ape::node.depth.edgelength(tree))
  # Get tree corrdinates for mig rates
  trace.mig$mig.out.x <- do.call("c",lapply(trace.mig$mig.out, function(x)  tree.plot$data$x[which(tree.plot$data$label==x)]-(tree.plot$data$branch.length[which(tree.plot$data$label==x)]/2)))
  trace.mig$mig.out.y <- do.call("c",lapply(trace.mig$mig.out, function(x)  tree.plot$data$y[which(tree.plot$data$label==x)]))
  trace.mig$mig.in.x <- do.call("c",lapply(trace.mig$mig.in, function(x)  tree.plot$data$x[which(tree.plot$data$label==x)]-(tree.plot$data$branch.length[which(tree.plot$data$label==x)]/2)))
  trace.mig$mig.in.y <- do.call("c",lapply(trace.mig$mig.in, function(x)  tree.plot$data$y[which(tree.plot$data$label==x)]))
  trace.mig$mig.out.x <-  tree.height-trace.mig$mig.out.x+rep(c(-0.1,0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.in.x <- tree.height-trace.mig$mig.in.x+rep(c(-0.1,0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.out.y <- trace.mig$mig.out.y+rep(c(0.1,-0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.in.y <- trace.mig$mig.in.y+rep(c(-0.1,0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.out.y.label <- trace.mig$mig.out.y-rep(c(-0.6,0.6),length.out = length(trace.mig$mig))
  trace.mig$mig.rate.perc <- (trace.mig$mig.rate*u)*100
  
  tree.plot$data$x.theta <- tree.height-tree.plot$data$x+(tree.plot$data$branch.length/2)
  tree.plot$data$y.theta <- tree.plot$data$y
  tree.plot$data$x.theta[tree.plot$data$label=="titia"|tree.plot$data$label=="amer.calv"] <- tree.height
  tree.plot$data$y.theta[tree.plot$data$label=="titia"|tree.plot$data$label=="amer.calv"] <- tree.plot$data$y[tree.plot$data$label=="titia"|tree.plot$data$label=="amer.calv"]-0.1
  
  trace.mig$mig.out.y.label
  # Plot tree
  
  p[[trace.name]] <- tree.plot +
    geom_tree(aes(linewidth = theta)) +
    #geom_rootedge(size = as.numeric(tree.plot$data$theta[4]*10),rootedge = 2) +
    #geom_segment(range = 'height_0.95_HPD', col = "deepskyblue", size = 4, alpha = .5) +
    geom_segment(aes(x = -lower, xend = -upper, y = , yend = y),col = "deepskyblue", linewidth = 10, alpha = .5,show.legend=F) +
    #geom_point(aes(x, y, size = theta)) +
    geom_tiplab(hjust = -0.2) +
    #coord_cartesian(clip = 'off') +
    # geom_nodelab(aes(x=x, label= paste(sprintf("%.1f", round(-x, digits = 1)))),size = 5, hjust= -0.4,vjust = 0.2,show.legend=F) +
    #geom_label(aes(x=x,y = y, label= paste0("Ne=",sprintf("%.1f", round(theta,1)))), size = 2.5, show.legend=F) +
    geom_label(aes(x=-x.theta,y = y.theta, label= paste0("Ne=",sprintf("%.1f", round(theta,1)))), size = 2.5, show.legend=F) +
    #geom_nodelab(aes(x=x, y=y, label=label, colour = "blue")) +
    #geom_treescale(fontsize=15, linesize=3, offset=2) +
    theme_tree2() +
    #coord_cartesian(clip = 'off') +
    theme_tree2(plot.margin = margin(1,1,1,1, "cm")) +
    theme(panel.grid.major=element_line(colour = "grey70"),
          panel.grid.major.y = element_blank()) +
    expand_limits(x = c(1,-10.2)) +
    xlab("Mya") +
    ggtitle(gsub(trace.name, pattern = "-run0", replacement = "")) +
    labs(size = "Effective population\nsize (millions)")
  
  if(dim(trace.mig)[1]>=1){
    p[[trace.name]] <- p[[trace.name]] + geom_segment(data = trace.mig, aes(x = -mig.out.x, y = mig.out.y, xend = -mig.in.x, yend = mig.in.y), arrow = arrow(length = unit(0.2,"cm"))) +
      geom_label(data = trace.mig, aes(x = -mig.out.x, y = mig.out.y.label, label = paste0(round(mig.rate.num.ind.per.gen, 3))), size = 2.5)
  }
  
  p[[trace.name]] <- revts(p[[trace.name]]) + scale_x_continuous(breaks = seq(-20, 0, 1),
                                                                 labels = abs(seq(20, 0, -1)))
  
  ggsave(plot = (p[[trace.name]]), filename = paste0(plot.dir,run,"-u-",
                                                     gsub(u, replacement = "-", pattern ="\\."),".png"),
         width =12, height = 6)
  
}

# Merge all tau
trace.df <- trace.df[row.names(trace.df)!=1,]
trace.df$model <- row.names(trace.df)

trace.df$species <- sapply(strsplit(trace.df$model, "-"), "[", 5)
trace.df$mig.type <- sapply(strsplit(trace.df$model, "-"), "[", 6)
trace.df$beta <- sapply(strsplit(trace.df$model, "-"), "[", 4)
names(trace.df)

for(i in 1:length(trace.df$model)){
  HPD.temp <- HPD.trace.list[[trace.df$model[i]]]
  trace.df$HPD[i] <- paste0(paste(HPD.temp[1,],"=", HPD.temp[2,]), collapse = " : ")
}

trace.df
write.table(trace.df, paste0("4_Manuscript/data/G-Phocs/","G-Phocs_results_",draft_genome,"_", gsub(u, replacement = "-", pattern ="\\."),".txt"),
            row.names = F, quote=F)

effect.size
names(p)
# Mean and range of divergence estiamtes across models
group_by(trace.df[trace.df$species=="titia",], mig.type) %>%
  summarise(mn.tau_titia = mean(tau_titia), min.tau_titia = min(tau_titia), max.tau_titia = max(tau_titia),
            mn.tau_titia.Atl = mean(tau_titia.Atl), min.tau_titia.Atl = min(tau_titia.Atl), max.tau_titia.Atl = max(tau_titia.Atl))

group_by(trace.df[trace.df$species=="americana",], mig.type) %>%
  summarise(mn.tau_amer.calv = mean(tau_amer.calv), min.tau_amer.calv = min(tau_amer.calv), max.tau_amer.calv = max(tau_amer.calv),
            mn.tau_amer = mean(tau_amer), min.tau_amer = min(tau_amer), max.tau_amer = max(tau_amer))

all_plot <- (p[[1]] + ggtitle(expression(paste(italic("H. titia"), " - No Migration")))) + 
  (p[[2]] + ggtitle(expression(paste(italic("H. americana/calverti"), " - No Migration")))) + 
  (p[[3]] + ggtitle(expression(paste(italic("H. titia"), " - Migration")))) +
  (p[[4]] + ggtitle(expression(paste(italic("H. americana/calverti"), " - Migration")))) +
  plot_annotation(title = "Heterospecific Loci calling", theme = theme(plot.title = element_text(size = 25)),
                  tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(ncol = 2)

man.plot <- ((p[[paste0("G-Phocs-a1-b20-titia-Iso-N3")]] + ggtitle(expression(paste(italic("H. titia"), " - No Migration")))) 
                      + xlim(x = c(-8,1)) + theme(legend.position="right") + 
             (p[[paste0("G-Phocs-a1-b20-americana-Iso-N3")]] + ggtitle(expression(paste(italic("H. americana/calverti"), " - No Migration"))))
                      + xlim(x = c(-8,1)) + theme(legend.position="right")) /
             ((p[[paste0("G-Phocs-a1-b20-titia-migR-N3")]] + ggtitle(expression(paste(italic("H. titia"), " - Migration"))))
                      + xlim(x = c(-8,1)) + theme(legend.position="right") + 
               (p[[paste0("G-Phocs-a1-b20-americana-migR-N3")]] + ggtitle(expression(paste(italic("H. americana/calverti"), " - Migration"))))
                      + xlim(x = c(-8,1)) + theme(legend.position="right")) + 
  plot_annotation(title = expression(paste(italic("denovo")," Loci calling")), theme = theme(plot.title = element_text(size = 25))
    ,tag_levels = 'a',tag_prefix = "(", tag_suffix = ")")

ggsave(plot = man.plot, filename = paste0(plot.dir,"Concat_plots_",draft_genome,"-u-",
                                          gsub(u, replacement = "-", pattern ="\\."),".png"),
       width =12, height = 10)




# # # # # # # # # # # # # # # # # # # # # # # #
#### Conspecific specific ####
# # # # # # # # # # # # # # # # # # # # # # # #

### Produce plots that use conspecific specific SNPs for g-phocs
input.dir.top <- "4_Manuscript/data/G-Phocs/"
model.names <- paste0(c("titia_dg/G-Phocs-a1-b20-titia-Iso-N","americana_dg/G-Phocs-a1-b20-americana-Iso-N",
                        "titia_dg/G-Phocs-a1-b20-titia-migR-N","americana_dg/G-Phocs-a1-b20-americana-migR-N"),
                      N_select)
u = 2.8*10^-9 
plot.dir <- "4_Manuscript/plots/G-Phocs/"
run <- "americana_dg/G-Phocs-a1-b20-titia-Iso-N3"
# Blank lists for filling in data
p <- list()
effect.size <- list()
run.df <- list()
p.map <- list()
HPD.trace.list <- list()
trace.df <- data.frame("tau_amer" = NA, "tau_amer.calv" = NA, "tau_titia.Atl" = NA,"tau_titia" = NA)
# Loop through every model
for(run in model.names){ 
  rm(trace.name,trace,trace.raw,trace.means.raw,trace.melt,trace.mig,HPD.trace,tree.text)
  trace.name <- run
  run <- stringr::str_split_i(run, pattern = "/", 2)
  print(trace.name)
  trace <- read.table(file = paste0(input.dir.top,
                                    trace.name, ".trace"), header = T)
  # Remove burn in
  burn.in <- round(length(trace$Sample)*0.1, 0)
  trace <- trace[burn.in:length(trace$Sample),]
  
  # Plot tracer
  # plot(trace[seq(1,length(trace$Sample),length.out = 5000),grep(names(trace), pattern = "tau")])
  
  effect.size[[trace.name]] <- effectiveSize(trace[,-1])
  effect.size[[trace.name]]
  
  # Allow visulasaion of trace file in ggplot
  trace.melt <- reshape2::melt(trace[,c(1,grep(names(trace), pattern = "tau"))], id.vars="Sample")
  table(trace.melt$variable)
  trace.melt <- reshape2::melt(trace[seq(1,length(trace$Sample),length.out = 5000),], id.vars="Sample")
  
  #ggplot(trace.melt) +
  #  geom_line(aes(y = value, x = Sample)) +
  #  facet_wrap(~variable, scales = "free", ncol = 2)
  
  # Raw means
  trace.raw <- trace
  trace.means.raw <- colMeans(trace.raw[,-1], na.rm = T)
  
  # Migration rates
  head(trace)
  trace.mig <- data.frame(mig = names(trace.raw[,grep(names(trace), pattern = "m_")]))
  if(dim(trace.mig)[1]>=1){
    trace.mig$mig.out <- gsub(sapply(strsplit(trace.mig$mig,"\\.\\."),"[[",1), pattern = "m_", replacement = "")
    trace.mig$mig.in <- sapply(strsplit(trace.mig$mig,"\\.\\."),"[[",2)
    if(dim(trace.mig)[1]>1){
      trace.mig$mig.rate <- colMeans(trace.raw[,grep(names(trace.raw), pattern = "m_")], na.rm = T)/0.001
    } else if(dim(trace.mig)[1]==1){trace.mig$mig.rate <- mean(trace.raw[,grep(names(trace), pattern = "m_")], na.rm = T)/0.001}
    
  }
  trace.mig$mig.in.theta <- do.call("c",lapply(trace.mig$mig.in, function(x) trace.means.raw[which(names(trace.means.raw)==paste0("theta_",x))]))
  trace.mig$mig.in.theta <- ((trace.mig$mig.in.theta/1000)/(4*u))
  trace.mig$mig.rate.prop <- trace.mig$mig.rate*u
  trace.mig$mig.rate.num.ind.per.gen <- trace.mig$mig.in.theta*trace.mig$mig.rate.prop
  trace.mig$mig.rate.num.ind.per.gen.pct <- trace.mig$mig.rate.num.ind.per.gen*100
  
  
  # Convert output to mya and number of individuals (millions)
  trace[,grep(names(trace), pattern = "tau")] <- ((trace[,grep(names(trace), pattern = "tau")]/1000)/u)*10^-6
  trace[,grep(names(trace), pattern = "m_")] <- ((trace[,grep(names(trace), pattern = "m_")]/0.001)*u)
  trace[,grep(names(trace), pattern = "theta")] <- ((trace[,grep(names(trace), pattern = "theta")]/1000)/(4*u))*10^-6
  
  
  for(i in colnames(trace.df)[colnames(trace.df)%in%names(trace)]){
    trace.df[run,i] <- mean(trace[,i])
  }
  
  #ggplot(trace.melt[]) +
  # geom_violin(aes(x = value, y = variable),orientation = "y") +
  # facet_wrap(~variable, scales = "free")
  
  # calculate HPD for all tau traces
  HPD.trace <- HPDinterval(as.mcmc(trace[,grep(names(trace), pattern = "tau")]))
  HPD.trace.list[[trace.name]] <- matrix(c(row.names(HPD.trace)[1:2],
                                           c(paste(HPD.trace[1,1:2], collapse = " - "), paste(HPD.trace[2,1:2], collapse = " - "))),
                                         nrow=2, byrow = T)
  # Add to trace.df
  trace.df
  
  trace.means <- colMeans(trace[,-1], na.rm = T)
  names(trace.means)
  trace.means
  # Read in tree file
  if(grepl(run, pattern = "titia")){
    tree.text <- paste0("((titia.NAtl:",trace.means["tau_titia.Atl"],",titia.SAtl:",trace.means["tau_titia.Atl"],")titia.Atl:",trace.means["tau_titia"]-trace.means["tau_titia.Atl"]
                        ,",titia.Pac:",trace.means["tau_titia"],")titia:",trace.means["tau_titia"],";")
  }
  if(grepl(run, pattern = "americana")){
    tree.text <- paste0("((amer.N:",trace.means["tau_amer"],",amer.S:",trace.means["tau_amer"],")amer:",trace.means["tau_amer.calv"]-trace.means["tau_amer"]
                        ,",calv.1:",trace.means["tau_amer.calv"],")amer.calv:",trace.means["tau_amer.calv"],";")
  }
  
  tree <- ape::read.tree(text = tree.text)
  #plot(tree)
  # Merge HPD and tree
  tree.plot <- ggtree(tree) %<+% data.frame(node.labes = gsub(row.names(HPD.trace),pattern = "tau_", replacement= ""),HPD.trace)
  
  # Add in HPD in format that ggtree understands (list in a tibble)
  height_0.95_HPD_list <- list()
  for(i in 1:length(tree.plot$data$parent)){
    if(tree.plot$data$lower[i]>0&!is.na(tree.plot$data$lower[i])){
      height_0.95_HPD_list[[i]] <- c(tree.plot$data$lower[i],tree.plot$data$upper[i])
    }else height_0.95_HPD_list[[i]] <- c(tree.plot$data$x[i],tree.plot$data$x[i])
  }
  
  ## Add effective population size
  theta.means <- trace.means[grep(pattern = "theta",names(trace.means))]
  theta.means <- theta.means[match(paste0("theta_", tree.plot$data$label),names(theta.means))]
  
  cbind(paste0("theta_", tree.plot$data$label),theta.means)
  
  tree.plot$data$theta <- theta.means
  
  # Add to tree data
  tree.plot$data$height_0.95_HPD <- height_0.95_HPD_list
  tree.height <- max(ape::node.depth.edgelength(tree))
  # Get tree corrdinates for mig rates
  trace.mig$mig.out.x <- do.call("c",lapply(trace.mig$mig.out, function(x)  tree.plot$data$x[which(tree.plot$data$label==x)]-(tree.plot$data$branch.length[which(tree.plot$data$label==x)]/2)))
  trace.mig$mig.out.y <- do.call("c",lapply(trace.mig$mig.out, function(x)  tree.plot$data$y[which(tree.plot$data$label==x)]))
  trace.mig$mig.in.x <- do.call("c",lapply(trace.mig$mig.in, function(x)  tree.plot$data$x[which(tree.plot$data$label==x)]-(tree.plot$data$branch.length[which(tree.plot$data$label==x)]/2)))
  trace.mig$mig.in.y <- do.call("c",lapply(trace.mig$mig.in, function(x)  tree.plot$data$y[which(tree.plot$data$label==x)]))
  trace.mig$mig.out.x <-  tree.height-trace.mig$mig.out.x+rep(c(-0.1,0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.in.x <- tree.height-trace.mig$mig.in.x+rep(c(-0.1,0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.out.y <- trace.mig$mig.out.y+rep(c(0.1,-0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.in.y <- trace.mig$mig.in.y+rep(c(-0.1,0.1),length.out = length(trace.mig$mig))
  trace.mig$mig.out.y.label <- trace.mig$mig.out.y-rep(c(-0.6,0.6),length.out = length(trace.mig$mig))
  trace.mig$mig.rate.perc <- (trace.mig$mig.rate*u)*100
  
  tree.plot$data$x.theta <- tree.height-tree.plot$data$x+(tree.plot$data$branch.length/2)
  tree.plot$data$y.theta <- tree.plot$data$y
  tree.plot$data$x.theta[tree.plot$data$label=="titia"|tree.plot$data$label=="amer.calv"] <- tree.height
  tree.plot$data$y.theta[tree.plot$data$label=="titia"|tree.plot$data$label=="amer.calv"] <- tree.plot$data$y[tree.plot$data$label=="titia"|tree.plot$data$label=="amer.calv"]-0.1
  
  trace.mig$mig.out.y.label
  # Plot tree
  
  p[[trace.name]] <- tree.plot +
    geom_tree(aes(linewidth = theta)) +
    #geom_rootedge(size = as.numeric(tree.plot$data$theta[4]*10),rootedge = 2) +
    #geom_segment(range = 'height_0.95_HPD', col = "deepskyblue", size = 4, alpha = .5) +
    geom_segment(aes(x = -lower, xend = -upper, y = , yend = y),col = "deepskyblue", linewidth = 10, alpha = .5,show.legend=F) +
    #geom_point(aes(x, y, size = theta)) +
    geom_tiplab(hjust = -0.2) +
    #coord_cartesian(clip = 'off') +
    # geom_nodelab(aes(x=x, label= paste(sprintf("%.1f", round(-x, digits = 1)))),size = 5, hjust= -0.4,vjust = 0.2,show.legend=F) +
    #geom_label(aes(x=x,y = y, label= paste0("Ne=",sprintf("%.1f", round(theta,1)))), size = 2.5, show.legend=F) +
    geom_label(aes(x=-x.theta,y = y.theta, label= paste0("Ne=",sprintf("%.1f", round(theta,1)))), size = 2.5, show.legend=F) +
    #geom_nodelab(aes(x=x, y=y, label=label, colour = "blue")) +
    #geom_treescale(fontsize=15, linesize=3, offset=2) +
    theme_tree2() +
    #coord_cartesian(clip = 'off') +
    theme_tree2(plot.margin = margin(1,1,1,1, "cm")) +
    theme(panel.grid.major=element_line(colour = "grey70"),
          panel.grid.major.y = element_blank()) +
    expand_limits(x = c(1,-5)) +
    xlab("Mya") +
    labs(size = "Effective population\nsize (millions)")
  
  if(dim(trace.mig)[1]>=1){
    p[[trace.name]] <- p[[trace.name]] + geom_segment(data = trace.mig, aes(x = -mig.out.x, y = mig.out.y, xend = -mig.in.x, yend = mig.in.y), arrow = arrow(length = unit(0.2,"cm"))) +
      geom_label(data = trace.mig, aes(x = -mig.out.x, y = mig.out.y.label, label = paste0(round(mig.rate.num.ind.per.gen, 3))), size = 2.5)
  }
  
  p[[trace.name]] <- revts(p[[trace.name]]) + scale_x_continuous(breaks = seq(-20, 0, 1),
                                                                 labels = abs(seq(20, 0, -1)))
}

all_plot <- (p[[1]] + ggtitle(expression(paste(italic("H. titia"), " - No Migration")))) + 
  (p[[2]] + ggtitle(expression(paste(italic("H. americana/calverti"), " - No Migration")))) + 
  (p[[3]] + ggtitle(expression(paste(italic("H. titia"), " - Migration")))) +
  (p[[4]] + ggtitle(expression(paste(italic("H. americana/calverti"), " - Migration")))) +
  plot_annotation(title = "Conspecific Loci calling", theme = theme(plot.title = element_text(size = 25)),
                  tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(ncol = 2)

ggsave(file = paste0(plot.dir, "conspecific_g-phocs-plots.png"), all_plot, height = 10, width = 12) 
ggsave(file = paste0(plot.dir, "conspecific_g-phocs-plots.pdf"), all_plot, height = 10, width = 12) 

# Merge all tau
trace.df <- trace.df[row.names(trace.df)!=1,]
trace.df$model <- row.names(trace.df)

trace.df$species <- sapply(strsplit(trace.df$model, "-"), "[", 5)
trace.df$mig.type <- sapply(strsplit(trace.df$model, "-"), "[", 6)
trace.df$beta <- sapply(strsplit(trace.df$model, "-"), "[", 4)
names(trace.df)

for(i in 1:length(trace.df$model)){
  HPD.temp <- HPD.trace.list[[i]]
  trace.df$HPD[i] <- paste0(paste(HPD.temp[1,],"=", HPD.temp[2,]), collapse = " : ")
}

trace.df
write.table(trace.df, paste0("4_Manuscript/data/G-Phocs/","G-Phocs_results_conspecific_genome_", gsub(u, replacement = "-", pattern ="\\."),".txt"),
            row.names = F, quote=F)
