## Whole genome sequence data for four samples including unknown BC1 backcross from CUAJ
#Hamilton LEA and PCR plots
.libPaths("/nobackup/tmjj24/apps/R/x86_64-pc-linux-gnu-library/4.2/")
library(vcfR)
library(LEA)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(introgress)
## Read in data

input_dir <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/VCF_chrom_r10000/"
output_dir <- "/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_manuscript/WGS_titia/VCF_chrom_r10000/results/"
plot.dir <- "/home/tmjj24/scripts/Github/Thesis-Phylogeographic-Hetaerina/4_Manuscript/plots/WGS_titia/"
dir.create(output_dir)
dir.create(plot.dir)
vcf <- read.vcfR(paste0(input_dir ,"WGS_titia_chr1-12.vcf"))

# Reformat names
gsub("^.*/","",colnames(vcf@gt))
colnames(vcf@gt) <- gsub("^.*/","",colnames(vcf@gt))
colnames(vcf@gt) <- substr(colnames(vcf@gt), 1, nchar(colnames(vcf@gt))-4)
colnames(vcf@gt) <- gsub(colnames(vcf@gt), pattern = "_R", replacement = "")
colnames(vcf@gt)[1] <- "FORMAT"
colnames(vcf@gt)

sample_map <- colnames(vcf@gt)[-1]

vcf@gt[1:10,]
vcf.bi <- vcf[is.biallelic(vcf)]
vcf.ply <- vcf.bi[is.polymorphic(vcf.bi, na.omit = T)]
vcf.ply@gt[1:10,1:6]
# Convert to geno
X_chrom <- "HetTit1.0.p.scaff-12-96647824"

table(vcf.ply@fix[,1]!=X_chrom)
# remove X chrom
geno <- extract.gt(vcf.ply[vcf.ply@fix[,1]!=X_chrom])

geno.mat <- as.matrix(geno)
# table(geno.mat)
geno.mat[geno.mat=="1/1"] <- 2
geno.mat[geno.mat=="0/1"] <- 1
geno.mat[geno.mat=="0/0"] <- 0
# remove samples that are entirely heterozgous
is.only.het <- apply(geno.mat, MARGIN = 2, function(x) gsub(paste0(unique(x), collapse = ""), pattern = "NA",replacement = "")=="1")
length(which(is.only.het))
if(any(is.only.het)){geno.mat <- geno.mat[,-which(is.only.het)]}

geno.mat[is.na(geno.mat)] <- 9

write.table(x = geno.mat, file = paste0(output_dir, "WGS_titia.geno"),
            col.names = F, row.names = F, quote = F, sep = "")

#Conduct PCA
geno2lfmm(paste0(output_dir, "WGS_titia.geno"),
          paste0(output_dir, "WGS_titia.geno.lfmm"))

#PCA
pc <- pca(paste0(output_dir, "WGS_titia.geno.lfmm"), scale = TRUE)

pc.sum <- summary(pc)

pca.data <- data.frame(pc$projections)
colnames(pca.data) <- paste0("pca", 1:dim(pc$projections)[2])
pca.labs <- paste("pca", 1:dim(pc$projections)[2], " (",round(as.numeric(pc.sum[2,1:dim(pc$projections)[2]])*100, 1), "%)", sep = "")
pca.data$samples <- colnames(vcf.bi@gt)[2:dim(vcf.bi@gt)[2]]

p <- ggplot(pca.data) +
  geom_point(aes(as.numeric(pca1), as.numeric(pca2)), size = 2) +
  geom_label_repel(aes(as.numeric(pca1), as.numeric(pca2), label = samples), size = 2) +
  xlab(pca.labs[1]) +
  ylab(pca.labs[2]) 

ggsave(paste0(plot.dir, "PCA_WGS_plot.pdf"), p)
## Conduct snmf
max.K <- 4
obj.at <- snmf(paste0(output_dir, "WGS_titia.geno"), K = 1:max.K, ploidy = 2, entropy = T,
              CPU = 5, project = "new", repetitions = 20, alpha = 100)
titia.snmf <- load.snmfProject(file = paste0(output_dir, "WGS_titia.snmfProject"))
titia.snmf.sum <- summary(titia.snmf)

# Cross entropy plot
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

ggsave(paste0(plot.dir,"WGS_cross_entropy.pdf"), plot = ce.plot)


# plot(titia.snmf, col = "blue4", cex = 1.4, pch = 19)
K = 2
best <- which.min(cross.entropy(titia.snmf, K = K))
qmatrix = Q(titia.snmf, K = K, run = best)

qtable <- cbind(rep(sample_map,K), rep(1:K, each = length(sample_map)), c(qmatrix[,1:K]))
qtable <-  data.frame(qtable)
colnames(qtable) <- c("sample","Qid", "Q")

p.bar <- ggplot(qtable)+
  geom_bar(stat="identity", aes(sample, as.numeric(Q), fill = Qid,), position = "stack", width = 1, col = "black") 
# theme(axis.text.x = element_text(size=6))
ggsave(paste0(plot.dir, "snmf_K", K, "_bar_WGS_plot.png"), p.bar)

## What is the coverage of the X chromosome of the hybrid individuals
X_chrom <- "HetTit1_0_p_scaff-12-96647824"

##### introgression
## install.packages("genetics")
## url <- "http://cran.nexr.com/src/contrib/introgress_1.2.3.tar.gz"
## pkgFile <- "introgress_1.2.3.tar.gz"
## download.file(url = url, destfile = pkgFile)
## install.packages(pkgs=pkgFile, type="source", repos=NULL)
library(introgress)

## Relevent sample sites
Atl.sites <- c("CT", "TXRS", "CUAJ", "MIXT")
Pac.sites <- c("PUMA", "RLPE", "ZANA", "TULI")

sites <- data.frame(samples = sample_map, site.sub = substr(sample_map, start = 1, stop = nchar(sample_map)-3))
#Subset sites
hybrid.sites <- sites[sites$site.sub%in%c(Atl.sites, Pac.sites),]

# Extract genotypes from vcf
mat <- extract.gt(vcf.ply)

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
atl.af<-conv.mat[,atl.samps]/2
pac.af<-(rowSums(conv.mat[,pac.samps], na.rm=T)/(rowSums(is.na(conv.mat[,pac.samps]) == FALSE)))/2
# atl.af<-(rowSums(conv.mat[,atl.samps], na.rm=T)/(rowSums(is.na(conv.mat[,atl.samps]) == FALSE)))/2
# pac.af<-(rowSums(conv.mat[,pac.samps], na.rm=T)/(rowSums(is.na(conv.mat[,pac.samps]) == FALSE)))/2
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
                       lg=vcf.ply@fix[,1][is.na(diff) == FALSE & diff > fixed_thres],
                       marker.pos=vcf.ply@fix[,2][is.na(diff) == FALSE & diff > fixed_thres])

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

#we now have a gt matrix in proper format for introgress
#convert genotype data into a matrix of allele counts
count.matrix<-prepare.data(admix.gen=gen.mat, loci.data=locus.info,
                           parental1="1",parental2="0", pop.id=F,
                           ind.id=F, fixed=T)

hi.index.sim<-est.h(introgress.data=count.matrix,loci.data=locus.info,
                    fixed=T, p1.allele="1", p2.allele="0")

png(paste0(plot.dir, "hybrid_geno_count_WGS_plot.png"), width = 2000, height = 1000, units = "px")
mk.image(introgress.data=count.matrix, loci.data=locus.info,
         hi.index=hi.index.sim, ylab.image="Individuals",
         xlab.h="population of Pacific ancestry", pdf=F,
         col.image=c(rgb(1,0,0,alpha=.5),rgb(0,0,0,alpha=.8),rgb(0,0,1,alpha=.5)))
dev.off()

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

locus.geno.type$sample <- as.factor(locus.geno.type$sample)

# Reorder samples by hybrid index
hybrid.sites$new.lat <- hybrid.sites$Lat
#Annoying that TXRS is plotted below CUAJ when using lat
hybrid.sites$new.lat[hybrid.sites$site.sub=="TXRS"] <- hybrid.sites$new.lat[hybrid.sites$site.sub=="TXRS"]+1
CUAJ.lat <- hybrid.sites$new.lat[hybrid.sites$samples=="CUAJa02"]
locus.geno.type$sample <- fct_relevel(locus.geno.type$sample, colnames(gen.mat)[order(as.numeric(hi.index.sim[,2]), decreasing = T)])
locus.geno.type$sample <- fct_relevel(locus.geno.type$sample, colnames(gen.mat)[order((hybrid.sites$new.lat-CUAJ.lat)+(-hi.index.sim[,2])*0.001, decreasing = F)])
hybrid.sites$samples <- fct_relevel(hybrid.sites$sample, hybrid.sites$samples[order((hybrid.sites$new.lat-CUAJ.lat)+(-hi.index.sim[,2])*0.001, decreasing = F)])

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

het.cols <- c("#AF0F09","#3E3C3A","#E5D9BA")

#horizontal
p.h <- ggplot(locus.geno.type) +
  geom_raster(aes(x = as.factor(BPcum), y = sample, fill = genotype)) +
  geom_vline(xintercept = vline_chrom$cumsum.chrom+1) +
  geom_text(data = vline_chrom, aes(x  = text_chrom_pos, y = -1, label = by), ) +
  scale_fill_manual(values = c(het.cols[c(3,2,1)],"white")) +
  coord_cartesian(ylim = c(1, length(unique(locus.geno.type$sample))), # This focuses the x-axis on the range of interest
                  clip = 'off') +
  theme(plot.margin = unit(c(3,1,3,1), "lines"),legend.position = c(0.5,1.05),legend.direction = "horizontal",
        axis.text.x = element_blank(), axis.title.x = element_blank()) +
  theme(legend.key=element_rect(colour="black")) +
  labs(fill = "Genotype")

ggsave(filename = paste0(plot.dir,"WGS_introgress_grid_titia_chr1-12_r10000.png"), plot = p.h)

#p.h

# Just X chromosome
p.x <- ggplot(locus.geno.type[locus.geno.type$lg=="X",]) +
  geom_raster(aes(x = as.factor(BPcum), y = sample, fill = genotype)) +
  scale_fill_manual(values = c(het.cols[c(3,2,1)],"white")) +
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

# Inbuilt hybrid plot that used base R plotting
introgress::triangle.plot(hi.index=hi.index.sim, int.het=het, pdf = F)

# attach het and hybrid to sample dataframe
hybrid.sites$het.fst <- het
hybrid.sites$hybrid.index <- hi.index.sim[,2]

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
  theme(legend.position = "bottom")

ggsave(file = paste0(plot.dir,"introgress_tri_plot_v2.png"), q, width = 6.5, height = 6)


#### CONVERTED CODE NOT TESTED FROM HERE


# Get depth of reads

# What is the coverage of samples across the genome
vcf.ply
# Extrat raw read coverage
depth_df <- extract.gt(vcf.ply, element = "DP", as.numeric = T)
# Melth into tidy format
depth_df <- reshape2::melt(depth_df, id = c("SNP", "SAMPLE"), value.name = "depth_var")

# Extract chromosome position and chromosome bnumber for each SNP
depth_df$chrom <- str_split_i(depth_df$Var1, pattern = "-", i = 2)
table(depth_df$chrom)
depth_df$chrom_length <- str_split_i(str_split_i(depth_df$Var1, pattern = "-", i = 3), pattern = "_", i = 1)
depth_df$position <- str_split_i(str_split_i(depth_df$Var1, pattern = "-", i = 3), pattern = "_", i = 2)
depth_df$chrom <- factor(depth_df$chrom, levels = sort(as.numeric(as.character(unique(depth_df$chrom)))))

# Extract depth for the hybrid individual
depth_CUAJ <- depth_df[depth_df$Var2=="CUAJa03",]
table(depth_CUAJ$chrom)
depth_CUAJ$chrom <- as.factor(depth_CUAJ$chrom)
depth_CUAJ$chrom <- factor(depth_CUAJ$chrom, levels = sort(as.numeric(as.character(unique(depth_CUAJ$chrom)))))

# Extract SNP genotype for hybrid
depth_het <- apply(sample_map, MARGIN = 1, function(x){
  print(x[2])
  depth_temp <- depth_df[depth_df$Var2==x[2],]
  vcf_temp <- vcf[sample = x[2]]
  depth_temp$is_het <- is.het(extract.gt(vcf_temp, element = "GT"), na_is_false = F)
  return(depth_temp)
  })

depth_het <- do.call("rbind", depth_het)

vcf_CUAJ <- vcf.ply[sample = "CUAJa02"]
# Get heterozgousity of each SNP
depth_CUAJ$is_het <- is.het(extract.gt(vcf_CUAJ, element = "GT"), na_is_false = F)

# Plot Read depth over heterozgousity for each SNP for each Chromosome
dp.plot <- ggplot(depth_CUAJ[depth_CUAJ$chrom%in%1:12&!is.na(depth_CUAJ$is_het),]) +
  geom_boxplot(aes(x = chrom, y = depth_var), fill = NA, outlier.alpha = 0, size = 2) +
  geom_jitter(aes(x = chrom, y = depth_var, fill = is_het), shape = 21, height = 0, width = 0.2, size = 3) +
  scale_fill_manual(values = c("#3E3C3A","#E5D9BA")) +
  labs(x = "Chromosome", y = "SNP raw read depth", fill = "Is heterozgous") +
  ylim(0,100) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(filename = paste0(plot.dir, "Depth_WGS_titia_chr1-12_r10000_CUAJa02.png"), plot = dp.plot)

dp_all <- ggplot(depth_het) +
  geom_jitter(aes(x = chrom, y = depth_var, col = is_het), shape = 19, height = 0, width = 0.2, size = 0.2) +
  geom_boxplot(aes(x = chrom, y = depth_var), col = "black", fill = NA, outlier.alpha = 0, size = 0.5, width = 0.5) +
  facet_wrap(~Var2) +
  scale_color_manual(values = c("deepskyblue","darkred")) +
  labs(x = "Chromosome", y = "SNP raw read depth", fill = "Is heterozgous") +
  ylim(0,100) +
  theme_bw() +
  theme(legend.position = "bottom")


ggsave(filename = "scripts/WGS_titia/WGS_all_samples_coverage_WGS_titia_chr1-12_r10000.png", plot = dp_all)


# What is the average hetero zgousity of each individual
chrom_sum_het <- depth_het %>%
  group_by(Var2) %>%
  group_by(chrom, .add = T) %>%
  summarise(pr_het = sum(is_het)/length(is_het), av_cov = mean(depth_var))

ggplot(chrom_sum_het) +
  geom_label(aes(x = chrom, label = Var2 , y = av_cov))

ggplot(chrom_sum_het) +
  geom_label(aes(x = chrom, label = Var2 , y = pr_het))

