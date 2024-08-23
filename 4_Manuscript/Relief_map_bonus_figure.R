library(ggplot2)
library(ggrepel)

# Outlines for each plot
SNP.library.name <- "Hetaerina_titia_ddRAD_titia_dg"
species <- "titia"
dir.path <- paste0("4_Manuscript/data/SNP_libraries/",SNP.library.name,"/")
analysis.name <- "H_titia_complete_snp_noX"

plot.dir <- paste0("4_Manuscript/plots/",SNP.library.name,"/")

world.e <- data.frame(Long = c(-120,-70), Lat = c(12,40))
mex.e <- data.frame( Long = c(-106,-92), Lat = c(22,14.5))
mex.e.zoom <- data.frame( Long = c(-96,-94), Lat = c(18,16))
CR.e <- data.frame(Long = c(-85,-82), Lat = c(8,11))

# World projection
crs_text <- "+proj=laea +x_0=0 +y_0=0 +lon_0=-74 +lat_0=40"
crs_text <- paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=",360+mean(sites$Long), " +lat_0=",mean(sites$Lat))

# Download terrrian and river data
source("4_Manuscript/hydrobasins_extract_code_sf.R")

# Create segments
num.segments <- 100
num.samples <- 200
longs <- seq(mex.e.zoom$Long[1], mex.e.zoom$Long[2], length.out = num.segments)
lats <- seq(mex.e.zoom$Lat[2], mex.e.zoom$Lat[1], length.out = num.samples)

# Combine into list
list.segments <- list()
for(i in 1:num.segments){
  list.segments[[i]] <- cbind(longs[i], lats, i)
}

# Convert list to df and then to sf
segments.df <- as.data.frame(do.call("rbind", list.segments))
segments.sf <- st_as_sf(segments.df, coords = c("V1","lats"), crs = projcrs)
segments.sf$Lon <- segments.df$V1
segments.sf$Lat <- segments.df$lats

# Plot to check
plot(elevation.Mex.z)
## Plot points takes a while if there are a lot
# points(segments.sf)

##  Function for long and elevation
relief_calculator <- function(Long, Elevation){
  temp_elevation_plot <- (Long)-((Elevation/max(Elevation, na.rm = T))*0.4)
  temp_elevation_plot[Elevation<=0] <- NA
  return(temp_elevation_plot)
}

# Extract elevation 
segments.sf$elevation <- extract(elevation.Mex.z, segments.sf)[,2]

ma.elv <- max(segments.sf$elevation)
segments.sf$elevation.plot <- relief_calculator(Long = segments.sf$Lon, Elevation = segments.sf$elevation)


sites <- read.table(paste0(dir.path, "coorrd_Titia_complete",analysis.name ,".txt"))

## Relevent sample sites
Atl.sites <- c("CT", "TXRS", "CUAJ", "MIXT")
Pac.sites <- c("PUMA", "RLPE", "ZANA", "TULI")

het.cols <- c("#AF0F09","#E5D9BA","#3E3C3A")
#Subset sites
hybrid.sites <- sites[sites$site.sub%in%c(Atl.sites, Pac.sites),]
hybrid.sites$elevation <- extract(elevation.Mex.z, st_as_sf(hybrid.sites, coords = c("Long", "Lat"), crs = projcrs))[,2]
hybrid.sites$elevation.plot <- relief_calculator(Long = hybrid.sites$Long, Elevation = hybrid.sites$elevation)

# get hyrdobasin (not currently working)
hydrobasins_geo_df <- as.data.frame(st_coordinates(hydrobasins_geo[23]))# , coords = c("X", "Y"), crs = projcrs))
hydrobasins_geo_df$elevation <- extract(elevation.Mex.z, hydrobasins_geo_df[,c("X","Y")])[,2]
hydrobasins_geo_df <- hydrobasins_geo_df[hydrobasins_geo_df$X>mex.e.zoom$Long[1]&hydrobasins_geo_df$X<mex.e.zoom$Long[2],]
hydrobasins_geo_df <- hydrobasins_geo_df[hydrobasins_geo_df$Y>mex.e.zoom$Lat[2]&hydrobasins_geo_df$Y<mex.e.zoom$Lat[1],]
hydrobasins_geo_df$elevation.plot <- relief_calculator(Long = hydrobasins_geo_df$X, Elevation = as.numeric(hydrobasins_geo_df$elevation))
hydrobasins_geo_df<-  hydrobasins_geo_df[order(hydrobasins_geo_df$Y),]
hydrobasins_geo_df<-  hydrobasins_geo_df[round(seq(1, dim(hydrobasins_geo_df)[1], length.out = num.segments)),]


summary(hydrobasins_geo_df$elevation.plot)
relief_plot <- ggplot(segments.sf) +
  #geom_ribbon(aes(x = Lat, ymax = elevation.plot, ymin=Lon*0.1, group = i), fill = "grey") +
  geom_line(aes(x = Lat, y = elevation.plot, group = i, col = elevation), size  = 1, lineend = "round") +
  # geom_point(data = hydrobasins_geo_df, aes(x = Y, y = elevation.plot)) +
  geom_label_repel(data = hybrid.sites[!duplicated(hybrid.sites$Site.ID),],
                   aes(x = Lat, y = elevation.plot, label = Site.ID),
             nudge_y = 0.25, show.legend = F) +
  geom_point(data = hybrid.sites, aes(x = Lat, y = elevation.plot, fill = Ocean.drainage), shape = 21, size  = 3) +
  # scale_color_gradientn(colours=c("#5C2700","#FFE3D4")) +
  scale_color_gradientn(colours=rev(c("#725636","#392b1b","#7c8485","#b68c5d"))) +
  # scale_color_gradientn(colours=c("#7c8d4c",terrain.colors(5)[1:4])) +
  # scale_color_gradientn(colours=c("#7c8d4c","#b5ba61","#725428","#725428","#725428","#725428","#e5d9c2")) +
  scale_fill_manual(values = het.cols[c(2,1)]) +
  ylim(max(segments.sf$elevation.plot, na.rm = T), min(segments.sf$elevation.plot, na.rm = T)) +
  xlim(min(segments.sf$Lat), max(segments.sf$Lat)) +
  labs(x = "Latitude", y = "Longitude (°E) + elevation", color = "Elevation (m)", fill = "Ocean\nDrainage") +
  theme_bw() +
  theme(text = element_text(size =20))
relief_plot

ggsave(paste0(plot.dir,"Relief_map.png"), relief_plot, width = 12, height = 10)
ggsave(paste0(plot.dir,"Relief_map.pdf"), relief_plot, width = 12, height = 10)

# "#7c8d4c","#b5ba61","#725428","#e5d9c2","#b6e3db"

# c("#725636","#392b1b","#7c8485","#f0f0b5","#b68c5d")