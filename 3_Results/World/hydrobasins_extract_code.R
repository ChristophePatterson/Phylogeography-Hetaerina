#Hydrobasins test map
#install.packages("ggplot2")
#install.packages("elevatr")
library(sf)
library(sp)
library(ggplot2)
library(elevatr)

library(raster)
worldmap <- rgdal::readOGR("/home/tmjj24/World_map_data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
#Cropping the extent of worldmap to reduce plotting time
e <- raster::extent(-130, -50, 0, 50)
worldmap <- raster::crop(worldmap, e)

#hydrobasins <- sf::st_read(/home/tmjj24/World_map_data/hydrosheds/hybas_na_lev01-06_v1c/hybas_na_lev03_v1c.shp")
#hydrobasins.sub <- hydrobasins[hydrobasins$ENDO==0,]
#hydrobasins_geo <- st_geometry(hydrobasins.sub)
#plot(hydrobasins_geo)
#hydrorivers <- sf::st_read("/home/tmjj24/World_map_data/hydrosheds/HydroRIVERS_v10_na_shp/HydroRIVERS_v10_na_shp/HydroRIVERS_v10_na.shp")
#hydrorivers.sub <- hydrorivers
#hydrorivers.sub <- hydrorivers.sub[hydrorivers.sub$ORD_CLAS<=1,]
#hydrorivers.sub <- hydrorivers.sub[hydrorivers.sub$ORD_FLOW<=5,]
#hydrorivers.sub <- hydrorivers.sub[hydrorivers.sub$ENDORHEIC==0,]
#plot(hydrorivers.sub)
#
#save(hydrorivers.sub, file = "/home/tmjj24/World_map_data/hydrosheds/hydroRIVERS/HydroRIVERS_v10_na_sub1.rda")
#save(hydrobasins.sub, file = "/home/tmjj24/World_map_data/hydrosheds/hydroBASIN/hybas_lake_na_lev03_v1c_sub1_rda")
load("/home/tmjj24/World_map_data/hydrosheds/hydroRIVERS/HydroRIVERS_v10_na_sub1.rda")
load("/home/tmjj24/World_map_data/hydrosheds/hydroBASIN/hybas_lake_na_lev03_v1c_sub1_rda")
hydrorivers_geo <- st_geometry(hydrorivers.sub)
#hydrobasins.sub <- hydrobasins[hydrobasins$AREA_SQKM>20,]
hydrobasins_geo <- st_geometry(hydrobasins.sub)

#Projection
#worldmap <- rgdal::readOGR("/home/tmjj24/World_map_data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")

#crs(worldmap)
projcrs <- "+proj=longlat +datum=WGS84 +no_defs"

#Setting extent of points that cover Mexico and CR site
World <- SpatialPoints(world.e)
e <- bbox(World)
Mex <- SpatialPoints(mex.e)
e.Mex <- bbox(Mex)
Mex.z <- SpatialPoints(mex.e.zoom)
CR <- SpatialPoints(CR.e)
e.CR <- bbox(CR)

#Downloading elevational data
elevation <- get_elev_raster(World, z = 4, prj = projcrs)
elevation <- crop(elevation, e)
elevation.Mex <- get_elev_raster(Mex, z = 6, prj = projcrs)
elevation.Mex.z <- get_elev_raster(Mex.z, z = 8, prj = projcrs)
elevation.CR <- get_elev_raster(CR, z = 7, prj = projcrs)
#Setting water height to 0
elevation[elevation<=0] <- 0
elevation.CR[elevation.CR<0] <- 0
elevation.Mex[elevation.Mex<0] <- 0
elevation.Mex.z[elevation.Mex.z<0] <- 0

#Calculating aspect for America
slope.raster <- terrain(elevation*10, opt='slope')
aspect.raster <- terrain(elevation*10, opt='aspect')
hill.raster <- hillShade(slope.raster, aspect.raster, 45, 90)
#plot(hill.raster)
hill.m <- rasterToPoints(hill.raster)
hill.df <-  data.frame(hill.m)
colnames(hill.df) <- c("lon", "lat", "hill")

#Calculating aspect for Mexico
slope.raster.Mex <- terrain(elevation.Mex*10, opt='slope')
aspect.raster.Mex <- terrain(elevation.Mex*10, opt='aspect')
hill.raster.Mex <- hillShade(slope.raster.Mex, aspect.raster.Mex, 45, 90)
#plot(hill.raster.Mex)
hill.m.Mex  <- rasterToPoints(hill.raster.Mex)
hill.df.Mex <-  data.frame(hill.m.Mex)
colnames(hill.df.Mex) <- c("lon", "lat", "hill")

#Calculating aspect for Mexico zoom
slope.raster.Mex.z <- terrain(elevation.Mex.z*10, opt='slope')
aspect.raster.Mex.z <- terrain(elevation.Mex.z*10, opt='aspect')
hill.raster.Mex.z <- hillShade(slope.raster.Mex.z, aspect.raster.Mex.z, 45, 90)
#plot(hill.raster.Mex)
hill.m.Mex.z  <- rasterToPoints(hill.raster.Mex.z)
hill.df.Mex.z <-  data.frame(hill.m.Mex.z)
colnames(hill.df.Mex.z) <- c("lon", "lat", "hill")

#Calculating aspect for CR
slope.raster.CR <- terrain(elevation.CR*10, opt='slope')
aspect.raster.CR <- terrain(elevation.CR*10, opt='aspect')
hill.raster.CR <- hillShade(slope.raster.CR, aspect.raster.CR, 45, 90)
#plot(hill.raster.CR)
hill.m.CR <- rasterToPoints(hill.raster.CR)
hill.df.CR <-  data.frame(hill.m.CR)
colnames(hill.df.CR) <- c("lon", "lat", "hill")


#pdf("test_hydrobasins.pdf")
#ggplot() +
#  geom_raster(data = hill.df, aes(lon, lat, fill = hill), alpha = 1) +
#  geom_sf(data = hydrobasins_geo, col = "black", fill = rgb(0,0,0,alpha = 0.3)) +
#  geom_sf(data = hydrorivers_geo, col = "black", lineend = "round") +
#  geom_polygon(data = worldmap, aes(long, lat, group = group), col = "black", fill = NA) +
#  theme(legend.position="none") +
#  scale_fill_gradientn(colours=c("black","white")) +
#  coord_sf(xlim = c(-97,-92), ylim = c(14,19))
#dev.off()
