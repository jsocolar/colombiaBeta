# This script generates a mask for the regions of colombia that we do not analyze.

library(sf)
library(raster)
library(ggplot2)

`%ni%` <- Negate(`%in%`)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/GIS_processing/hydrosheds_extraction.R")

raster_elev_AEA <- raster("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/elev_raster/raster_elev_AEA.grd")

pac_raster <- fasterize::fasterize(st_transform(pacific, AEAstring), raster_elev_AEA)
pac_raster[raster_elev_AEA > 1100] <- NA
mask_raster <- pac_raster
tac_raster <- fasterize::fasterize(st_transform(tacarcuna, AEAstring), raster_elev_AEA)
mask_raster[tac_raster == 1] <- 1
snsm_raster <- fasterize::fasterize(st_transform(snsm, AEAstring), raster_elev_AEA)
snsm_raster[raster_elev_AEA < 3000] <- NA
snsm_raster[is.na(raster_elev_AEA)] <- NA
mask_raster[snsm_raster == 1] <- 1

devtools::install_github("nebulae-co/colmaps")
dpts <- colmaps::departamentos

mask_dpts <- dpts[dpts$depto %in% c("Cesar", "La Guajira", "Magdalena", "Bolívar", "Atlántico", "Córdoba", "Sucre", "Antioquia", "Chocó"), ]

md <- st_buffer(st_transform(st_as_sf(mask_dpts), AEAstring),3000)
md <- st_cast(md, "MULTIPOLYGON")

md_raster <- fasterize::fasterize(md, raster_elev_AEA)
md_raster[raster_elev_AEA > 1100] <- NA

mask_raster[md_raster == 1] <- 1
raster::writeRaster(mask_raster, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/mask_raster/mask.grd", overwrite=T)
