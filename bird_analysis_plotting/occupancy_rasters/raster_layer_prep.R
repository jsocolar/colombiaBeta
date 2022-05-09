# This script prepares raster layers necessary to simulate the posterior occupancy probability for each 
# species across Colombia at 2 km resolution for one posterior iteration

library(sf)
library(reticulate)
library(raster)

`%ni%` <- Negate(`%in%`)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

##### Get elevation raster across Colombia #####
# Set up GEE session
use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()             # Trigger the authentication
np <- import("numpy")       # Import Numpy        needed for converting gee raster to R raster object
pd <- import("pandas")      # Import Pandas       ditto the above
# Get elevations for Colombia
countries <- ee$FeatureCollection('USDOS/LSIB_SIMPLE/2017')
roi <- countries$filterMetadata('country_na', 'equals', 'Colombia')  
DEM_full <- ee$Image("JAXA/ALOS/AW3D30/V2_2")
DEM <- DEM_full$select(list("AVE_DSM"),list("elevation"))
latlng <- ee$Image$pixelLonLat()$addBands(DEM)
latlng <- latlng$reduceRegion(reducer = ee$Reducer$toList(),
                              geometry = roi,
                              maxPixels = 10^9,
                              scale=2000)
# Convert to arrays
lats <- np$array((ee$Array(latlng$get("latitude"))$getInfo()))
lngs <- np$array((ee$Array(latlng$get("longitude"))$getInfo()))
elevs <- np$array((ee$Array(latlng$get("elevation"))$getInfo()))
# Convert to elevation raster
elevation <- data.frame(x = lngs, y = lats, elevation = elevs)
raster_elev <- raster::rasterFromXYZ(elevation,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# Reproject
raster_elev_AEA <- raster::projectRaster(raster_elev, crs = sp::CRS(AEAstring))
# Crop
source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/GIS_processing/get_mainland.R")
mainland <- st_transform(mainland, AEAstring)
raster_elev_AEA <- crop(raster_elev_AEA, extent(mainland))
# Save raster
raster::writeRaster(raster_elev, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/elev_raster/raster_elev.grd", overwrite = T)
raster::writeRaster(raster_elev_AEA, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/elev_raster/raster_elev_AEA.grd", overwrite = T)

##### Rasterize the Ayerbe maps #####
ayerbe_list_updated <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/ayerbe_maps/ayerbe_list_updated.RDS')
ayerbe_buffered_ranges_updated <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/ayerbe_maps/ayerbe_buffered_ranges_updated.RDS")

pb <- txtProgressBar(min = 0, max = length(ayerbe_buffered_ranges_updated), initial = 0, style = 3) 
for(i in 1:length(ayerbe_buffered_ranges_updated)){
  setTxtProgressBar(pb,i)
  ayerbe_raster <- fasterize::fasterize(st_sf(st_union(ayerbe_list_updated[[i]])), raster_elev_AEA)
  raster::writeRaster(ayerbe_raster, paste0('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Ayerbe_rasters/regular/', names(ayerbe_list_updated)[i], '.grd'))
  ayerbe_buffer_raster <- fasterize::fasterize(st_sf(ayerbe_buffered_ranges_updated[[i]]), raster_elev_AEA)
  raster::writeRaster(ayerbe_buffer_raster, paste0('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Ayerbe_rasters/buffered/', names(ayerbe_list_updated)[i], '_buffered.grd'))
}

##### Distance-to-range raster for each species #####
library(sf)
library(raster)
source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/GIS_processing/get_mainland.R")
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
mainland <- st_transform(mainland, AEAstring)
mainland_inward <- st_buffer(mainland, -7000)
ayerbe_list_updated <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/ayerbe_maps/ayerbe_list_updated.RDS')
bird_data <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data6_package.RDS")
birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/birds.RDS")
dtr_offset <-bird_data$means_and_sds$distance_to_range_offset
dtr_sd <- bird_data$means_and_sds$distance_to_range_sd
dtr_rescale <- bird_data$means_and_sds$distance_to_range_logit_rescale

# get all points
raster_elev_AEA <- raster::raster("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/elev_raster/raster_elev_AEA.grd")
pfd <- raster::as.data.frame(raster_elev_AEA, xy=T)

# now keep track of point indices that don't require distances
mask_raster <- raster::raster("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/mask_raster/mask.grd")
mdf <- raster::as.data.frame(mask_raster, xy=T)
all.equal(mdf[,c("x","y")], pfd[,c("x","y")])
points <- st_as_sf(raster::coordinates(raster_elev_AEA, spatial = T))

indices_include <- is.na(mdf$layer) & st_intersects(points, mainland, sparse = F)

coords <- st_coordinates(points)
sp_list <- unique(birds$species)

pb <- txtProgressBar(min = 0, max = length(sp_list), initial = 0, style = 3) 
for(i in 1:length(sp_list)){
  setTxtProgressBar(pb,i)
  sp <- sp_list[i]
  ayerbe_buffer_polygon <- st_union(ayerbe_buffered_ranges_updated[[gsub("_", " ", sp)]])
  indices_sp <- st_intersects(points, ayerbe_buffer_polygon, sparse = F)
  indices_m <- which(indices_sp & indices_include)
  
  points_m <- points[indices_m,]

  ayerbe_polygon <- st_union(ayerbe_list_updated[[gsub("_", " ", sp)]])
  range_linestring <- st_cast(ayerbe_polygon, "MULTILINESTRING")
  range_linestring_cropped <- st_intersection(mainland_inward, range_linestring)  
  inside <- as.numeric(as.numeric(st_distance(points_m, ayerbe_polygon)) == 0)
  if(nrow(range_linestring_cropped) > 0){
    distance_inside <- -1*st_distance(points_m[as.logical(inside),], range_linestring_cropped)
  }else{
    distance_inside <- rep(-2e+06, sum(inside))
  }
  distance_outside <- st_distance(points_m[!inside,], range_linestring)
  distances <- rep(NA, nrow(points))
  distances[indices_m][as.logical(inside)] <- distance_inside
  distances[indices_m][!inside] <- distance_outside
  transformed_distances <- boot::inv.logit(dtr_rescale*(distances - dtr_offset)/dtr_sd)
  td_df <- cbind(coords, transformed_distances)
  ayerbe_dist_raster_i <- raster::rasterFromXYZ(td_df,crs=AEAstring)
  raster::writeRaster(ayerbe_dist_raster_i, paste0('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Ayerbe_rasters/transformed_distance/', sp_list[i], '.grd'), overwrite = T)
}


pb <- txtProgressBar(min = 0, max = length(sp_list), initial = 0, style = 3) 
for(i in 1:length(sp_list)){
  setTxtProgressBar(pb,i)
  sp <- sp_list[i]
  ayerbe_buffer_polygon <- st_union(ayerbe_buffered_ranges_updated[[gsub("_", " ", sp)]])
  indices_sp <- st_intersects(points, ayerbe_buffer_polygon, sparse = F)
  indices_m <- which(indices_sp)
  
  points_m <- points[indices_m,]
  
  ayerbe_polygon <- st_union(ayerbe_list_updated[[gsub("_", " ", sp)]])
  range_linestring <- st_cast(ayerbe_polygon, "MULTILINESTRING")
  range_linestring_cropped <- st_intersection(mainland_inward, range_linestring)  
  inside <- as.numeric(as.numeric(st_distance(points_m, ayerbe_polygon)) == 0)
  if(nrow(range_linestring_cropped) > 0){
    distance_inside <- -1*st_distance(points_m[as.logical(inside),], range_linestring_cropped)
  }else{
    distance_inside <- rep(-2e+06, sum(inside))
  }
  distance_outside <- st_distance(points_m[!inside,], range_linestring)
  distances <- rep(NA, nrow(points))
  distances[indices_m][as.logical(inside)] <- distance_inside
  distances[indices_m][!inside] <- distance_outside
  transformed_distances <- boot::inv.logit(dtr_rescale*(distances - dtr_offset)/dtr_sd)
  td_df <- cbind(coords, transformed_distances)
  ayerbe_dist_raster_i <- raster::rasterFromXYZ(td_df,crs=AEAstring)
  raster::writeRaster(ayerbe_dist_raster_i, paste0('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Ayerbe_rasters/transformed_distance_unmasked/', sp_list[i], '.grd'), overwrite = T)
}





