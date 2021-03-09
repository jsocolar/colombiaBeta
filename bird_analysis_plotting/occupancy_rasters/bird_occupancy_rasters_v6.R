# This script simulates the posterior occupancy probability for each species-point across Colombia at 2 km resolution
# for one posterior iteration

library(sf)
library(reticulate)

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
  raster::writeRaster(ayerbe_raster, paste0('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Ayerbe_rasters/regular/', names(ayerbe_list_updated)[i], '.grd'), overwrite = T)
  ayerbe_buffer_raster <- fasterize::fasterize(st_sf(ayerbe_buffered_ranges_updated[[i]]), raster_elev_AEA)
  raster::writeRaster(ayerbe_buffer_raster, paste0('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Ayerbe_rasters/buffered/', names(ayerbe_list_updated)[i], '_buffered.grd'), overwrite = T)
}

##### Distance-to-range raster for each species #####
library(sf)
library(raster)
source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/GIS_processing/get_mainland.R")
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
mainland <- st_transform(mainland, AEAstring)
mainland_inward <- st_buffer(mainland, -7000)
ayerbe_list_updated <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/ayerbe_maps/ayerbe_list_updated.RDS')

raster_elev_AEA <- crop(raster("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/elev_raster/raster_elev_AEA.grd"), extent(mainland))
bird_data <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data6_package.RDS")
birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/birds.RDS")
dtr_offset <-bird_data$means_and_sds$distance_to_range_offset
dtr_sd <- bird_data$means_and_sds$distance_to_range_sd
dtr_rescale <- bird_data$means_and_sds$distance_to_range_logit_rescale
points <- st_as_sf(raster::coordinates(raster_elev_AEA, spatial = T)) # For reasons that I don't understand, it is more than an order of magnitude 
                                                                      # faster to get the distances under the AEA projection than under epsg 4326.
indices_m <- st_intersects(points, mainland, sparse = F)
indices_n <- 1 - indices_m
points_m <- points[indices_m,]
points_n <- points[indices_n,]
coords <- st_coordinates(points)
sp_list <- unique(birds$species)

pb <- txtProgressBar(min = 0, max = length(sp_list), initial = 0, style = 3) 
for(i in 1:length(sp_list)){
  setTxtProgressBar(pb,i)
  sp <- sp_list[i]
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

##### Predict occupancy across Colombia #####
raster_elev_AEA <- raster::raster("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/elev_raster/raster_elev_AEA.grd")
elev_df <- raster::as.data.frame(raster_elev_AEA, xy = T) 
elev_df$cell_id <- 1:nrow(elev_df)
names(elev_df)[3] <- "elevation"

source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/bird_analysis_plotting/get_posterior/get_posterior_z_v6.R")
bird_data <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data6_package.RDS")
birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/birds.RDS")
draws <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/v6_draws/draws.RDS")
z_info <- data.frame(bird_data$data[8:41])
z_info$point <- birds$point
z_info$species <- birds$species

iter <- 1
pc <- get_prediction_components(draws, iter, z_info)

forest_probs <- pasture_probs <- elev_df[!is.na(elev_df$elevation),c("x","y","cell_id")]
pb <- txtProgressBar(min = 0, max = 1614, initial = 0, style = 3) 
for(i in 1:1614){
  setTxtProgressBar(pb,i)
  sp <- unique(birds$species[bird_data$data$id_sp == i])
  sp2 <- gsub("_", " ", sp)
  if(length(sp)!= 1){stop('sp does not have length 1')}
  sp_df <- elev_df
  # Load raster of transformed distances and add to sp_df
  dist_raster <- raster::raster(paste0("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Ayerbe_rasters/transformed_distance/", sp, ".grd"))
  sp_df$transformed_distance <- raster::as.data.frame(dist_raster)
  # Load raster of buffered range and set out-of-range pixels to zero.
  buffer_raster <- raster::raster(paste0('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Ayerbe_rasters/buffered/', sp2, '_buffered.grd'))
  buffer_df <- raster::as.data.frame(buffer_raster)
  buffer_df$layer[is.na(buffer_df$layer)] <- 0
  sp_df$in_range <- buffer_df$layer
  # Retain only rows that are in Colombia and in range
  sp_df <- sp_df[!is.na(sp_df$elevation) & sp_df$in_range == 1, ]
  # Get elevational min/max
  sp_lower <- unique(birds$lower[birds$species == sp])
  sp_upper <- unique(birds$upper[birds$species == sp])
  sp_breadth <- sp_upper - sp_lower
  # Remove rows outside of buffered elevation
  sp_df <- sp_df[(sp_df$elevation > (sp_lower - sp_breadth)) & (sp_df$elevation < (sp_upper + sp_breadth)), ]
  # Compute relev
  sp_df$relev <- ((sp_df$elevation - sp_lower)/sp_breadth - bird_data$means_and_sds$relev_offset)/bird_data$means_and_sds$relev_sd
  sp_df$relev2 <- sp_df$relev^2
  # Compute occupancy probabilities
  sp_pc <- pc[i, ]
  sp_forest_logit <- as.numeric((sp_pc$logit_psi_forest + sp_df$relev * sp_pc$b1_relev_sp + sp_df$relev2 * sp_pc$b1_relev2_sp +
    sp_df$relev * sp_pc$lowland * sp_pc$b1_x_lowland_relev + sp_df$relev2 * sp_pc$lowland * sp_pc$b1_x_lowland_relev2 +
    sp_df$transformed_distance * sp_pc$b5_distance_to_range_sp)[,1])
  sp_pasture_logit <- sp_forest_logit + sp_pc$logit_psi_pasture_offset
  
  integrate_forest <- function(k){integrate(function(x){dnorm(x, sp_forest_logit[k], sp_pc$sigma_sp_cl)*boot::inv.logit(x)},
                                          sp_forest_logit[k]-6*sp_pc$sigma_sp_cl, sp_forest_logit[k]+6*sp_pc$sigma_sp_cl)$value}
  sp_df$forest_prob <- unlist(lapply(1:nrow(sp_df), integrate_forest))
  
  integrate_pasture <- function(k){integrate(function(x){dnorm(x, sp_pasture_logit[k], sp_pc$sigma_sp_cl)*boot::inv.logit(x)},
                                            sp_pasture_logit[k]-6*sp_pc$sigma_sp_cl, sp_pasture_logit[k]+6*sp_pc$sigma_sp_cl)$value}
  sp_df$pasture_prob <- unlist(lapply(1:nrow(sp_df), integrate_pasture))
  
  forest_probs <- merge(forest_probs, sp_df[,c("cell_id","forest_prob")], by = "cell_id", all = T)
  forest_probs[,i+3][is.na(forest_probs[,i+3])] <- 0
  pasture_probs <- merge(pasture_probs, sp_df[,c("cell_id", "pasture_prob")], by = "cell_id", all = T)
  pasture_probs[,i+3][is.na(pasture_probs[,i+3])] <- 0
  names(forest_probs)[i + 3] <- names(pasture_probs)[i + 3] <- sp
}

saveRDS(forest_probs, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v6_predictions/iteration_1/forest_probs.RDS")
saveRDS(pasture_probs, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v6_predictions/iteration_1/pasture_probs.RDS")

