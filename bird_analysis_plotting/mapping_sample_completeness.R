library(sf)
library(reticulate)

`%ni%` <- Negate(`%in%`)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"


##### Set up GEE session #####
use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()             # Trigger the authentication
np <- import("numpy")       # Import Numpy        needed for converting gee raster to R raster object
pd <- import("pandas")      # Import Pandas       ditto the above
countries <- ee$FeatureCollection('USDOS/LSIB_SIMPLE/2017')
roi <- countries$filterMetadata('country_na', 'equals', 'Colombia')  
DEM_full <- ee$Image("JAXA/ALOS/AW3D30/V2_2")
DEM <- DEM_full$select(list("AVE_DSM"),list("elevation"))
latlng <- ee$Image$pixelLonLat()$addBands(DEM)
latlng <- latlng$reduceRegion(reducer = ee$Reducer$toList(),
                                geometry = roi,
                                maxPixels = 10^9,
                                scale=500)
  
# Convert to arrays
lats <- np$array((ee$Array(latlng$get("latitude"))$getInfo()))
lngs <- np$array((ee$Array(latlng$get("longitude"))$getInfo()))
elevs <- np$array((ee$Array(latlng$get("elevation"))$getInfo()))
  
# Convert to elevation raster
elevation <- data.frame(x = lngs, y = lats, elevation = elevs)
raster_elev <- raster::rasterFromXYZ(elevation,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

raster_elev_AEA <- raster::projectRaster(raster_elev, crs = sp::CRS(AEAstring))

raster::writeRaster(raster_elev, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/elev_raster/raster_elev.grd")
raster::writeRaster(raster_elev_AEA, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/elev_raster/raster_elev_AEA.grd")

###### Rasterize the Ayerbe maps
ayerbe_list_updated <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/ayerbe_maps/ayerbe_list_updated.RDS')
ayerbe_buffered_ranges_updated <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/ayerbe_maps/ayerbe_buffered_ranges_updated.RDS")

for(i in 1:length(ayerbe_buffered_ranges_updated)){
  if(floor(i/10) == i/10){print(i)}
  ayerbe_raster <- fasterize::fasterize(st_transform(st_sf(st_union(ayerbe_list_updated[[i]])), 4326), raster_elev)
  raster::writeRaster(ayerbe_raster, paste0('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Ayerbe_rasters/regular/', names(ayerbe_list_updated)[i], '.grd'))
  ayerbe_buffer_raster <- fasterize::fasterize(st_transform(st_sf(ayerbe_buffered_ranges_updated[[i]]), 4326), raster_elev)
  raster::writeRaster(ayerbe_buffer_raster, paste0('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Ayerbe_rasters/buffered/', names(ayerbe_list_updated)[i], '_buffered.grd'))
}


completeness_prop <- raster::as.data.frame(raster_elev, xy = T)
completeness_prop$recorded_1elev <- completeness_prop$present_1elev <- completeness_prop$recorded_buffer_1elev <- completeness_prop$present_buffer_1elev <- 
  completeness_prop$recorded_2elev <- completeness_prop$present_2elev <- completeness_prop$recorded_buffer_2elev <- completeness_prop$present_buffer_2elev <- 
  completeness_prop$fullsp_noelev <- completeness_prop$allsp_noelev <- completeness_prop$recordedsp_noelev <- 0


species_names <- gsub(" ", "_", names(ayerbe_list_updated))

birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_data_trimmed.RDS")
traitdata <- birds[!duplicated(birds$species),]
traitdata$elev2_lower <- traitdata$lower - (traitdata$upper - traitdata$lower)
traitdata$elev2_upper <- traitdata$upper + (traitdata$upper - traitdata$lower)
all_species <- traitdata$species
recorded_species <- unique(birds$species[birds$Q == 1])
elevations <- completeness_prop$elevation
elevations[is.na(elevations)] <- -99999

for(i in 1:length(species_names)){
  print(i)
  species <- species_names[i]
  ayerbe_raster <- raster::raster(paste0('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Ayerbe_rasters/regular/', gsub("_", " ", species), '.grd'))
  raster_df <- raster::as.data.frame(ayerbe_raster)
  raster_df$layer[is.na(raster_df$layer)] <- 0
  completeness_prop$fullsp_noelev <- completeness_prop$fullsp_noelev + raster_df$layer
  if(species %in% all_species){completeness_prop$allsp_noelev <- completeness_prop$allsp_noelev + raster_df$layer}
  if(species %in% recorded_species){completeness_prop$recordedsp_noelev <- completeness_prop$recordedsp_noelev + raster_df$layer}
  
  
  if(species %in% all_species){
    td <- traitdata[traitdata$species == species, ]
    raster_df$layer[(elevations < td$elev2_lower) | (td$elev2_upper < elevations)] <- 0
    completeness_prop$present_2elev <- completeness_prop$present_2elev + raster_df$layer
    if(species %in% recorded_species){
      completeness_prop$recorded_2elev <- completeness_prop$recorded_2elev + raster_df$layer
    }
    raster_df$layer[(elevations < td$lower) | (td$upper < elevations)] <- 0
    completeness_prop$present_1elev <- completeness_prop$present_1elev + raster_df$layer
    if(species %in% recorded_species){
      completeness_prop$recorded_1elev <- completeness_prop$recorded_1elev + raster_df$layer
    }
    
    ayerbe_raster_buffer <- raster::raster(paste0('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Ayerbe_rasters/buffered/', gsub("_", " ", species), '_buffered.grd'))
    raster_df <- raster::as.data.frame(ayerbe_raster_buffer)
    raster_df$layer[is.na(raster_df$layer)] <- 0
    raster_df$layer[(elevations < td$elev2_lower) | (td$elev2_upper < elevations)] <- 0
    completeness_prop$present_buffer_2elev <- completeness_prop$present_buffer_2elev + raster_df$layer
    if(species %in% recorded_species){
      completeness_prop$recorded_buffer_2elev <- completeness_prop$recorded_buffer_2elev + raster_df$layer
    }
    raster_df$layer[(elevations < td$lower) | (td$upper < elevations)] <- 0
    completeness_prop$present_buffer_1elev <- completeness_prop$present_buffer_1elev + raster_df$layer
    if(species %in% recorded_species){
      completeness_prop$recorded_buffer_1elev <- completeness_prop$recorded_buffer_1elev + raster_df$layer
    }
  }
}

saveRDS(completeness_prop, file = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/sample_completeness/completeness_prop.RDS")

recorded_over_full <- cbind(completeness_prop$x, completeness_prop$y, completeness_prop$recordedsp_noelev/completeness_prop$fullsp_noelev)
recorded_over_full_raster <- raster::rasterFromXYZ(recorded_over_full)
raster::plot(recorded_over_full_raster)

recorded_over_all <- cbind(completeness_prop$x, completeness_prop$y, completeness_prop$recordedsp_noelev/completeness_prop$allsp_noelev)
recorded_over_all_raster <- raster::rasterFromXYZ(recorded_over_all)
raster::plot(recorded_over_all_raster)

all_over_full <- cbind(completeness_prop$x, completeness_prop$y, completeness_prop$allsp_noelev/completeness_prop$fullsp_noelev)
all_over_full_raster <- raster::rasterFromXYZ(all_over_full)
raster::plot(all_over_full_raster)



recorded1_over_present1 <- cbind(completeness_prop$x, completeness_prop$y, completeness_prop$recorded_1elev/completeness_prop$present_1elev)
recorded1_over_present1_raster <- raster::rasterFromXYZ(recorded1_over_present1)
raster::plot(recorded1_over_present1_raster)


recorded2_over_present2 <- cbind(completeness_prop$x, completeness_prop$y, completeness_prop$recorded_2elev/completeness_prop$present_2elev)
recorded2_over_present2_raster <- raster::rasterFromXYZ(recorded2_over_present2)
raster::plot(recorded2_over_present2_raster)

recorded2_over_present2_buffer <- cbind(completeness_prop$x, completeness_prop$y, completeness_prop$recorded_buffer_2elev/completeness_prop$present_buffer_2elev)
recorded2_over_present2_buffer_raster <- raster::rasterFromXYZ(recorded2_over_present2_buffer)
raster::plot(recorded2_over_present2_buffer_raster)

recorded1_over_present1_buffer <- cbind(completeness_prop$x, completeness_prop$y, completeness_prop$recorded_buffer_1elev/completeness_prop$present_buffer_1elev)
recorded1_over_present1_buffer_raster <- raster::rasterFromXYZ(recorded1_over_present1_buffer)
raster::plot(recorded1_over_present1_buffer_raster)

species_results_prelim <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/species_results_prelim.RDS')
birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_data_trimmed.RDS")
birds$sp_id <- as.integer(as.factor(birds$species))
traitdata <- birds[!duplicated(birds$species), ]
sd_elev <- sd(birds$elev_sp_standard)

elevations <- completeness_prop$elevation
elevations[is.na(elevations)] <- -99999

local_pres_forest <- local_pres_pasture <- matrix(0, nrow = length(elevations), ncol = 2)
total_pres_forest <- total_pres_pasture <- matrix(0, nrow = nrow(traitdata), ncol = 2)

for(j in 1:2){
  for(i in 1:nrow(traitdata)){
    print(c(j, i))
    sp_id <- traitdata$sp_id[i]
    species <- traitdata$species[i]
    td <- traitdata[i,]
    ayerbe_raster_buffer <- raster::raster(paste0('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Ayerbe_rasters/buffered/', gsub("_", " ", species), '_buffered.grd'))
    raster_df <- raster::as.data.frame(ayerbe_raster_buffer)
    raster_df$layer[is.na(raster_df$layer)] <- 0
    elevations_sp_scaled <- ((elevations - td$lower)/(td$upper - td$lower) - .5)/sd_elev
    elevations_sp_scaled2 <- elevations_sp_scaled^2
    raster_df$layer[(elevations_sp_scaled < -1) | (elevations_sp_scaled > 3)] <- 0

    forest_logit <- species_results_prelim$forest_int[sp_id,j] + species_results_prelim$elev1_coef[sp_id,j]*elevations_sp_scaled + 
        species_results_prelim$elev2_coef[sp_id,j]*elevations_sp_scaled2
    pasture_logit <- forest_logit + species_results_prelim$pasture_offset[sp_id,j]
    species_probs_forest <- boot::inv.logit(forest_logit)*raster_df$layer
    species_probs_pasture <- boot::inv.logit(pasture_logit)*raster_df$layer
    
    total_pres_forest[i,j] <- sum(species_probs_forest)
    total_pres_pasture[i,j] <- sum(species_probs_pasture)
    
    local_pres_forest[,j] <- local_pres_forest[,j] + species_probs_forest
    local_pres_pasture[,j] <- local_pres_pasture[,j] + species_probs_pasture
  }
}



fp_total_mean_pres_ratio <- colSums(total_pres_pasture/total_pres_forest)/1614
fp_total_mean_pres_logratio <- colSums(log(total_pres_pasture/total_pres_forest))/1614


fp_local_mean_pres_ratio <- local_pres_pasture[local_pres_pasture[,1] != 0,]/local_pres_forest[local_pres_pasture[,1] != 0,]
summary(fp_local_mean_pres_ratio[,2])

all_pts <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Points/all_pts.RDS")


birds_q <- birds[birds$Q == 1, ]
quantile(birds_q$elev_sp_standard, .0155)
quantile(birds_q$elev_sp_standard, .95562)

# 94% of species-point combinations with detections are within Ayerbe elevational ranges