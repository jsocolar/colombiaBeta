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
                                scale=1000)
  
# Convert to arrays
lats <- np$array((ee$Array(latlng$get("latitude"))$getInfo()))
lngs <- np$array((ee$Array(latlng$get("longitude"))$getInfo()))
elevs <- np$array((ee$Array(latlng$get("elevation"))$getInfo()))
  
# Convert to elevation raster
elevation <- data.frame(x = lngs, y = lats, elevation = elevs)
raster_elev <- raster::rasterFromXYZ(elevation,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

raster_elev_AEA <- raster::projectRaster(raster_elev, crs = sp::CRS(AEAstring))

raster::writeRaster(raster_elev, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/elev_raster/raster_elev.grd", overwrite = T)
raster::writeRaster(raster_elev_AEA, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/elev_raster/raster_elev_AEA.grd", overwrite = T)

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

completeness_prop <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/sample_completeness/completeness_prop.RDS")

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

breaks <- seq(0,1,.01)
col=c(rep("gray90", 63), viridis::viridis(37))

raster::plot(recorded1_over_present1_raster, col=col, breaks = breaks)



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

logit_normal_mean_1 <- function(logit.mu, logit.sigma){
  integrate(function(x){boot::inv.logit(x)*dnorm(x,logit.mu, logit.sigma)}, lower = logit.mu - 7*logit.sigma, upper = logit.mu + 7*logit.sigma, subdivisions = 500L)
}

logit_normal_mean <- function(l.mu, l.sigma){
  integrals <- sapply(l.mu, logit_normal_mean_1, logit.sigma = l.sigma)
  if(!all(integrals[4,]=="OK")){stop('some integrals not OK!')}
  return(as.numeric(integrals[1,]))
}

jmax <- 1

local_ratios01 <- local_ratios05 <- local_ratios10 <- 
  use01_sum <- use05_sum <- use10_sum <-
  matrix(0, nrow = length(elevations), ncol = jmax)

total_pres_forest01 <- total_pres_pasture01 <- total_pres_forest05 <- total_pres_pasture05 <- 
  total_pres_forest10 <- total_pres_pasture10 <- matrix(0, nrow = nrow(traitdata), ncol = jmax)


raster_probs <- raster_elev
raster_probs$pasture <- raster_probs$forest <- 0
raster_probs <- raster::dropLayer(raster_probs, 1)

for(j in 1:jmax){
  dir.create(paste0("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/species_predictions/iteration_",j))
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
    raster_df$layer[(elevations_sp_scaled < -1.5) | (elevations_sp_scaled > 1.5)] <- 0
    
    elevations_sp_scaled_forcalc <- ((seq(0,6000,20) - td$lower)/(td$upper - td$lower) - .5)/sd_elev
    elevations_sp_scaled2_forcalc <- elevations_sp_scaled_forcalc^2

    forest_logits_abstract <- species_results_prelim$forest_int[sp_id,j] + species_results_prelim$elev1_coef[sp_id,j]*elevations_sp_scaled_forcalc + 
        species_results_prelim$elev2_coef[sp_id,j]*elevations_sp_scaled2_forcalc
    pasture_logits_abstract <- forest_logits_abstract + species_results_prelim$pasture_offset[sp_id,j]
    
    species_probs_forest_abstract <- logit_normal_mean(forest_logits_abstract, species_results_prelim$cluster_sigma[j])
    species_probs_pasture_abstract <- logit_normal_mean(pasture_logits_abstract, species_results_prelim$cluster_sigma[j])
    
    relevant_elevs <- elevations[raster_df$layer!=0]
    relevant_elevs[relevant_elevs < 0] <- 0
    relevant_elevs_index <- floor(relevant_elevs*(301/6000)) + 1

    species_probs_forest <- species_probs_pasture <- raster_df$layer
    species_probs_forest[species_probs_forest!=0] <- species_probs_forest_abstract[relevant_elevs_index]
    species_probs_pasture[species_probs_pasture!=0] <- species_probs_pasture_abstract[relevant_elevs_index]
    
    raster_probs$forest <- species_probs_forest
    raster_probs$pasture <- species_probs_pasture
    
    raster::writeRaster(raster_probs, filename = paste0("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/species_predictions/iteration_",
                                                        j, "/", species, ".grd"))
    

    use01 <- species_probs_forest > .01 | species_probs_pasture > .01
    use05 <- species_probs_forest > .05 | species_probs_pasture > .05
    use10 <- species_probs_forest > .10 | species_probs_pasture > .10
    
    use01_sum[,j] <- use01_sum[,j] + use01
    use05_sum[,j] <- use05_sum[,j] + use05
    use10_sum[,j] <- use10_sum[,j] + use10
    
    lr01 <- lr05 <- lr10 <- log(species_probs_forest/species_probs_pasture)
    lr01[!use01] <- 0
    lr05[!use05] <- 0
    lr10[!use10] <- 0
    
    local_ratios01[,j] <- local_ratios01[,j] + lr01
    local_ratios05[,j] <- local_ratios05[,j] + lr05
    local_ratios10[,j] <- local_ratios10[,j] + lr10
    
    
    total_pres_forest01[i,j] <- sum(species_probs_forest[use01])
    total_pres_pasture01[i,j] <- sum(species_probs_pasture[use01])
    
    total_pres_forest05[i,j] <- sum(species_probs_forest[use05])
    total_pres_pasture05[i,j] <- sum(species_probs_pasture[use05])
    
    total_pres_forest10[i,j] <- sum(species_probs_forest[use10])
    total_pres_pasture10[i,j] <- sum(species_probs_pasture[use10])
  }
}


tp05 <- total_pres_pasture05[total_pres_pasture05 > 0]
tf05 <- total_pres_forest05[total_pres_forest05 > 0]
mean(log(tf05/tp05))
exp(mean(log(tf05/tp05)))


summary(local_ratios05[use05_sum > 0]/use05_sum[use05_sum > 0])

raster_local <- raster_elev
raster_local$local <- local_ratios05/use05_sum

fp_local_mean_pres_ratio <- local_pres_pasture[local_pres_pasture[,1] != 0,]/local_pres_forest[local_pres_pasture[,1] != 0,]
summary(fp_local_mean_pres_ratio[,2])

all_pts <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Points/all_pts.RDS")


birds_q <- birds[birds$Q == 1, ]
quantile(birds_q$elev_sp_standard, .0155)
quantile(birds_q$elev_sp_standard, .95562)

# 94% of species-point combinations with detections are within Ayerbe elevational ranges