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
bird_data <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data4_package.RDS")
birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/birds.RDS")
dtr_mean <- mean(birds$distance_from_range)
dtr_sd <- sd(birds$distance_from_range)

points <- st_as_sf(raster::coordinates(raster_elev_AEA, spatial = T)) # For reasons that I don't understand, it is more than an order of magnitude 
# faster to get the distances under the AEA projection than under epsg 4326.
coords <- st_coordinates(points)
sp_list <- unique(birds$species)

pb <- txtProgressBar(min = 0, max = length(sp_list), initial = 0, style = 3) 
for(i in 1:length(sp_list)){
  setTxtProgressBar(pb,i)
  sp <- sp_list[i]
  ayerbe_polygon <- st_union(ayerbe_list_updated[[gsub("_", " ", sp)]])
  distances <- as.numeric(st_distance(points, st_cast(ayerbe_polygon, to = "MULTILINESTRING")))*(2*(as.numeric(st_distance(points, ayerbe_polygon))>0) - 1) # The second part gives positive distances for outside-of-range and negative distances for in-range.  Turns out that as.numeric(st_distance(points, ayerbe))>0) is much faster than !st_within(points, ayerbe) 
  transformed_distances <- boot::inv.logit(4.7*(distances - dtr_mean)/dtr_sd)
  td_df <- cbind(coords, transformed_distances)
  ayerbe_dist_raster_i <- raster::rasterFromXYZ(td_df,crs=AEAstring)
  raster::writeRaster(ayerbe_dist_raster_i, paste0('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Ayerbe_rasters/transformed_distance/', sp_list[i], '.grd'), overwrite = F)
}


##### Predict occupancy across Colombia #####
raster_elev_AEA <- raster::raster("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/elev_raster/raster_elev_AEA.grd")
elev_df <- raster::as.data.frame(raster_elev_AEA, xy = T) 
elev_df$cell_id <- 1:nrow(elev_df)
names(elev_df)[3] <- "elevation"

source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/bird_analysis_plotting/get_posterior/get_posterior_z.R")
v5 <- cmdstanr::read_cmdstan_csv("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/v5_first_run/occupancy_v5_threads-202012282018-1-261afe.csv")
draws <- posterior::as_draws_df(v5$post_warmup_draws[1:2000,,])
# Long format: 
bird_data <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data4_package.RDS")
birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/birds.RDS")

z_info <- data.frame(bird_data$data[8:43])
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

saveRDS(forest_probs, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v5_predictions/iteration_1/forest_probs.RDS")
saveRDS(pasture_probs, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v5_predictions/iteration_1/pasture_probs.RDS")



###############################

forest_probs <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v5_predictions/iteration_1/forest_probs.RDS")
pasture_probs <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v5_predictions/iteration_1/pasture_probs.RDS")


cell_logratios <- get_avg_cell_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05)
colombia_logratio <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = NULL)

mean(cell_logratios$avg_logratio, na.rm = T)
weighted.mean(cell_logratios$avg_logratio, w=cell_logratios$n, na.rm=T)
spatstat::weighted.median(cell_logratios$med_logratio, w=cell_logratios$n, na.rm = T)

cell_logratios_df <- data.frame(x=forest_probs$x, y = forest_probs$y, logratio=cell_logratios$avg_logratio)
cell_logratios_raster <- raster::rasterFromXYZ(cell_logratios_df)
raster::plot(cell_logratios_raster)

cell_richness_df <- data.frame(x=forest_probs$x, y=forest_probs$y, richness = cell_logratios$n)
cell_richness_raster <- raster::rasterFromXYZ(cell_richness_df)
raster::plot(cell_richness_raster)


colombia_logratio_51_100 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 50 & cell_logratios$n < 101))
colombia_logratio_101_200 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 100 & cell_logratios$n < 201))
colombia_logratio_201_300 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 200 & cell_logratios$n < 301))
colombia_logratio_301_max <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 300))

mean(cell_logratios$avg_logratio[cell_logratios$n > 50 & cell_logratios$n < 101])
colombia_logratio_51_100

mean(cell_logratios$avg_logratio[cell_logratios$n > 100 & cell_logratios$n < 201])
colombia_logratio_101_200

mean(cell_logratios$avg_logratio[cell_logratios$n > 200 & cell_logratios$n < 301])
colombia_logratio_201_300

mean(cell_logratios$avg_logratio[cell_logratios$n > 300])
colombia_logratio_301_max



sum(!is.nan(overall))
overall2 <- overall[1:614]

mean(overall2[!is.nan(overall2)])

pointwise <- numerator/denominator

mean(pointwise[!is.nan(pointwise)])
min(pointwise[!is.nan(pointwise)])



maxforest <- maxpasture <- mfp <- rep(0,39)
for(i in 1:39){
  print(i)
  sp <- species_results_prelim$species[i]
  sp_df <- readRDS(paste0('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/species_predictions/', sp, '.RDS'))
  maxforest[i] <- max(sp_df$sp_forest)
  maxpasture[i] <- max(sp_df$sp_pasture)
  mfp[i] <- max(sp_df$sp_forest + sp_df$sp_pasture)/2
}





completeness_prop <- raster::as.data.frame(raster_elev, xy = T)
completeness_prop$recorded_1elev <- completeness_prop$present_1elev <- completeness_prop$recorded_buffer_1elev <- completeness_prop$present_buffer_1elev <- 
  completeness_prop$recorded_2elev <- completeness_prop$present_2elev <- completeness_prop$recorded_buffer_2elev <- completeness_prop$present_buffer_2elev <- 
  completeness_prop$fullsp_noelev <- completeness_prop$allsp_noelev <- completeness_prop$recordedsp_noelev <- 0


species_names <- gsub(" ", "_", names(ayerbe_list_updated))

bird_stan_data1_package <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data1_package.RDS")
birds <- bird_stan_data1_package$data
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