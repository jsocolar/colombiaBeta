
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