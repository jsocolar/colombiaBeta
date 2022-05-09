##### Predict occupancy across Colombia #####

raster_elev_AEA <- raster::raster("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/elev_raster/raster_elev_AEA.grd")
raster_agg <- raster::aggregate(raster_elev_AEA, 10)
raster::values(raster_agg) <- 1:(raster::ncell(raster_agg))
raster_disagg <- raster::disaggregate(raster_agg, 10)
elev_df <- raster::as.data.frame(raster_elev_AEA, xy = T) 
elev_df$cell_id <- 1:nrow(elev_df)
sr_df <- raster::as.data.frame(raster::crop(raster_disagg, raster_elev_AEA), xy = T)
all.equal(elev_df[,c("x","y")], sr_df[,c("x","y")])
elev_df$sr_id <- sr_df$elevation
names(elev_df)[3] <- "elevation"

source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/bird_analysis_plotting/get_posterior/get_posterior_z_v6.R")
bird_data <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data6_package.RDS")
birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/birds.RDS")
draws <- posterior::as_draws_df(readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/v9_final/draws.RDS"))

# draws <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/v6_draws/draws.RDS")
z_info <- data.frame(bird_data$data[8:41])
z_info$point <- birds$point
z_info$species <- birds$species

iter <- 4000
pc <- get_prediction_components(draws, iter, z_info)

rm(draws)
gc()

forest_probs_integrated <- forest_probs_semi <- forest_probs_rep <- 
  pasture_probs_integrated <- pasture_probs_semi <- pasture_probs_rep <- 
  elev_df[!is.na(elev_df$elevation),c("x","y","cell_id", "sr_id")]


# Split to multiple loops do deal with memory issues
##### Integrated #####
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
  sp_df <- sp_df[(sp_df$elevation > (sp_lower - sp_breadth)) & (sp_df$elevation < (sp_upper + sp_breadth)) & !is.na(sp_df$transformed_distance), ]
  # Compute relev
  sp_df$relev <- ((sp_df$elevation - sp_lower)/sp_breadth - bird_data$means_and_sds$relev_offset)/bird_data$means_and_sds$relev_sd
  sp_df$relev2 <- sp_df$relev^2
  # Compute occupancy probabilities
  sp_pc <- pc[i, ]
  sp_forest_logit <- as.numeric((sp_pc$logit_psi_0 - sp_pc$logit_psi_pasture_offset + sp_df$relev * sp_pc$b1_relev_sp + sp_df$relev2 * sp_pc$b1_relev2_sp +
    sp_df$relev * sp_pc$lowland * sp_pc$b1_x_lowland_relev + sp_df$relev2 * sp_pc$lowland * sp_pc$b1_x_lowland_relev2 +
    sp_df$transformed_distance * sp_pc$b5_distance_to_range_sp)[,1])
  sp_pasture_logit <- sp_forest_logit + 2*sp_pc$logit_psi_pasture_offset
  
  sigma_spatial <- sqrt(sp_pc$sigma_sp_cl^2 + sp_pc$sigma_sp_sr^2)
  
  ##### Integrating over cluster and subregion effects #####
  integrate_forest <- function(k){integrate(function(x){dnorm(x, sp_forest_logit[k], sigma_spatial)*boot::inv.logit(x)},
                                          lower = sp_forest_logit[k] - 6*sigma_spatial, 
                                          upper = sp_forest_logit[k] + 6*sigma_spatial,
                                          rel.tol = .001)$value}
  sp_df$forest_prob_integrated <- unlist(lapply(1:nrow(sp_df), integrate_forest))

  integrate_pasture <- function(k){integrate(function(x){dnorm(x, sp_pasture_logit[k], sigma_spatial)*boot::inv.logit(x)},
                                            lower = sp_pasture_logit[k] - 6*sigma_spatial, 
                                            upper = sp_pasture_logit[k] + 6*sigma_spatial,
                                            rel.tol = .001)$value}
  sp_df$pasture_prob_integrated <- unlist(lapply(1:nrow(sp_df), integrate_pasture))
  
  forest_probs_integrated <- merge(forest_probs_integrated, sp_df[,c("cell_id","forest_prob_integrated")], by = "cell_id", all = T)
  forest_probs_integrated[,i+4][is.na(forest_probs_integrated[,i+4])] <- 0
  pasture_probs_integrated <- merge(pasture_probs_integrated, sp_df[,c("cell_id", "pasture_prob_integrated")], by = "cell_id", all = T)
  pasture_probs_integrated[,i+4][is.na(pasture_probs_integrated[,i+4])] <- 0
  names(forest_probs_integrated)[i + 4] <- names(pasture_probs_integrated)[i + 4] <- sp
  
}

saveRDS(forest_probs_integrated, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v6_predictions/iteration_4000/forest_probs_integrated.RDS")
saveRDS(pasture_probs_integrated, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v6_predictions/iteration_4000/pasture_probs_integrated.RDS")

rm(forest_probs_integrated)
rm(pasture_probs_integrated)
gc()

##### Semi-integrated #####
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
  sp_df <- sp_df[(sp_df$elevation > (sp_lower - sp_breadth)) & (sp_df$elevation < (sp_upper + sp_breadth)) & !is.na(sp_df$transformed_distance), ]
  # Compute relev
  sp_df$relev <- ((sp_df$elevation - sp_lower)/sp_breadth - bird_data$means_and_sds$relev_offset)/bird_data$means_and_sds$relev_sd
  sp_df$relev2 <- sp_df$relev^2
  # Compute occupancy probabilities
  sp_pc <- pc[i, ]
  sp_forest_logit <- as.numeric((sp_pc$logit_psi_0 - sp_pc$logit_psi_pasture_offset + sp_df$relev * sp_pc$b1_relev_sp + sp_df$relev2 * sp_pc$b1_relev2_sp +
                                   sp_df$relev * sp_pc$lowland * sp_pc$b1_x_lowland_relev + sp_df$relev2 * sp_pc$lowland * sp_pc$b1_x_lowland_relev2 +
                                   sp_df$transformed_distance * sp_pc$b5_distance_to_range_sp)[,1])
  sp_pasture_logit <- sp_forest_logit + 2*sp_pc$logit_psi_pasture_offset
  
  ##### Integrating over clusters but sampling subregions #####
  sr_effects <- rnorm(max(elev_df$sr_id), 0, sp_pc$sigma_sp_sr)
  integrate_forest <- function(k){integrate(function(x){dnorm(x, sp_forest_logit[k] + sr_effects[sp_df$sr_id[k]], sp_pc$sigma_sp_cl)*boot::inv.logit(x)},
                                            lower = sp_forest_logit[k] + sr_effects[sp_df$sr_id[k]] - 6*sp_pc$sigma_sp_cl, 
                                            upper = sp_forest_logit[k] + sr_effects[sp_df$sr_id[k]] + 6*sp_pc$sigma_sp_cl,
                                            rel.tol = .001)$value}
  sp_df$forest_prob_semi <- unlist(lapply(1:nrow(sp_df), integrate_forest))
  Sys.time() - a
  integrate_pasture <- function(k){integrate(function(x){dnorm(x, sp_pasture_logit[k] + sr_effects[sp_df$sr_id[k]], sp_pc$sigma_sp_cl)*boot::inv.logit(x)},
                                             lower = sp_pasture_logit[k] + sr_effects[sp_df$sr_id[k]] - 6*sp_pc$sigma_sp_cl, 
                                             upper = sp_pasture_logit[k] + sr_effects[sp_df$sr_id[k]] + 6*sp_pc$sigma_sp_cl,
                                             rel.tol = .001)$value}
  sp_df$pasture_prob_semi <- unlist(lapply(1:nrow(sp_df), integrate_pasture))
  
  forest_probs_semi <- merge(forest_probs_semi, sp_df[,c("cell_id","forest_prob_semi")], by = "cell_id", all = T)
  forest_probs_semi[,i+4][is.na(forest_probs_semi[,i+4])] <- 0
  pasture_probs_semi <- merge(pasture_probs_semi, sp_df[,c("cell_id", "pasture_prob_semi")], by = "cell_id", all = T)
  pasture_probs_semi[,i+4][is.na(pasture_probs_semi[,i+4])] <- 0
  names(forest_probs_semi)[i + 4] <- names(pasture_probs_semi)[i + 4] <- sp
}

saveRDS(forest_probs_semi, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v6_predictions/iteration_4000/forest_probs_semi.RDS")
saveRDS(pasture_probs_semi, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v6_predictions/iteration_4000/pasture_probs_semi.RDS")

rm(forest_probs_semi)
rm(pasture_probs_semi)
gc()

##### Reps ######
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
  sp_df <- sp_df[(sp_df$elevation > (sp_lower - sp_breadth)) & (sp_df$elevation < (sp_upper + sp_breadth)) & !is.na(sp_df$transformed_distance), ]
  # Compute relev
  sp_df$relev <- ((sp_df$elevation - sp_lower)/sp_breadth - bird_data$means_and_sds$relev_offset)/bird_data$means_and_sds$relev_sd
  sp_df$relev2 <- sp_df$relev^2
  # Compute occupancy probabilities
  sp_pc <- pc[i, ]
  sp_forest_logit <- as.numeric((sp_pc$logit_psi_0 - sp_pc$logit_psi_pasture_offset + sp_df$relev * sp_pc$b1_relev_sp + sp_df$relev2 * sp_pc$b1_relev2_sp +
                                   sp_df$relev * sp_pc$lowland * sp_pc$b1_x_lowland_relev + sp_df$relev2 * sp_pc$lowland * sp_pc$b1_x_lowland_relev2 +
                                   sp_df$transformed_distance * sp_pc$b5_distance_to_range_sp)[,1])
  sp_pasture_logit <- sp_forest_logit + 2*sp_pc$logit_psi_pasture_offset
  
  sr_effects <- rnorm(max(elev_df$sr_id), 0, sp_pc$sigma_sp_sr)
  cell_effects <- rnorm(max(elev_df$cell_id), 0, sp_pc$sigma_sp_cl)
  
  sp_df$forest_prob_rep <- boot::inv.logit(sp_forest_logit + sr_effects[sp_df$sr_id] + cell_effects[sp_df$cell_id])
  sp_df$pasture_prob_rep <- boot::inv.logit(sp_pasture_logit + sr_effects[sp_df$sr_id] + cell_effects[sp_df$cell_id])
  
  
  forest_probs_rep <- merge(forest_probs_rep, sp_df[,c("cell_id","forest_prob_rep")], by = "cell_id", all = T)
  forest_probs_rep[,i+4][is.na(forest_probs_rep[,i+4])] <- 0
  pasture_probs_rep <- merge(pasture_probs_rep, sp_df[,c("cell_id", "pasture_prob_rep")], by = "cell_id", all = T)
  pasture_probs_rep[,i+4][is.na(pasture_probs_rep[,i+4])] <- 0
  names(forest_probs_rep)[i + 4] <- names(pasture_probs_rep)[i + 4] <- sp
}

saveRDS(forest_probs_rep, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v6_predictions/iteration_4000/forest_probs_rep.RDS")
saveRDS(pasture_probs_rep, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v6_predictions/iteration_4000/pasture_probs_rep.RDS")

rm(forest_probs_rep)
rm(pasture_probs_rep)

gc()
##### Sampled #####


forest_sample <- pasture_sample <- elev_df[!is.na(elev_df$elevation),c("x","y","cell_id", "sr_id")]

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
  sp_df <- sp_df[(sp_df$elevation > (sp_lower - sp_breadth)) & (sp_df$elevation < (sp_upper + sp_breadth)) & !is.na(sp_df$transformed_distance), ]
  # Compute relev
  sp_df$relev <- ((sp_df$elevation - sp_lower)/sp_breadth - bird_data$means_and_sds$relev_offset)/bird_data$means_and_sds$relev_sd
  sp_df$relev2 <- sp_df$relev^2
  # Compute occupancy probabilities
  sp_pc <- pc[i, ]
  sp_forest_logit <- as.numeric((sp_pc$logit_psi_0 - sp_pc$logit_psi_pasture_offset + sp_df$relev * sp_pc$b1_relev_sp + sp_df$relev2 * sp_pc$b1_relev2_sp +
                                   sp_df$relev * sp_pc$lowland * sp_pc$b1_x_lowland_relev + sp_df$relev2 * sp_pc$lowland * sp_pc$b1_x_lowland_relev2 +
                                   sp_df$transformed_distance * sp_pc$b5_distance_to_range_sp)[,1])
  sp_pasture_logit <- sp_forest_logit + 2*sp_pc$logit_psi_pasture_offset
  
  sr_effects <- rnorm(max(elev_df$sr_id), 0, sp_pc$sigma_sp_sr)
  cluster_effects <- matrix(rnorm(16*max(elev_df$cell_id), 0, sp_pc$sigma_sp_cl), ncol = 16)
  sp_df$forest_sample <- sp_df$pasture_sample <- 0
  for(jj in 1:16){
    fpr <- boot::inv.logit(sp_forest_logit + sr_effects[sp_df$sr_id] + cluster_effects[sp_df$cell_id,jj])
    ppr <- boot::inv.logit(sp_pasture_logit + sr_effects[sp_df$sr_id] + cluster_effects[sp_df$cell_id,jj])
    sp_df$forest_sample <- sp_df$forest_sample + (rbinom(length(fpr), 5, fpr) > 0)
    sp_df$pasture_sample <- sp_df$pasture_sample + (rbinom(length(ppr), 5, ppr) > 0)
  }
  forest_sample <- merge(forest_sample, sp_df[,c("cell_id","forest_sample")], by = "cell_id", all = T)
  forest_sample[,i+4][is.na(forest_sample[,i+4])] <- 0
  pasture_sample <- merge(pasture_sample, sp_df[,c("cell_id", "pasture_sample")], by = "cell_id", all = T)
  pasture_sample[,i+4][is.na(pasture_sample[,i+4])] <- 0
  names(forest_sample)[i + 4] <- names(pasture_sample)[i + 4] <- sp
}


saveRDS(forest_sample, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v6_predictions/iteration_4000/forest_sample.RDS")
saveRDS(pasture_sample, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v6_predictions/iteration_4000/pasture_sample.RDS")
