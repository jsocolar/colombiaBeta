##### Predict occupancy across Colombia #####

# Still need to do: assign each cell a subregion value on a 20-km grid and simulate within this and/or
# integrate over this.



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

