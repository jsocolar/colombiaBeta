# Script to fit gdms to raw data and posterior Z matrix
library(reticulate)
library(gdm)
library(sf)
library(gtools)
library(raster)
library(betapart)

##### Formating covariates #####
# load bird dataset and restrict to points with a detection
birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/birds.RDS")
birds2 <- birds[birds$Q == 1,]

# Extract points and some covariates
birdpts <- birds[!duplicated(birds$point), c("point", "lat", "lon", "elev_ALOS", "pasture")]

# # use gee to get precipitation data, and save result
# use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
# ee <- import("ee")          # Import the Earth Engine library
# ee$Initialize()             # Trigger the authentication
# BioClim <- ee$Image('WORLDCLIM/V1/BIO')
# BioClim_precip <- BioClim$select('bio12')
# geompts <- sapply(1:nrow(birdpts),function(x)ee$Geometry$Point(c(birdpts$lon[x],birdpts$lat[x])))
# geompts <- ee$FeatureCollection(c(unlist(geompts)))
# pts_precip <- BioClim$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
# BioClim_precip <- sapply(c(1:length(pts_precip$features)),function(x)pts_precip$features[[x]]$properties$bio12)
# birdpts$precip <- BioClim_precip
# CHIRP <- raster("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/weather/SA_MAP_1km_CHIRP/SA_MAP_1km_CHIRP.tif")
# birdpts$precip_ceccherini <- extract(CHIRP, cbind(birdpts$lon, birdpts$lat))
# IDEAM <- raster("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/weather/Escenario_Precipitacion_1976_2005/ECC_Prcp_GeoTiff_2011_2040/ECC_Prcp_1976_2005_100K_2015.tif")
# birdpts$precip_IDEAM <- extract(IDEAM, cbind(birdpts$lon, birdpts$lat))
# birdpts[birdpts$precip_ceccherini<1000,]
# birdpts[birdpts$precip<1000,]
# birdpts[grep("CC", birdpts$point),]
# saveRDS(birdpts, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/birdpts_w_precip.RDS")

birdpts <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/birdpts_w_precip.RDS")

# get biogeographic regions for each point
source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/GIS_processing/hydrosheds_extraction.R")
bp2 <- st_as_sf(birdpts, coords = c("lon", "lat"))
st_crs(bp2) <- 4326
east <- st_union(catatumbo, amazon_orinoco)
eAndes <- st_union(east, magdalena_east)
cAndes <- st_union(magdalena_west, cauca_east)
wAndes <- st_union(pacific, cauca_west)
all <- st_union(eAndes, st_union(cAndes, st_union(wAndes, st_union(snsm, pasto))))
if(sum(!st_covered_by(bp2, all, sparse=F)) == 0){print("All points assigned to a region :)")}
birdpts$east <- st_covered_by(bp2, east, sparse = F)
birdpts$range <- "0"
birdpts$range[st_covered_by(bp2, eAndes, sparse = F) | st_covered_by(bp2, pasto, sparse = F)] <- "east_south"
birdpts$range[st_covered_by(bp2, cAndes, sparse = F)] <- "central"
birdpts$range[st_covered_by(bp2, wAndes, sparse = F)] <- "west"
birdpts$range[st_covered_by(bp2, snsm, sparse = F)] <- "snsm"

# Separate to forest points and pasture points
birdpts_f <- birdpts_f_stable <- birdpts[birdpts$pasture == 0, names(birdpts) != "pasture"] # the "stable" version will be a version to come back to when we start modifying birdpts_f later on
birdpts_p <- birdpts_p_stable <- birdpts[birdpts$pasture == 1, names(birdpts) != "pasture"]
birds_f <- birds2[birds2$pasture == 0, c("species", "point", "lat", "lon")]
birds_p <- birds2[birds2$pasture == 1, c("species", "point", "lat", "lon")]
birdpts_p <- birdpts_p[birdpts_p$point %in% birds_p$point,] # exclude points with no detections in the raw pasture data (there are no such points in the forest data)
points_f <- birdpts_f$point
points_p <- birdpts_p$point

# For forest points:
#   distance matrix for montane barriers
    montane_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
    for(i in 1:length(points_f)){
      for(j in 1:length(points_f)){
        if(birdpts_f$east[i] != birdpts_f$east[j]){
          montane_matrix_f[i,j] <- (4100 - max(birdpts_f$elev_ALOS[i], birdpts_f$elev_ALOS[j]))
        }
      }
    }
    montane_barrier_f <- cbind(points_f, data.frame(montane_matrix_f))
    names(montane_barrier_f)[1] <- "point"
#   distance matrix for valley barriers
    valley_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
    for(i in 1:length(points_f)){
      for(j in 1:length(points_f)){
        if(birdpts_f$range[i] != birdpts_f$range[j]){
          valley_matrix_f[i,j] <- (min(birdpts_f$elev_ALOS[i], birdpts_f$elev_ALOS[j]))
        }
      }
    }
    valley_barrier_f <- cbind(points_f, data.frame(valley_matrix_f))
    names(valley_barrier_f)[1] <- "point"

# Pasture points
#   distance matrix for montane barriers
    montane_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
    for(i in 1:length(points_p)){
      for(j in 1:length(points_p)){
        if(birdpts_p$east[i] != birdpts_p$east[j]){
          montane_matrix_p[i,j] <- (4100 - max(birdpts_p$elev_ALOS[i], birdpts_p$elev_ALOS[j]))
        }
      }
    }
    montane_barrier_p <- cbind(points_p, data.frame(montane_matrix_p))
    names(montane_barrier_p)[1] <- "point"
#   distance matrix for valley barriers
    valley_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
    for(i in 1:length(points_p)){
      for(j in 1:length(points_p)){
        if(birdpts_p$range[i] != birdpts_p$range[j]){
          valley_matrix_p[i,j] <- (min(birdpts_p$elev_ALOS[i], birdpts_p$elev_ALOS[j]))
        }
      }
    }
    valley_barrier_p <- cbind(points_p, data.frame(valley_matrix_p))
    names(valley_barrier_p)[1] <- "point"

    
##### Raw data gdm: sorensen #####
# Format site-pair tables for forest and pasture
gdmTab_f <- formatsitepair(birds_f, bioFormat=2, XColumn="lon", YColumn="lat",
                           sppColumn="species", siteColumn="point", 
                           predData = birdpts_f[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")],
                           distPreds = list(montane_barrier_f=montane_barrier_f, valley_barrier_f=valley_barrier_f))
gdmTab_p <- formatsitepair(birds_p, bioFormat=2, XColumn="lon", YColumn="lat",
                           sppColumn="species", siteColumn="point", 
                           predData = birdpts_p[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")],
                           distPreds = list(montane_barrier_p=montane_barrier_p, valley_barrier_p=valley_barrier_p))

# Fit GDM for forest points
#   Observed data
    forest_gdm_obs <- gdm(gdmTab_f, geo = T)
#   bayesian bootstrap resample
    n1 <- paste0(gdmTab_f$s1.xCoord, "_", gdmTab_f$s1.yCoord)
    n2 <- paste0(gdmTab_f$s2.xCoord, "_", gdmTab_f$s2.yCoord)
    s_all <- unique(c(n1,n2))
    s1 <- match(n1, s_all)
    s2 <- match(n2, s_all)
    withweights_f <- gdmTab_f
    forest_gdm_bb <- list()
    for(i in 1:100){
      print(i)
      site_weights_f <- rdirichlet(1, rep(1, length(points_f)))
      pair_weights_f <- site_weights_f[s1] * site_weights_f[s2]
      withweights_f$weights <- pair_weights_f
      forest_gdm_bb[[i]] <- gdm(withweights_f, geo = T)
    }

# Fit GDM for pasture points
#   Observed data
    pasture_gdm_obs <- gdm(gdmTab_p, geo = T)
#   bayesian bootstrap resample
    n1 <- paste0(gdmTab_p$s1.xCoord, "_", gdmTab_p$s1.yCoord)
    n2 <- paste0(gdmTab_p$s2.xCoord, "_", gdmTab_p$s2.yCoord)
    s_all <- unique(c(n1,n2))
    s1 <- match(n1, s_all)
    s2 <- match(n2, s_all)
    withweights_p <- gdmTab_p
    pasture_gdm_bb <- list()
    for(i in 1:100){
      print(i)
      site_weights_p <- rdirichlet(1, rep(1, length(points_p)))
      pair_weights_p <- site_weights_p[s1] * site_weights_p[s2]
      withweights_p$weights <- pair_weights_p
      pasture_gdm_bb[[i]] <- gdm(withweights_p, geo = T)
    }

gdms_raw <- list(forest_obs = forest_gdm_obs, pasture_obs = pasture_gdm_obs, forest_bb = forest_gdm_bb, pasture_bb = pasture_gdm_bb)
saveRDS(gdms_raw, file = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/gdm_outputs/birds/gdms_raw.RDS")
    
##### Modeled data GDM #####
draws <- posterior::as_draws_df(readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/v9_final/draws_thinned_500.RDS"))

bird_data <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data6_package.RDS")
z_info <- data.frame(bird_data$data[8:41])
z_info$point <- birds$point
z_info$species <- birds$species
    
source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/bird_analysis_plotting/get_posterior/get_posterior_z_v6.R")
source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/bird_analysis_plotting/posterior_predictive_checks/discrepancy_functions.R")
Z_rep <- list()
for(i in 1:100){
  print(i)
  iter = 5*i
  Z_probs_include <- get_Z_probs(draws, iter, z_info, spatial_effect = "include")
  Z_rep[[i]] <- data.frame(point = Z_probs_include$point, species = Z_probs_include$species,
                           Z = rbinom(nrow(Z_probs_include), 1, Z_probs_include$Z_prob))
}

saveRDS(Z_rep, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/temporary/v6_Zreps.RDS")
Z_rep <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/temporary/v6_Zreps.RDS")
    
forest_gdm_rep <- pasture_gdm_rep <- list()
for(k in 1:100){ # k indexes the posterior iteration
  print(k)
  birds3 <- birds[Z_rep[[k]]$Z == 1, ]
  birds_f <- birds3[birds3$pasture == 0, c("species", "point", "lat", "lon")]
  birds_p <- birds3[birds3$pasture == 1, c("species", "point", "lat", "lon")]
  birdpts_f <- birdpts_f_stable[birdpts_f_stable$point %in% birds_f$point,] # exclude points with no detections
  birdpts_p <- birdpts_p_stable[birdpts_p_stable$point %in% birds_p$point,] # exclude points with no detections
  # forest
  points_f <- birdpts_f$point
  montane_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
  for(i in 1:length(points_f)){
    for(j in 1:length(points_f)){
      if(birdpts_f$east[i] != birdpts_f$east[j]){
        montane_matrix_f[i,j] <- (4100 - max(birdpts_f$elev_ALOS[i], birdpts_f$elev_ALOS[j]))
      }
    }
  }
  montane_barrier_f <- cbind(points_f, data.frame(montane_matrix_f))
  names(montane_barrier_f)[1] <- "point"
  valley_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
  for(i in 1:length(points_f)){
    for(j in 1:length(points_f)){
      if(birdpts_f$range[i] != birdpts_f$range[j]){
        valley_matrix_f[i,j] <- (min(birdpts_f$elev_ALOS[i], birdpts_f$elev_ALOS[j]))
      }
    }
  }
  valley_barrier_f <- cbind(points_f, data.frame(valley_matrix_f))
  names(valley_barrier_f)[1] <- "point"
  # pasture
  points_p <- birdpts_p$point
  montane_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
  for(i in 1:length(points_p)){
    for(j in 1:length(points_p)){
      if(birdpts_p$east[i] != birdpts_p$east[j]){
        montane_matrix_p[i,j] <- (4100 - max(birdpts_p$elev_ALOS[i], birdpts_p$elev_ALOS[j]))
      }
    }
  }
  montane_barrier_p <- cbind(points_p, data.frame(montane_matrix_p))
  names(montane_barrier_p)[1] <- "point"
  valley_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
  for(i in 1:length(points_p)){
    for(j in 1:length(points_p)){
      if(birdpts_p$range[i] != birdpts_p$range[j]){
        valley_matrix_p[i,j] <- (min(birdpts_p$elev_ALOS[i], birdpts_p$elev_ALOS[j]))
      }
    }
  }
  valley_barrier_p <- cbind(points_p, data.frame(valley_matrix_p))
  names(valley_barrier_p)[1] <- "point"
  
  gdmTab_f <- formatsitepair(birds_f, bioFormat=2, XColumn="lon", YColumn="lat",
                             sppColumn="species", siteColumn="point", 
                             predData = birdpts_f[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")],
                             distPreds = list(montane_barrier_f=montane_barrier_f, valley_barrier_f=valley_barrier_f))
  gdmTab_p <- formatsitepair(birds_p, bioFormat=2, XColumn="lon", YColumn="lat",
                             sppColumn="species", siteColumn="point", 
                             predData = birdpts_p[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")],
                             distPreds = list(montane_barrier_p=montane_barrier_p, valley_barrier_p=valley_barrier_p))
  
  # Bootstrap replicate (one per posterior iteration)
  n1 <- paste0(gdmTab_f$s1.xCoord, "_", gdmTab_f$s1.yCoord)
  n2 <- paste0(gdmTab_f$s2.xCoord, "_", gdmTab_f$s2.yCoord)
  s_all <- unique(c(n1,n2))
  s1 <- match(n1, s_all)
  s2 <- match(n2, s_all)
  withweights_f <- gdmTab_f
  site_weights_f <- rdirichlet(1, rep(1, length(points_f)))
  pair_weights_f <- site_weights_f[s1] * site_weights_f[s2]
  withweights_f$weights <- pair_weights_f
  
  n1 <- paste0(gdmTab_p$s1.xCoord, "_", gdmTab_p$s1.yCoord)
  n2 <- paste0(gdmTab_p$s2.xCoord, "_", gdmTab_p$s2.yCoord)
  s_all <- unique(c(n1,n2))
  s1 <- match(n1, s_all)
  s2 <- match(n2, s_all)
  withweights_p <- gdmTab_p
  site_weights_p <- rdirichlet(1, rep(1, length(points_p)))
  pair_weights_p <- site_weights_p[s1] * site_weights_p[s2]
  withweights_p$weights <- pair_weights_p
  
  forest_gdm_rep[[k]] <- gdm(withweights_f, geo = T)
  pasture_gdm_rep[[k]] <- gdm(withweights_p, geo = T)
}

gdms_modeled <- list(forest_gdm_rep_bb = forest_gdm_rep, pasture_gdm_rep_bb = pasture_gdm_rep)
saveRDS(gdms_modeled, file = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/gdm_outputs/birds/gdms_modeled_v6.RDS")





##### Raw data gdm: turnover (Simpsons) #####
birdpts_f <- birdpts_f_stable <- birdpts[birdpts$pasture == 0, names(birdpts) != "pasture"] # the "stable" version will be a version to come back to when we start modifying birdpts_f later on
birdpts_p <- birdpts_p_stable <- birdpts[birdpts$pasture == 1, names(birdpts) != "pasture"]
birds_f <- birds2[birds2$pasture == 0, c("species", "point", "lat", "lon")]
birds_p <- birds2[birds2$pasture == 1, c("species", "point", "lat", "lon")]
birdpts_p <- birdpts_p[birdpts_p$point %in% birds_p$point,] # exclude points with no detections in the raw pasture data (there are no such points in the forest data)
points_f <- birdpts_f$point
points_p <- birdpts_p$point

# Format site-pair tables for forest and pasture
bird_spp_f <- unique(birds_f$species)
birds_matrix_f <- matrix(data=0, nrow=nrow(birdpts_f), ncol = length(bird_spp_f))
birds_matrix_f[cbind(match(birds_f$point, birdpts_f$point), match(birds_f$species, bird_spp_f))] <- 1
forest_betapart <- beta.pair(birds_matrix_f)
bird_format3_f <- cbind(birdpts_f$point, as.matrix(forest_betapart$beta.sim))
colnames(bird_format3_f)[1] <- "point"

bird_spp_p <- unique(birds_p$species)
birds_matrix_p <- matrix(data=0, nrow=nrow(birdpts_p), ncol = length(bird_spp_p))
birds_matrix_p[cbind(match(birds_p$point, birdpts_p$point), match(birds_p$species, bird_spp_p))] <- 1
pasture_betapart <- beta.pair(birds_matrix_p)
bird_format3_p <- cbind(birdpts_p$point, as.matrix(pasture_betapart$beta.sim))
colnames(bird_format3_p)[1] <- "point"


gdmTab_f <- formatsitepair(bird_format3_f, bioFormat=3, XColumn="lon", YColumn="lat",
                           sppColumn="species", siteColumn="point", 
                           predData = birdpts_f[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")],
                           distPreds = list(montane_barrier_f=montane_barrier_f, valley_barrier_f=valley_barrier_f))
gdmTab_p <- formatsitepair(bird_format3_p, bioFormat=3, XColumn="lon", YColumn="lat",
                           sppColumn="species", siteColumn="point", 
                           predData = birdpts_p[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")],
                           distPreds = list(montane_barrier_p=montane_barrier_p, valley_barrier_p=valley_barrier_p))

# Fit GDM for forest points
#   Observed data
forest_gdm_obs <- gdm(gdmTab_f, geo = T)
#   bayesian bootstrap resample
n1 <- paste0(gdmTab_f$s1.xCoord, "_", gdmTab_f$s1.yCoord)
n2 <- paste0(gdmTab_f$s2.xCoord, "_", gdmTab_f$s2.yCoord)
s_all <- unique(c(n1,n2))
s1 <- match(n1, s_all)
s2 <- match(n2, s_all)
withweights_f <- gdmTab_f
forest_gdm_bb <- list()
for(i in 1:100){
  print(i)
  site_weights_f <- rdirichlet(1, rep(1, length(points_f)))
  pair_weights_f <- site_weights_f[s1] * site_weights_f[s2]
  withweights_f$weights <- pair_weights_f
  forest_gdm_bb[[i]] <- gdm(withweights_f, geo = T)
}

# Fit GDM for pasture points
#   Observed data
pasture_gdm_obs <- gdm(gdmTab_p, geo = T)
#   bayesian bootstrap resample
n1 <- paste0(gdmTab_p$s1.xCoord, "_", gdmTab_p$s1.yCoord)
n2 <- paste0(gdmTab_p$s2.xCoord, "_", gdmTab_p$s2.yCoord)
s_all <- unique(c(n1,n2))
s1 <- match(n1, s_all)
s2 <- match(n2, s_all)
withweights_p <- gdmTab_p
pasture_gdm_bb <- list()
for(i in 1:100){
  print(i)
  site_weights_p <- rdirichlet(1, rep(1, length(points_p)))
  pair_weights_p <- site_weights_p[s1] * site_weights_p[s2]
  withweights_p$weights <- pair_weights_p
  pasture_gdm_bb[[i]] <- gdm(withweights_p, geo = T)
}

gdms_raw_turnover <- list(forest_obs = forest_gdm_obs, pasture_obs = pasture_gdm_obs, forest_bb = forest_gdm_bb, pasture_bb = pasture_gdm_bb)
saveRDS(gdms_raw_turnover, file = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/gdm_outputs/birds/gdms_raw_turnover.RDS")





#################
##### Modeled data: Turnover #####

forest_gdm_bb <-  pasture_gdm_bb <- list()
for(k in 1:100){ # k indexes the posterior iteration
  print(k)
  
  
  birds3 <- birds[Z_rep[[k]]$Z == 1, ]
  birds_f <- birds3[birds3$pasture == 0, c("species", "point", "lat", "lon")]
  birds_p <- birds3[birds3$pasture == 1, c("species", "point", "lat", "lon")]
  birdpts_f <- birdpts_f_stable[birdpts_f_stable$point %in% birds_f$point,] # exclude points with no detections
  birdpts_p <- birdpts_p_stable[birdpts_p_stable$point %in% birds_p$point,] # exclude points with no detections
  # forest
  points_f <- birdpts_f$point
  montane_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
  for(i in 1:length(points_f)){
    for(j in 1:length(points_f)){
      if(birdpts_f$east[i] != birdpts_f$east[j]){
        montane_matrix_f[i,j] <- (4100 - max(birdpts_f$elev_ALOS[i], birdpts_f$elev_ALOS[j]))
      }
    }
  }
  montane_barrier_f <- cbind(points_f, data.frame(montane_matrix_f))
  names(montane_barrier_f)[1] <- "point"
  valley_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
  for(i in 1:length(points_f)){
    for(j in 1:length(points_f)){
      if(birdpts_f$range[i] != birdpts_f$range[j]){
        valley_matrix_f[i,j] <- (min(birdpts_f$elev_ALOS[i], birdpts_f$elev_ALOS[j]))
      }
    }
  }
  valley_barrier_f <- cbind(points_f, data.frame(valley_matrix_f))
  names(valley_barrier_f)[1] <- "point"
  # pasture
  points_p <- birdpts_p$point
  montane_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
  for(i in 1:length(points_p)){
    for(j in 1:length(points_p)){
      if(birdpts_p$east[i] != birdpts_p$east[j]){
        montane_matrix_p[i,j] <- (4100 - max(birdpts_p$elev_ALOS[i], birdpts_p$elev_ALOS[j]))
      }
    }
  }
  montane_barrier_p <- cbind(points_p, data.frame(montane_matrix_p))
  names(montane_barrier_p)[1] <- "point"
  valley_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
  for(i in 1:length(points_p)){
    for(j in 1:length(points_p)){
      if(birdpts_p$range[i] != birdpts_p$range[j]){
        valley_matrix_p[i,j] <- (min(birdpts_p$elev_ALOS[i], birdpts_p$elev_ALOS[j]))
      }
    }
  }
  valley_barrier_p <- cbind(points_p, data.frame(valley_matrix_p))
  names(valley_barrier_p)[1] <- "point"
  
  
  
  # Format site-pair tables for forest and pasture
  bird_spp_f <- unique(birds_f$species)
  birds_matrix_f <- matrix(data=0, nrow=nrow(birdpts_f), ncol = length(bird_spp_f))
  birds_matrix_f[cbind(match(birds_f$point, birdpts_f$point), match(birds_f$species, bird_spp_f))] <- 1
  forest_betapart <- beta.pair(birds_matrix_f)
  bird_format3_f <- cbind(birdpts_f$point, as.matrix(forest_betapart$beta.sim))
  colnames(bird_format3_f)[1] <- "point"
  
  bird_spp_p <- unique(birds_p$species)
  birds_matrix_p <- matrix(data=0, nrow=nrow(birdpts_p), ncol = length(bird_spp_p))
  birds_matrix_p[cbind(match(birds_p$point, birdpts_p$point), match(birds_p$species, bird_spp_p))] <- 1
  pasture_betapart <- beta.pair(birds_matrix_p)
  bird_format3_p <- cbind(birdpts_p$point, as.matrix(pasture_betapart$beta.sim))
  colnames(bird_format3_p)[1] <- "point"
  
  
  gdmTab_f <- formatsitepair(bird_format3_f, bioFormat=3, XColumn="lon", YColumn="lat",
                             sppColumn="species", siteColumn="point", 
                             predData = birdpts_f[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")],
                             distPreds = list(montane_barrier_f=montane_barrier_f, valley_barrier_f=valley_barrier_f))
  gdmTab_p <- formatsitepair(bird_format3_p, bioFormat=3, XColumn="lon", YColumn="lat",
                             sppColumn="species", siteColumn="point", 
                             predData = birdpts_p[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")],
                             distPreds = list(montane_barrier_p=montane_barrier_p, valley_barrier_p=valley_barrier_p))
  
  
  
  # Fit GDM for forest points
  #   Observed data
  forest_gdm_obs <- gdm(gdmTab_f, geo = T)
  #   bayesian bootstrap resample
  n1 <- paste0(gdmTab_f$s1.xCoord, "_", gdmTab_f$s1.yCoord)
  n2 <- paste0(gdmTab_f$s2.xCoord, "_", gdmTab_f$s2.yCoord)
  s_all <- unique(c(n1,n2))
  s1 <- match(n1, s_all)
  s2 <- match(n2, s_all)
  withweights_f <- gdmTab_f


    site_weights_f <- rdirichlet(1, rep(1, length(points_f)))
    pair_weights_f <- site_weights_f[s1] * site_weights_f[s2]
    withweights_f$weights <- pair_weights_f
    forest_gdm_bb[[k]] <- gdm(withweights_f, geo = T)
  
  # Fit GDM for pasture points
  #   Observed data
  pasture_gdm_obs <- gdm(gdmTab_p, geo = T)
  #   bayesian bootstrap resample
  n1 <- paste0(gdmTab_p$s1.xCoord, "_", gdmTab_p$s1.yCoord)
  n2 <- paste0(gdmTab_p$s2.xCoord, "_", gdmTab_p$s2.yCoord)
  s_all <- unique(c(n1,n2))
  s1 <- match(n1, s_all)
  s2 <- match(n2, s_all)
  withweights_p <- gdmTab_p
  

    site_weights_p <- rdirichlet(1, rep(1, length(points_p)))
    pair_weights_p <- site_weights_p[s1] * site_weights_p[s2]
    withweights_p$weights <- pair_weights_p
    pasture_gdm_bb[[k]] <- gdm(withweights_p, geo = T)
  
}

gdms_modeled_simpsons <- list(forest_gdm_rep_bb = forest_gdm_bb, pasture_gdm_rep_bb = pasture_gdm_bb)
saveRDS(gdms_modeled_simpsons, file = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/gdm_outputs/birds/gdms_modeled_simpsons_v6.RDS")







##########
##### Raw data gdm: Raup-crick #####
birdpts_f <- birdpts_f_stable <- birdpts[birdpts$pasture == 0, names(birdpts) != "pasture"] # the "stable" version will be a version to come back to when we start modifying birdpts_f later on
birdpts_p <- birdpts_p_stable <- birdpts[birdpts$pasture == 1, names(birdpts) != "pasture"]
birds_f <- birds2[birds2$pasture == 0, c("species", "point", "lat", "lon")]
birds_p <- birds2[birds2$pasture == 1, c("species", "point", "lat", "lon")]
birdpts_p <- birdpts_p[birdpts_p$point %in% birds_p$point,] # exclude points with no detections in the raw pasture data (there are no such points in the forest data)
points_f <- birdpts_f$point
points_p <- birdpts_p$point

# Format site-pair tables for forest and pasture
bird_spp_f <- unique(birds_f$species)
birds_matrix_f <- matrix(data=0, nrow=nrow(birdpts_f), ncol = length(bird_spp_f))
birds_matrix_f[cbind(match(birds_f$point, birdpts_f$point), match(birds_f$species, bird_spp_f))] <- 1
forest_raup <- vegan::raupcrick(birds_matrix_f)
bird_format3_f <- cbind(birdpts_f$point, as.matrix(forest_raup))
colnames(bird_format3_f)[1] <- "point"

bird_spp_p <- unique(birds_p$species)
birds_matrix_p <- matrix(data=0, nrow=nrow(birdpts_p), ncol = length(bird_spp_p))
birds_matrix_p[cbind(match(birds_p$point, birdpts_p$point), match(birds_p$species, bird_spp_p))] <- 1
pasture_raup <- vegan::raupcrick(birds_matrix_p)
bird_format3_p <- cbind(birdpts_p$point, as.matrix(pasture_raup))
colnames(bird_format3_p)[1] <- "point"


gdmTab_f <- formatsitepair(bird_format3_f, bioFormat=3, XColumn="lon", YColumn="lat",
                           sppColumn="species", siteColumn="point", 
                           predData = birdpts_f[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")],
                           distPreds = list(montane_barrier_f=montane_barrier_f, valley_barrier_f=valley_barrier_f))
gdmTab_p <- formatsitepair(bird_format3_p, bioFormat=3, XColumn="lon", YColumn="lat",
                           sppColumn="species", siteColumn="point", 
                           predData = birdpts_p[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")],
                           distPreds = list(montane_barrier_p=montane_barrier_p, valley_barrier_p=valley_barrier_p))

# Fit GDM for forest points
#   Observed data
forest_gdm_obs <- gdm(gdmTab_f, geo = T)
#   bayesian bootstrap resample
n1 <- paste0(gdmTab_f$s1.xCoord, "_", gdmTab_f$s1.yCoord)
n2 <- paste0(gdmTab_f$s2.xCoord, "_", gdmTab_f$s2.yCoord)
s_all <- unique(c(n1,n2))
s1 <- match(n1, s_all)
s2 <- match(n2, s_all)
withweights_f <- gdmTab_f
forest_gdm_bb <- list()
for(i in 1:100){
  print(i)
  site_weights_f <- rdirichlet(1, rep(1, length(points_f)))
  pair_weights_f <- site_weights_f[s1] * site_weights_f[s2]
  withweights_f$weights <- pair_weights_f
  forest_gdm_bb[[i]] <- gdm(withweights_f, geo = T)
}

# Fit GDM for pasture points
#   Observed data
pasture_gdm_obs <- gdm(gdmTab_p, geo = T)
#   bayesian bootstrap resample
n1 <- paste0(gdmTab_p$s1.xCoord, "_", gdmTab_p$s1.yCoord)
n2 <- paste0(gdmTab_p$s2.xCoord, "_", gdmTab_p$s2.yCoord)
s_all <- unique(c(n1,n2))
s1 <- match(n1, s_all)
s2 <- match(n2, s_all)
withweights_p <- gdmTab_p
pasture_gdm_bb <- list()
for(i in 1:100){
  print(i)
  site_weights_p <- rdirichlet(1, rep(1, length(points_p)))
  pair_weights_p <- site_weights_p[s1] * site_weights_p[s2]
  withweights_p$weights <- pair_weights_p
  pasture_gdm_bb[[i]] <- gdm(withweights_p, geo = T)
}

gdms_raw_raup <- list(forest_obs = forest_gdm_obs, pasture_obs = pasture_gdm_obs, forest_bb = forest_gdm_bb, pasture_bb = pasture_gdm_bb)
saveRDS(gdms_raw_raup, file = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/gdm_outputs/birds/gdms_raw_raup.RDS")


#############
##########
##### Raw data gdm: Raup-crick; no dist #####
# Fit GDM for forest points
#   Observed data
forest_gdm_obs <- gdm(gdmTab_f, geo = F)
#   bayesian bootstrap resample
n1 <- paste0(gdmTab_f$s1.xCoord, "_", gdmTab_f$s1.yCoord)
n2 <- paste0(gdmTab_f$s2.xCoord, "_", gdmTab_f$s2.yCoord)
s_all <- unique(c(n1,n2))
s1 <- match(n1, s_all)
s2 <- match(n2, s_all)
withweights_f <- gdmTab_f
forest_gdm_bb <- list()
for(i in 1:100){
  print(i)
  site_weights_f <- rdirichlet(1, rep(1, length(points_f)))
  pair_weights_f <- site_weights_f[s1] * site_weights_f[s2]
  withweights_f$weights <- pair_weights_f
  forest_gdm_bb[[i]] <- gdm(withweights_f, geo = F)
}

# Fit GDM for pasture points
#   Observed data
pasture_gdm_obs <- gdm(gdmTab_p, geo = F)
#   bayesian bootstrap resample
n1 <- paste0(gdmTab_p$s1.xCoord, "_", gdmTab_p$s1.yCoord)
n2 <- paste0(gdmTab_p$s2.xCoord, "_", gdmTab_p$s2.yCoord)
s_all <- unique(c(n1,n2))
s1 <- match(n1, s_all)
s2 <- match(n2, s_all)
withweights_p <- gdmTab_p
pasture_gdm_bb <- list()
for(i in 1:100){
  print(i)
  site_weights_p <- rdirichlet(1, rep(1, length(points_p)))
  pair_weights_p <- site_weights_p[s1] * site_weights_p[s2]
  withweights_p$weights <- pair_weights_p
  pasture_gdm_bb[[i]] <- gdm(withweights_p, geo = F)
}

gdms_raw_raup <- list(forest_obs = forest_gdm_obs, pasture_obs = pasture_gdm_obs, forest_bb = forest_gdm_bb, pasture_bb = pasture_gdm_bb)
saveRDS(gdms_raw_raup, file = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/gdm_outputs/birds/gdms_raw_raup.RDS")