# Script to fit gdms to raw data and posterior Z matrix
library(reticulate)
library(gdm)
library(sf)
library(gtools)

##### Raw data GDM ####
# load bird dataset and restrict to points with a detection
birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/birds.RDS")
birds2 <- birds[birds$Q == 1,]

# Extract points and some covariates
birdpts <- birds[!duplicated(birds$point), c("point", "lat", "lon", "elev_ALOS", "pasture")]

# use gee to get precipitation data
use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()             # Trigger the authentication
BioClim <- ee$Image('WORLDCLIM/V1/BIO')
BioClim_precip <- BioClim$select('bio12')
geompts <- sapply(1:nrow(birdpts),function(x)ee$Geometry$Point(c(birdpts$lon[x],birdpts$lat[x])))
geompts <- ee$FeatureCollection(c(unlist(geompts)))
pts_precip <- BioClim$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
BioClim_precip <- sapply(c(1:length(pts_precip$features)),function(x)pts_precip$features[[x]]$properties$bio12)
birdpts$precip <- BioClim_precip

# get biogeographic regions for mountain and valley barriers
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
birdpts_f <- birdpts_f_stable <- birdpts[birdpts$pasture == 0, names(birdpts) != "pasture"] # the "stable" version will be a version to come back to when we start modifying birdpts_f later on
birdpts_p <- birdpts_p_stable <- birdpts[birdpts$pasture == 1, names(birdpts) != "pasture"]
birds_f <- birds2[birds2$pasture == 0, c("species", "point", "lat", "lon")]
birds_p <- birds2[birds2$pasture == 1, c("species", "point", "lat", "lon")]
birdpts_p <- birdpts_p[birdpts_p$point %in% birds_p$point,] # exclude points with no detections in the raw pasture data (there are no such points in the forest data)

  # get biogeographic predictors for forest points
points_f <- birdpts_f$point
montane_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
for(i in 1:length(points_f)){
  for(j in 1:length(points_f)){
    if(birdpts_f$east[i] != birdpts_f$east[j]){
      montane_matrix_f[i,j] <- (4100 - max(birdpts_f$elev_ALOS[i], birdpts_f$elev_ALOS[j]))/2000
    }
  }
}
montane_barrier_f <- cbind(points_f, data.frame(montane_matrix_f))
names(montane_barrier_f)[1] <- "point"

valley_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
for(i in 1:length(points_f)){
  for(j in 1:length(points_f)){
    if(birdpts_f$range[i] != birdpts_f$range[j]){
      valley_matrix_f[i,j] <- (min(birdpts_f$elev_ALOS[i], birdpts_f$elev_ALOS[j]))/2000
    }
  }
}
valley_barrier_f <- cbind(points_f, data.frame(valley_matrix_f))
names(valley_barrier_f)[1] <- "point"

  # get biogeographic predictors for pasture
points_p <- birdpts_p$point

montane_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
for(i in 1:length(points_p)){
  for(j in 1:length(points_p)){
    if(birdpts_p$east[i] != birdpts_p$east[j]){
      montane_matrix_p[i,j] <- (4100 - max(birdpts_p$elev_ALOS[i], birdpts_p$elev_ALOS[j]))/2000
    }
  }
}
montane_barrier_p <- cbind(points_p, data.frame(montane_matrix_p))
names(montane_barrier_p)[1] <- "point"

valley_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
for(i in 1:length(points_p)){
  for(j in 1:length(points_p)){
    if(birdpts_p$range[i] != birdpts_p$range[j]){
      valley_matrix_p[i,j] <- (min(birdpts_p$elev_ALOS[i], birdpts_p$elev_ALOS[j]))/2000
    }
  }
}
valley_barrier_p <- cbind(points_p, data.frame(valley_matrix_p))
names(valley_barrier_p)[1] <- "point"

# Format site-pair tables for forest and pasture
gdmTab_f <- formatsitepair(birds_f, bioFormat=2, XColumn="lon", YColumn="lat",
                           sppColumn="species", siteColumn="point", 
                           predData = birdpts_f[,c("point", "lat", "lon", "elev_ALOS", "precip")],
                           distPreds = list(montane_barrier_f=montane_barrier_f, valley_barrier_f=valley_barrier_f))
gdmTab_p <- formatsitepair(birds_p, bioFormat=2, XColumn="lon", YColumn="lat",
                           sppColumn="species", siteColumn="point", 
                           predData = birdpts_p[,c("point", "lat", "lon", "elev_ALOS", "precip")],
                           distPreds = list(montane_barrier_p=montane_barrier_p, valley_barrier_p=valley_barrier_p))

# Fit GDM for forest points
forest_gdm <- gdm(gdmTab_f, geo = T)
forest_splines <- isplineExtract(forest_gdm)
  # bayesian bootstrap resample
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
pasture_gdm <- gdm(gdmTab_p, geo = T)
pasture_splines <- isplineExtract(pasture_gdm)
  # bayesian bootstrap resample
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

# Plot results
# spatial distance
pred <- "Geographic"
ymax <- 0
for(k in 1:100){
  forest_splines <- isplineExtract(forest_gdm_bb[[k]])
  pasture_splines <- isplineExtract(pasture_gdm_bb[[k]])
  ymax <- max(ymax, max(forest_splines$y[,pred]))
  ymax <- max(ymax, max(pasture_splines$y[,pred]))
}
plot(forest_splines$x[,pred], forest_splines$y[,pred],
     type="n", xlab="Geographic distance", ylab="Partial ecological distance",
     ylim = c(0,ymax))
for(i in 1:100){
  f_spline_rep <- isplineExtract(forest_gdm_bb[[i]])
  lines(f_spline_rep$x[,pred], f_spline_rep$y[,pred],
        type="l", col = "slateblue3", lwd=.5)
  p_spline_rep <- isplineExtract(pasture_gdm_bb[[i]])
  lines(p_spline_rep$x[,pred], p_spline_rep$y[,pred],
        type="l", col = "lightsalmon3", lwd=.5)
}
lines(pasture_splines$x[,pred], pasture_splines$y[,pred], 
      type="l", col = "black", lwd=3)
lines(forest_splines$x[,pred], forest_splines$y[,pred],
      type="l", col = "slateblue3", lwd=3)
# elevation
pred <- "elev_ALOS"
ymax <- 0
for(k in 1:100){
  forest_splines <- isplineExtract(forest_gdm_bb[[k]])
  pasture_splines <- isplineExtract(pasture_gdm_bb[[k]])
  ymax <- max(ymax, max(forest_splines$y[,pred]))
  ymax <- max(ymax, max(pasture_splines$y[,pred]))
}
plot(forest_splines$x[,pred], forest_splines$y[,pred],
     type="n", xlab="Elevation", ylab="Partial ecological distance", 
     ylim = c(0,ymax))
for(i in 1:100){
  f_spline_rep <- isplineExtract(forest_gdm_bb[[i]])
  lines(f_spline_rep$x[,pred], f_spline_rep$y[,pred],
        type="l", col = "slateblue3", lwd=.5)
  p_spline_rep <- isplineExtract(pasture_gdm_bb[[i]])
  lines(p_spline_rep$x[,pred], p_spline_rep$y[,pred],
        type="l", col = "lightsalmon3", lwd=.5)
}
lines(pasture_splines$x[,pred], pasture_splines$y[,pred], 
      type="l", col = "black", lwd=3)
lines(forest_splines$x[,pred], forest_splines$y[,pred],
      type="l", col = "black", lwd=3)
# precipitation
pred <- "precip"
ymax <- 0
for(k in 1:100){
  forest_splines <- isplineExtract(forest_gdm_bb[[k]])
  pasture_splines <- isplineExtract(pasture_gdm_bb[[k]])
  ymax <- max(ymax, max(forest_splines$y[,pred]))
  ymax <- max(ymax, max(pasture_splines$y[,pred]))
}
plot(forest_splines$x[,pred], forest_splines$y[,pred],
     type="n", xlab="Annual precipitation", ylab="Partial ecological distance",
     ylim = c(0,ymax))
for(i in 1:100){
  f_spline_rep <- isplineExtract(forest_gdm_bb[[i]])
  lines(f_spline_rep$x[,pred], f_spline_rep$y[,pred],
        type="l", col = "slateblue3", lwd=.5)
  p_spline_rep <- isplineExtract(pasture_gdm_bb[[i]])
  lines(p_spline_rep$x[,pred], p_spline_rep$y[,pred],
        type="l", col = "lightsalmon3", lwd=.5)
}
lines(pasture_splines$x[,pred], pasture_splines$y[,pred], 
      type="l", col = "black", lwd=3)
lines(forest_splines$x[,pred], forest_splines$y[,pred],
      type="l", col = "slateblue3", lwd=3)
# mountain barrier
pred <- "matrix_1"
ymax <- 0
for(k in 1:100){
  forest_splines <- isplineExtract(forest_gdm_bb[[k]])
  pasture_splines <- isplineExtract(pasture_gdm_bb[[k]])
  ymax <- max(ymax, max(forest_splines$y[,pred]))
  ymax <- max(ymax, max(pasture_splines$y[,pred]))
}
plot(forest_splines$x[,pred], forest_splines$y[,pred],
     type="n", xlab="Mountain barrier", ylab="Partial ecological distance",
     ylim = c(0,ymax))
for(i in 1:100){
  f_spline_rep <- isplineExtract(forest_gdm_bb[[i]])
  lines(f_spline_rep$x[,pred], f_spline_rep$y[,pred],
        type="l", col = "slateblue3", lwd=.5)
  p_spline_rep <- isplineExtract(pasture_gdm_bb[[i]])
  lines(p_spline_rep$x[,pred], p_spline_rep$y[,pred],
        type="l", col = "lightsalmon3", lwd=.5)
}
lines(pasture_splines$x[,pred], pasture_splines$y[,pred], 
      type="l", col = "black", lwd=3)
lines(forest_splines$x[,pred], forest_splines$y[,pred],
      type="l", col = "slateblue3", lwd=3)
# valley barrier
pred <- "matrix_2"
ymax <- 0
for(k in 1:100){
  forest_splines <- isplineExtract(forest_gdm_bb[[k]])
  pasture_splines <- isplineExtract(pasture_gdm_bb[[k]])
  ymax <- max(ymax, max(forest_splines$y[,pred]))
  ymax <- max(ymax, max(pasture_splines$y[,pred]))
}
plot(forest_splines$x[,pred], forest_splines$y[,pred],
     type="n", xlab="Valley barrier", ylab="Partial ecological distance",
     ylim = c(0,ymax))
for(i in 1:100){
  f_spline_rep <- isplineExtract(forest_gdm_bb[[i]])
  lines(f_spline_rep$x[,pred], f_spline_rep$y[,pred],
        type="l", col = "slateblue3", lwd=.5)
  p_spline_rep <- isplineExtract(pasture_gdm_bb[[i]])
  lines(p_spline_rep$x[,pred], p_spline_rep$y[,pred],
        type="l", col = "lightsalmon3", lwd=.5)
}
lines(pasture_splines$x[,pred], pasture_splines$y[,pred], 
      type="l", col = "black", lwd=3)
lines(forest_splines$x[,pred], forest_splines$y[,pred],
      type="l", col = "slateblue3", lwd=3)


##### Posterior distribution of Z #####
# v5 <- cmdstanr::read_cmdstan_csv("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/v5_first_run/occupancy_v5_threads-202012282018-1-261afe.csv")
# draws <- posterior::as_draws_df(v5$post_warmup_draws[1:2000,,])
# saveRDS(draws, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/temporary/v5_draws.RDS")
draws <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/temporary/v5_draws.RDS")

bird_data <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data4_package.RDS")
z_info <- data.frame(bird_data$data[8:43])
z_info$point <- birds$point
z_info$species <- birds$species

source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/bird_analysis_plotting/get_posterior_z.R")
source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/bird_analysis_plotting/posterior_predictive_checks/discrepancy_functions.R")
Z_rep <- list()
for(i in 1:100){
  print(i)
  iter = 20*i
  Z_probs_include <- get_Z_probs(draws, iter, z_info, cluster_effect = "include")
  Z_rep[[i]] <- data.frame(point = Z_probs_include$point, species = Z_probs_include$species, 
                                      Z = rbinom(nrow(Z_probs_include), 1, Z_probs_include$Z_prob))
}

saveRDS(Z_rep, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/temporary/v5_Zreps.RDS")
Z_rep <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/temporary/v5_Zreps.RDS")

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
        montane_matrix_f[i,j] <- (4100 - max(birdpts_f$elev_ALOS[i], birdpts_f$elev_ALOS[j]))/2000
      }
    }
  }
  montane_barrier_f <- cbind(points_f, data.frame(montane_matrix_f))
  names(montane_barrier_f)[1] <- "point"
  valley_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
  for(i in 1:length(points_f)){
    for(j in 1:length(points_f)){
      if(birdpts_f$range[i] != birdpts_f$range[j]){
        valley_matrix_f[i,j] <- (min(birdpts_f$elev_ALOS[i], birdpts_f$elev_ALOS[j]))/2000
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
        montane_matrix_p[i,j] <- (4100 - max(birdpts_p$elev_ALOS[i], birdpts_p$elev_ALOS[j]))/2000
      }
    }
  }
  montane_barrier_p <- cbind(points_p, data.frame(montane_matrix_p))
  names(montane_barrier_p)[1] <- "point"
  valley_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
  for(i in 1:length(points_p)){
    for(j in 1:length(points_p)){
      if(birdpts_p$range[i] != birdpts_p$range[j]){
        valley_matrix_p[i,j] <- (min(birdpts_p$elev_ALOS[i], birdpts_p$elev_ALOS[j]))/2000
      }
    }
  }
  valley_barrier_p <- cbind(points_p, data.frame(valley_matrix_p))
  names(valley_barrier_p)[1] <- "point"
  
  gdmTab_f <- formatsitepair(birds_f, bioFormat=2, XColumn="lon", YColumn="lat",
                             sppColumn="species", siteColumn="point", 
                             predData = birdpts_f[,c("point", "lat", "lon", "elev_ALOS", "precip")],
                             distPreds = list(montane_barrier_f=montane_barrier_f, valley_barrier_f=valley_barrier_f))
  gdmTab_p <- formatsitepair(birds_p, bioFormat=2, XColumn="lon", YColumn="lat",
                             sppColumn="species", siteColumn="point", 
                             predData = birdpts_p[,c("point", "lat", "lon", "elev_ALOS", "precip")],
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


pred <- "Geographic"
ymax <- 0
for(k in 1:100){
  forest_splines <- isplineExtract(forest_gdm_rep[[k]])
  pasture_splines <- isplineExtract(pasture_gdm_rep[[k]])
  ymax <- max(ymax, max(forest_splines$y[,pred]))
  ymax <- max(ymax, max(pasture_splines$y[,pred]))
}
plot(forest_splines$x[,pred], forest_splines$y[,pred],
     type="n", xlab="Geographic distance", ylab="Partial ecological distance",
     ylim = c(0, ymax))
for(k in 1:100){
  forest_splines <- isplineExtract(forest_gdm_rep[[k]])
  pasture_splines <- isplineExtract(pasture_gdm_rep[[k]])
  lines(pasture_splines$x[,pred], pasture_splines$y[,pred], 
        type="l", col = "lightsalmon3", lwd=.5)
  lines(forest_splines$x[,pred], forest_splines$y[,pred],
        type="l", col = "slateblue3", lwd=.5)
}

pred <- "elev_ALOS"
ymax <- 0
for(k in 1:100){
  forest_splines <- isplineExtract(forest_gdm_rep[[k]])
  pasture_splines <- isplineExtract(pasture_gdm_rep[[k]])
  ymax <- max(ymax, max(forest_splines$y[,pred]))
  ymax <- max(ymax, max(pasture_splines$y[,pred]))
}
plot(forest_splines$x[,pred], forest_splines$y[,pred],
     type="n", xlab="Elevation", ylab="Partial ecological distance",
     ylim = c(0, ymax))
for(k in 1:100){
  forest_splines <- isplineExtract(forest_gdm_rep[[k]])
  pasture_splines <- isplineExtract(pasture_gdm_rep[[k]])
  lines(pasture_splines$x[,pred], pasture_splines$y[,pred], 
        type="l", col = "lightsalmon3", lwd=.5)
  lines(forest_splines$x[,pred], forest_splines$y[,pred],
        type="l", col = "slateblue3", lwd=.5)
}

pred <- "precip"
ymax <- 0
for(k in 1:100){
  forest_splines <- isplineExtract(forest_gdm_rep[[k]])
  pasture_splines <- isplineExtract(pasture_gdm_rep[[k]])
  ymax <- max(ymax, max(forest_splines$y[,pred]))
  ymax <- max(ymax, max(pasture_splines$y[,pred]))
}
plot(forest_splines$x[,pred], forest_splines$y[,pred],
     type="n", xlab="Annual precipitation", ylab="Partial ecological distance",
     ylim = c(0, ymax))
for(k in 1:100){
  forest_splines <- isplineExtract(forest_gdm_rep[[k]])
  pasture_splines <- isplineExtract(pasture_gdm_rep[[k]])
  lines(pasture_splines$x[,pred], pasture_splines$y[,pred], 
        type="l", col = "lightsalmon3", lwd=.5)
  lines(forest_splines$x[,pred], forest_splines$y[,pred],
        type="l", col = "slateblue3", lwd=.5)
}

pred <- "matrix_1"
ymax <- 0
for(k in 1:100){
  forest_splines <- isplineExtract(forest_gdm_rep[[k]])
  pasture_splines <- isplineExtract(pasture_gdm_rep[[k]])
  ymax <- max(ymax, max(forest_splines$y[,pred]))
  ymax <- max(ymax, max(pasture_splines$y[,pred]))
}
plot(forest_splines$x[,pred], forest_splines$y[,pred],
     type="n", xlab="Mountain barrier", ylab="Partial ecological distance",
     ylim = c(0, ymax))
for(k in 1:100){
  forest_splines <- isplineExtract(forest_gdm_rep[[k]])
  pasture_splines <- isplineExtract(pasture_gdm_rep[[k]])
  lines(pasture_splines$x[,pred], pasture_splines$y[,pred], 
        type="l", col = "lightsalmon3", lwd=.5)
  lines(forest_splines$x[,pred], forest_splines$y[,pred],
        type="l", col = "slateblue3", lwd=.5)
}

pred <- "matrix_2"
ymax <- 0
for(k in 1:100){
  forest_splines <- isplineExtract(forest_gdm_rep[[k]])
  pasture_splines <- isplineExtract(pasture_gdm_rep[[k]])
  ymax <- max(ymax, max(forest_splines$y[,pred]))
  ymax <- max(ymax, max(pasture_splines$y[,pred]))
}
plot(forest_splines$x[,pred], forest_splines$y[,pred],
     type="n", xlab="Valley barrier", ylab="Partial ecological distance",
     ylim = c(0, ymax))
for(k in 1:100){
  forest_splines <- isplineExtract(forest_gdm_rep[[k]])
  pasture_splines <- isplineExtract(pasture_gdm_rep[[k]])
  lines(pasture_splines$x[,pred], pasture_splines$y[,pred], 
        type="l", col = "lightsalmon3", lwd=.5)
  lines(forest_splines$x[,pred], forest_splines$y[,pred],
        type="l", col = "slateblue3", lwd=.5)
}
