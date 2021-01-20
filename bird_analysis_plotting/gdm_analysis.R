library(reticulate)
library(gdm)
library(sf)

birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/birds.RDS")
birds2 <- birds[!duplicated(paste0(birds$species, "__", birds$point)) & birds$Q == 1,]

birdpts <- birds[!duplicated(birds$point), c("point", "lat", "lon", "elev_ALOS", "pasture")]

# get ALOS elevations
use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()             # Trigger the authentication

BioClim <- ee$Image('WORLDCLIM/V1/BIO')
BioClim_precip <- BioClim$select('bio12')
# Featurecollection of point geometries
geompts <- sapply(1:nrow(birdpts),function(x)ee$Geometry$Point(c(birdpts$lon[x],birdpts$lat[x])))
geompts <- ee$FeatureCollection(c(unlist(geompts)))
# Extract ALOS elevations for all points - combine into dataframe
pts_precip <- BioClim$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
BioClim_precip <- sapply(c(1:length(pts_precip$features)),function(x)pts_precip$features[[x]]$properties$bio12)
birdpts$precip <- BioClim_precip

source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/GIS_processing/hydrosheds_extraction.R")
bp2 <- st_as_sf(birdpts, coords = c("lon", "lat"))
st_crs(bp2) <- 4326

east <- st_union(catatumbo, amazon_orinoco)
eAndes <- st_union(east, magdalena_east)
cAndes <- st_union(magdalena_west, cauca_east)
wAndes <- st_union(pacific, cauca_west)

all <- st_union(eAndes, st_union(cAndes, st_union(wAndes, st_union(snsm, pasto))))

bp2[!st_covered_by(bp2, all, sparse=F),]

birdpts$east <- st_covered_by(bp2, east, sparse = F)
birdpts$range <- "0"
birdpts$range[st_covered_by(bp2, eAndes, sparse = F) | st_covered_by(bp2, pasto, sparse = F)] <- "east_south"
birdpts$range[st_covered_by(bp2, cAndes, sparse = F)] <- "central"
birdpts$range[st_covered_by(bp2, wAndes, sparse = F)] <- "west"
birdpts$range[st_covered_by(bp2, snsm, sparse = F)] <- "snsm"

birdpts_f <- birdpts[birdpts$pasture == 0, names(birdpts) != "pasture"]
birdpts_p <- birdpts[birdpts$pasture == 1, names(birdpts) != "pasture"]

birds_f <- birds2[birds2$pasture == 0, c("species", "point", "lat", "lon")]
birds_p <- birds2[birds2$pasture == 1, c("species", "point", "lat", "lon")]

birdpts_p <- birdpts_p[birdpts_p$point %in% birds_p$point,] # exclude points with no detections


# distPreds
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



######

gdmTab_f <- formatsitepair(birds_f, bioFormat=2, XColumn="lon", YColumn="lat",
                           sppColumn="species", siteColumn="point", 
                           predData = birdpts_f[,c("point", "lat", "lon", "elev_ALOS", "precip")],
                           distPreds = list(montane_barrier_f=montane_barrier_f, valley_barrier_f=valley_barrier_f))

gdmTab_p <- formatsitepair(birds_p, bioFormat=2, XColumn="lon", YColumn="lat",
                           sppColumn="species", siteColumn="point", 
                           predData = birdpts_p[,c("point", "lat", "lon", "elev_ALOS", "precip")],
                           distPreds = list(montane_barrier_p=montane_barrier_p, valley_barrier_p=valley_barrier_p))

# distPreds will be a list of distance matricies:
# montane barrier separation, valley barrier separation, interactions



forest_gdm <- gdm(gdmTab_f, geo = T)
pasture_gdm <- gdm(gdmTab_p, geo = T)
plot(pasture_gdm)
plot(forest_gdm)

forest_splines <- isplineExtract(forest_gdm)
pasture_splines <- isplineExtract(pasture_gdm)

plot(forest_splines$x[,"Geographic"], forest_splines$y[,"Geographic"],
     type="n", xlab="Geographic distance", ylab="Partial ecological distance")
lines(pasture_splines$x[,"Geographic"], pasture_splines$y[,"Geographic"], 
      type="l", col = "lightsalmon3", lwd=3)
lines(forest_splines$x[,"Geographic"], forest_splines$y[,"Geographic"],
      type="l", col = "slateblue3", lwd=3)

plot(forest_splines$x[,"elev_ALOS"], forest_splines$y[,"elev_ALOS"],
     type="n", xlab="Elevation", ylab="Partial ecological distance")
lines(pasture_splines$x[,"elev_ALOS"], pasture_splines$y[,"elev_ALOS"], 
      type="l", col = "lightsalmon3", lwd=3)
lines(forest_splines$x[,"elev_ALOS"], forest_splines$y[,"elev_ALOS"],
      type="l", col = "slateblue3", lwd=3)

plot(forest_splines$x[,"precip"], forest_splines$y[,"precip"],
     type="n", xlab="Annual precip", ylab="Partial ecological distance")
lines(pasture_splines$x[,"precip"], pasture_splines$y[,"precip"], 
      type="l", col = "lightsalmon3", lwd=3)
lines(forest_splines$x[,"precip"], forest_splines$y[,"precip"],
      type="l", col = "slateblue3", lwd=3)

plot(forest_splines$x[,"matrix_1"], forest_splines$y[,"matrix_1"],
     type="n", xlab="Mountain barrier", ylab="Partial ecological distance")
lines(pasture_splines$x[,"matrix_1"], pasture_splines$y[,"matrix_1"], 
      type="l", col = "lightsalmon3", lwd=3)
lines(forest_splines$x[,"matrix_1"], forest_splines$y[,"matrix_1"],
      type="l", col = "slateblue3", lwd=3)

plot(forest_splines$x[,"matrix_2"], forest_splines$y[,"matrix_2"],
     type="n", xlab="Valley barrier", ylab="Partial ecological distance")
lines(pasture_splines$x[,"matrix_2"], pasture_splines$y[,"matrix_2"], 
      type="l", col = "lightsalmon3", lwd=3)
lines(forest_splines$x[,"matrix_2"], forest_splines$y[,"matrix_2"],
      type="l", col = "slateblue3", lwd=3)


######### Posterior distribution of Z #########
# v5 <- cmdstanr::read_cmdstan_csv("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/v5_first_run/occupancy_v5_threads-202012282018-1-261afe.csv")
# draws <- posterior::as_draws_df(v5$post_warmup_draws[1:2000,,])
# saveRDS(draws, "/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/temporary/v5_draws.RDS")
draws <- readRDS("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/temporary/v5_draws.RDS")

bird_data <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data4_package.RDS")
birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/birds.RDS")
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
