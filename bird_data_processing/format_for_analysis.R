# This script takes previously ingested and standardized data on birds, traits, ranges, points, and visits and produces a unified data object for analysis

##### Script dependencies: combined_bird_maps.R, bird_import_and_cleaning.R, elevations_prep_and_exploration.R, points_formatting.R, migratory_dates.R, species_covariate_formatting.R #####

`%ni%` <- Negate(`%in%`)

# Get formatted bird surveys object
bird_surveys <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_surveys_current.RDS')

# Get a matrix (actually a dataframe) of the distance from each species range to each sampling point.
# This is just for checking whether to include a species in analysis, not for creating the distance-
# to-range covariate.
point_distances <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/point_distances/point_distances_biogeographic_clip_ayerbe.RDS")

# Get a list of species with ranges overlapping our points
include_species <- vector()
for(i in 1:nrow(point_distances)){
  include_species[i] <- sum(point_distances[i,] > 0) != ncol(point_distances)
}

# subset the point distances to just the species that overlap at least one point
point_distances_include <- point_distances[include_species, ]

# Get a list of the species that overlap at least one point in underscore format
species_list <- gsub(" ", "_", row.names(point_distances_include))

# Confirm that all detected species are in the species list
which(bird_surveys$species_names %ni% species_list)

# Extract detection array and pad with zeros for all never-detected species in species_list
det_array <- bird_surveys$detection_array[,1:4,]
det_array_padded <- abind::abind(det_array, array(data = 0, dim = c(848, 4, length(species_list) - dim(det_array)[3])), along = 3)

# Species names for det_array_padded
species_names <- c(bird_surveys$species_names, species_list[species_list %ni% bird_surveys$species_names])

# Create flattened data object, where each species-point gets its own row
nrow_flat <- sum(point_distances == 0)

flattened_data <- as.data.frame(matrix(data = 0, nrow = nrow_flat, ncol = 6))
names(flattened_data) <- c('species', 'point', 'v1', 'v2', 'v3', 'v4')

counter <- 0
for(i in 1:length(species_list)){
  print(i)
  species <- species_list[i]
  det_array_ind <- which(species_names == species)
  for(j in 1:length(bird_surveys$point_names)){
    point <- bird_surveys$point_names[j]
    if(point_distances_include[i, which(names(point_distances_include) == point)] == 0){
      counter <- counter + 1
      flattened_data$species[counter] <- species
      flattened_data$point[counter] <- point
      if(det_array_ind <= dim(det_array)[3]){
        flattened_data[counter, 3:6] <- det_array[j, , det_array_ind]
      }
    }
  }
}
# saveRDS(flattened_data, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/flattened_data.RDS")
flattened_data <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/flattened_data.RDS")

# Column for whether the species is ever detected at the point
flattened_data$Q <- as.numeric(rowSums(flattened_data[,3:6], na.rm = T) > 0)

# Read in point covariate information and merge with flattened_data
all_pts <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Points/all_pts.RDS")
fd <- merge(flattened_data, all_pts, by.x = "point", by.y = "point", all.x = T)
fd$v4[fd$nv %in% c(2,3)] <- NA
fd$v3[fd$nv == 2] <- NA

# Read in species-trait covariate information and merge with fd
traits <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/traits/traits.RDS")
# Confirm that we have trait covariates for every species of interest
all(flattened_data$species %in% traits$latin_underscore)
fd <- merge(fd, traits, by.x = "species", by.y = "latin_underscore", all.x = T)
# Compute species-standardized elevations
fd$elev_sp_standard <- (fd$elev_ALOS - fd$lower)/(fd$upper - fd$lower)


# Read in migratory date information and merge with fd
mig_dates <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/traits/mig_dates.RDS")
mig_dates$species <- gsub(" ", "_", mig_dates$latin)
fd <- merge(fd, mig_dates, by.x = "species", by.y = "species", all.x = T)

# For each row, determine if the species is within the appropriate migratory period
fd$in_date_range <- 0
for(i in 1:nrow(fd)){
  oday1 <- fd$oday1[i]
  if(is.na(fd$oday2[i])){oday2 <- oday1}else{oday2 <- fd$oday2[i]}
  if(is.na(fd$oday3[i])){oday3 <- oday2}else{oday3 <- fd$oday3[i]}
  if(is.na(fd$oday4[i])){oday4 <- oday3}else{oday4 <- fd$oday4[i]}
  
  if(is.na(fd$start1[i])){
    fd$in_date_range[i] <- 1
  }else if((fd$start1[i] < fd$end1[i]) & oday1 >= fd$start1[i] & oday1 <= fd$end1[i]){
    fd$in_date_range[i] <- 1
  }else if((fd$start1[i] < fd$end1[i]) & oday4 >= fd$start1[i] & oday4 <= fd$end1[i]){
    fd$in_date_range[i] <- 1
  }else if((fd$start1[i] > fd$end1[i]) & (oday1 >= fd$start1[i] | oday1 <= fd$end1[i])){
    fd$in_date_range[i] <- 1
  }else if((fd$start1[i] > fd$end1[i]) & (oday4 >= fd$start1[i] | oday4 <= fd$end1[i])){
    fd$in_date_range[i] <- 1
  }else if(fd$start2[i] != fd$start1[i]){
    if(fd$start2[i] > fd$end2[i]){
      stop()
    }else if(oday1 >= fd$start2[i] & oday1 <= fd$end2[i]){
      fd$in_date_range[i] <- 1
    }else if(oday4 >= fd$start2[i] & oday4 <= fd$end2[i]){
      fd$in_date_range[i] <- 1
    }
  }
}

sum(fd$in_date_range)
# Confirm that no detections fall outside the migratory date range.
sum(fd$Q == 1 & fd$in_date_range == 0)

# Remove superfluous columns from fd
flattened_data_full <- fd[, names(fd) %ni% c("birds", "beetles", "habitat", "other", "latin.x", "latin.y")]
# saveRDS(flattened_data_full, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/flattened_data_full.RDS")

flattened_data_full <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/flattened_data_full.RDS")
# Look at statistics of species-standardized elevations at species-points with a detection  
fdq <- flattened_data_full[flattened_data_full$Q == 1,]
max(fdq$elev_sp_standard)
min(fdq$elev_sp_standard)
hist(fdq$elev_sp_standard)

a <- seq(-1,2,.2)
nq <- vector()
nall <- vector()
for(i in 2:length(a)){
  nq[i-1] <- sum(fdq$elev_sp_standard > a[i-1] & fdq$elev_sp_standard <= a[i])
  nall[i-1] <- sum(flattened_data_full$elev_sp_standard > a[i-1] & flattened_data_full$elev_sp_standard <= a[i])
}

plot(boot::logit(nq/nall) ~ seq(-.9, 1.9, .2))
plot((nq/nall) ~ seq(-.9, 1.9, .2))

plot(nall ~ seq(-.9, 1.9, .2))
min(nall)


# Subset to include only species-points that are in the date range and in an elevational range of (-1, 2)
bird_data_trimmed <- flattened_data_full[flattened_data_full$in_date_range == 1 & flattened_data_full$elev_sp_standard > -1 & flattened_data_full$elev_sp_standard < 2, ]
# saveRDS(bird_data_trimmed, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_data_trimmed.RDS")
bird_data_trimmed <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_data_trimmed.RDS")
# Examine statistics of final dataset
nrow(bird_data_trimmed)
sum(bird_data_trimmed$Q)
mean(bird_data_trimmed$Q)


############
vscale <- function(x){return(as.vector(scale(x)))}

birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_data_trimmed.RDS")
birds$sp_cl <- paste(birds$species, birds$cluster, sep = "__")
birds$sp_sr <- paste(birds$species, birds$cluster, sep = "__")
birds$elev_median <- rowMeans(cbind(birds$lower, birds$upper))
birds$det_data <- as.matrix(birds[,c("v1", "v2", "v3", "v4")])
birds$det_data[is.na(birds$det_data)] <- -1
birds$obsSM <- matrix(as.numeric(c(birds$obs1 == "SCM", birds$obs2 == "SCM", birds$obs3 == "SCM", birds$obs4 == "SCM")), ncol = 4)
birds$obsSM[is.na(birds$obsSM)] <- 0
birds$obsDE <- matrix(as.numeric(c(birds$obs1 == "DPE", birds$obs2 == "DPE", birds$obs3 == "DPE", birds$obs4 == "DPE")), ncol = 4)
birds$obsDE[is.na(birds$obsDE)] <- 0
birds$obsJG <- matrix(as.numeric(c(birds$obs1 == "JJG", birds$obs2 == "JJG", birds$obs3 == "JJG", birds$obs4 == "JJG")), ncol = 4)
birds$obsJG[is.na(birds$obsJG)] <- 0
birds$time <- matrix((scale(c(birds$hps1, birds$hps2, birds$hps3, birds$hps4))), ncol = 4)
birds$time[is.na(birds$time)] <- 0  # these zeros correspond to visits that don't actually exist
birds$relev <- (birds$elev_sp_standard -.5)/sd(birds$elev_sp_standard)
birds$relev2 <- birds$relev^2
birds$elev_median_scaled <- vscale(birds$elev_median)
birds$elev_breadth_scaled <- vscale(birds$elev_breadth)
birds$log_mass_scaled <- vscale(log(birds$BodyMass.Value))
birds$lowland <- as.numeric(birds$lower == 0)
birds$mountain_limited <- birds$east_only | birds$west_only
birds$valley_limited <- birds$snsm_only | birds$wandes_absent | birds$eandes_absent

#### Add distance-from-range covariate ####
library(sf)
library(ggplot2)
source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/GIS_processing/get_mainland.R")
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
mainland <- st_transform(mainland, AEAstring)
mainland_inward <- st_buffer(mainland, -7000)

ayerbe_list_updated <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/ayerbe_maps/ayerbe_list_updated.RDS')
birds$distance_from_range <- 0
sp_list <- unique(birds$species)
for(i in 1:length(sp_list)){  # this takes ~ 5 minutes
  print(i)
  sp <- sp_list[i]
  ayerbe <- st_union(ayerbe_list_updated[[gsub("_", " ", sp)]])
  points <- st_as_sf(birds[birds$species == sp, c("lon", "lat")], coords = c("lon", "lat"))
  st_crs(points) <- st_crs("WGS84")
  points <- st_transform(points, st_crs(ayerbe))
  range_linestring <- st_cast(ayerbe, "MULTILINESTRING")
  range_linestring_cropped <- st_intersection(mainland_inward, range_linestring)
  if(nrow(range_linestring_cropped) > 0){
    inside <- as.numeric(as.numeric(st_distance(points, ayerbe)) == 0)
    distance_inside <- -1 * st_distance(points, range_linestring_cropped)
    distance_outside <- st_distance(points, range_linestring)
    distances <- as.numeric(inside * distance_inside + (1 - inside) * distance_outside)
  }else{
    birds$distance_from_range[birds$species == sp] <- -2e-06
  }
}

birds[which(birds$distance_from_range > 160000),]
birds$distance_from_range_scaled <- birds$distance_from_range/sd(birds$distance_from_range)

# Explore good functional form for distance covariate:
hist(birds$distance_from_range)
hist(birds$distance_from_range[birds$Q==1])
n0 <- n1 <- vector()
for(i in 1:40){
  n0[i] <- sum(birds$distance_from_range > 16000*(i-31) & birds$distance_from_range > 16000*(i-30) & birds$Q == 0)
  n1[i] <- sum(birds$distance_from_range > 16000*(i-31) & birds$distance_from_range > 16000*(i-30) & birds$Q == 1)
}
range_data <- data.frame(prop_det = n1/(n1+n0), distance = 16000*(-30:9)+8000)
plot(prop_det ~ distance, data = range_data[range_data$distance>0,])
plot(prop_det ~ distance, data = range_data)

range_data$logit_prop <- boot::logit(range_data$prop_det)
range_data$logit_prop[range_data$logit_prop == -Inf] <- boot::logit(min(range_data$prop_det[range_data$prop_det > 0])*.5)
plot(logit_prop ~ distance, data = range_data)

range_data$distance_scaled <- range_data$distance/sd(birds$distance_from_range)
plot(logit_prop ~ distance_scaled, data = range_data)

range_data$dist_trans <- boot::inv.logit(range_data$distance_scaled*4)
plot(logit_prop ~ dist_trans, data = range_data)

fit_qual <- vector()
for(i in 1:400){
  range_data$dist_trans <- boot::inv.logit(range_data$distance_scaled*(1 + i/10))
  fit_qual[i] <- summary(lm(logit_prop ~ dist_trans, data = range_data))$adj.r.squared
}
plot(fit_qual)
1 + which(fit_qual == max(fit_qual))/10
range_data$dist_trans <- boot::inv.logit(range_data$distance_scaled*11.3)
plot(logit_prop ~ dist_trans, data = range_data)

birds$distance_from_range_scaled2 <- boot::inv.logit(birds$distance_from_range_scaled*11.3)

birds$sp_sr <- paste0(birds$species, "__", birds$subregion)
birds$sp_obs1 <- paste0(birds$species, "__", birds$obs1)
birds$sp_obs2 <- paste0(birds$species, "__", birds$obs2)
birds$sp_obs3 <- paste0(birds$species, "__", birds$obs3)
birds$sp_obs4 <- paste0(birds$species, "__", birds$obs4)

sp_obs_columns <- c(birds$sp_obs1, birds$sp_obs2, birds$sp_obs3, birds$sp_obs4)
sp_obs_columns[grep("__NA", sp_obs_columns)] <- NA
sp_obs_matrix <- matrix(as.integer(as.factor(sp_obs_columns)), ncol=4)
sp_obs_matrix[is.na(sp_obs_matrix)] <- 0
birds$sp_obs_matrix <- sp_obs_matrix

# saveRDS(birds, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/birds.RDS")
birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/birds.RDS")
ec <- c(-1,1)
###########
bird_stan_data1 <- list(
  # Grainsize for reduce_sum
  grainsize = 1,
  
  # Dimensions
  n_spCl = length(unique(birds$sp_cl)),
  n_sp = length(unique(birds$species)),
  n_fam = length(unique(birds$Family)),
  n_tot = nrow(birds),
  n_visit_max = max(birds$nv),
  
  # Detection matrix
  det_data = det_data,
  
  # Q and nv
  Q = birds$Q,
  nv = birds$nv,
  
  # Random effect IDs
  id_spCl = as.numeric(as.factor(birds$sp_cl)),
  id_sp = as.numeric(as.factor(birds$species)),
  id_fam = as.numeric(as.factor(birds$Family)),
  
  # Covariates
  relev = birds$relev,
  relev2 = birds$relev2,
  pasture = birds$pasture,
  eastOnly = birds$east_only,
  westOnly = birds$west_only,
  snsmOnly = birds$snsm_only,
  notWandes = birds$wandes_absent,
  notEandes = birds$eandes_absent,
  lowland = birds$lowland,
  elevMedian = birds$elev_median_scaled,
  elevBreadth = birds$elev_breadth_scaled,
  forestPresent = birds$forest_present,
  forestSpecialist = birds$forest_specialist,
  tfSpecialist = birds$tf_specialist,
  dryForestPresent = birds$dry_forest_present,
  floodDrySpecialist = birds$flood_dry_specialist,
  floodSpecialist = birds$floodplain_specialist,
  aridPresent = birds$arid_present,
  migratory = as.numeric(!is.na(birds$start1)),
  mass = birds$log_mass_scaled,
  dietInvert = as.numeric(birds$Diet.5Cat == "Invertebrate"),
  dietCarn = as.numeric(birds$Diet.5Cat == "VertFishScav"),
  dietFruitNect = as.numeric(birds$Diet.5Cat == "FruiNect"),
  dietGran = as.numeric(birds$Diet.5Cat == "PlantSeed"),
  time = time,
  obsSM = obsSM,
  obsJG = obsJG,
  obsDE = obsDE)

bird_standata1_means_and_sds <- list(time_mean = time_mean, time_sd = time_sd,
                      relev_offset = relev_offset, relev_sd = relev_sd,
                      elev_median_mean = elev_median_mean, elev_median_sd = elev_median_sd,
                      elev_breadth_mean = elev_breadth_mean, elev_breadth_sd = elev_breadth_sd,
                      log_mass_mean = log_mass_mean, log_mass_sd = log_mass_sd)

bird_stan_data1_package <- list(data = bird_stan_data1,
                           means_and_sds = bird_standata1_means_and_sds)

saveRDS(bird_stan_data1_package, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data1_package.RDS")





#####
bird_stan_data4 <- list(
  # Grainsize for reduce_sum
  grainsize = 1,
  
  # Dimensions
  n_spCl = length(unique(birds$sp_cl)),
  n_sp = length(unique(birds$species)),
  n_fam = length(unique(birds$Family)),
  n_tot = nrow(birds),
  n_visit_max = max(birds$nv),
  
  # Detection matrix
  det_data = det_data,
  
  # Q and nv
  Q = birds$Q,
  nv = birds$nv,
  
  # Random effect IDs
  id_spCl = as.numeric(as.factor(birds$sp_cl)),
  id_sp = as.numeric(as.factor(birds$species)),
  id_fam = as.numeric(as.factor(birds$Family)),
  
  # Covariates
  relev = birds$relev,
  relev2 = birds$relev2,
  pasture = birds$pasture,
  lowland = birds$lowland,
  mountain_barrier = birds$mountain_limited,
  valley_barrier = birds$valley_limited,
  eastOnly = birds$east_only,
  westOnly = birds$west_only,
  snsmOnly = birds$snsm_only,
  notEandes = birds$eandes_absent,
  notWandes = birds$wandes_absent,
  elevMedian = birds$elev_median_scaled,
  elevBreadth = birds$elev_breadth_scaled,
  forestPresent = birds$forest_present,
  forestSpecialist = birds$forest_specialist,
  tfSpecialist = birds$tf_specialist,
  dryForestPresent = birds$dry_forest_present,
  floodDrySpecialist = birds$flood_dry_specialist,
  floodSpecialist = birds$floodplain_specialist,
  aridPresent = birds$arid_present,
  migratory = as.numeric(!is.na(birds$start1)),
  mass = birds$log_mass_scaled,
  dietInvert = as.numeric(birds$Diet.5Cat == "Invertebrate"),
  dietCarn = as.numeric(birds$Diet.5Cat == "VertFishScav"),
  dietFruitNect = as.numeric(birds$Diet.5Cat == "FruiNect"),
  dietGran = as.numeric(birds$Diet.5Cat == "PlantSeed"),
  distance_to_range = as.vector(birds$distance_from_range_scaled2),
  time = time,
  obsSM = obsSM,
  obsJG = obsJG,
  obsDE = obsDE)
bird_standata4_means_and_sds <- list(time_mean = time_mean, time_sd = time_sd,
                                     relev_offset = relev_offset, relev_sd = relev_sd,
                                     elev_median_mean = elev_median_mean, elev_median_sd = elev_median_sd,
                                     elev_breadth_mean = elev_breadth_mean, elev_breadth_sd = elev_breadth_sd,
                                     log_mass_mean = log_mass_mean, log_mass_sd = log_mass_sd, distance_to_range_offset = 0,
                                     distance_to_range_sd = sd(birds$distance_from_range), distance_to_range_logit_rescale = 11.3)
bird_stan_data4_package_prelim <- list(data = bird_stan_data4,
                                means_and_sds = bird_standata4_means_and_sds)

saveRDS(bird_stan_data4_package_prelim, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data4_package_prelim.RDS")


#####
bird_stan_data6 <- list(
  # Grainsize for reduce_sum
    grainsize = 1,
  # Dimensions
    n_spCl = length(unique(birds$sp_cl)),
    n_spSr = length(unique(birds$sp_sr)),
    n_sp = length(unique(birds$species)),
    n_fam = length(unique(birds$Family)),
    n_spObs = length(unique(as.vector(birds$sp_obs_matrix[!is.na(as.vector(birds$sp_obs_matrix))]))),
    n_tot = nrow(birds),
    n_visit_max = max(birds$nv),
  # Detection matrix
    det_data = birds$det_data,
  # Q and nv
    Q = birds$Q,
    nv = birds$nv,
  # Random effect IDs
    id_spCl = as.numeric(as.factor(birds$sp_cl)),
    id_spSr = as.numeric(as.factor(birds$sp_sr)),
    id_sp = as.numeric(as.factor(birds$species)),
    id_fam = as.numeric(as.factor(birds$Family)),
    id_spObs = birds$sp_obs_matrix,
  # Covariates
    relev = birds$relev,
    relev2 = birds$relev2,
    lowland = ec[birds$lowland+1],
    pasture = ec[birds$pasture+1],
    mountain_barrier = ec[birds$mountain_limited+1],
    valley_barrier = ec[birds$valley_limited+1],
    elevMedian = birds$elev_median_scaled,
    elevBreadth = birds$elev_breadth_scaled,
    forestPresent = ec[birds$forest_present+1],
    forestSpecialist = ec[birds$forest_specialist+1],
    tfSpecialist = ec[birds$tf_specialist+1],
    dryForestPresent = ec[birds$dry_forest_present+1],
    floodDrySpecialist = ec[birds$flood_dry_specialist+1],
    aridPresent = ec[birds$arid_present+1],
    migratory = ec[as.numeric(!is.na(birds$start1))+1],
    mass = birds$log_mass_scaled,
    dietInvert = ec[as.numeric(birds$Diet.5Cat == "Invertebrate")+1],
    dietCarn = ec[as.numeric(birds$Diet.5Cat == "VertFishScav")+1],
    dietFruitNect = ec[as.numeric(birds$Diet.5Cat == "FruiNect")+1],
    dietGran = ec[as.numeric(birds$Diet.5Cat == "PlantSeed")+1],
    distance_to_range = as.vector(birds$distance_from_range_scaled2),
    time = birds$time,
    obsSM = matrix(ec[birds$obsSM+1], ncol=4),
    obsJG = matrix(ec[birds$obsJG+1], ncol=4),
    obsDE = matrix(ec[birds$obsDE+1], ncol=4)
)

bird_standata6_means_and_sds <- list(time_mean = mean(c(birds$hps1, birds$hps2, birds$hps3, birds$hps4), na.rm = T), 
                                     time_sd =  sd(c(birds$hps1, birds$hps2, birds$hps3, birds$hps4), na.rm = T),
                                     relev_offset = .5, relev_sd = sd(birds$elev_sp_standard),
                                     elev_median_mean =  mean(birds$elev_median), elev_median_sd = sd(birds$elev_median),
                                     elev_breadth_mean = mean(birds$elev_breadth), elev_breadth_sd = sd(birds$elev_breadth),
                                     log_mass_mean = mean(log(birds$BodyMass.Value)), log_mass_sd = sd(log(birds$BodyMass.Value)), 
                                     distance_to_range_offset = 0, distance_to_range_sd = sd(birds$distance_from_range), 
                                     distance_to_range_logit_rescale = 11.3)
bird_stan_data6_package <- list(data = bird_stan_data6,
                                       means_and_sds = bird_standata6_means_and_sds)
saveRDS(bird_stan_data6_package, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data6_package.RDS")




##### Stan data 7 #####
ec <- c(-1,1)

grt <- function(predictor){
  return(list(u = unique(predictor), id = match(predictor, unique(predictor))))
}
grt_b <- function(predictor){
  if(sum(predictor != 1 & predictor != -1) != 0){stop('predictor must be effects coded')}
  predictor[predictor == 1] <- 2
  predictor[predictor == -1] <- 1
  return(predictor)
}

grt_m <- function(predictor){
  v <- as.vector(predictor)
  u <- unique(v)
  out <- matrix(match(v, u), ncol=ncol(predictor))
  return(out)
}


bird_stan_data7 <- list(
  # Grainsize for reduce_sum
  grainsize = 1,
  # Dimensions
    # Random effects
      n_spCl = length(unique(birds$sp_cl)),
      n_spSr = length(unique(birds$sp_sr)),
      n_sp = length(unique(birds$species)),
      n_fam = length(unique(birds$Family)),
      n_spObs = length(unique(as.vector(birds$sp_obs_matrix[!is.na(as.vector(birds$sp_obs_matrix))]))),
    # Dataset size
      n_tot = nrow(birds),
      n_visit_max = max(birds$nv),
    # Unique values sizes
      n_relev = length(unique(birds$relev)),
      n_relev2 = length(unique(birds$relev2)),
      n_lowland_x_relev = length(unique(ec[birds$lowland+1] * birds$relev)),
      n_lowland_x_relev2 = length(unique(ec[birds$lowland+1] * birds$relev2)),
      n_elevMedian = length(unique(birds$elev_median_scaled)),
      n_elevBreadth = length(unique(birds$elev_breadth_scaled)),
      n_mass = length(unique(birds$log_mass_scaled)),
      n_time = length(unique(as.vector(birds$time))),
              # skip distance to range because it doesn't compress
      n_elevMedian_x_forestPresent = length(unique(birds$elev_median_scaled * ec[birds$forest_present + 1])),
      n_elevMedian_x_forestSpecialist = length(unique(birds$elev_median_scaled * ec[birds$forest_specialist + 1])),
      n_elevMedian_x_pasture = length(unique(birds$elev_median_scaled * ec[birds$pasture + 1])),
      n_elevBreadth_x_pasture = length(unique(birds$elev_breadth_scaled * ec[birds$pasture + 1])),
      n_mass_x_pasture = length(unique(birds$log_mass_scaled * ec[birds$pasture + 1])),
      n_elevMedian_x_forestPresent_x_pasture = length(unique(birds$elev_median_scaled * ec[birds$forest_present + 1] * ec[birds$pasture + 1])),
      n_elevMedian_x_forestSpecialist_x_pasture = length(unique(birds$elev_median_scaled * ec[birds$forest_specialist + 1] * ec[birds$pasture + 1])),
      n_time_x_elev = length(unique(as.vector(sweep(birds$time, MARGIN = 1, birds$elev_median_scaled, FUN = `*`)))),
  # Integer data matrix
    main_data = data.frame(
      # Detection data
        det_data = birds$det_data,  # 1-4
      # Q and nv
        Q = birds$Q,               # 5
        nv = birds$nv,             # 6
      # Random effect IDs
        id_spCl = as.numeric(as.factor(birds$sp_cl)),     # 7
        id_spSr = as.numeric(as.factor(birds$sp_sr)),     # 8
        id_sp = as.numeric(as.factor(birds$species)),     # 9
        id_fam = as.numeric(as.factor(birds$Family)),     # 10
        id_spObs = birds$sp_obs_matrix,                   # 11-14
      # Covariate indices
        lowland = grt_b(ec[birds$lowland+1]),             # 15
        pasture = grt_b(ec[birds$pasture+1]),                          # 16
        mountain_barrier = grt_b(ec[birds$mountain_limited+1]),        # 17
        valley_barrier = grt_b(ec[birds$valley_limited+1]),            # 18
        forestPresent = grt_b(ec[birds$forest_present+1]),             # 19
        forestSpecialist = grt_b(ec[birds$forest_specialist+1]),       # 20
        tfSpecialist = grt_b(ec[birds$tf_specialist+1]),               # 21
        dryForestPresent = grt_b(ec[birds$dry_forest_present+1]),      # 22
        floodDrySpecialist = grt_b(ec[birds$flood_dry_specialist+1]),  # 23
        aridPresent = grt_b(ec[birds$arid_present+1]),                 # 24
        migratory = grt_b(ec[as.numeric(!is.na(birds$start1))+1]),     # 25
        dietInvert = grt_b(ec[as.numeric(birds$Diet.5Cat == "Invertebrate")+1]),    # 26
        dietCarn = grt_b(ec[as.numeric(birds$Diet.5Cat == "VertFishScav")+1]),      # 27
        dietFruitNect = grt_b(ec[as.numeric(birds$Diet.5Cat == "FruiNect")+1]),     # 28
        dietGran = grt_b(ec[as.numeric(birds$Diet.5Cat == "PlantSeed")+1]),         # 29
        mountainBarrier_x_pasture = grt_b(ec[birds$mountain_limited + 1] * ec[birds$pasture + 1]),        # 30
        valleyBarrier_x_pasture = grt_b(ec[birds$valley_limited + 1] * ec[birds$pasture + 1]),            # 31
        forestPresent_x_pasture = grt_b(ec[birds$forest_present + 1] * ec[birds$pasture + 1]),            # 32
        forestSpecialist_x_pasture = grt_b(ec[birds$forest_specialist + 1] * ec[birds$pasture + 1]),      # 33
        tfSpecialist_x_pasture = grt_b(ec[birds$tf_specialist + 1] * ec[birds$pasture + 1]),              # 34
        dryForestPresent_x_pasture = grt_b(ec[birds$dry_forest_present + 1] * ec[birds$pasture + 1]),     # 35
        floodDrySpecialist_x_pasture = grt_b(ec[birds$flood_dry_specialist + 1] * ec[birds$pasture + 1]), # 36
        aridPresent_x_pasture = grt_b(ec[birds$arid_present + 1] * ec[birds$pasture + 1]),                # 37
        migratory_x_pasture = grt_b(ec[as.numeric(!is.na(birds$start1))+1] * ec[birds$pasture + 1]),      # 38
        dietInvert_x_pasture = grt_b(ec[as.numeric(birds$Diet.5Cat == "Invertebrate")+1] * ec[birds$pasture + 1]),   # 39
        dietCarn_x_pasture = grt_b(ec[as.numeric(birds$Diet.5Cat == "VertFishScav")+1] * ec[birds$pasture + 1]),     # 40
        dietFruitNect_x_pasture = grt_b(ec[as.numeric(birds$Diet.5Cat == "FruiNect")+1] * ec[birds$pasture + 1]),    # 41
        dietGran_x_pasture = grt_b(ec[as.numeric(birds$Diet.5Cat == "PlantSeed")+1] * ec[birds$pasture + 1]),        # 42
        obsSM = matrix(grt_b(ec[birds$obsSM+1]), ncol=4),   # 43-46
        obsJG = matrix(grt_b(ec[birds$obsJG+1]), ncol=4),   # 47-50
        obsDE = matrix(grt_b(ec[birds$obsDE+1]), ncol=4)    # 51-54
  ),
  # Continuous distance-to-range (does not compress well)
    distance_to_range = as.vector(birds$distance_from_range_scaled2),
  # continuous covariates
    relev = birds$relev,
    relev2 = birds$relev2,
    lowland_x_relev = ec[birds$lowland+1] * birds$relev,
    lowland_x_relev2 = ec[birds$lowland+1] * birds$relev2,
    elevMedian = birds$elev_median_scaled,
    elevBreadth = birds$elev_breadth_scaled,
    mass = birds$log_mass_scaled,
    elevMedian_x_forestPresent = birds$elev_median_scaled * ec[birds$forest_present + 1],
    elevMedian_x_forestSpecialist = birds$elev_median_scaled * ec[birds$forest_specialist + 1],
    elevMedian_x_pasture = birds$elev_median_scaled * ec[birds$pasture + 1],
    elevBreadth_x_pasture = birds$elev_breadth_scaled * ec[birds$pasture + 1],
    mass_x_pasture = birds$log_mass_scaled * ec[birds$pasture + 1],
    elevMedian_x_forestPresent_x_pasture = birds$elev_median_scaled * ec[birds$forest_present + 1] * ec[birds$pasture + 1],
    elevMedian_x_forestSpecialist_x_pasture = birds$elev_median_scaled * ec[birds$forest_specialist + 1] * ec[birds$pasture + 1],
    time = birds$time,
    time_x_elev = sweep(birds$time, MARGIN = 1, birds$elev_median_scaled, FUN = `*`)
)

bird_standata7_means_and_sds <- list(time_mean = mean(c(birds$hps1, birds$hps2, birds$hps3, birds$hps4), na.rm = T), 
                                     time_sd =  sd(c(birds$hps1, birds$hps2, birds$hps3, birds$hps4), na.rm = T),
                                     relev_offset = .5, relev_sd = sd(birds$elev_sp_standard),
                                     elev_median_mean =  mean(birds$elev_median), elev_median_sd = sd(birds$elev_median),
                                     elev_breadth_mean = mean(birds$elev_breadth), elev_breadth_sd = sd(birds$elev_breadth),
                                     log_mass_mean = mean(log(birds$BodyMass.Value)), log_mass_sd = sd(log(birds$BodyMass.Value)), 
                                     distance_to_range_offset = 0, distance_to_range_sd = sd(birds$distance_from_range), 
                                     distance_to_range_logit_rescale = 11.3)
bird_stan_data7_package <- list(data = bird_stan_data7,
                                means_and_sds = bird_standata7_means_and_sds)
saveRDS(bird_stan_data7_package, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data7_package.RDS")



###### Stan data 8 #####
ec <- c(-1,1)

grt <- function(predictor){
  return(list(u = unique(predictor), id = match(predictor, unique(predictor))))
}
grt_b <- function(predictor){
  if(sum(predictor != 1 & predictor != -1) != 0){stop('predictor must be effects coded')}
  predictor[predictor == 1] <- 2
  predictor[predictor == -1] <- 1
  return(predictor)
}

bird_stan_data8 <- list(
  # Grainsize for reduce_sum
  grainsize = 1,
  # Dimensions
  # Random effects
  n_spCl = length(unique(birds$sp_cl)),
  n_spSr = length(unique(birds$sp_sr)),
  n_sp = length(unique(birds$species)),
  n_fam = length(unique(birds$Family)),
  n_spObs = length(unique(as.vector(birds$sp_obs_matrix[!is.na(as.vector(birds$sp_obs_matrix))]))),
  n_tot = nrow(birds),
  n_visit_max = max(birds$nv),
  # Unique values for continuous covariates
  n_relev = length(unique(birds$relev)),
  n_relev2 = length(unique(birds$relev2)),
  n_lowland_x_relev = length(unique(ec[birds$lowland+1] * birds$relev)),
  n_lowland_x_relev2 = length(unique(ec[birds$lowland+1] * birds$relev2)),
  n_elevMedian = length(unique(birds$elev_median_scaled)),
  n_elevBreadth = length(unique(birds$elev_breadth_scaled)),
  n_mass = length(unique(birds$log_mass_scaled)),
  n_time = length(unique(birds$time)),
  n_elevMedian_x_forestPresent = length(unique(birds$elev_median_scaled * ec[birds$forest_present + 1])),
  n_elevMedian_x_forestSpecialist = length(unique(birds$elev_median_scaled * ec[birds$forest_specialist + 1])),
  n_elevMedian_x_pasture = length(unique(birds$elev_median_scaled * ec[birds$pasture + 1])),
  n_elevBreadth_x_pasture = length(unique(birds$elev_breadth_scaled * ec[birds$pasture + 1])),
  n_mass_x_pasture = length(unique(birds$log_mass_scaled * ec[birds$pasture + 1])),
  n_elevMedian_x_forestPresent_x_pasture = length(unique(birds$elev_median_scaled * ec[birds$forest_present + 1] * ec[birds$pasture + 1])),
  n_elevMedian_x_forestSpecialist_x_pasture = length(unique(birds$elev_median_scaled * ec[birds$forest_specialist + 1] * ec[birds$pasture + 1])),
  n_time_x_elev = length(unique(sweep(birds$time, MARGIN = 1, birds$elev_median_scaled, FUN = `*`))),
  

  det_data = birds$det_data,
  # Q and nv
  Q = birds$Q,
  nv = birds$nv,
  # Random effect IDs
  id_spCl = as.numeric(as.factor(birds$sp_cl)),
  id_spSr = as.numeric(as.factor(birds$sp_sr)),
  id_sp = as.numeric(as.factor(birds$species)),
  id_fam = as.numeric(as.factor(birds$Family)),
  id_spObs = birds$sp_obs_matrix,

  # Covariates
  relev = birds$relev,
  relev2 = birds$relev2,
  lowland = ec[birds$lowland+1],
  pasture = ec[birds$pasture+1],
  mountain_barrier = grt_b(ec[birds$mountain_limited+1]),
  valley_barrier = grt_b(ec[birds$valley_limited+1]),
  
  elevMedian_u = grt(birds$elev_median_scaled)$u,
  elevMedian_id = grt(birds$elev_median_scaled)$id,
  
  elevBreadth_u = grt(birds$elev_breadth_scaled)$u,
  elevBreadth_id = grt(birds$elev_breadth_scaled)$id,
  
  forestPresent = grt_b(ec[birds$forest_present+1]),
  forestSpecialist = grt_b(ec[birds$forest_specialist+1]),
  tfSpecialist = grt_b(ec[birds$tf_specialist+1]),
  dryForestPresent = grt_b(ec[birds$dry_forest_present+1]),
  floodDrySpecialist = grt_b(ec[birds$flood_dry_specialist+1]),
  aridPresent = grt_b(ec[birds$arid_present+1]),
  migratory = grt_b(ec[as.numeric(!is.na(birds$start1))+1]),
  
  mass_u = grt(birds$log_mass_scaled)$u,
  mass_id = grt(birds$log_mass_scaled)$id,
  
  dietInvert = grt_b(ec[as.numeric(birds$Diet.5Cat == "Invertebrate")+1]),
  dietCarn = grt_b(ec[as.numeric(birds$Diet.5Cat == "VertFishScav")+1]),
  dietFruitNect = grt_b(ec[as.numeric(birds$Diet.5Cat == "FruiNect")+1]),
  dietGran = grt_b(ec[as.numeric(birds$Diet.5Cat == "PlantSeed")+1]),
  distance_to_range = as.vector(birds$distance_from_range_scaled2),
  time = birds$time,
  obsSM = matrix(ec[birds$obsSM+1], ncol=4),
  obsJG = matrix(ec[birds$obsJG+1], ncol=4),
  obsDE = matrix(ec[birds$obsDE+1], ncol=4),
  # pre-computed interactions
  mountainBarrier_x_pasture = grt_b(ec[birds$mountain_limited + 1] * ec[birds$pasture + 1]),
  valleyBarrier_x_pasture = grt_b(ec[birds$valley_limited + 1] * ec[birds$pasture + 1]),
  elevMedian_x_forestPresent = birds$elev_median_scaled * ec[birds$forest_present + 1],
  elevMedian_x_forestSpecialist = birds$elev_median_scaled * ec[birds$forest_specialist + 1],
  forestPresent_x_pasture = grt_b(ec[birds$forest_present + 1] * ec[birds$pasture + 1]),
  forestSpecialist_x_pasture = grt_b(ec[birds$forest_specialist + 1] * ec[birds$pasture + 1]),
  tfSpecialist_x_pasture = grt_b(ec[birds$tf_specialist + 1] * ec[birds$pasture + 1]),
  dryForestPresent_x_pasture = grt_b(ec[birds$dry_forest_present + 1] * ec[birds$pasture + 1]),
  floodDrySpecialist_x_pasture = grt_b(ec[birds$flood_dry_specialist + 1] * ec[birds$pasture + 1]),
  aridPresent_x_pasture = grt_b(ec[birds$arid_present + 1] * ec[birds$pasture + 1]),
  migratory_x_pasture = grt_b(ec[as.numeric(!is.na(birds$start1))+1] * ec[birds$pasture + 1]),
  dietInvert_x_pasture = grt_b(ec[as.numeric(birds$Diet.5Cat == "Invertebrate")+1] * ec[birds$pasture + 1]),
  dietCarn_x_pasture = grt_b(ec[as.numeric(birds$Diet.5Cat == "VertFishScav")+1] * ec[birds$pasture + 1]),
  dietFruitNect_x_pasture = grt_b(ec[as.numeric(birds$Diet.5Cat == "FruiNect")+1] * ec[birds$pasture + 1]),
  dietGran_x_pasture = grt_b(ec[as.numeric(birds$Diet.5Cat == "PlantSeed")+1] * ec[birds$pasture + 1]),
  elevMedian_x_forestPresent_x_pasture = birds$elev_median_scaled * ec[birds$forest_present + 1] * ec[birds$pasture + 1],
  elevMedian_x_forestSpecialist_x_pasture = birds$elev_median_scaled * ec[birds$forest_specialist + 1] * ec[birds$pasture + 1],
  
  elevMedian_x_pasture = birds$elev_median_scaled * ec[birds$pasture + 1],
  elevBreadth_x_pasture = birds$elev_breadth_scaled * ec[birds$pasture + 1],
  mass_x_pasture = birds$log_mass_scaled * ec[birds$pasture + 1],
  time_x_elev = sweep(birds$time, MARGIN = 1, birds$elev_median_scaled, FUN = `*`),
  
  lowland_x_relev = ec[birds$lowland+1] * birds$relev,
  lowland_x_relev2 = ec[birds$lowland+1] * birds$relev2
  
)

bird_standata8_means_and_sds <- list(time_mean = mean(c(birds$hps1, birds$hps2, birds$hps3, birds$hps4), na.rm = T), 
                                     time_sd =  sd(c(birds$hps1, birds$hps2, birds$hps3, birds$hps4), na.rm = T),
                                     relev_offset = .5, relev_sd = sd(birds$elev_sp_standard),
                                     elev_median_mean =  mean(birds$elev_median), elev_median_sd = sd(birds$elev_median),
                                     elev_breadth_mean = mean(birds$elev_breadth), elev_breadth_sd = sd(birds$elev_breadth),
                                     log_mass_mean = mean(log(birds$BodyMass.Value)), log_mass_sd = sd(log(birds$BodyMass.Value)), 
                                     distance_to_range_offset = 0, distance_to_range_sd = sd(birds$distance_from_range), 
                                     distance_to_range_logit_rescale = 11.3)
bird_stan_data8_package <- list(data = bird_stan_data8,
                                means_and_sds = bird_standata8_means_and_sds)
saveRDS(bird_stan_data8_package, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data8_package.RDS")



##### Stan data 9 #####
ec <- c(-1,1)

grt <- function(predictor){
  return(list(u = unique(predictor), id = match(predictor, unique(predictor))))
}
grt_b <- function(predictor){
  if(sum(predictor != 1 & predictor != -1) != 0){stop('predictor must be effects coded')}
  predictor[predictor == 1] <- 2
  predictor[predictor == -1] <- 1
  return(predictor)
}

grt_m <- function(predictor){
  v <- as.vector(predictor)
  u <- unique(v)
  out <- matrix(match(v, u), ncol=ncol(predictor))
  return(out)
}


bird_stan_data9 <- list(
  # Grainsize for reduce_sum
  grainsize = 1,
  # Dimensions
  # Random effects
  n_spCl = length(unique(birds$sp_cl)),
  n_spSr = length(unique(birds$sp_sr)),
  n_sp = length(unique(birds$species)),
  n_fam = length(unique(birds$Family)),
  n_spObs = length(unique(as.vector(birds$sp_obs_matrix[!is.na(as.vector(birds$sp_obs_matrix))]))),
  # Dataset size
  n_tot = nrow(birds),
  n_visit_max = max(birds$nv),
  # Unique values sizes
  n_relev = length(unique(birds$relev)),
  n_relev2 = length(unique(birds$relev2)),
  n_lowland_x_relev = length(unique(ec[birds$lowland+1] * birds$relev)),
  n_lowland_x_relev2 = length(unique(ec[birds$lowland+1] * birds$relev2)),
  n_elevMedian = length(unique(birds$elev_median_scaled)),
  n_elevBreadth = length(unique(birds$elev_breadth_scaled)),
  n_mass = length(unique(birds$log_mass_scaled)),
  n_time = length(unique(as.vector(birds$time))),
  # skip distance to range because it doesn't compress
  n_elevMedian_x_forestPresent = length(unique(birds$elev_median_scaled * ec[birds$forest_present + 1])),
  n_elevMedian_x_forestSpecialist = length(unique(birds$elev_median_scaled * ec[birds$forest_specialist + 1])),
  n_elevMedian_x_pasture = length(unique(birds$elev_median_scaled * ec[birds$pasture + 1])),
  n_elevBreadth_x_pasture = length(unique(birds$elev_breadth_scaled * ec[birds$pasture + 1])),
  n_mass_x_pasture = length(unique(birds$log_mass_scaled * ec[birds$pasture + 1])),
  n_elevMedian_x_forestPresent_x_pasture = length(unique(birds$elev_median_scaled * ec[birds$forest_present + 1] * ec[birds$pasture + 1])),
  n_elevMedian_x_forestSpecialist_x_pasture = length(unique(birds$elev_median_scaled * ec[birds$forest_specialist + 1] * ec[birds$pasture + 1])),
  n_time_x_elev = length(unique(as.vector(sweep(birds$time, MARGIN = 1, birds$elev_median_scaled, FUN = `*`)))),
  # Integer data matrix
  integer_data = data.frame(
    # Detection data
    det_data = birds$det_data,  # 1-4
    # Q and nv
    Q = birds$Q,               # 5
    nv = birds$nv,             # 6
    # Random effect IDs
    id_spCl = as.numeric(as.factor(birds$sp_cl)),     # 7
    id_spSr = as.numeric(as.factor(birds$sp_sr)),     # 8
    id_sp = as.numeric(as.factor(birds$species)),     # 9
    id_fam = as.numeric(as.factor(birds$Family)),     # 10
    id_spObs = birds$sp_obs_matrix,                   # 11-14
    # Covariate indices
    lowland = grt_b(ec[birds$lowland+1]),             # 15
    pasture = grt_b(ec[birds$pasture+1]),                          # 16
    mountain_barrier = grt_b(ec[birds$mountain_limited+1]),        # 17
    valley_barrier = grt_b(ec[birds$valley_limited+1]),            # 18
    forestPresent = grt_b(ec[birds$forest_present+1]),             # 19
    forestSpecialist = grt_b(ec[birds$forest_specialist+1]),       # 20
    tfSpecialist = grt_b(ec[birds$tf_specialist+1]),               # 21
    dryForestPresent = grt_b(ec[birds$dry_forest_present+1]),      # 22
    floodDrySpecialist = grt_b(ec[birds$flood_dry_specialist+1]),  # 23
    aridPresent = grt_b(ec[birds$arid_present+1]),                 # 24
    migratory = grt_b(ec[as.numeric(!is.na(birds$start1))+1]),     # 25
    dietInvert = grt_b(ec[as.numeric(birds$Diet.5Cat == "Invertebrate")+1]),    # 26
    dietCarn = grt_b(ec[as.numeric(birds$Diet.5Cat == "VertFishScav")+1]),      # 27
    dietFruitNect = grt_b(ec[as.numeric(birds$Diet.5Cat == "FruiNect")+1]),     # 28
    dietGran = grt_b(ec[as.numeric(birds$Diet.5Cat == "PlantSeed")+1]),         # 29
    mountainBarrier_x_pasture = grt_b(ec[birds$mountain_limited + 1] * ec[birds$pasture + 1]),        # 30
    valleyBarrier_x_pasture = grt_b(ec[birds$valley_limited + 1] * ec[birds$pasture + 1]),            # 31
    forestPresent_x_pasture = grt_b(ec[birds$forest_present + 1] * ec[birds$pasture + 1]),            # 32
    forestSpecialist_x_pasture = grt_b(ec[birds$forest_specialist + 1] * ec[birds$pasture + 1]),      # 33
    tfSpecialist_x_pasture = grt_b(ec[birds$tf_specialist + 1] * ec[birds$pasture + 1]),              # 34
    dryForestPresent_x_pasture = grt_b(ec[birds$dry_forest_present + 1] * ec[birds$pasture + 1]),     # 35
    floodDrySpecialist_x_pasture = grt_b(ec[birds$flood_dry_specialist + 1] * ec[birds$pasture + 1]), # 36
    aridPresent_x_pasture = grt_b(ec[birds$arid_present + 1] * ec[birds$pasture + 1]),                # 37
    migratory_x_pasture = grt_b(ec[as.numeric(!is.na(birds$start1))+1] * ec[birds$pasture + 1]),      # 38
    dietInvert_x_pasture = grt_b(ec[as.numeric(birds$Diet.5Cat == "Invertebrate")+1] * ec[birds$pasture + 1]),   # 39
    dietCarn_x_pasture = grt_b(ec[as.numeric(birds$Diet.5Cat == "VertFishScav")+1] * ec[birds$pasture + 1]),     # 40
    dietFruitNect_x_pasture = grt_b(ec[as.numeric(birds$Diet.5Cat == "FruiNect")+1] * ec[birds$pasture + 1]),    # 41
    dietGran_x_pasture = grt_b(ec[as.numeric(birds$Diet.5Cat == "PlantSeed")+1] * ec[birds$pasture + 1]),        # 42
    obsSM = matrix(grt_b(ec[birds$obsSM+1]), ncol=4),   # 43-46
    obsJG = matrix(grt_b(ec[birds$obsJG+1]), ncol=4),   # 47-50
    obsDE = matrix(grt_b(ec[birds$obsDE+1]), ncol=4)    # 51-54
  ),
  # Continuous distance-to-range (does not compress well)
  distance_to_range = as.vector(birds$distance_from_range_scaled2),
  # continuous covariates
  relev = birds$relev,
  relev2 = birds$relev2,
  lowland_x_relev = ec[birds$lowland+1] * birds$relev,
  lowland_x_relev2 = ec[birds$lowland+1] * birds$relev2,
  elevMedian = birds$elev_median_scaled,
  elevBreadth = birds$elev_breadth_scaled,
  mass = birds$log_mass_scaled,
  elevMedian_x_forestPresent = birds$elev_median_scaled * ec[birds$forest_present + 1],
  elevMedian_x_forestSpecialist = birds$elev_median_scaled * ec[birds$forest_specialist + 1],
  elevMedian_x_pasture = birds$elev_median_scaled * ec[birds$pasture + 1],
  elevBreadth_x_pasture = birds$elev_breadth_scaled * ec[birds$pasture + 1],
  mass_x_pasture = birds$log_mass_scaled * ec[birds$pasture + 1],
  elevMedian_x_forestPresent_x_pasture = birds$elev_median_scaled * ec[birds$forest_present + 1] * ec[birds$pasture + 1],
  elevMedian_x_forestSpecialist_x_pasture = birds$elev_median_scaled * ec[birds$forest_specialist + 1] * ec[birds$pasture + 1],
  time = birds$time,
  time_x_elev = sweep(birds$time, MARGIN = 1, birds$elev_median_scaled, FUN = `*`)
)

bird_standata9_means_and_sds <- list(time_mean = mean(c(birds$hps1, birds$hps2, birds$hps3, birds$hps4), na.rm = T), 
                                     time_sd =  sd(c(birds$hps1, birds$hps2, birds$hps3, birds$hps4), na.rm = T),
                                     relev_offset = .5, relev_sd = sd(birds$elev_sp_standard),
                                     elev_median_mean =  mean(birds$elev_median), elev_median_sd = sd(birds$elev_median),
                                     elev_breadth_mean = mean(birds$elev_breadth), elev_breadth_sd = sd(birds$elev_breadth),
                                     log_mass_mean = mean(log(birds$BodyMass.Value)), log_mass_sd = sd(log(birds$BodyMass.Value)), 
                                     distance_to_range_offset = 0, distance_to_range_sd = sd(birds$distance_from_range), 
                                     distance_to_range_logit_rescale = 11.3)
bird_stan_data9_package <- list(data = bird_stan_data9,
                                means_and_sds = bird_standata9_means_and_sds)
saveRDS(bird_stan_data9_package, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data9_package.RDS")

