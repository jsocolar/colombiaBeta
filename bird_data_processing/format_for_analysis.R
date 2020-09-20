# This script takes previously ingested and standardized data on birds, traits, ranges, points, and visits and produces a unified data object for analysis

##### Script dependencies: combined_bird_maps.R, bird_import_and_cleaning.R, elevations_prep_and_exploration.R, points_formatting.R, migratory_dates.R, species_covariate_formatting.R #####

`%ni%` <- Negate(`%in%`)

# Get formatted bird surveys object
bird_surveys <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_surveys_current.RDS')

# Get a matrix (actually a dataframe) of the distance from each species range to each sampling point.
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
saveRDS(flattened_data, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/flattened_data.RDS")
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
saveRDS(flattened_data_full, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/flattened_data_full.RDS")

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

plot(nq/nall ~ seq(-.9, 1.9, .2))
plot(nall ~ seq(-.9, 1.9, .2))
min(nall)


# Subset to include only species-points that are in the date range and in an elevational range of (-1, 2)
bird_data_trimmed <- flattened_data_full[flattened_data_full$in_date_range == 1 & flattened_data_full$elev_sp_standard > -1 & flattened_data_full$elev_sp_standard < 2, ]
saveRDS(bird_data_trimmed, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_data_trimmed.RDS")

# Examine statistics of final dataset
nrow(bird_data_trimmed)
sum(bird_data_trimmed$Q)
mean(bird_data_trimmed$Q)



