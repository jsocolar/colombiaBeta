# This script takes previously ingested and standardized data on birds, traits, ranges, points, and visits and produces a unified data object for analysis

##### Script dependencies: combined_bird_maps.R, bird_import_and_cleaning.R, elevations_prep_and_exploration.R, eventually the points formatting script #####

`%ni%` <- Negate(`%in%`)

bird_surveys <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_surveys_current.RDS')
point_distances <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/point_distances/point_distances_biogeographic_clip_ayerbe.RDS")

include_species <- vector()
for(i in 1:nrow(point_distances)){
  include_species[i] <- sum(point_distances[i,] > 0) != ncol(point_distances)
}

point_distances_include <- point_distances[include_species, ]

species_list <- gsub(" ", "_", row.names(point_distances_include))

which(bird_surveys$species_names %ni% species_list)

det_array <- bird_surveys$detection_array[,1:4,]
det_array_padded <- abind::abind(det_array, array(data = 0, dim = c(848, 4, length(species_list) - dim(det_array)[3])), along = 3)

species_names <- c(bird_surveys$species_names, species_list[species_list %ni% bird_surveys$species_names])

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


flattened_data$nv <- 4
flattened_data$nv[is.na(flattened_data$v4)] <- 3
flattened_data$nv[is.na(flattened_data$v3)] <- 2

flattened_data$Q <- as.numeric(rowSums(flattened_data[,3:6], na.rm = T) > 0)


all_pts <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Points/all_pts.RDS")


fd <- merge(flattened_data, all_pts, by.x = "point", by.y = "point", all.x = T)




traits <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/traits/traits.RDS")

fd <- merge(fd, traits, by.x = "species", by.y = "latin_underscore", all.x = T)
fd$elev_sp_standard <- (fd$elev_ALOS - fd$lower)/(fd$upper - fd$lower)

all(flattened_data$species %in% traits$latin_underscore)

####
mig_dates <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/traits/mig_dates.RDS")
mig_dates$species <- gsub(" ", "_", mig_dates$latin)

fd <- merge(fd, mig_dates, by.x = "species", by.y = "species", all.x = T)

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
sum(fd$Q == 1 & fd$in_date_range == 0)

flattened_data_full <- fd[, names(fd) %ni% c("birds", "beetles", "habitat", "other", "latin.x", "latin.y")]
saveRDS(flattened_data_full, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/flattened_data_full.RDS")
  
  
  
final_data <- fd[fd$in_date_range == 1 & fd$elev_sp_standard > -1 & fd$elev_sp_standard < 2, ]

nrow(final_data)
sum(final_data$Q)
nrow(final_data)


fdq <- fd[fd$Q == 1,]

View(fdq[fdq$elev_sp_standard > 1.9 | fdq$elev_sp_standard < -.9,])

length(unique(fd$species[fd$elev_sp_standard > -1 & fd$elev_sp_standard < 2]))
sum(fd$elev_sp_standard > -1 & fd$elev_sp_standard < 2, na.rm = T)

sum(is.na(fd$elev_sp_standard > -1 & fd$elev_sp_standard < 2))

a <- seq(-1,2,.2)
nq <- vector()
nall <- vector()
for(i in 2:length(a)){
  nq[i-1] <- sum(fdq$elev_sp_standard > a[i-1] & fdq$elev_sp_standard <= a[i])
  nall[i-1] <- sum(fd$elev_sp_standard > a[i-1] & fd$elev_sp_standard <= a[i])
}

plot(nq/nall ~ seq(-.9, 1.9, .2))
plot(nall ~ seq(-.9, 1.9, .2))
min(nall)


cutoff <- 2
sum(fd$elev_sp_standard > cutoff, na.rm = T) + sum(fd$elev_sp_standard < (1 - cutoff), na.rm = T)


