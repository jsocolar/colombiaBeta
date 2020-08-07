`%ni%` <- Negate(`%in%`)

bird_surveys <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_surveys_current.RDS')
point_distances <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/point_distances/point_distances_biogeographic_clip.RDS")


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

Q <- apply(bird_surveys$detection_array, MARGIN = c(1,3), FUN = function(x){return(sum(x, na.rm = T) > 0)})
