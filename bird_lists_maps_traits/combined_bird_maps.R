# This script produces analysis-ready updated Ayerbe/BirdLife range maps and regionally buffered maps.

library(sf)
library(ggplot2)
`%ni%` <- Negate(`%in%`)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

##### Get colombia mainland shapefile #####
colombia <- st_read('/Users/JacobSocolar/Dropbox/Work/Colombia/Data/GIS/colombia_maps/gadm36_COL_shp/gadm36_COL_0.shp')
# file consists of many disjoint polygons, representing the mainland and numerous islands. Figure out which is the mainland and extract it.
npoly <- length(colombia$geometry[[1]])
size <- rep(0,npoly)
for(i in 1:npoly){
  size[i] <- dim(colombia$geometry[[1]][[i]][[1]])[1]
}
mainland <- colombia$geometry[[1]][[which(size == max(size))]] %>%
  st_polygon() %>% st_sfc() %>% st_sf()
st_crs(mainland) <- st_crs(colombia)
# Transform to AEA conic centered on Colombia
mainland <- st_transform(mainland, AEAstring)
buffered_mainland <- st_buffer(mainland, 100000)
buffered_mainland2 <- st_buffer(mainland, 200000)
m2 <- st_simplify(mainland, dTolerance = 10000)


##### Combined BirdLife + Ayerbe maps #####
initial_species_list <- read.csv("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/species_list_creation/initial_species_list.csv")
load('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/birdlife_maps/recast_range_maps.Rdata')  # from Species_lists.R
ayerbe_maps <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/ayerbe_maps/ayerbe_maps.RDS') # from ayerbe_maps.R

# extract only the relevant species
rrm <- recast_range_maps[recast_range_maps$SCINAME %in% initial_species_list$HBW,]

# under WGS84, crop to a 200-km buffer of mainland Colombia (I think that switching to AEA before cropping perhaps introduces some weird artifacts, perhaps from distant parts of the globe)
birdlife_valid <- st_make_valid(rrm)
save(birdlife_valid, file = '/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/combined_maps/birdlife_valid.Rdata')

birdlife_cropped <- st_intersection(birdlife_valid, st_transform(buffered_mainland2, st_crs(birdlife_valid)))

birdlife_transformed <- st_transform(birdlife_cropped, AEAstring)

birdlife_colombia <- st_intersection(birdlife_transformed, buffered_mainland)

saveRDS(birdlife_colombia, file = '/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/combined_maps/birdlife_colombia.RDS')

birdlife_colombia <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/combined_maps/birdlife_colombia.RDS')

blc <- list()
for(i in 1:nrow(initial_species_list)){
  print(i)
  blc[[i]] <- st_union(birdlife_colombia[birdlife_colombia$SCINAME == initial_species_list$HBW[i], ])
}

combined_maps <- list()
for(i in 1:nrow(initial_species_list)){
  print(i)
  combined_maps[[i]] <- st_union(blc[[i]], ayerbe_maps[ayerbe_maps$Species == initial_species_list$HBW[i], ])
}
names(combined_maps) <- initial_species_list$HBW

##### Deal with Grallaria rufula #####
rufula_lato <- combined_maps[["Grallaria rufula"]]
rufula_ig <- st_transform(st_read("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/rufula_updates/Iguaque.kml"), AEAstring)
rufula_west <- st_transform(st_read("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/rufula_updates/westAndes.kml"), AEAstring)
rufula_central <- st_transform(st_read("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/rufula_updates/centralSouthernAndes.kml"), AEAstring)
rufula_east <- st_transform(st_read("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/rufula_updates/eastAndes.kml"), AEAstring)
rufula_snsm <- st_transform(st_read("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/rufula_updates/snsm.kml"), AEAstring)

G_alvarezi <- st_difference(rufula_lato, st_union(rufula_central, st_union(rufula_east, rufula_snsm)))
G_spatiator <- st_difference(rufula_lato, st_union(rufula_west, st_union(rufula_central, rufula_east)))
G_rufula <- st_difference(rufula_lato, st_union(rufula_ig, st_union(rufula_west, st_union(rufula_central, rufula_snsm))))
G_saturata <- st_union(rufula_ig, st_difference(rufula_lato, st_union(rufula_west, st_union(rufula_east, rufula_snsm))))

combined_maps[["Grallaria rufula"]] <- G_rufula
combined_maps[["Grallaria spatiator"]] <- G_spatiator
combined_maps[["Grallaria alvarezi"]] <- G_alvarezi
combined_maps[["Grallaria saturata"]] <- G_saturata

saveRDS(combined_maps, file = '/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/combined_maps/combined_maps.RDS')
combined_maps <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/combined_maps/combined_maps.RDS')

##### Check distances from occurrences in dataset to combined maps #####
bird_surveys <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_surveys_current.RDS')
colombia_points <- st_as_sf(data.frame(readxl::read_excel('/Users/jacobsocolar/Downloads/RAWdata_samplingpointsCO_MASTER.xlsx')), 
                            coords = c('long', 'lat'), crs = 4326)
for(i in 10:12){
  colombia_points$point_id <- gsub(paste0('CHA', i, 'd'), paste0('CHA', i, 'D'), colombia_points$point_id)
}
bird_points <- st_transform(colombia_points[colombia_points$point_id %in% bird_surveys$point_names, ], AEAstring)

distances <- list()
dsp <- vector()
sp_pts <- list()
for(i in 1:length(bird_surveys$species_names)){
    species <- bird_surveys$species_names[i]
    species_range <- combined_maps[[gsub('_', ' ', species)]]
    surveys <- bird_surveys$detection_array[,,i]
    point_names <- bird_surveys$point_names[rowSums(surveys, na.rm = T) > 0]
    points <- bird_points[bird_points$point_id %in% point_names, ]
    distance_matrix <- st_distance(points, species_range)
    distances[[i]] <- as.vector(distance_matrix)
    dsp[i] <- species
    sp_pts[[i]] <- points
}

mdist <- vector()
for(i in 1:length(dsp)){
  mdist[i] <- max(distances[[i]])
}

d <- 100000
bad_species <- dsp[mdist > d]

bad_species



for(i in 1:length(bad_species)){
  print(ggplot(m2) + geom_sf() +
    geom_sf(data = birdlife_colombia[birdlife_colombia$SCINAME == gsub("_", " ", bad_species[i]),], inherit.aes = F, fill = 'cornflowerblue') +
    geom_sf(data = ayerbe_maps[ayerbe_maps$Species == gsub("_", " ", bad_species[i]),], inherit.aes = F, fill = 'red', alpha = .2) +
    geom_sf(data = sp_pts[[which(dsp == bad_species[i])]], inherit.aes = F) +
    ggtitle(bad_species[i]))
}

source('/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/GIS_processing/hydrosheds_extraction.R')

clipping_polygons <- list(amazon_orinoco = st_transform(amazon_orinoco, AEAstring), magdalena_east = st_transform(magdalena_east, AEAstring),
                          magdalena_west = st_transform(magdalena_west, AEAstring), cauca_east = st_transform(cauca_east, AEAstring),
                          cauca_west = st_transform(cauca_west, AEAstring), pasto = st_transform(pasto, AEAstring), 
                          pacific = st_transform(pacific, AEAstring), snsm = st_transform(snsm, AEAstring), 
                          guajira_perija = st_transform(guajira_perija, AEAstring), catatumbo = st_transform(catatumbo, AEAstring))

buffered_ranges <- list()
for(i in 1:nrow(initial_species_list)){
  print(i)
  species <- gsub("_", " ", initial_species_list$HBW[i])
  range <- combined_maps[[species]]
  buffer_list <- list()
  counter <- 0
  for(j in 1:length(clipping_polygons)){
    cp <- clipping_polygons[[j]]
    clipped_range <- st_intersection(cp, range)
    if(nrow(clipped_range) > 0){
      if(nrow(clipped_range) != 1){clipped_range <- st_union(clipped_range)}   # this is an issue for Ayerbe splits/HBW lumps, e.g. Frederickena
      counter <- counter + 1
      cr <- clipped_range[, 1]
      buffered_cr <- st_buffer(cr, 100000)
      buffer_list[[counter]] <- st_intersection(buffered_cr, cp)
    }
  }
  buffered_output <- st_union(buffer_list[[1]])
  if(length(buffer_list) > 1){
    for(k in 2:length(buffer_list)){
      buffered_output <- st_union(buffered_output, buffer_list[[k]])
    }
  }
  buffered_ranges[[i]] <- buffered_output
  names(buffered_ranges)[i] <- species
}

saveRDS(buffered_ranges, file = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/combined_maps/buffered_ranges.RDS")
buffered_ranges <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/combined_maps/buffered_ranges.RDS")

distances <- list()
dsp <- vector()
sp_pts <- list()
for(i in 1:length(bird_surveys$species_names)){
  species <- bird_surveys$species_names[i]
  species_range <- buffered_ranges[[gsub('_', ' ', species)]]
  surveys <- bird_surveys$detection_array[,,i]
  point_names <- bird_surveys$point_names[rowSums(surveys, na.rm = T) > 0]
  points <- bird_points[bird_points$point_id %in% point_names, ]
  distance_matrix <- st_distance(points, species_range)
  distances[[i]] <- as.vector(distance_matrix)
  dsp[i] <- species
  sp_pts[[i]] <- points
}

mdist <- vector()
for(i in 1:length(dsp)){
  mdist[i] <- max(distances[[i]])
}

d <- 0
bad_species <- dsp[mdist > d]

bad_species



for(i in 1:length(bad_species)){
  print(ggplot(m2) + geom_sf() +
          geom_sf(data = buffered_ranges[[gsub("_", " ", bad_species[i])]], fill = 'green') +
          geom_sf(data = birdlife_colombia[birdlife_colombia$SCINAME == gsub("_", " ", bad_species[i]),], inherit.aes = F, fill = 'cornflowerblue', alpha = .2) +
          geom_sf(data = ayerbe_maps[ayerbe_maps$Species == gsub("_", " ", bad_species[i]),], inherit.aes = F, fill = 'red', alpha = .2) +
          geom_sf(data = sp_pts[[which(dsp == bad_species[i])]], inherit.aes = F) +
          ggtitle(bad_species[i]))
}


# Needs update; new extension:
# Crypturellus brevirostris, Hemitriccus striaticollis, Leistes militaris, Thamnophilus unicolor

# known map error:
# Thlypopsis superciliaris, Momotus subrufescens, Endomodestes coracinus, Myrmotherula ambigua


# Suspected dataset errors
# Elaenia pallatangae, Myiozetetes similis, Phyllomyias uropygialis


##### Combined map errors #####
combined_maps <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/combined_maps/combined_maps.RDS')

update_files <- list.files('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/combined_map_updates', recursive = T, full.names = T)
update_species <- gsub("_", " ", list.files('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/combined_map_updates'))

for(i in 1:length(update_files)){
  update <- st_transform(st_zm(st_read(update_files[i])), AEAstring)
  combined_maps[[update_species[i]]] <- st_union(combined_maps[[update_species[i]]], update)
}

combined_maps_updated <- combined_maps
saveRDS(combined_maps_updated, file = '/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/combined_maps/combined_maps_updated.RDS')
combined_maps_updated <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/combined_maps/combined_maps_updated.RDS')

# Redo the buffering for the updated ranges
buffered_ranges_updated <- buffered_ranges
for(i in which(names(buffered_ranges) %in% update_species)){
  print(i)
  us <- names(buffered_ranges)[i]
  print(us)
  range <- combined_maps[[us]]
  buffer_list <- list()
  counter <- 0
  for(j in 1:length(clipping_polygons)){
    cp <- clipping_polygons[[j]]
    clipped_range <- st_intersection(cp, range)
    if(nrow(clipped_range) > 0){
      if(nrow(clipped_range) != 1){clipped_range <- st_union(clipped_range)}   # this is an issue for Ayerbe splits/HBW lumps, e.g. Frederickena
      counter <- counter + 1
      cr <- clipped_range[, 1]
      buffered_cr <- st_buffer(cr, 100000)
      buffer_list[[counter]] <- st_intersection(buffered_cr, cp)
    }
  }
  buffered_output <- st_union(buffer_list[[1]])
  if(length(buffer_list) > 1){
    for(k in 2:length(buffer_list)){
      buffered_output <- st_union(buffered_output, buffer_list[[k]])
    }
  }
  buffered_ranges_updated[[i]] <- buffered_output
}

saveRDS(buffered_ranges_updated, file = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/combined_maps/buffered_ranges_updated.RDS")
buffered_ranges_updated <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/combined_maps/buffered_ranges_updated.RDS")


distances <- list()
dsp <- vector()
sp_pts <- list()
for(i in 1:length(bird_surveys$species_names)){
  species <- bird_surveys$species_names[i]
  species_range <- buffered_ranges_updated[[gsub('_', ' ', species)]]
  surveys <- bird_surveys$detection_array[,,i]
  point_names <- bird_surveys$point_names[rowSums(surveys, na.rm = T) > 0]
  points <- bird_points[bird_points$point_id %in% point_names, ]
  distance_matrix <- st_distance(points, species_range)
  distances[[i]] <- as.vector(distance_matrix)
  dsp[i] <- species
  sp_pts[[i]] <- points
}

mdist <- vector()
for(i in 1:length(dsp)){
  mdist[i] <- max(distances[[i]])
}

d <- 0
bad_species <- dsp[mdist > d]

bad_species


for(i in 1:length(bad_species)){
  print(ggplot(m2) + geom_sf() +
          geom_sf(data = buffered_ranges_updated[[gsub("_", " ", bad_species[i])]], fill = 'green') +
          geom_sf(data = buffered_ranges[[gsub("_", " ", bad_species[i])]], fill = 'limegreen') +
          geom_sf(data = birdlife_colombia[birdlife_colombia$SCINAME == gsub("_", " ", bad_species[i]),], inherit.aes = F, fill = 'cornflowerblue', alpha = .2) +
          geom_sf(data = ayerbe_maps[ayerbe_maps$Species == gsub("_", " ", bad_species[i]),], inherit.aes = F, fill = 'red', alpha = .2) +
          geom_sf(data = sp_pts[[which(dsp == bad_species[i])]], inherit.aes = F) +
          ggtitle(bad_species[i]))
}


point_distances <- as.data.frame(matrix(data = NA, nrow = nrow(initial_species_list), ncol = nrow(bird_points)))
names(point_distances) <- bird_points$point_id
row.names(point_distances) <- initial_species_list$HBW

st_distance(bird_points, buffered_ranges_updated[[1]])

for(i in 1:length(buffered_ranges_updated)){
  print(i)
  point_distances[i, ] <- st_distance(bird_points, buffered_ranges_updated[[i]])
}

saveRDS(point_distances, file = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/point_distances/point_distances_biogeographic_clip.RDS")

include <- exclude <- vector()
for(i in 1:nrow(point_distances)){include[i] <- sum(point_distances[i,] > 0) != ncol(point_distances)}
for(i in 1:nrow(point_distances)){exclude[i] <- sum(point_distances[i,] > 0) == ncol(point_distances)}

initial_species_list$HBW[exclude]

sum(point_distances == 0)

