# This script produces analysis-ready updated Ayerbe range maps and buffered 
# Ayerbe range maps

###### Script dependencies: 
#   species_lists.R
#   ayerbe_maps.R
#   bird_import_and_cleaning.R
#   hydrosheds_extraction.R

# housekeeping ----
library(sf); library(ggplot2)
`%ni%` <- Negate(`%in%`)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 
    +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

##### Get colombia mainland shapefile #####
colombia <- st_read('inputs/gadm36_COL_shp/gadm36_COL_0.shp')
# file consists of many disjoint polygons, representing the mainland and 
# numerous islands. Figure out which is the mainland and extract it.
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


##### Buffered Ayerbe maps #####
initial_species_list <- read.csv("outputs/initial_species_list.csv")
ayerbe_maps <- readRDS('outputs/ayerbe_maps.RDS') # from ayerbe_maps.R

ayerbe_maps$Species[ayerbe_maps$Species %ni% initial_species_list$HBW]
# Pseudastur occidentalis; occurrence in Colombia is sketchy, and in any case 
# this is not relevant to our points
# Trogon caligatus; updated lump
# Corapipo altera; updated lump
# Chlorostilbon melanorhynchus; updated lump
# Basileuterus tacarcunae; updated lump
# Turdus assimilis; updated lump
initial_species_list$HBW[initial_species_list$HBW %ni% ayerbe_maps$Species]

ayerbe_list_fname <- 'outputs/ayerbe_list.RDS'
if(!file.exists(ayerbe_list_fname)) {
    ayerbe_list <- list()
    for(i in 1:nrow(initial_species_list)){
        print(i)
        ayerbe_list[[i]] <- st_make_valid(ayerbe_maps[ayerbe_maps$Species == initial_species_list$HBW[i], ])
    }
    names(ayerbe_list) <- initial_species_list$HBW
    ayerbe_list[["Trogon violaceus"]] <- st_union(ayerbe_list[["Trogon violaceus"]], 
                                                  ayerbe_maps[ayerbe_maps$Species == "Trogon caligatus", ])[,1]
    ayerbe_list[["Corapipo leucorrhoa"]] <- st_union(ayerbe_list[["Corapipo leucorrhoa"]], 
                                                     ayerbe_maps[ayerbe_maps$Species == "Corapipo altera", ])[,1]
    ayerbe_list[["Chlorostilbon mellisugus"]] <- st_union(ayerbe_list[["Chlorostilbon mellisugus"]], 
                                                          ayerbe_maps[ayerbe_maps$Species == "Chlorostilbon melanorhynchus", ])[,1]
    ayerbe_list[["Basileuterus tristriatus"]] <- st_union(ayerbe_list[["Basileuterus tristriatus"]], 
                                                          ayerbe_maps[ayerbe_maps$Species == "Basileuterus tacarcunae", ])[,1]
    ayerbe_list[["Turdus albicollis"]] <- st_union(ayerbe_list[["Turdus albicollis"]], 
                                                   ayerbe_maps[ayerbe_maps$Species == "Turdus assimilis", ])[,1]
    
    
    ayerbe_list[["Vireo olivaceus"]] <- mainland
    T_ochr <- st_transform(st_read("inputs/ayerbe_updates/ayerbe_missing2/Troglodytes_ochraceus/Troglodytes_ochraceus.kml"), 
                           AEAstring)
    ayerbe_list[["Troglodytes ochraceus"]] <- T_ochr
    
    ##### Deal with Grallaria rufula #####
    
    rufula_lato <- ayerbe_list[["Grallaria rufula"]]
    rufula_ig <- st_transform(st_read("inputs/rufula_updates/Iguaque.kml"), AEAstring)
    rufula_west <- st_transform(st_read("inputs/rufula_updates/westAndes.kml"), AEAstring)
    rufula_central <- st_transform(st_read("inputs/rufula_updates/centralSouthernAndes.kml"), AEAstring)
    rufula_east <- st_transform(st_read("inputs/rufula_updates/eastAndes.kml"), AEAstring)
    rufula_snsm <- st_transform(st_read("inputs/rufula_updates/snsm.kml"), AEAstring)
    
    G_alvarezi <- st_difference(rufula_lato, st_union(rufula_central, st_union(rufula_east, rufula_snsm)))
    G_spatiator <- st_difference(rufula_lato, st_union(rufula_west, st_union(rufula_central, rufula_east)))
    G_rufula <- st_difference(rufula_lato, st_union(rufula_ig, st_union(rufula_west, st_union(rufula_central, rufula_snsm))))
    G_saturata <- st_union(rufula_ig, st_difference(rufula_lato, st_union(rufula_west, st_union(rufula_east, rufula_snsm))))
    
    ayerbe_list[["Grallaria rufula"]] <- G_rufula
    ayerbe_list[["Grallaria spatiator"]] <- G_spatiator
    ayerbe_list[["Grallaria alvarezi"]] <- G_alvarezi
    ayerbe_list[["Grallaria saturata"]] <- G_saturata
    
    saveRDS(ayerbe_list, file = ayerbe_list_fname)
} else {
    ayerbe_list <- readRDS(ayerbe_list_fname)   
}

##### Buffered ayerbe maps #####
ayerbe_fname <- "outputs/ayerbe_buffered_ranges.RDS"
if(!file.exists(ayerbe_fname)) {
    source('code/GIS_processing/hydrosheds_extraction.R')
    
    clipping_polygons <- list(amazon_orinoco = st_transform(amazon_orinoco, AEAstring), 
                              magdalena_east = st_transform(magdalena_east, AEAstring),
                              magdalena_west = st_transform(magdalena_west, AEAstring), 
                              cauca_east = st_transform(cauca_east, AEAstring),
                              cauca_west = st_transform(cauca_west, AEAstring), 
                              pasto = st_transform(pasto, AEAstring),
                              pacific = st_transform(pacific, AEAstring), 
                              snsm = st_transform(snsm, AEAstring),
                              guajira_perija = st_transform(guajira_perija, AEAstring), 
                              catatumbo = st_transform(catatumbo, AEAstring),
                              tacarcuna = st_transform(tacarcuna, AEAstring))
    
    
    ayerbe_buffered_ranges <- list()
    for(i in 1:nrow(initial_species_list)){
        print(i)
        species <- initial_species_list$HBW[i]
        sp_range <- ayerbe_list[[species]]
        if(nrow(sp_range) == 0){stop()}
        buffer_list <- list()
        counter <- 0
        for(j in 1:length(clipping_polygons)){
            cp <- clipping_polygons[[j]]
            clipped_range <- st_intersection(cp, sp_range)
            if(nrow(clipped_range) > 0){
                # this is an issue for Ayerbe splits/HBW lumps, e.g. Frederickena
                if(nrow(clipped_range) != 1){clipped_range <- st_union(clipped_range)}  
                counter <- counter + 1
                if(i != 1045){cr <- clipped_range[, 1]}else{cr <- clipped_range}
                buffered_cr <- st_buffer(cr, 160000)
                buffer_list[[counter]] <- st_intersection(buffered_cr, cp)
            }
        }
        buffered_output <- st_union(buffer_list[[1]])
        if(length(buffer_list) > 1){
            for(k in 2:length(buffer_list)){
                buffered_output <- st_union(buffered_output, buffer_list[[k]])
            }
        }
        ayerbe_buffered_ranges[[i]] <- buffered_output
        names(ayerbe_buffered_ranges)[i] <- species
    }
    
    saveRDS(ayerbe_buffered_ranges, file = ayerbe_fname)
} else {
    ayerbe_buffered_ranges <- readRDS(ayerbe_fname)
}

bird_surveys <- readRDS('outputs/bird_surveys_current.RDS')
colombia_points <- st_as_sf(readRDS("outputs/all_pts.RDS"), 
                            coords = c('lon', 'lat'), crs = 4326)
bird_points <- st_transform(colombia_points[colombia_points$point %in% bird_surveys$point_names, ], AEAstring)

ayerbe_distances <- list()
dsp <- vector()
sp_pts <- list()
for(i in 1:length(bird_surveys$species_names)){
  species <- bird_surveys$species_names[i]
  species_range <- ayerbe_buffered_ranges[[gsub('_', ' ', species)]]
  surveys <- bird_surveys$detection_array[,,i]
  point_names <- bird_surveys$point_names[rowSums(surveys, na.rm = T) > 0]
  points <- bird_points[bird_points$point %in% point_names, ]
  if(!is.na(species_range)){
    distance_matrix <- st_distance(points, species_range)
    ayerbe_distances[[i]] <- as.vector(distance_matrix)
    dsp[i] <- species
    sp_pts[[i]] <- points
  }else{
    dsp[i] <- species
    ayerbe_distances[[i]] <- NA
    sp_pts[[i]] <- points
  }
}

mdist <- vector()
for(i in 1:length(dsp)){
  mdist[i] <- max(ayerbe_distances[[i]])
}

d <- 0
bad_species <- dsp[mdist > d]
bad_species <- bad_species[!is.na(bad_species)]
bad_species


for(i in 1:length(bad_species)){
  print(ggplot(m2) + geom_sf() +
          geom_sf(data = ayerbe_buffered_ranges[[gsub("_", " ", bad_species[i])]], fill = 'green') +
          geom_sf(data = sp_pts[[which(dsp == bad_species[i])]], inherit.aes = F) +
          ggtitle(bad_species[i]))
}

for(i in 1:length(bad_species)){
  print(ggplot(m2) + geom_sf() +
          geom_sf(data = ayerbe_list[[gsub("_", " ", bad_species[i])]], fill = 'blue') +
          geom_sf(data = sp_pts[[which(dsp == bad_species[i])]], inherit.aes = F) +
          ggtitle(bad_species[i]))
}


##### Combined map errors #####
ayerbe_upd_list_fname <- 'outputs/ayerbe_list_updated.RDS'
if(!exists(ayerbe_upd_list_fname)) {
    update_files <- list.files('inputs/ayerbe_map_updates', recursive = T, full.names = T)
    update_species <- gsub("_", " ", list.files('inputs/ayerbe_map_updates'))
    
    for(i in 1:length(update_files)){
        update <- st_transform(st_zm(st_read(update_files[i])), AEAstring)
        ayerbe_list[[update_species[i]]] <- st_union(ayerbe_list[[update_species[i]]], update)
    }
    
    ayerbe_list_updated <- ayerbe_list
    saveRDS(ayerbe_list_updated, file = ayerbe_upd_list_fname)
} else {
    ayerbe_list_updated <- readRDS(ayerbe_upd_list_fname)  
}

# Redo the buffering for the updated ranges
update_buff_ranges_fname <- "outputs/ayerbe_buffered_ranges_updated.RDS"
if(!exists(update_buff_ranges_fname)) {
    ayerbe_buffered_ranges_updated <- ayerbe_buffered_ranges
    for(i in which(names(ayerbe_buffered_ranges) %in% update_species)){
        print(i)
        us <- names(ayerbe_buffered_ranges)[i]
        print(us)
        range <- ayerbe_list_updated[[us]]
        buffer_list <- list()
        counter <- 0
        for(j in 1:length(clipping_polygons)){
            cp <- clipping_polygons[[j]]
            clipped_range <- st_intersection(cp, range)
            if(nrow(clipped_range) > 0){
                # this is an issue for Ayerbe splits/HBW lumps, e.g. Frederickena
                if(nrow(clipped_range) != 1){clipped_range <- st_union(clipped_range)}
                counter <- counter + 1
                cr <- clipped_range[, 1]
                buffered_cr <- st_buffer(cr, 160000)
                buffer_list[[counter]] <- st_intersection(buffered_cr, cp)
            }
        }
        buffered_output <- st_union(buffer_list[[1]])
        if(length(buffer_list) > 1){
            for(k in 2:length(buffer_list)){
                buffered_output <- st_union(buffered_output, buffer_list[[k]])
            }
        }
        ayerbe_buffered_ranges_updated[[i]] <- buffered_output
    }
    
saveRDS(ayerbe_buffered_ranges_updated, file = update_buff_ranges_fname)
} else {
    ayerbe_buffered_ranges_updated <- readRDS(update_buff_ranges_fname)
}

distances <- list()
dsp <- vector()
sp_pts <- list()
for(i in 1:length(bird_surveys$species_names)){
  species <- bird_surveys$species_names[i]
  species_range <- ayerbe_buffered_ranges_updated[[gsub('_', ' ', species)]]
  surveys <- bird_surveys$detection_array[,,i]
  point_names <- bird_surveys$point_names[rowSums(surveys, na.rm = T) > 0]
  points <- bird_points[bird_points$point %in% point_names, ]
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

point_distances <- as.data.frame(matrix(data = NA, 
                                        nrow = nrow(initial_species_list), 
                                        ncol = nrow(bird_points)))
names(point_distances) <- bird_points$point
row.names(point_distances) <- initial_species_list$HBW

# This is just for checking whether anything is out-of-range at 160km scale, 
# not for actually assigning distance covariates
for(i in 1:length(ayerbe_buffered_ranges_updated)){
  print(i)
  point_distances[i, ] <- st_distance(bird_points, ayerbe_buffered_ranges_updated[[i]])
}

saveRDS(point_distances, file = "outputs/point_distances_biogeographic_clip_ayerbe.RDS")
