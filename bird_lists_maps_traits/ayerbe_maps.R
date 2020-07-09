library(sf)
`%ni%` <- Negate(`%in%`)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"


# ayerbe_files <- list.files('/Users/jacobsocolar/Google Drive/Simon_data/data/rangemaps_Quinones')
# ayerbe_shpfiles <- ayerbe_files[grep('.shp$', ayerbe_files)]
# 
# initial_map_list <- list()
# for(i in 1:length(ayerbe_shpfiles)){
#   initial_map_list[[i]] <- st_make_valid(st_read(paste0('/Users/jacobsocolar/Google Drive/Simon_data/data/rangemaps_Quinones/', ayerbe_shpfiles[i])))
# }
# initial_map_sf <- do.call(rbind, initial_map_list)
# rm(initial_map_list)
# 
# names(initial_map_sf)[1] <- 'Species'
# initial_map_sf$Species <- gsub('[[:space:]]', ' ', initial_map_sf$Species) # At least one entry contains some kind of crazy whitespace character that gets printed as " " in the R console but which does not match " " for the purposes of string matching
# 
# saveRDS(initial_map_sf, file = '/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/ayerbe_maps/initial_map_sf.RDS')

initial_map_sf <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/ayerbe_maps/initial_map_sf.RDS')
initial_species_list <- read.csv("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/species_list_creation/initial_species_list.csv")
initial_species_list$HBW_underscore <- gsub(" ", "_", initial_species_list$HBW)
# Handle clear-cut one-to-one (from a Colombian perspective) synonymy that is resolved by HBW/eBird/EltonTraits lookup
for(i in 1:nrow(initial_map_sf)){
  if((initial_map_sf$Species[i] %ni% initial_species_list$HBW) &             #  if the Ayerbe name is missing from HBW    
        (sum(initial_species_list$eBird == initial_map_sf$Species[i]) == 1)){  # but is present exactly once in the eBird synonymy
    initial_map_sf$Species[i] <- initial_species_list$HBW[initial_species_list$eBird == initial_map_sf$Species[i]]   # Replace Ayerbe name with corresponding HBW name
  }
}

for(i in 1:nrow(initial_map_sf)){
  if((initial_map_sf$Species[i] %ni% c(initial_species_list$HBW, initial_species_list$eBird)) &   #  if the Ayerbe name is missing from HBW & eBird
     (sum(initial_species_list$eltontraits == initial_map_sf$Species[i]) == 1)){  # but is present exactly once in the EltonTraits synonymy
    initial_map_sf$Species[i] <- initial_species_list$HBW[initial_species_list$eltontraits == initial_map_sf$Species[i]]   # Replace Ayerbe name with corresponding HBW name
  }
}

# Handle additional clear-cut one-to-one synonymy
initial_map_sf$Species[initial_map_sf$Species == "Oxypogon stubelii"] <- "Oxypogon stuebelii"
initial_map_sf$Species[initial_map_sf$Species == "Geospizopis unicolor"] <- "Geospizopsis unicolor"
initial_map_sf$Species[initial_map_sf$Species == "Picumnus olivaceus_M"] <- "Picumnus olivaceus"
initial_map_sf$Species[initial_map_sf$Species == "Tinamus tao_1"] <- "Tinamus tao"
initial_map_sf$Species[initial_map_sf$Species == "Scytalopus latrans_M"] <- "Scytalopus latrans"
initial_map_sf$Species[initial_map_sf$Species == "Picumnus spilogaster orinocensis VU"] <- "Picumnus spilogaster"
initial_map_sf$Species[initial_map_sf$Species == "Parkesia noveboracensis_M"] <- "Parkesia noveboracensis"
initial_map_sf$Species[initial_map_sf$Species == "Crypturellus soui_1"] <- "Crypturellus soui"
initial_map_sf$Species[initial_map_sf$Species == "Hylopezus macularius_"] <- "Hylopezus macularius"

# A weird case: Zimmerius vilissimus gets split to parvus and improbus, but because Elton, eBird, and HBW all split improbus but
# keep parvus lumped with vilissimus, the autoconversion changes vilissimus directly to parvus.  Here we change it back so that
# we can handle the splitting later on.
initial_map_sf$Species[initial_map_sf$Species == "Zimmerius parvus"] <- "Zimmerius vilissimus"


initial_species_list$HBW[initial_species_list$HBW %ni% initial_map_sf$Species]

# Missing species
# Rallus limicola
# Fulica ardesiaca
# Gallinago imperialis
# Hoploxypterus (=Vanellus) cayanus
# Myiothlypis conspicillata                 *
# Cacicus sclateri
# Dacnis berlepschi                         
# Entomodestes coracinus                    *
# Grallaria_bangsi                          *
# Larus serranus
# Leptotila verreauxi                       *
# Megascops clarkii
# Metallura iracunda
# Metriopela melanoptera
# Nothocercus julius                        *
# Nothoprocta curvirostris
# Veniliornis callonotus
# Pyrrhura viridicata
# Trogon melanurus                          *
# Drymophila hellmayri                      *
# Aprositornis disjuncta
# Muscisaxicola fulviatilis
# Pheugopedius columbianus (= sclateri)
# Snowornis subalaris
# Thripadectes virgaticeps                  *
# Cyanolyca turcosa
# Scytalopus sanctaemartae                  *
# Thripadectes ignobilis                    *
# Troglodytes ochraceus
# Hemitriccus inornatus
# Heterocercus aurantiivertex
# Cyphorhinus phaeocephalus


# Splits
# Pyrrhura chapmani/pacifica/melanura                                 Masks Done *
# Vireo chivi/olivaceus                                               
# Coeligena conradii/torquata                                         Masks Done *
# Coeligena consita/bonaparteii                                       Masks Done
# Forpus spengeli/xanthopterygius                                     Masks Done
# Galbula chalcocephala/albirostris                                   Masks Done * Based range of albirostris on ML photos from Mitu and HBW mention of E. Meta; assumed that chalcocephala must have contiguous range within Amazonian forest around the edge of albirostris
# Hypnelus bicinctus/ruficollis                                       Masks Done * HBW mention of hybrids in Catatumbo, but everyone treats ssp coloratus with the nominate group
# Sittasomus griseus/griseicapillus                                   Masks Done
# Cyanoloxia (=Cyanocompsa) cyanoides/rothschildii                    Masks Done
# Dacnis egregia/lineata                                              Masks Done
# Pyrrhura subandina/caeruleiceps (both from P. picta)                Masks Done
# Notharchus subtectus/tectus                                         Masks Done
# Cyanolyca quindiuna/armillata                                       Masks Done
# Arremon axillaris/taciturnus                                        Masks Done
# Grallaria alticola/quitensis                                        Masks Done
# Pteroglossus sanguineus/torquatus                                   Masks Done *
# Ramphastos culminatus/citreolaemus (both from vitellinus)           Masks Done
# Schistes albogularis/geoffroyi                                      Masks Done *
# Atlapetes nigrifrons/latinuchus                                     Masks Done
# Urochroa leucura/bougueri                                           Masks Done
# Vireolanius leucotis/mikettae                                       Masks Done
# Deconychura pallida/typica (both from longicauda)                   Masks Done
# Dendrocolaptes punctipectus/sanctithomae                            Masks Done *
# Formicivora intermedia/grisea                                       Masks Done
# Furnarius longirostris/cinnamomeus/leucopus                         Masks Done
# Grallaria salutensis/rufula                                         Masks Done
# Automolus virgatus/subulatus                                        Masks Done
# Myiobius sulphureipygius/barbatus                                   Masks Done * semiflavus debacle
# Myiodynastes solitarius/maculatus                                   Masks Done *
# Myrmornis stictoptera/torquata                                      Masks Done
# Myiopagis parambae/cinerea (both from caniceps)                     Masks Done
# Myiophobus crypterythrus/fasciatus                                  Masks Done
# Myiotriccus phoenicurus/ornatus                                     Masks Done
# Sakesphorus pulchellus/canadensis                                   Masks Done
# Tolmomyias viridiceps/flaviventris                                  Masks Done *
# Thamnistes aequatorialis/anabatinus                                 Masks Done
# Tolmomyias flavotectus/assimilis                                    Masks Done
# Turdus debilis/arthuri/ignobilis                                    Masks Done * Stiles & Avedaño for arthuri
# Zimmerius improbus/parvus (both from vilissimus, which gets autocoverted to parvus because improbus is consistent in eBird and Elton)  Masks Done
# Xiphorhynchus chunchotambo/beauperthuysii (both from ocellatus)     Masks Done * beauperthuysii needs a sliver of buffer
# Cryptopipo litae/holochlora                                         Masks Done
# Cacicus flavicrissus/cela                                           Masks Done * contra maps in Birdlife and Donegan, but consistent with HBW range descriptions, Jaramillo, and common sense, here treated Catatumbo populations as cela
# Cacicus uropygialis/pacificus                                       Masks Done
# Dubusia carrikeri/taeniata                                          Masks Done
# Habia (=Chlorothraupis) frenata/carmioli                            Masks Done
# Myioborus chrysops/ornatus                                          Masks Done
# Tangara lunigera/parzudakii                                         Masks Done
# Xiphorhynchus guttatoides/guttatus                                  Masks Done * eBird photos look like guttatoides near Inirida and Arauca; HBW indicates that guttatus is know from eastern Vichada
# Amazona bodini/festiva                                              Masks Done
# Manacus vitellinus/manacus                                          Masks Done * Birdlife shows zone of ambiguous overlap; possibly based on photo-supported eBird record near Puerto Liberatdores that appears to be a hybrid.  Here drew the line a bit west of that location.  Buffering will create 60 km overlap band anyway.
# Onychorhynchus mexicanus/coronatus                                  Masks Done
# Zimmerius albigularis/chrysops                                      Masks Done * eBird, XC, and ML indicate that these species are not allopatric replacements, and range of albigularis is much wider than shown by BirdLife or Donegan; range of chrysops also extends in to areas mapped exclusively for albigularis by these sources.  Masking on this basis.
# Aulacorhynchus albivitta/caeruleogularis (both from prasinus)       Masks Done
# Campephilus splendens/hameatogaster                                 Masks Done
# Saltator olivascens/coerulescens                                    Masks Done * line based on Boseman's HBW alive note plus eBird photos from Casanare, continuous cluster of eBird records from Venezuela down the white-sand savannas of NE Colombia, and XC material from Venezuela
# Ramphocelus icteronotus/flammigerus                                 Masks Done *

mask_update <- data.frame(new_species = list.files('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/ayerbe_updates/ayerbe_mask'),
                          path = list.files('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/ayerbe_updates/ayerbe_mask', full.names = T),
                          old_species = NA)

for(i in 1:nrow(mask_update)){
  nf <- list.files(mask_update$path[i])
  mask_update$old_species[i] <- gsub('.txt', '', nf[grep('.txt', nf)])
}

new_map_list <- list()
for(i in 1:nrow(mask_update)){
  old_map <- initial_map_sf[initial_map_sf$Species == mask_update$old_species[i],]
  files <- list.files(mask_update$path[i])
  mps <- files[grep('.kml', files)]
  masking_poly <- st_zm(st_read(paste0(mask_update$path[i], '/', mps[1])))
  new_map <- st_difference(old_map, masking_poly)[,c(1,4)]
  if(length(mps) > 1){
    for(j in 2:length(mps)){
      masking_poly <- st_zm(st_read(paste0(mask_update$path[i], '/', mps[j])))
      new_map <- st_difference(new_map, masking_poly)[,c(1,4)]
    }
  }
  new_map$Species <- gsub('_', ' ', mask_update$new_species[i])
#  plot(new_map, main = new_map$Species)
#  plot(st_buffer(new_map, dist = .03), main = new_map$Species)
  new_map_list[[i]] <- new_map
}

ayerbe_splits <- do.call(rbind, new_map_list)
  
ayerbe_maps_prelim <- initial_map_sf[initial_map_sf$Species %ni% mask_update$old_species, ]

ayerbe_maps <- st_transform(rbind(ayerbe_maps_prelim, ayerbe_splits), AEAstring)

bird_surveys <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_surveys_current.RDS')
colombia_points <- st_as_sf(data.frame(readxl::read_excel('/Users/jacobsocolar/Downloads/RAWdata_samplingpointsCO_MASTER.xlsx')), 
                            coords = c('long', 'lat'), crs = 4326)
for(i in 10:12){
  colombia_points$point_id <- gsub(paste0('CHA', i, 'd'), paste0('CHA', i, 'D'), colombia_points$point_id)
}
bird_points <- st_transform(colombia_points[colombia_points$point_id %in% bird_surveys$point_names, ], AEAstring)

missing <- c('Myiothlypis_conspicillata', 'Entomodestes_coracinus', 'Grallaria_bangsi', 'Leptotila_verreauxi',
             'Nothocercus_julius', 'Trogon_melanurus', 'Drymophila_hellmayri', 'Thripadectes_virgaticeps',
             'Scytalopus_sanctaemartae', 'Thripadectes_ignobilis', 'Vireo_olivaceus')


load('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/birdlife_maps/recast_range_maps.Rdata')
birdlife_relevant <- st_transform(recast_range_maps[recast_range_maps$SCINAME %in% ayerbe_maps$Species,], AEAstring)

distances <- distances_birdlife <- list()
dsp <- vector()
sp_pts <- list()
for(i in 1:length(bird_surveys$species_names)){
  if(bird_surveys$species_names[i] %ni% missing){
    species <- bird_surveys$species_names[i]
    species_range <- ayerbe_maps[ayerbe_maps$Species == gsub('_', ' ', species), ]
    species_range_birdlife <- birdlife_relevant[birdlife_relevant$SCINAME == gsub('_', ' ', species), ]
    surveys <- bird_surveys$detection_array[,,i]
    point_names <- bird_surveys$point_names[rowSums(surveys, na.rm = T) > 0]
    points <- bird_points[bird_points$point_id %in% point_names, ]
    distance_matrix <- st_distance(points, species_range)
    distance_matrix_birdlife <- st_distance(points, species_range_birdlife)
    distances[[i]] <- as.vector(distance_matrix)
    distances_birdlife[[i]] <-  apply(distance_matrix_birdlife, 1, min)
    dsp[i] <- species
    sp_pts[[i]] <- points
  }
}

mdist <- mdist2 <- vector()
for(i in 1:length(dsp)){
  mdist[i] <- max(distances[[i]])
  mdist2[i] <- max(pmin(distances[[i]], distances_birdlife[[i]]))
}

d <- 50000
bad_species <- dsp[mdist > d]
bad_species2 <- dsp[mdist2 > d]


colombia <- st_read('/Users/JacobSocolar/Dropbox/Work/Colombia/Data/GIS/colombia_maps/gadm36_COL_shp/gadm36_COL_0.shp')

# file consists of many disjoint polygons, representing the mainland and numerous islands. Figure out which is the mainland
# and extract it.
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
m2 <- st_simplify(mainland, dTolerance = 10000)

# for(i in 1:length(bad_species2)){
#   print(ggplot(m2) + geom_sf() + 
#     geom_sf(data = birdlife_relevant[birdlife_relevant$SCINAME == gsub("_", " ", bad_species2[i]),], inherit.aes = F, fill = 'cornflowerblue') +
#     geom_sf(data = ayerbe_maps[ayerbe_maps$Species == gsub("_", " ", bad_species2[i]),], inherit.aes = F, fill = 'red', alpha = .2) +
#     geom_sf(data = sp_pts[[which(dsp == bad_species2[i])]], inherit.aes = F) +
#     ggtitle(bad_species2[i]))
# }



clipping_polygons <- list(amazon_orinoco = st_transform(amazon_orinoco, AEAstring), magdalena_east = st_transform(magdalena_east, AEAstring),
                          magdalena_west = st_transform(magdalena_west, AEAstring), cauca_east = st_transform(cauca_east, AEAstring),
                          cauca_west = st_transform(cauca_west, AEAstring), pasto = st_transform(pasto, AEAstring), 
                          pacific = st_transform(pacific, AEAstring), snsm = st_transform(snsm, AEAstring), 
                          guajira_perija = st_transform(guajira_perija, AEAstring), catatumbo = st_transform(catatumbo, AEAstring))

rrm <- recast_range_maps[recast_range_maps$SCINAME %in% initial_species_list$HBW,]

birdlife_relevant <- st_intersection(st_make_valid(st_transform(rrm, AEAstring)), buffered_mainland)
saveRDS(birdlife_relevant, file = '/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/combined_maps/birdlife_relevant.RDS')

rangeclips <- list()
bufferclips <- list()
for(i in 1:length(clipping_polygons)){
  cp <- clipping_polygons[[i]]
  clip_list <- buffer_list <- list()
  counter <- 0
  for(k in 1:length(bird_surveys$species_names)){
    print(c(i, k))
    species <- bird_surveys$species_names[k]
    range <- birdlife_range <-  st_union(st_make_valid(birdlife_relevant[birdlife_relevant$SCINAME == gsub('_', ' ', species), ]))
    if(bird_surveys$species_names[k] %ni% missing){
      ayerbe_range <- ayerbe_maps[ayerbe_maps$Species == gsub('_', ' ', species), ]
      range <- st_union(st_union(ayerbe_range), birdlife_range)
    }
    clipped_range <- st_intersection(cp, range)
    if(nrow(clipped_range) > 0){
      if(nrow(clipped_range) != 1){stop()}
      counter <- counter + 1
      clipped_range$species <- species
      cr <- clipped_range[, 'species']
      names(cr) <- c('species', 'geometry')
      st_geometry(cr) <- 'geometry'
      clip_list[[counter]] <- cr
      
      buffered_cr <- st_buffer(cr, 100000)
      clipped_buffered_cr <- st_intersection(buffered_cr, cp)
      clipped_buffered_cr$species <- species
      cbc <- clipped_buffered_cr[, 'species']
      names(cbc) <- c('species', 'geometry')
      st_geometry(cbc) <- 'geometry'
      buffer_list[[counter]] <- cbc
    }
  }
  rangeclips[[i]] <- do.call(rbind, clip_list)
  bufferclips[[i]] <- do.call(rbind, buffer_list)
}
names(rangeclips) <- names(bufferclips) <- names(clipping_polygons)

bird_points$bioregion <- NA
for(i in 1:nrow(bird_points)){
  v <- rep(NA, length(clipping_polygons))
  for(j in 1:length(clipping_polygons)){
    v[j] <- st_intersects(clipping_polygons[[j]], bird_points[i, ], sparse = F)
  }
  if(sum(v) != 1){stop()}
  bird_points$bioregion[i] <- names(clipping_polygons)[which(v)]
}
head(bird_points)


point_distances <- as.data.frame(matrix(data = NA, nrow = nrow(initial_species_list), ncol = nrow(bird_points)))
names(point_distances) <- bird_points$point_id
row.names(point_distances) <- initial_species_list$HBW

for(i in 1:nrow(bird_points)){
  print(i)
  region_range <- rangeclips[[bird_points$bioregion[i]]]
  for(k in 1:nrow(initial_species_list)){
    if(initial_species_list$HBW_underscore[k] %in% region_range$species){
      point_distances[k, i] <- st_distance(bird_points[i, ], region_range[region_range$species == initial_species_list$HBW_underscore[k], ])
    }
  }
}

saveRDS(point_distances, file = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/point_distances/point_distances_biogeographic_clip.RDS")

Q <- apply(bird_surveys$detection_array, MARGIN = c(1,3), FUN = function(x){return(sum(x, na.rm = T) > 0)})

pt_distances <- data.frame(species = rep(NA, sum(Q)), point = rep(NA, sum(Q)), distance = rep(NA, sum(Q)))
counter <- 0
for(i in 1:length(bird_surveys$point_names)){
  pt <- bird_surveys$point_names[i]
  for(k in 1:length(bird_surveys$species_names)){
    if(Q[i, k] == 1){
      counter <- counter + 1
      sp <- gsub("_", " ", bird_surveys$species_names[k])
      pt_distances$species[counter] <- sp
      pt_distances$point[counter] <- pt
      pt_distances$distance[counter] <- point_distances[sp, pt]
    }
  }
}

pt_distances[is.na(pt_distances$distance), ]

pd2 <- pt_distances[!is.na(pt_distances$distance), ]

pd2[pd2$distance > 50000,]

ggplot(birdlife_relevant[birdlife_relevant$SCINAME == "Aburria aburri", ]) + geom_sf() + 
  geom_sf(data = ayerbe_maps[ayerbe_maps$Species == "Aburria aburri", ], inherit.aes = F) + 
  geom_sf(data = region_range[region_range$species == initial_species_list$HBW_underscore[k], ], fill = 'red', inherit.aes = F)
