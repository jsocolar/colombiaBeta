# This script reads in the Ayerbe maps and resolves taxonomy with HBW

###### Script dependencies: Species_lists.R #####

# housekeeping ----
library(sf)

`%ni%` <- Negate(`%in%`)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 
    +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

##### Get colombia mainland shapefile #####
colombia <- st_read('inputs/gadm36_COL_shp/gadm36_COL_0.shp')
# file consists of many disjoint polygons, representing the mainland and numerous 
# islands. Figure out which is the mainland and extract it.
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


##### Read Ayerbe shapefiles provided by Humboldt, and bind to sf frame #####
ayerbe_files <- list.files('inputs/rangemaps_Quinones')
ayerbe_shpfiles <- ayerbe_files[grep('.shp$', ayerbe_files)]

ayerbe_sf_fname <- 'outputs/initial_ayerbe_map_sf.RDS'
if(!exists(ayerbe_sf_fname)) {
    initial_map_list <- list()
    for(i in 1:length(ayerbe_shpfiles)){
        initial_map_list[[i]] <- st_make_valid(
            st_read(paste0('inputs/rangemaps_Quinones/', ayerbe_shpfiles[i]))
        )
    }
    initial_map_sf <- do.call(rbind, initial_map_list)
    rm(initial_map_list)
    names(initial_map_sf)[1] <- 'Species'
    # At least one entry contains some kind of crazy whitespace character that 
    # gets printed as " " in the R console but which does not match " " for the 
    # purposes of string matching
    initial_map_sf$Species <- gsub('[[:space:]]', ' ', initial_map_sf$Species) 
    saveRDS(initial_map_sf, file = ayerbe_sf_fname)   
}

##### Standardize taxonomy to HBW species list #####
initial_map_sf <- readRDS(ayerbe_sf_fname)
# ?? directory issue
initial_species_list <- read.csv("inputs/initial_species_list.csv")
initial_species_list$HBW_underscore <- gsub(" ", "_", initial_species_list$HBW)

# One-to-one (from a Colombian perspective) synonymy that is resolved by 
# HBW/eBird/EltonTraits lookup
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

# Additional one-to-one synonymy
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

##### What HBW species are now missing a map? #####
initial_species_list$HBW[initial_species_list$HBW %ni% initial_map_sf$Species]
# Missing species
# Rallus limicola
# Fulica ardesiaca
# Gallinago imperialis
# Hoploxypterus (=Vanellus) cayanus
# Myiothlypis conspicillata                 
# Cacicus sclateri
# Dacnis berlepschi                         
# Entomodestes coracinus                    
# Grallaria_bangsi                          
# Larus serranus
# Leptotila verreauxi                       
# Megascops clarkii
# Metallura iracunda
# Metriopela melanoptera
# Nothocercus julius                        
# Nothoprocta curvirostris
# Veniliornis callonotus
# Pyrrhura viridicata
# Trogon melanurus                          
# Drymophila hellmayri                      
# Aprositornis disjuncta
# Muscisaxicola fulviatilis
# Pheugopedius columbianus (= sclateri)
# Snowornis subalaris
# Thripadectes virgaticeps                  
# Cyanolyca turcosa
# Scytalopus sanctaemartae                  
# Thripadectes ignobilis                    
# Troglodytes ochraceus
# Cyphorhinus phaeocephalus

# In the below, * indicates splits that are mapped as parapatric by Ayerbe, such that the location of the division
# was determined from literature rather than self-evident from the map.
# Splits
# Pyrrhura chapmani/pacifica/melanura                                 Masks Done *
# Vireo chivi/olivaceus                                               Ayerbe doesn't explicitly map olivaceus; we'll just rely on the BirdLife map to pick it up
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
# Myiobius sulphureipygius/barbatus                                   Masks Done * semiflavus is an issue; here followed Todd and others in treating it as part of barbatus, and therfore parapatric with sulphureipygus
# Myiodynastes solitarius/maculatus                                   Masks Done *
# Myrmornis stictoptera/torquata                                      Masks Done
# Myiopagis parambae/cinerea (both from caniceps)                     Masks Done
# Myiophobus crypterythrus/fasciatus                                  Masks Done
# Myiotriccus phoenicurus/ornatus                                     Masks Done
# Sakesphorus pulchellus/canadensis                                   Masks Done
# Tolmomyias viridiceps/flaviventris                                  Masks Done *
# Thamnistes aequatorialis/anabatinus                                 Masks Done
# Tolmomyias flavotectus/assimilis                                    Masks Done
# Turdus debilis/arthuri/ignobilis                                    Masks Done * Stiles & Avedaño is consulted to delimit arthuri range
# Zimmerius improbus/parvus (both from vilissimus, which gets autocoverted to parvus because improbus is consistent in eBird and Elton)  Masks Done
# Xiphorhynchus chunchotambo/beauperthuysii (both from ocellatus)     Masks Done
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


##### Splits #####
# Splits are handled by creating masking polygons that cover up superflous parts 
# of the parent taxon's range.  These polygons were drawn by hand in Google Earth 
# with the relevant Ayerbe maps simultaneously imported.

# Get dataframe giving the HBW species, the file path for the masking polygons, and the Ayerbe species to which the mask must be applied
mask_update <- data.frame(new_species = list.files('inputs/ayerbe_updates/ayerbe_mask'),
                          path = list.files('inputs/ayerbe_updates/ayerbe_mask', full.names = T),
                          old_species = NA)
for(i in 1:nrow(mask_update)){
  nf <- list.files(mask_update$path[i])
  mask_update$old_species[i] <- gsub('.txt', '', nf[grep('.txt', nf)])
}

# Do the masking
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
  new_map_list[[i]] <- new_map
}

ayerbe_splits <- do.call(rbind, new_map_list)

##### Missing species #####
missing_species <- data.frame(species = list.files('ayerbe_updates/ayerbe_missing'),
                              path = list.files('ayerbe_updates/ayerbe_missing', full.names = T))
missing_map_list <- list()
counter <- 0
for(i in 1:nrow(missing_species)){
  files <- list.files(missing_species$path[[i]], full.names = T)
  if(sum(grepl('empty.txt', files)) == 0){
    counter <- counter + 1
    range_files <- files[grep('.kml', files)]
    rangemap <- st_zm(st_read(range_files[1]))
    if(length(range_files) > 1){
      for(j in 2:length(range_files)){
        rangemap2 <- st_zm(st_read(range_files[j]))
        rangemap <- st_union(rangemap, rangemap2)
      }
    }
    if(sum(grepl('HOLE', files)) != 0){
      holefiles <- list.files(files[grep('HOLE', files)], full.names = T)
      hole <- st_zm(st_read(holefiles[1]))
      rangemap <- st_difference(rangemap, hole)
      if(length(holefiles) > 1){
        for(j in 2:length(holefiles)){
          hole <- st_zm(st_read(holefiles[j]))
          rangemap <- st_difference(rangemap, hole)
        }
      }
    }
    rangemap$Species <- gsub('_', ' ', missing_species$species[i])
    missing_map_list[[counter]] <- rangemap[, c('Species', 'geometry')]
  }
}

ayerbe_missing_prelim <- do.call(rbind, missing_map_list)
# Crop to mainland, as hand-drawn ranges were not clipped to Colombian border
ayerbe_missing <- st_intersection(ayerbe_missing_prelim, 
                                  st_transform(mainland, st_crs(ayerbe_missing_prelim))
                                  ) 

##### Combine for complete Ayerbe map set #####
# note that there is no map, and no line in the sf frame, for Vireo olivaceus or 
# Troglodytes ochraceus
ayerbe_maps_prelim <- initial_map_sf[initial_map_sf$Species %ni% 
                                         mask_update$old_species, ]
ayerbe_maps_prelim2 <- rbind(ayerbe_maps_prelim, ayerbe_splits)
ayerbe_maps <- st_transform(rbind(ayerbe_maps_prelim2, ayerbe_missing), AEAstring)
ayerbe_maps[ayerbe_]
initial_species_list$HBW[initial_species_list$HBW %ni% ayerbe_maps$Species]

saveRDS(ayerbe_maps, file = 'outputs/ayerbe_maps/ayerbe_maps.RDS')
