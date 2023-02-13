# This script formats unified table of species covariates from the elevations file, 
# range maps, eltonTraits, & Parker

##### Script dependencies: parker_standardization.R, species_lists.R, 
# combined_bird_maps.R, elevations_prep_and_exploration.R

library(sf); library(ggplot2)
`%ni%` <- Negate(`%in%`)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 
+y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

##### Begin with elevations file #####
ef <- read.csv("inputs/elevations_final.csv")
ef$latin_underscore <- gsub(" ", "_", ef$latin)
elevations_final <- ef[is.na(ef$omit & !is.na(ef$is_na)), c("latin", "latin_underscore", "lower", "upper")]
elevations_final$elev_breadth <- elevations_final$upper - elevations_final$lower

##### Add basic info about biogeographic restrictions from Ayerbe maps #####
ayerbe_maps_updated <- readRDS('outputs/ayerbe_buffered_ranges_updated.RDS')
# Create traits file
traits <- elevations_final[elevations_final$latin %in% names(ayerbe_maps_updated), ]
traits$east_only <- traits$west_only <- traits$eandes_absent <- traits$wandes_absent <- traits$snsm_only <- NA

# Get biogeographic polygons
source('code/GIS_processing/hydrosheds_extraction.R')
cps <- list(amazon_orinoco = st_transform(amazon_orinoco, AEAstring), magdalena_east = st_transform(magdalena_east, AEAstring),
            magdalena_west = st_transform(magdalena_west, AEAstring), cauca_east = st_transform(cauca_east, AEAstring),
            cauca_west = st_transform(cauca_west, AEAstring), pasto = st_transform(pasto, AEAstring), 
            pacific = st_transform(pacific, AEAstring), snsm = st_transform(snsm, AEAstring), 
            guajira_perija = st_transform(guajira_perija, AEAstring), catatumbo = st_transform(catatumbo, AEAstring))

# Everything east of E Andes
east_polygon <- st_union(cps$amazon_orinoco, cps$catatumbo)
# Everything west of E Andes
west_polygon <- st_union(cps$magdalena_east, st_union(cps$magdalena_west, st_union(cps$cauca_east, st_union(cps$cauca_west,
                         st_union(cps$pacific, st_union(cps$snsm, st_union(cps$guajira_perija, cps$pasto)))))))
# Everything not on SN Santa Marta
not_snsm <- st_union(east_polygon, st_union(cps$magdalena_east, st_union(cps$magdalena_west, st_union(cps$cauca_east, st_union(cps$cauca_west,
                     st_union(cps$pacific, st_union(cps$guajira_perija, cps$pasto)))))))

# West Andes (both slopes)
wandes_polygon <- st_union(cps$cauca_west, cps$pacific)
# East Andes (both slopes)
eandes_polygon <- st_union(cps$amazon_orinoco, cps$magdalena_east)
# Central Andes (both slopes)
candes_polygon <- st_union(cps$magdalena_west, cps$cauca_east)

# read in GADM colombia shapefile
colombia <- st_read('inputs/gadm36_COL_shp/gadm36_COL_0.shp')

# File consists of many disjoint polygons, representing the mainland and numerous 
# islands. Figure out which is the mainland and extract it.
npoly <- length(colombia$geometry[[1]])
size <- rep(0,npoly)
for(i in 1:npoly){
  size[i] <- dim(colombia$geometry[[1]][[i]][[1]])[1]
}

mainland <- colombia$geometry[[1]][[which(size == max(size))]] %>%
  st_polygon() %>% st_sfc() %>% st_sf()
st_crs(mainland) <- st_crs(colombia)

# Transform to AEA conic centered on Colombia, for accurate 100 km buffering
mainland <- st_transform(mainland, AEAstring)

# Populate biogeographic restriction columns
for(i in 1:nrow(traits)){
  print(i)
  species <- traits$latin[i]
  map <- st_buffer(st_intersection(st_union(ayerbe_maps_updated[[species]]), mainland), -100) # inward buffer of 100 meters to make sure that adjacent polygons don't appear to overlap
  traits$east_only[i] <- as.numeric(st_intersects(east_polygon, map, sparse = F)[1,1] & !st_intersects(west_polygon, map, sparse = F)[1,1])   # Species east but not west of E Andes
  traits$west_only[i] <- as.numeric(!st_intersects(east_polygon, map, sparse = F)[1,1] & st_intersects(west_polygon, map, sparse = F)[1,1])   # Species west but not east of E Andes
  traits$wandes_absent[i] <- as.numeric(!st_intersects(wandes_polygon, map, sparse = F)[1,1] & (st_intersects(eandes_polygon, map, sparse = F)[1,1] | st_intersects(candes_polygon, map, sparse = F)[1,1]) & !traits$east_only[i])  # Species not on either slope of W Andes
  traits$eandes_absent[i] <- as.numeric(!traits$wandes_absent[i] & (st_intersects(wandes_polygon, map, sparse = F)[1,1] | st_intersects(candes_polygon, map, sparse = F)[1,1]) & !st_intersects(eandes_polygon, map, sparse = F)[1,1])  # Species not on either slope of E Andes
  traits$snsm_only[i] <- as.numeric(st_intersects(cps$snsm, map, sparse = F)[1,1] & !st_intersects(not_snsm, map, sparse = F)[1,1])    # Santa Marta endemics (in a Colombian context)
}

saveRDS(traits, "outputs/traits_prelim.RDS")

##### Add Family and Genus #####
birdlife_list <- readxl::read_xlsx("inputs/HBW-BirdLife_Checklist_v4_Dec19/Handbook of the Birds of the World and BirdLife International digital checklist of the birds of the world_Version_4.xlsx",
                                 skip = 1)
birdlife_v4 <- birdlife_list[!is.na(birdlife_list$Sequence),]
families <- birdlife_v4[,c("Scientific name", "Order", "Family name")]
traits <- merge(traits, families, by.x = "latin", by.y = "Scientific name", all.x = T)
traits$Genus <- sapply(strsplit(as.character(traits$latin), " "), "[[", 1)
names(traits)[names(traits) == "Family name"] <- "Family"
traits$Family[is.na(traits$Family) & traits$Genus == "Grallaria"] <- "Grallariidae"
traits$Order[is.na(traits$Order) & traits$Genus == "Grallaria"] <- "PASSERIFORMES"

##### Add Parker traits #####
parker <- read.csv('inputs/new_parker.csv')
for(j in 5:ncol(parker)){
  for(i in 1:nrow(parker)){
    if(parker[i,j] %in% c("Y", "Q", "E")){
      parker[i,j] <- 1
    }else{
      parker[i,j] <- 0
    }
  }
  parker[,j] <- as.numeric(parker[,j])
}

F_ind <- which(names(parker) == "F1") - 1
N_ind <- which(names(parker) == "N1") - 1
A_ind <- which(names(parker) == "A1") - 1

parker$forest_present <- as.numeric(rowSums(parker[ , (F_ind+1):(F_ind+13)]) > 0)
parker$forest_specialist <- as.numeric(parker$forest_present == 1 & rowSums(parker[ , (N_ind+1):(A_ind+12)]) == 0)
parker$tf_specialist <- as.numeric(parker$F1 == 1 & ((parker$F2 + parker$F3 + parker$F8) == 0))
parker$dry_forest_present <- parker$F7
parker$flood_dry_specialist <- as.numeric(parker$forest_specialist == 1 & (parker$F1 + parker$F4 + parker$F5 + parker$F6 + parker$F12) == 0 & ((parker$F2 + parker$F3 + parker$F7 + parker$F8) != 0))
parker$floodplain_specialist <- as.numeric(parker$flood_dry_specialist == 1 & (parker$F7 == 0))
parker$arid_present <- as.numeric((parker$N1 + parker$N2) != 0)

parker_pairstest <- parker[,46:52]
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(parker_pairstest, lower.panel = panel.cor, upper.panel = panel.smooth)

parker_traits <- parker[,names(parker) %ni% c("X", "parker", "parker2")]

traits <- merge(traits, parker_traits, by.x = "latin", by.y = "HBW", all.x = T)

##### Add BirdLife traits #####
load('outputs/birdlife_traits.Rdata')

bts <- data.frame(latin = birdlife_traits$names, generation_birdlife = 0, 
                  mass_birdlife = unlist(birdlife_traits$body_mass), 
                  migrat_birdlife = unlist(birdlife_traits$migratory_status),
                  forestDep_birdlife = unlist(birdlife_traits$forest_dep))
for(i in 1:nrow(bts)){
  if(length(birdlife_traits$generation[[i]]) > 0){
    bts$generation_birdlife[i] <- birdlife_traits$generation[[i]]
  }else{
    bts$generation_birdlife[i] <- NA
  }
}


traits <- merge(traits, bts, by.x = "latin", by.y = "latin", all.x = T)
traits$migrat_birdlife[is.na(traits$migrat_birdlife)] <- "not a migrant"

# Add habitats from birdlife
birdlife_habitats <- data.frame(latin = birdlife_traits$names, 
                                birdlife_hab_forest_trop_moist_lowland = NA, 
                                birdlife_hab_forest_trop_moist_montane = NA, 
                                birdlife_hab_forest_trop_dry = NA, 
                                birdlife_hab_forest_tropical_swamp = NA,
                                birdlife_hab_savanna = NA,
                                birdlife_hab_shrubland_trop_dry = NA, 
                                birdlife_hab_shrubland_trop_high = NA, 
                                birdlife_hab_shrubland_trop_moist = NA,
                                birdlife_hab_grassland_trop_dry = NA, 
                                birdlife_hab_grassland_trop_high = NA, 
                                birdlife_hab_grassland_trop_wet = NA, 
                                birdlife_hab_desert_hot = NA, 
                                birdlife_hab_wetlands = NA)
for(i in 1:nrow(birdlife_habitats)){
  habs <- birdlife_traits$habitats[[i]]
  birdlife_habitats$birdlife_hab_forest_trop_moist_lowland[i] <- as.numeric(sum(grepl("Forest; Subtropical/Tropical Moist Lowland;", habs)) > 0)
  birdlife_habitats$birdlife_hab_forest_trop_moist_montane[i] <- as.numeric(sum(grepl("Forest; Subtropical/Tropical Moist Montane;", habs)) > 0)
  birdlife_habitats$birdlife_hab_forest_trop_dry[i] <- as.numeric(sum(grepl("Forest; Subtropical/Tropical Dry;", habs)) > 0)
  birdlife_habitats$birdlife_hab_forest_tropical_swamp[i] <- as.numeric(sum(grepl("Forest; Subtropical/Tropical Swamp;", habs)) > 0)
  birdlife_habitats$birdlife_hab_savanna[i] <- as.numeric(sum(grepl("Savanna; ", habs)) > 0)
  birdlife_habitats$birdlife_hab_shrubland_trop_dry[i] <- as.numeric(sum(grepl("Shrubland; Subtropical/Tropical Dry;", habs)) > 0)
  birdlife_habitats$birdlife_hab_shrubland_trop_high[i] <- as.numeric(sum(grepl("Shrubland; Subtropical/Tropical High Altitude;", habs)) > 0)
  birdlife_habitats$birdlife_hab_shrubland_trop_moist[i] <- as.numeric(sum(grepl("Shrubland; Subtropical/Tropical Dry;", habs)) > 0)
  birdlife_habitats$birdlife_hab_grassland_trop_dry[i] <- as.numeric(sum(grepl("Grassland; Subtropical/Tropical Moise;", habs)) > 0)
  birdlife_habitats$birdlife_hab_grassland_trop_high[i] <- as.numeric(sum(grepl("Grassland; Subtropical/Tropical High Altitude;", habs)) > 0)
  birdlife_habitats$birdlife_hab_grassland_trop_wet[i] <- as.numeric(sum(grepl("Grassland; Subtropical/Tropical Seasonally Wet/Flooded;", habs)) > 0)
  birdlife_habitats$birdlife_hab_desert_hot[i] <- as.numeric(sum(grepl("Desert; Hot;", habs)) > 0)
  birdlife_habitats$birdlife_hab_wetlands[i] <- as.numeric(sum(grepl("Wetlands (inland);", habs)) > 0)
}

traits <- merge(traits, birdlife_habitats, by.x = "latin", by.y = "latin", all.x = T)


##### Add EltonTraits traits #####
initial_species_list <- read.csv("outputs/initial_species_list.csv")
elton <- read.delim("inputs/elton.txt")

elton_traits <- merge(initial_species_list, elton, by.x = "eltontraits", 
                      by.y = "Scientific")[,c(2, 16:26, 30:36, 42)]

traits <- merge(traits, elton_traits, by.x = "latin", by.y = "HBW", all.x = T)

# Confirm high correlation between EltonTraits masses and Birdlife masses
summary(lm(traits$mass_birdlife ~ traits$BodyMass.Value))

saveRDS(traits, "outputs/traits.RDS")
