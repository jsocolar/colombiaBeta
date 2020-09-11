# This script formats unified table of species covariates from the elevations file, range maps, eltonTraits, & Parker

##### Script dependencies: parker_standardization.R, species_lists.R, combined_bird_maps.R, elevations_prep_and_exploration.R

library(sf)
library(ggplot2)
`%ni%` <- Negate(`%in%`)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

ef <- read.csv("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/traits/elevations/elevations_final.csv")
ef$latin_underscore <- gsub(" ", "_", ef$latin)
elevations_final <- ef[is.na(ef$omit & !is.na(ef$is_na)), c("latin", "latin_underscore", "lower", "upper")]



parker <- read.csv('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/species_list_creation/new_parker.csv')
head(elevations_final)
head(parker)
elevations_final$latin[elevations_final$latin %ni% parker$HBW]



source('/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/GIS_processing/hydrosheds_extraction.R')
cps <- list(amazon_orinoco = st_transform(amazon_orinoco, AEAstring), magdalena_east = st_transform(magdalena_east, AEAstring),
            magdalena_west = st_transform(magdalena_west, AEAstring), cauca_east = st_transform(cauca_east, AEAstring),
            cauca_west = st_transform(cauca_west, AEAstring), pasto = st_transform(pasto, AEAstring), 
            pacific = st_transform(pacific, AEAstring), snsm = st_transform(snsm, AEAstring), 
            guajira_perija = st_transform(guajira_perija, AEAstring), catatumbo = st_transform(catatumbo, AEAstring))

east_polygon <- st_union(cps$amazon_orinoco, cps$catatumbo)
west_polygon <- st_union(cps$magdalena_east, st_union(cps$magdalena_west, st_union(cps$cauca_east, st_union(cps$cauca_west,
                         st_union(cps$pacific, st_union(cps$snsm, st_union(cps$guajira_perija, cps$pasto)))))))
not_snsm <- st_union(east_polygon, st_union(cps$magdalena_east, st_union(cps$magdalena_west, st_union(cps$cauca_east, st_union(cps$cauca_west,
                     st_union(cps$pacific, st_union(cps$guajira_perija, cps$pasto)))))))


wandes_polygon <- st_union(cps$cauca_west, cps$pacific)
eandes_polygon <- st_union(cps$amazon_orinoco, cps$magdalena_east)
candes_polygon <- st_union(cps$magdalena_west, cps$cauca_east)


combined_maps_updated <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/combined_maps/combined_maps_updated.RDS')

traits <- elevations_final[elevations_final$latin %in% names(combined_maps_updated), ]
traits$east_only <- traits$west_only <- traits$eandes_absent <- traits$wandes_absent <- traits$snsm_only <- NA

for(i in 1:nrow(traits)){
  print(i)
  species <- traits$latin[i]
  map <- st_union(combined_maps_updated[[species]])
  traits$east_only[i] <- as.numeric(st_intersects(east_polygon, map, sparse = F)[1,1] & !st_intersects(west_polygon, map, sparse = F)[1,1])
  traits$west_only[i] <- as.numeric(!st_intersects(east_polygon, map, sparse = F)[1,1] & st_intersects(west_polygon, map, sparse = F)[1,1])
  traits$eandes_absent[i] <- as.numeric((st_intersects(wandes_polygon, map, sparse = F)[1,1] | st_intersects(candes_polygon, map, sparse = F)[1,1]) & !st_intersects(eandes_polygon, map, sparse = F)[1,1])
  traits$wandes_absent[i] <- as.numeric(!st_intersects(wandes_polygon, map, sparse = F)[1,1] & (st_intersects(eandes_polygon, map, sparse = F)[1,1] | st_intersects(candes_polygon, map, sparse = F)[1,1]) & !traits$east_only[i])
  traits$snsm_only[i] <- as.numeric(st_intersects(cps$snsm, map, sparse = F)[1,1] & !st_intersects(not_snsm, map, sparse = F)[1,1])
}

traits$elev_breadth <- traits$upper - traits$lower

saveRDS(traits, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/traits/traits.RDS")
