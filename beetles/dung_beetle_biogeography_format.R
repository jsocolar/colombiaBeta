dist <- read.csv("/Users/jacobsocolar/Downloads/DB_distri-20-07-21.csv", sep = ";")
points <- readRDS("/Users/jacobsocolar/Downloads/all_pts.rds")

db <- read.csv("/Users/jacobsocolar/Downloads/all_datasets.csv", sep = ";")

`%ni%` <- Negate(`%in%`)

sum(db$point %ni% points$point)
unique(db$point[db$point %ni% points$point])

sum(points$beetles == 1)
unique(points$point[points$point %ni% db$point])

source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/GIS_processing/hydrosheds_extraction.R")

points_sf <- st_as_sf(points, coords = c("lon", "lat"), crs = 4326)
points_sf$dist_slope <- points_sf$dist_region <- NA

points_sf$dist_region[st_intersects(points_sf, snsm, sparse = F)] <- "SNS Marta"
points_sf$dist_region[!st_intersects(points_sf, snsm, sparse = F) &
                        ((!st_intersects(points_sf, amazon_orinoco, sparse = F)) |
                           (points_sf$elev_ALOS > 500))] <- "Andean"
points_sf$dist_region[grepl("leguizamo", points_sf$subregion)] <- "Amazon"
points_sf$dist_region[grepl("chiribiquete", points_sf$subregion)] <- "Amazon"
points_sf$dist_region[grepl("guaviare", points_sf$subregion)] <- "Amazon"
points_sf$dist_region[grepl("llanos", points_sf$subregion)] <- "Llanos"



points_sf$dist_slope[st_intersects(points_sf, amazon_orinoco, sparse = F) |
                       st_intersects(points_sf, catatumbo, sparse = F)] <- "EC: Eastern"
points_sf$dist_slope[st_intersects(points_sf, magdalena_east, sparse = F)] <- "EC: Western"
points_sf$dist_slope[st_intersects(points_sf, magdalena_west, sparse = F)] <- "CC: Eastern"
points_sf$dist_slope[st_intersects(points_sf, cauca_east, sparse = F)] <- "CC: Western"
points_sf$dist_slope[st_intersects(points_sf, cauca_west, sparse = F)] <- "WC: Eastern"
points_sf$dist_slope[st_intersects(points_sf, pacific, sparse = F)] <- "WC: Western"
points_sf$dist_slope[grep("PUP", points_sf$point)] <- "CC: Western slope"

db1 <- db[!is.na(db$abundance), ]
db1 <- db1[db$abundance > 0, ]

counter1 <- 0
counter2 <- 0
region_mismatch <- slope_mismatch <- vector()
for(i in 1:nrow(db1)){
  point <- db1$point[i]
  species <- gsub(" ", "_", db1$scientificName[i])
  if(point %ni% c("BEP1 ", paste0("BF", 1:6), paste0("BP", 1:6), NA, "sep-01", "sep-02", "sep-03",
                  "TAP1",
                  "FOREST_NEAR1", "FOREST_NEAR1A","FOREST_NEAR2", "FOREST_NEAR2A", "FOREST_NEAR3") &
     species %ni% c("", "Uroxys_sp._3", "Coprophanaeus_gr_pluto", "Canthon_colombianus",
                    "Ateuchus_sp._01H", "Canthon_fulgidus", "Ateuchus_sp._09H", "Scybalocanthon_arnaudi_",
                    "Dichotomius_boreus", "Ateuchus_sp._2", "Ateuchus_sp._08H", "Ateuchus_sp._10H",
                    "Gromphas_lemoinei_", "Onthophagus_sp._2J", "Onthophagus_landolti_", "Oxysternon_silenus_",
                    "non_spp_record", "Anisocanthon_villosus_", "Sylvicanthon_sp._03H", "Ateuchus_sp._05H",
                    "Onthophagus_curivicornis", "Onthophagus_sp.08H", "Canthidium_sp._4FE")){
    if(point %ni% points_sf$point){stop("point issue")}
    if(species %ni% dist$scientificName){stop("name issue")}
    distRegion <- dist$Region[dist$scientificName == species]
    pointsRegion <- points_sf$dist_region[points$point == point]
    if(!grepl(pointsRegion, distRegion)){
      counter1 <- counter1+1
      region_mismatch[counter1] <- i
    }
    distSlope <- dist$Slope[dist$scientificName == species]
    pointsSlope <- points_sf$dist_slope[points$point == point]
    
    if(grepl(pointsRegion, "Andean")){
      pointsSlope <- gsub(" ", "", pointsSlope)
      distSlope <- gsub(" ", "", distSlope)
      pointsSlope <- unlist(strsplit(pointsSlope, ":"))
      distSlope <- unlist(strsplit(distSlope, "(:|;)"))
      
      if (pointsSlope[1] %ni% distSlope){
        counter2 <- counter2+1
        slope_mismatch[counter2] <- i
      } else {
        n <- which(distSlope == pointsSlope[1])
        l <- length(distSlope)
        if(l < (n+2)) {
          if(!grepl(pointsSlope[2], distSlope[l])) {
            counter2 <- counter2+1
            slope_mismatch[counter2] <- i
          }
        } else if (!grepl(pointsSlope[2], distSlope[(n+1):(n+2)])) {
          counter2 <- counter2+1
          slope_mismatch[counter2] <- i
        }
      }
    }
  }
}

region_mismatch_df <- data.frame(species = db1$scientificName[region_mismatch], 
                                 point = db1$point[region_mismatch],
                                 species_dist = NA, point_dist = NA)
for(i in 1:nrow(region_mismatch_df)){
  region_mismatch_df$species_dist[i] <- dist$Region[dist$scientificName == gsub(" ", "_", region_mismatch_df$species[i])]
  region_mismatch_df$point_dist[i] <- points_sf$dist_region[points_sf$point == region_mismatch_df$point[i]]
}


slope_mismatch_df <- data.frame(species = db1$scientificName[slope_mismatch], 
                                 point = db1$point[slope_mismatch],
                                 species_dist = NA, point_dist = NA)
for(i in 1:nrow(slope_mismatch_df)){
  slope_mismatch_df$species_dist[i] <- dist$Slope[dist$scientificName == gsub(" ", "_", slope_mismatch_df$species[i])]
  slope_mismatch_df$point_dist[i] <- points_sf$dist_slope[points_sf$point == slope_mismatch_df$point[i]]
}
slope_mismatch_df <- slope_mismatch_df[order(slope_mismatch_df$point_dist),]

mismatch_species <- unique(slope_mismatch_df$species)

slope_mismatch_df[slope_mismatch_df$species == mismatch_species[41], ]
