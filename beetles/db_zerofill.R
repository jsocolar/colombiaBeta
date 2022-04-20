socolar.desktop <- file.exists('/Users/jacobsocolar/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
socolar.laptop <- file.exists('/Users/jacob/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
if(socolar.desktop){
  dir.path <- "/Users/JacobSocolar/Dropbox"
}else if(socolar.laptop){
  dir.path <- "/Users/jacob/Dropbox"
}

##### Ingest and clean beetle data #####
db <- read.csv(paste0(dir.path, "/CO_DBdata/DB_all_datasets.csv"), sep = ";")
db2 <- db[, c("point", "day", "scientificName", "abundance")]

# Read point data
j_points <- readRDS(paste0(dir.path, "/Work/Colombia/Data/GIS/Points/all_pts.RDS"))

### Identify and fix points in db2 that don't exist in j_points
unique(db2$point[!(db2$point %in% j_points$point)])

# Naming errors
db2$point[db2$point == "CHA10d"] <- "CHA10D"
db2$point[db2$point == "CHA11d"] <- "CHA11D"
db2$point[db2$point == "CHA12d"] <- "CHA12D"
db2$point[db2$point == "SEP_1"] <- "SEP1"
db2$point[db2$point == "SEP_2"] <- "SEP2"
db2$point[db2$point == "SEP_3"] <- "SEP3"
db2$point[db2$point == "BEP1 "] <- "BEP1"

# Superfluous points: TAP1 (located too close to forest) and CHA7-9 (not properly collected)
db2 <- db2[!(db2$point %in% c("TAP1", "CHA7", "CHA8", "CHA9")), ]

unique(db2$point[!(db2$point %in% j_points$point)])

### Inspect and fix points in j_point that don't exist in db2
missing_pts <- unique(j_points$point[!(j_points$point %in% db2$point) & (j_points$habitat != "PALM") & (j_points$habitat != "Sy")])
missing_pts

# Confirm that nothing unexpected is happening here
missing_pts2 <- j_points$point[!(j_points$point %in% db2$point) & (j_points$natural | j_points$pasture)]
missing_pts[!(missing_pts %in% missing_pts2)]

#View(j_points[j_points$point %in% missing_pts, ])
# Most of these are llanos or wandes points for which we lack beetle collections. 
# Some are MOF points whose collections were improperly and unrecoverably cataloged
# The remainder are PLP13-24, SAP10, and BEA1-21.
pt_additions <- c(paste0("PLP", 13:24), "SAP10", paste0("BEA", 1:21))
extra_frame <- data.frame(point = rep(pt_additions, each = 4),
                          day = paste0("R", 1:4),
                          scientificName = "",
                          abundance = NA)
db2 <- rbind(db2, extra_frame)

### Remove points with insufficient collection day info
no_day_points <- unique(db2$point[db2$day == ""])
no_day_points
db2 <- db2[!(db2$point %in% no_day_points), ]

### Get all point-days and species-point-days that exist in the data
point_day <- paste(db2$point, db2$day, sep = "__")
species_point_day <- paste(db2$scientificName, point_day, sep = "___")

# Remove duplicated rows by summing.  In all cases, either (1) at least one
# of the duplicate rows had abundance 0 or (2) the duplicate rows are known
# cases confirmed by Diego of beetles recorded separately when pinned vs 
# preserved in alcohol
unique(db2$scientificName[is.na(db2$abundance)]) # all na abundances correspond to missing species
db2_duplicates <- duplicated(species_point_day)
for(i in rev(which(db2_duplicates))) {
  i1 <- min(which(species_point_day == species_point_day[i]))
  db2$abundance[i1] <- db2$abundance[i1] + db2$abundance[i]
}
db2 <- db2[!db2_duplicates, ]

# remove rows without a species (Important that we already generated point_day)
db2 <- db2[db2$scientificName != "", ]

##### Zero-fill #####
# Get all species; all point-days
species <- unique(db2$scientificName)
point_days <- unique(point_day)

# Get all species-point-days that exist in db2
species_point_day_db2 <- paste(db2$scientificName, db2$point, db2$day, sep = "__")

# Assemble data-frame with all possible species-point-days
all_spd <- data.frame(point = NA, day = NA, 
                      scientificName = rep(species, length(point_days)), 
                      abundance = 0, 
                      point_day = rep(point_days, each = length(species)))
species_point_day_all <- paste(all_spd$scientificName, all_spd$point_day, sep = "__")
first_cols <- do.call(rbind, strsplit(all_spd$point_day, "__"))
all_spd$point <- first_cols[,1]
all_spd$day <- first_cols[,2]

all_spd <- all_spd[, c("point", "day", "scientificName", "abundance")]

head(all_spd)
head(db2)

db2_additions <- all_spd[!(species_point_day_all %in% species_point_day_db2), ]

# Confirm that the dimension makes sense
nrow(all_spd) == nrow(db2_additions) + nrow(db2)

# rbind
db3 <- rbind(db2, db2_additions)

##### Add point data #####
j_points <- j_points[,c("point", "lat", "lon", "site", "cluster", 
                        "habitat", "natural", "paramo", "pasture", 
                        "other", "mixed_cluster", "elev_ALOS", "subregion")]

source(paste0(dir.path, "/Work/Code/colombiaBeta/GIS_processing/hydrosheds_extraction.R"))

points_sf <- st_as_sf(j_points, coords = c("lon", "lat"), crs = 4326)
j_points$pt_slope <- j_points$pt_region <- NA

j_points$pt_region[st_intersects(points_sf, snsm, sparse = F)] <- "SNS Marta"
j_points$pt_region[!st_intersects(points_sf, snsm, sparse = F) &
                        ((!st_intersects(points_sf, amazon_orinoco, sparse = F)) |
                           (points_sf$elev_ALOS > 500))] <- "Andean"
j_points$pt_region[grepl("leguizamo", points_sf$subregion)] <- "Amazon"
j_points$pt_region[grepl("chiribiquete", points_sf$subregion)] <- "Amazon"
j_points$pt_region[grepl("guaviare", points_sf$subregion)] <- "Amazon"
j_points$pt_region[grepl("llanos", points_sf$subregion)] <- "Llanos"

j_points$pt_slope[st_intersects(points_sf, amazon_orinoco, sparse = F) |
                       st_intersects(points_sf, catatumbo, sparse = F)] <- "EC: Eastern"
j_points$pt_slope[st_intersects(points_sf, magdalena_east, sparse = F)] <- "EC: Western"
j_points$pt_slope[st_intersects(points_sf, magdalena_west, sparse = F)] <- "CC: Eastern"
j_points$pt_slope[st_intersects(points_sf, cauca_east, sparse = F)] <- "CC: Western"
j_points$pt_slope[st_intersects(points_sf, cauca_west, sparse = F)] <- "WC: Eastern"
j_points$pt_slope[st_intersects(points_sf, pacific, sparse = F)] <- "WC: Western"
j_points$pt_slope[j_points$pt_region == "SNS Marta"] <- "SNSM"

# Deal with minor details of points very near to the triple-point in Purace.
# Treat points on east side of the Turbio/Patia valley as west slope of central cordillera. 
# Also treat points near valencia as west slope of central
j_points$pt_slope[points_sf$point %in% paste0("PUF", 4:12)] <- "EC: Western"
j_points$pt_slope[points_sf$point %in% paste0("PUP", 1:12)] <- "CC: Western"

db4 <- merge(db3, j_points, by.x = "point", by.y = "point", all = FALSE)

##### Format biogeography data #####
biogeography <- read.csv(paste0(dir.path, "/CO_DBdata/DB_Distributions_traits.csv"), sep = ";")
biogeography$Slope <- gsub(" ", "", biogeography$Slope)

slope_split <- strsplit(biogeography$Slope, split = ";")
unique(unlist(slope_split))

biogeography2 <- data.frame(scientificName = biogeography$scientificName, 
                            sp_elev_lower = biogeography$Lower, sp_elev_upper = biogeography$Upper,
                            sp_region_amazon = NA, sp_region_llanos = NA,
                            sp_region_caribbean = NA, sp_region_snsm = NA,
                            sp_region_andes = NA,
                            sp_slope_ECe = NA, sp_slope_ECw = NA,
                            sp_slope_CCe = NA, sp_slope_CCw = NA,
                            sp_slope_WCe = NA, sp_slope_WCw = NA,
                            sp_slope_SNSM = NA)


as.integer2 <- function (x) {
  if (length(x) == 0) {
    return(0)
  } else {
    return(max(as.integer(x)))
  }
}

for (i in 1:nrow(biogeography)) {
  biogeography2$sp_region_eastern[i] <- as.integer2(grepl("Eastern lowlands", biogeography$Region[i]))
  biogeography2$sp_region_caribbean[i] <- as.integer2(grepl("Caribbean", biogeography$Region[i]))
  biogeography2$sp_region_snsm[i] <- as.integer2(grepl("SNS Marta", biogeography$Region[i]))
  biogeography2$sp_region_andes[i] <- as.integer2(grepl("Andean", biogeography$Region[i]))
  
  slope_data <- slope_split[[i]]
  
  ec <- slope_data[grep("EC:", slope_data)]
  biogeography2$sp_slope_ECe[i] <- as.integer2(grepl("Eastern", ec))
  biogeography2$sp_slope_ECw[i] <- as.integer2(grepl("Western", ec))
  
  cc <- slope_data[grep("CC:", slope_data)]
  biogeography2$sp_slope_CCe[i] <- as.integer2(grepl("Eastern", cc))
  biogeography2$sp_slope_CCw[i] <- as.integer2(grepl("Western", cc))
  
  wc <- slope_data[grep("WC:", slope_data)]
  biogeography2$sp_slope_WCe[i] <- as.integer2(grepl("Eastern", wc))
  biogeography2$sp_slope_WCw[i] <- as.integer2(grepl("Western", wc))
  
  if(biogeography2$sp_region_eastern[i] == 1) {
    biogeography2$sp_slope_ECe[i] <- 1
  }
  
  # lowland species always cross valleys
  if ((biogeography2$sp_slope_ECw[i] != biogeography2$sp_slope_CCe[i]) & (biogeography2$sp_elev_lower[i] < 300)) {
    biogeography2$sp_slope_ECw[i] <- biogeography2$sp_slope_CCe[i] <- 1
  }

  if ((biogeography2$sp_slope_CCw[i] != biogeography2$sp_slope_WCe[i]) & (biogeography2$sp_elev_lower[i] < 300)) {
    biogeography2$sp_slope_CCw[i] <- biogeography2$sp_slope_WCe[i] <- 1
  }
  # 
  # # Caribbean species all potentially enter the Magdalena and reach the SNSM (max lower for these species is 102)
  # if (biogeography2$sp_region_caribbean[i] == 1) {
  #   biogeography2$sp_slope_ECw[i] <- biogeography2$sp_slope_CCe[i] <- biogeography2$sp_slope_SNSM[i] <- 1
  # }
  # 
  # # no gaps
  # ci <- which(colnames(biogeography2) == "sp_slope_ECe")
  # cf <- which(colnames(biogeography2) == "sp_slope_WCw")
  # the_numbers <- biogeography2[i, ci:cf]
  # if(sum(the_numbers) > 0){
  #   first_n <- min(which(the_numbers == 1))
  #   last_n <- max(which(the_numbers == 1))
  #   the_numbers[first_n:last_n] <- 1
  #   biogeography2[i, ci:cf] <- the_numbers
  # }
}
biogeography2$sp_slope_SNSM <- biogeography2$sp_region_snsm

# update problems: EC E
biogeography2$sp_slope_ECe[which(biogeography2$scientificName == "Deltochilum_sp._18H")] <- 1
biogeography2$sp_slope_ECe[which(biogeography2$scientificName == "Deltochilum_sp._38H")] <- 1

# update problems: EC W
biogeography2$sp_slope_ECw[which(biogeography2$scientificName == "Uroxys_sp._06H")] <- 1

# update problems: WC W
biogeography2$sp_slope_WCw[which(biogeography2$scientificName == "Canthidium_sp._07H")] <- 1
biogeography2$sp_slope_WCw[which(biogeography2$scientificName == "Canthidium_sp._19H")] <- 1

db5 <- merge(db4, biogeography2, by = "scientificName", all = TRUE)

db5$sp_elev_lower2 <- pmin(db5$sp_elev_lower, (db5$sp_elev_lower + db5$sp_elev_upper)/2 - 250)
db5$sp_elev_upper2 <- pmax(db5$sp_elev_upper, (db5$sp_elev_lower + db5$sp_elev_upper)/2 + 250)

##### Check biogeography #####
db_obs <- db5[db5$abundance > 0, ]
# View(db_obs[db_obs$pt_region == "Andean" & db_obs$sp_region_andes == 0, ])
# View(db_obs[db_obs$pt_region == "Amazon" & db_obs$sp_region_amazon == 0, ])
# View(db_obs[db_obs$pt_region == "Llanos" & db_obs$sp_region_llanos == 0, ])
# View(db_obs[db_obs$pt_region == "SNS Marta" & db_obs$sp_region_snsm == 0, ])


sum(db_obs$pt_slope == "SNSM" & db_obs$sp_slope_SNSM == 0)
View(db_obs[db_obs$pt_slope == "EC: Eastern" & db_obs$sp_slope_ECe == 0,])
View(db_obs[db_obs$pt_slope == "EC: Western" & db_obs$sp_slope_ECw == 0,])
sum(db_obs$pt_slope == "CC: Eastern" & db_obs$sp_slope_CCe == 0)
sum(db_obs$pt_slope == "CC: Western" & db_obs$sp_slope_CCw == 0)
sum(db_obs$pt_slope == "WC: Eastern" & db_obs$sp_slope_WCe == 0)
sum(db_obs$pt_slope == "WC: Western" & db_obs$sp_slope_WCw == 0)

View(db_obs[db_obs$pt_slope == "WC: Western" & db_obs$sp_slope_WCw == 0,])


View(db_obs[db_obs$elev_ALOS > db_obs$sp_elev_upper2 + (db_obs$sp_elev_upper2 - db_obs$sp_elev_lower2),])
View(db_obs[db_obs$elev_ALOS < db_obs$sp_elev_lower2 - (db_obs$sp_elev_upper2 - db_obs$sp_elev_lower2),])


##### Biogeographic clipping #####
db5$elev_standard <- 2*(db5$elev_ALOS - db5$sp_elev_lower2)/(db5$sp_elev_upper2 - db5$sp_elev_lower2) - 1
hist(db5$elev_standard[db5$abundance > 0])

View(db5[db5$elev_standard > 3 & db5$abundance > 0, ])

db_ocana <- db5[db5$subregion == "subregion_Ocana", ]
db_other <- db5[db5$subregion != "subregion_Ocana", ]


db_ocana_clipped <- db_ocana[(db_ocana$sp_slope_ECe + db_ocana$sp_slope_ECw) > 0, ]
db_other_clipped <- db_other[((db_other$pt_slope == "EC: Eastern") & (db_other$sp_slope_ECe == 1)) |
                               ((db_other$pt_slope == "EC: Western") & (db_other$sp_slope_ECw == 1)) |
                               ((db_other$pt_slope == "CC: Eastern") & (db_other$sp_slope_CCe == 1)) |
                               ((db_other$pt_slope == "CC: Western") & (db_other$sp_slope_CCw == 1)) |
                               ((db_other$pt_slope == "WC: Eastern") & (db_other$sp_slope_WCe == 1)) |
                               ((db_other$pt_slope == "WC: Western") & (db_other$sp_slope_WCw == 1)) |
                               ((db_other$pt_slope == "SNSM") & (db_other$sp_slope_SNSM == 1)), 
                             ]


db_clipped <- rbind(db_ocana_clipped, db_other_clipped)
saveRDS(db_clipped, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Beetles/db_clipped.RDS")


sum(db_clipped$elev_standard < -3)



db_clipped2 <- db_clipped[db_clipped$elev_standard < 3 & db_clipped$elev_standard > -3, ]

##### brms model #####
library(brms)
db_clipped$elev_standard_squared <- db_clipped$elev_standard^2
db_clipped$subregion_species <- paste0(db_clipped$subregion, "__", db_clipped$scientificName)
db_clipped$cluster_species <- paste0(db_clipped$cluster, "__", db_clipped$scientificName)

db_mod1 <- 
  brm(abundance ~ pasture + elev_standard + elev_standard_squared + day + 
      (1 + pasture + elev_standard + elev_standard_squared + day | scientificName) +
      (1 | subregion_species) + (1 | cluster_species), 
    family = "negbinomial", data = db_clipped, 
    chains = 3, cores = 3, backend = 'cmdstanr',
    refresh = 10)
