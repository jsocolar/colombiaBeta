library(data.table)
library(brms)

setwd("/Users/jacobsocolar/Dropbox/Work/Colombia/Data")

##### Read in data #####
# Sampling data
db_data <- read.csv("Beetles/DB_final/all_datasets_22-06-21.csv", sep = ";")

# View(db_data)
unique(db_data$habitat)
unique(db_data$habitat_all_points)
db_data$forest <- db_data$habitat_all_points %in% c("Forest", "FOREST", "Paramo", "P", "Sm")
db_data$pasture <- db_data$habitat_all_points %in% c("Pasture", "PASTURE", "G")

db_data$include <- (db_data$forest + db_data$pasture) > 0

db_data <- db_data[db_data$include, ]

unique(db_data$habitat[db_data$habitat_all_points == "P"])
sum(db_data$habitat[db_data$habitat_all_points == "Forest"] == "Pasture")

# Trait data
db_traits <- read.csv("Beetles/DB_final/traits-distri-22-06-21.csv", sep = ";")

sum(is.na(db_data$abundance))
unique(db_data$abundance)
unique(db_data$day)
db_data[db_data$day == "",]


unique(db_data$taxonRank)


leguizamo_data <- db_data[db_data$county == "Legizamo",]

a <- unique(paste(soata_data$point, soata_data$day))
a[order(a)]
a


db_data <- db_data[db_data$scientificName != "", ]

species <- unique(db_data$scientificName)
points <- unique(db_data$point)

zero_filled_data <- as.data.frame(matrix(NA, nrow = length(species) * length(points), ncol = ncol(db_dt)))
names(zero_filled_data) <- names(db_dt)
zero_filled_data$scientificName <- rep(species, length(points))
zero_filled_data$point <- rep(points, times = length(species))

zero_filled_data <- rbind(db_data, zero_filled_data)
db_dt <- as.data.table(zero_filled_data)


flattened_data <- db_dt[,.(abun = sum(abundance, na.rm = T)), by = .(scientificName, point)]


db_data_final1 <- merge(flattened_data, db_traits, by = "scientificName", all = T)

points <- read.csv("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Points/CO_sampling_points_metafile.csv")

db_data_final <- merge(db_data_final, points, by.x = "point", by.y = "point_id")

forest_point_data <- db_data[!duplicated(db_data$point), c("point", "forest")]

db_data_final <- merge(db_data_final, forest_point_data, by = "point")


detections <- db_data_final[db_data_final$abun > 0, ]
nrow(detections)

rel_elev <- (detections$elev_ALOS30m - (detections$Lower))/((detections$Upper) - (detections$Lower))

rel_elev <- (detections$elev_ALOS30m - (detections$Lower - 200))/((detections$Upper + 200) - (detections$Lower - 200))

hist(rel_elev)


nrow(detections[rel_elev < 1.5, ])
View(detections[rel_elev < -1, c("point", "scientificName", "Lower", "Upper", "elev_ALOS30m")])

db_data_final$rel_elev <- (db_data_final$elev_ALOS30m - (db_data_final$Lower - 200))/((db_data_final$Upper + 200) - (db_data_final$Lower - 200)) - 0.5
db_data_final$rel_elev2 <- db_data_final$rel_elev^2

db_data_final$nest_guild[db_data_final$nest_guild == "Endocoprid"] <- "dweller"

brm(abun ~ rel_elev + rel_elev2 + (1 + rel_elev + rel_elev2 | scientificName), 
    data = db_data_final, family = zero_inflated_negbinomial())

a <- lme4::glmer(abun ~ rel_elev + rel_elev2 + forest*nest_guild + (1|scientificName), data = db_data_final, family = poisson)

summary(a)
