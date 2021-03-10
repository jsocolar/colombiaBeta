source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/bird_analysis_plotting/biodiversity_loss/compute_loss.R")

forest_probs <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v5_predictions/iteration_1/forest_probs.RDS")
pasture_probs <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v5_predictions/iteration_1/pasture_probs.RDS")


cell_logratios <- get_avg_cell_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05)
colombia_logratio <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = NULL)

mean(cell_logratios$avg_logratio, na.rm = T)
weighted.mean(cell_logratios$avg_logratio, w=cell_logratios$n, na.rm=T)
spatstat::weighted.median(cell_logratios$med_logratio, w=cell_logratios$n, na.rm = T)
colombia_logratio

cell_logratios_df <- data.frame(x=forest_probs$x, y = forest_probs$y, logratio=cell_logratios$avg_logratio)
cell_logratios_raster <- raster::rasterFromXYZ(cell_logratios_df)
raster::plot(cell_logratios_raster)

cell_richness_df <- data.frame(x=forest_probs$x, y=forest_probs$y, richness = cell_logratios$n)
cell_richness_raster <- raster::rasterFromXYZ(cell_richness_df)
raster::plot(cell_richness_raster)

cell_logratios_raster_51_100 <- cell_logratios_raster_101_150 <- cell_logratios_raster_151_200 <-
  cell_logratios_raster_201_250 <- cell_logratios_raster_251_300 <- cell_logratios_raster_301_350 <-
  cell_logratios_raster_351_400 <- cell_logratios_raster_401_max <- cell_logratios_raster
cell_logratios_raster_51_100[cell_richness_raster <= 50 | cell_richness_raster > 100] <- NA
cell_logratios_raster_101_150[cell_richness_raster <= 100 | cell_richness_raster > 150] <- NA
cell_logratios_raster_151_200[cell_richness_raster <= 150 | cell_richness_raster > 200] <- NA
cell_logratios_raster_201_250[cell_richness_raster <= 200 | cell_richness_raster > 250] <- NA
cell_logratios_raster_251_300[cell_richness_raster <= 250 | cell_richness_raster > 300] <- NA
cell_logratios_raster_301_350[cell_richness_raster <= 300 | cell_richness_raster > 350] <- NA
cell_logratios_raster_351_400[cell_richness_raster <= 350 | cell_richness_raster > 400] <- NA
cell_logratios_raster_401_max[cell_richness_raster <= 400] <- NA

plot(cell_logratios_raster_51_100)
plot(cell_logratios_raster_101_150)
plot(cell_logratios_raster_151_200)
plot(cell_logratios_raster_201_250)
plot(cell_logratios_raster_251_300)
plot(cell_logratios_raster_301_350)
plot(cell_logratios_raster_351_400)
plot(cell_logratios_raster_401_max)

colombia_logratio_51_100 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 50 & cell_logratios$n < 101))
colombia_logratio_101_150 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 100 & cell_logratios$n < 151))
colombia_logratio_151_200 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 150 & cell_logratios$n < 201))
colombia_logratio_201_250 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 200 & cell_logratios$n < 251))
colombia_logratio_251_300 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 250 & cell_logratios$n < 301))
colombia_logratio_301_350 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 300 & cell_logratios$n < 351))
colombia_logratio_351_400 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 350 & cell_logratios$n < 401))
colombia_logratio_401_max <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 400))

mean(cell_logratios$avg_logratio[cell_logratios$n > 50 & cell_logratios$n < 101])
colombia_logratio_51_100

mean(cell_logratios$avg_logratio[cell_logratios$n > 100 & cell_logratios$n < 151])
colombia_logratio_101_150

mean(cell_logratios$avg_logratio[cell_logratios$n > 150 & cell_logratios$n < 201])
colombia_logratio_151_200

mean(cell_logratios$avg_logratio[cell_logratios$n > 200 & cell_logratios$n < 251])
colombia_logratio_201_250

mean(cell_logratios$avg_logratio[cell_logratios$n > 250 & cell_logratios$n < 301])
colombia_logratio_251_300

mean(cell_logratios$avg_logratio[cell_logratios$n > 300 & cell_logratios$n < 351])
colombia_logratio_301_350

mean(cell_logratios$avg_logratio[cell_logratios$n > 350  & cell_logratios$n < 401])
colombia_logratio_351_400

mean(cell_logratios$avg_logratio[cell_logratios$n > 400])
colombia_logratio_401_max

#######
library(dggridR)
library(sf)
library(raster)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
points <- st_as_sf(forest_probs[,1:3], coords = c("x", "y"), crs = AEAstring)
points_latlong <- st_transform(points, 4326)
points_coords <- as.data.frame(st_coordinates(points_latlong))

dggs_4 <- dgconstruct(res=4)
dggs_7 <- dgconstruct(res=7)
dggs_10 <- dgconstruct(res=10)
dggs_13 <- dgconstruct(res=13)

cells_4 <- dgGEO_to_SEQNUM(dggs_4, points_coords$X, points_coords$Y)$seqnum
cells_7 <- dgGEO_to_SEQNUM(dggs_7, points_coords$X, points_coords$Y)$seqnum
cells_10 <- dgGEO_to_SEQNUM(dggs_10, points_coords$X, points_coords$Y)$seqnum
cells_13 <- dgGEO_to_SEQNUM(dggs_13, points_coords$X, points_coords$Y)$seqnum

uc_4 <- unique(cells_4)
uc_7 <- unique(cells_7)
uc_10 <- unique(cells_10)
uc_13 <- unique(cells_13)

cell_logratios_4 <- cell_numbers_4 <- vector()
cell_logratios_4_expanded <- rep(NA, nrow(points))
for(i in 1:length(uc_4)){
  print(i)
  cell_logratios_4[i] <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cells_4 == uc_4[i]))$avg_logratio
  cell_logratios_4_expanded[which(cells_4 == uc_4[i])] <- cell_logratios_4[i]
  cell_numbers_4[i] <- sum(cells_4 == uc_4[i])
}
mean(cell_logratios_4)
mean(cell_logratios_4[c(3,4,6)])
raster_4 <- rasterFromXYZ(cbind(forest_probs[,2:3], cell_logratios_4_expanded))
plot(raster_4)


cell_logratios_7 <- cell_numbers_7 <- vector()
cell_logratios_7_expanded <- rep(NA, nrow(points))
for(i in 1:length(uc_7)){
  print(i)
  cell_logratios_7[i] <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cells_7 == uc_7[i]))$avg_logratio
  cell_logratios_7_expanded[which(cells_7 == uc_7[i])] <- cell_logratios_7[i]
  cell_numbers_7[i] <- sum(cells_7 == uc_7[i])
}
mean(cell_logratios_7, na.rm=T)
cell_numbers_7
raster_7 <- rasterFromXYZ(cbind(forest_probs[,2:3], cell_logratios_7_expanded))
plot(raster_7)

cell_logratios_10 <- cell_numbers_10 <- vector()
cell_logratios_10_expanded <- rep(NA, nrow(points))
for(i in 1:length(uc_10)){
  print(i)
  cell_logratios_10[i] <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cells_10 == uc_10[i]))$avg_logratio
  cell_logratios_10_expanded[which(cells_10 == uc_10[i])] <- cell_logratios_10[i]
  cell_numbers_10[i] <- sum(cells_10 == uc_10[i])
}
mean(cell_logratios_10, na.rm=T)
cell_numbers_10
raster_10 <- rasterFromXYZ(cbind(forest_probs[,2:3], cell_logratios_10_expanded))
plot(raster_10)


