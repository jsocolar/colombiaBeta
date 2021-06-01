# This script analyzes patterns of biodiversity loss across different collections of points
library(raster)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
llanos_amazon <- sf::st_transform(sf::st_read("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/masking_polygons/E_Amazon_Llanos.kml"),
                                  AEAstring)
# Functions to compute biodiveristy loss over arbitrary collections of points
source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/bird_analysis_plotting/biodiversity_loss/compute_loss.R")

# Point-scale occupancy probabilites for each species in forest and pasture
forest_probs_integrated <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v6_predictions/iteration_4000/forest_probs_integrated.RDS")
pasture_probs_integrated <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v6_predictions/iteration_4000/pasture_probs_integrated.RDS")

forest_probs_semi <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v6_predictions/iteration_4000/forest_probs_semi.RDS")
pasture_probs_semi <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v6_predictions/iteration_4000/pasture_probs_semi.RDS")

forest_probs_rep <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v6_predictions/iteration_4000/forest_probs_rep.RDS")
pasture_probs_rep <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/v6_predictions/iteration_4000/pasture_probs_rep.RDS")

forest_probs <- forest_probs_rep[, 1:1618]
pasture_probs <- pasture_probs_rep[, 1:1618]

forest_probs_points <- sf::st_as_sf(forest_probs [,1:3], coords = c("x", "y"), crs = AEAstring)
use_points <- !sf::st_intersects(forest_probs_points, llanos_amazon, sparse = F)

forest_probs <- forest_probs[use_points,]
pasture_probs <- pasture_probs[use_points,]

# Get pointwise average logratios (average across species) for all points
cell_logratios <- get_avg_cell_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.5)
# Get regional (Colombia-wide) logratio (average across species)
colombia_logratio <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.5, cell_positions = NULL)

##### Colombia-wide comparison #####
mean(cell_logratios$med_logratio, na.rm = T)
colombia_logratio

# Map the local averages
cell_logratios_df <- data.frame(x=forest_probs$x, y = forest_probs$y, logratio=cell_logratios$avg_logratio)
cell_logratios_raster <- rasterFromXYZ(cell_logratios_df)
plot(cell_logratios_raster)
# Map species richness
cell_richness_df <- data.frame(x=forest_probs$x, y=forest_probs$y, richness = cell_logratios$n)
cell_richness_raster <- raster::rasterFromXYZ(cell_richness_df)
raster::plot(cell_richness_raster)


# ##### Comparison within bands of constant richness #####
# # Relationship between local richness and local average
# plot(cell_logratios_df$logratio ~ cell_richness_df$richness, pch = ".")
# 
# # rasters for different richness bands
# cell_logratios_raster_51_100 <- cell_logratios_raster_101_150 <- cell_logratios_raster_151_200 <-
#   cell_logratios_raster_201_250 <- cell_logratios_raster_251_300 <- cell_logratios_raster_301_350 <-
#   cell_logratios_raster_351_400 <- cell_logratios_raster_401_max <- cell_logratios_raster_1_250 <- cell_logratios_raster
# cell_logratios_raster_51_100[cell_richness_raster <= 50 | cell_richness_raster > 100] <- NA
# cell_logratios_raster_101_150[cell_richness_raster <= 100 | cell_richness_raster > 150] <- NA
# cell_logratios_raster_151_200[cell_richness_raster <= 150 | cell_richness_raster > 200] <- NA
# cell_logratios_raster_201_250[cell_richness_raster <= 200 | cell_richness_raster > 250] <- NA
# cell_logratios_raster_251_300[cell_richness_raster <= 250 | cell_richness_raster > 300] <- NA
# cell_logratios_raster_301_350[cell_richness_raster <= 300 | cell_richness_raster > 350] <- NA
# cell_logratios_raster_351_400[cell_richness_raster <= 350 | cell_richness_raster > 400] <- NA
# cell_logratios_raster_401_max[cell_richness_raster <= 400] <- NA
# 
# cell_logratios_raster_1_250[cell_richness_raster >= 250] <- NA
# 
# 
# # Plotting
# plot(cell_logratios_raster_51_100)
# plot(cell_logratios_raster_101_150)
# plot(cell_logratios_raster_151_200)
# plot(cell_logratios_raster_201_250)
# plot(cell_logratios_raster_251_300)
# plot(cell_logratios_raster_301_350)
# plot(cell_logratios_raster_351_400)
# plot(cell_logratios_raster_401_max)
# plot(cell_logratios_raster_1_250)
# 
# 
# # Get the regional logratios within each richness band
# colombia_logratio_51_100 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 50 & cell_logratios$n < 101))
# colombia_logratio_101_150 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 100 & cell_logratios$n < 151))
# colombia_logratio_151_200 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 150 & cell_logratios$n < 201))
# colombia_logratio_201_250 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 200 & cell_logratios$n < 251))
# colombia_logratio_251_300 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 250 & cell_logratios$n < 301))
# colombia_logratio_301_350 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 300 & cell_logratios$n < 351))
# colombia_logratio_351_400 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 350 & cell_logratios$n < 401))
# colombia_logratio_401_max <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 400))
# 
# colombia_logratio_1_250 <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.05, cell_positions = which(cell_logratios$n > 0 & cell_logratios$n < 251))
# 
# 
# # Compare pointwise averages to regional average within bands
# median(cell_logratios$med_logratio[cell_logratios$n > 50 & cell_logratios$n < 101])
# colombia_logratio_51_100
# 
# median(cell_logratios$med_logratio[cell_logratios$n > 100 & cell_logratios$n < 151])
# colombia_logratio_101_150
# 
# median(cell_logratios$med_logratio[cell_logratios$n > 150 & cell_logratios$n < 201])
# colombia_logratio_151_200
# 
# median(cell_logratios$med_logratio[cell_logratios$n > 200 & cell_logratios$n < 251])
# colombia_logratio_201_250
# 
# median(cell_logratios$med_logratio[cell_logratios$n > 250 & cell_logratios$n < 301])
# colombia_logratio_251_300
# 
# median(cell_logratios$med_logratio[cell_logratios$n > 300 & cell_logratios$n < 351])
# colombia_logratio_301_350
# 
# median(cell_logratios$med_logratio[cell_logratios$n > 350  & cell_logratios$n < 401])
# colombia_logratio_351_400
# 
# median(cell_logratios$med_logratio[cell_logratios$n > 400])
# colombia_logratio_401_max
# 
# 
# median(cell_logratios$med_logratio[cell_logratios$n > 0 & cell_logratios$n < 350])
# colombia_logratio_1_350
# 
# median(cell_logratios$med_logratio[cell_logratios$n > 0 & cell_logratios$n < 250])
# colombia_logratio_1_250

##### Grids with varying pixel size #####
library(dggridR)
library(sf)
library(raster)

# Assign each point to a hexagonal cell on a grid of variable resolution
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
points <- st_as_sf(forest_probs[,1:3], coords = c("x", "y"), crs = AEAstring)
points_latlong <- st_transform(points, 4326)
points_coords <- as.data.frame(st_coordinates(points_latlong))

dggs_4 <- dgconstruct(res=4)
dggs_7 <- dgconstruct(res=7)
dggs_10 <- dgconstruct(res=10)
dggs_8 <- dgconstruct(res=8)
dggs_9 <- dgconstruct(res=9)


cells_4 <- dgGEO_to_SEQNUM(dggs_4, points_coords$X, points_coords$Y)$seqnum
cells_7 <- dgGEO_to_SEQNUM(dggs_7, points_coords$X, points_coords$Y)$seqnum
cells_10 <- dgGEO_to_SEQNUM(dggs_10, points_coords$X, points_coords$Y)$seqnum
cells_8 <- dgGEO_to_SEQNUM(dggs_8, points_coords$X, points_coords$Y)$seqnum
cells_9 <- dgGEO_to_SEQNUM(dggs_9, points_coords$X, points_coords$Y)$seqnum

uc_4 <- unique(cells_4)
uc_7 <- unique(cells_7)
uc_10 <- unique(cells_10)
uc_8 <- unique(cells_8)
uc_9 <- unique(cells_9)

# Get average (mean) pointwise medians and the regional medians for each grid cell
cell_pointwise_logratios_4 <- cell_logratios_4 <- cell_richness_4 <- cell_pointwise_richness_4 <- cell_numbers_4 <- vector()
cell_pointwise_logratios_4_expanded <- cell_logratios_4_expanded <- rep(NA, nrow(points))
for(i in 1:length(uc_4)){
  print(i)
  grl <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.2, cell_positions = which(cells_4 == uc_4[i]))
  cell_logratios_4[i] <- grl$med_logratio
  cell_richness_4[i] <- grl$n
  cell_logratios_4_expanded[which(cells_4 == uc_4[i])] <- cell_logratios_4[i]
  cell_pointwise_logratios_4[i] <- mean(cell_logratios$med_logratio[which(cells_4 == uc_4[i])], na.rm=T)
  cell_pointwise_richness_4[i] <- mean(cell_logratios$n[cells_4 == uc_4[i] & !is.na(cell_logratios$avg_logratio)])
  
  cell_pointwise_logratios_4_expanded[which(cells_4 == uc_4[i])] <- cell_pointwise_logratios_4[i]
  cell_numbers_4[i] <- sum(cells_4 == uc_4[i])
}
raster_regional_4 <- rasterFromXYZ(cbind(forest_probs[,2:3], cell_logratios_4_expanded))
raster_pointwise_4 <- rasterFromXYZ(cbind(forest_probs[,2:3], cell_pointwise_logratios_4_expanded))
raster_difference_4 <- raster_regional_4 - raster_pointwise_4

cell_pointwise_logratios_7 <- cell_logratios_7 <- cell_richness_7 <- cell_pointwise_richness_7 <- cell_numbers_7 <- vector()
cell_pointwise_logratios_7_expanded <- cell_logratios_7_expanded <-  rep(NA, nrow(points))
for(i in 1:length(uc_7)){
  print(i)
  grl <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.2, cell_positions = which(cells_7 == uc_7[i]))
  cell_logratios_7[i] <- grl$med_logratio
  cell_richness_7[i] <- grl$n
  cell_logratios_7_expanded[which(cells_7 == uc_7[i])] <- cell_logratios_7[i]
  cell_pointwise_logratios_7[i] <- mean(cell_logratios$med_logratio[which(cells_7 == uc_7[i])], na.rm=T)
  cell_pointwise_richness_7[i] <- mean(cell_logratios$n[cells_7 == uc_7[i] & !is.na(cell_logratios$avg_logratio)])
  
  cell_pointwise_logratios_7_expanded[which(cells_7 == uc_7[i])] <- cell_pointwise_logratios_7[i]
  cell_numbers_7[i] <- sum(cells_7 == uc_7[i])
}
raster_regional_7 <- rasterFromXYZ(cbind(forest_probs[,2:3], cell_logratios_7_expanded))
raster_pointwise_7 <- rasterFromXYZ(cbind(forest_probs[,2:3], cell_pointwise_logratios_7_expanded))
raster_difference_7 <- raster_regional_7 - raster_pointwise_7

cell_pointwise_logratios_10 <- cell_logratios_10 <- cell_richness_10 <- cell_pointwise_richness_10 <- cell_numbers_10 <- vector()
cell_pointwise_logratios_10_expanded <- cell_logratios_10_expanded <-  rep(NA, nrow(points))
for(i in 1:length(uc_10)){
  print(i)
  grl <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.2, cell_positions = which(cells_10 == uc_10[i]))
  cell_logratios_10[i] <- grl$med_logratio
  cell_richness_10[i] <- grl$n
  cell_logratios_10_expanded[which(cells_10 == uc_10[i])] <- cell_logratios_10[i]
  cell_pointwise_logratios_10[i] <- mean(cell_logratios$med_logratio[which(cells_10 == uc_10[i])], na.rm=T)
  cell_pointwise_richness_10[i] <- mean(cell_logratios$n[cells_10 == uc_10[i] & !is.na(cell_logratios$avg_logratio)])
  
  cell_pointwise_logratios_10_expanded[which(cells_10 == uc_10[i])] <- cell_pointwise_logratios_10[i]
  cell_numbers_10[i] <- sum(cells_10 == uc_10[i])
}
raster_regional_10 <- rasterFromXYZ(cbind(forest_probs[,2:3], cell_logratios_10_expanded))
raster_pointwise_10 <- rasterFromXYZ(cbind(forest_probs[,2:3], cell_pointwise_logratios_10_expanded))
raster_difference_10 <- raster_regional_10 - raster_pointwise_10




# cell_pointwise_logratios_8 <- cell_logratios_8 <- cell_richness_8 <- cell_pointwise_richness_8 <- cell_numbers_8 <- vector()
# cell_pointwise_logratios_8_expanded <- cell_logratios_8_expanded <-  rep(NA, nrow(points))
# for(i in 1:length(uc_8)){
#   print(i)
#   grl <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.2, cell_positions = which(cells_8 == uc_8[i]))
#   cell_logratios_8[i] <- grl$med_logratio
#   cell_richness_8[i] <- grl$n
#   cell_logratios_8_expanded[which(cells_8 == uc_8[i])] <- cell_logratios_8[i]
#   cell_pointwise_logratios_8[i] <- mean(cell_logratios$med_logratio[which(cells_8 == uc_8[i])], na.rm=T)
#   cell_pointwise_richness_8[i] <- mean(cell_logratios$n[cells_8 == uc_8[i] & !is.na(cell_logratios$avg_logratio)])
#   
#   cell_pointwise_logratios_8_expanded[which(cells_8 == uc_8[i])] <- cell_pointwise_logratios_8[i]
#   cell_numbers_8[i] <- sum(cells_8 == uc_8[i])
# }
# raster_regional_8 <- rasterFromXYZ(cbind(forest_probs[,2:3], cell_logratios_8_expanded))
# raster_pointwise_8 <- rasterFromXYZ(cbind(forest_probs[,2:3], cell_pointwise_logratios_8_expanded))
# raster_difference_8 <- raster_regional_8 - raster_pointwise_8
# 
# 
# 
# cell_pointwise_logratios_9 <- cell_logratios_9 <- cell_richness_9 <- cell_pointwise_richness_9 <- cell_numbers_9 <- vector()
# cell_pointwise_logratios_9_expanded <- cell_logratios_9_expanded <-  rep(NA, nrow(points))
# for(i in 1:length(uc_9)){
#   print(i)
#   grl <- get_regional_logratios(forest_probs, pasture_probs, cutoff_type="absolute", cutoff=.2, cell_positions = which(cells_9 == uc_9[i]))
#   cell_logratios_9[i] <- grl$med_logratio
#   cell_richness_9[i] <- grl$n
#   cell_logratios_9_expanded[which(cells_9 == uc_9[i])] <- cell_logratios_9[i]
#   cell_pointwise_logratios_9[i] <- mean(cell_logratios$med_logratio[which(cells_9 == uc_9[i])], na.rm=T)
#   cell_pointwise_richness_9[i] <- mean(cell_logratios$n[cells_9 == uc_9[i] & !is.na(cell_logratios$avg_logratio)])
#   
#   cell_pointwise_logratios_9_expanded[which(cells_9 == uc_9[i])] <- cell_pointwise_logratios_9[i]
#   cell_numbers_9[i] <- sum(cells_9 == uc_9[i])
# }
# raster_regional_9 <- rasterFromXYZ(cbind(forest_probs[,2:3], cell_logratios_9_expanded))
# raster_pointwise_9 <- rasterFromXYZ(cbind(forest_probs[,2:3], cell_pointwise_logratios_9_expanded))
# raster_difference_9 <- raster_regional_9 - raster_pointwise_9


# Plot the regional losses by cell
plot(raster_regional_10)
plot(raster_regional_7)
plot(raster_regional_4)

# Plot the difference between the regional and the mean pointwise losses by cell
plot(raster_difference_10)
# plot(raster_difference_9)
# plot(raster_difference_8)
plot(raster_difference_7)
plot(raster_difference_4)

# Check how the mean (across cells) of the regional losses decreases with decreasing cell size
colombia_logratio$avg_logratio
mean(cell_logratios_4[c(3,4,6)], na.rm=T)
mean(cell_logratios_7, na.rm=T)
mean(cell_logratios_10, na.rm=T)
mean(cell_logratios$med_logratio, na.rm = T)

a <- (cell_logratios_10 - cell_pointwise_logratios_10)
b <- (cell_richness_10/cell_pointwise_richness_10)

c <- (cell_logratios_7 - cell_pointwise_logratios_7)
d <- (cell_richness_7/cell_pointwise_richness_7)

e <- (cell_logratios_4 - cell_pointwise_logratios_4)
f <- (cell_richness_4/cell_pointwise_richness_4)


# a <- (cell_logratios_8 - cell_pointwise_logratios_8)
# b <- (cell_richness_8/cell_pointwise_richness_8)

# a <- (cell_logratios_9 - cell_pointwise_logratios_9)
# b <- (cell_richness_9/cell_pointwise_richness_9)


plot(a ~ b, xlim = c(1, 15), ylim = c(-2,2))
points(c ~ d, col = "red")
points(e ~ f, col = "blue")
# points(4.3, .38, col = "green")


points(6.5, .59, col = "green")


summary(lm(a ~ b))
summary(lm(c ~ d))
summary(lm(e ~ f))

##### Elevational bands #####











sample_declines <- get_sample_percent_decline(forest_sample[use_points, ], pasture_sample[use_points, ], 10000)
sample_declines$nsp
sample_declines$total
hist(sample_declines$pointwise, breaks = 50)
mean(sample_declines$pointwise, na.rm = T)
mean(sample_declines$pointwise <= sample_declines$total, na.rm = T)
