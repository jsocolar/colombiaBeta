# Generate prediction info data objects in varying formats (stars, data.table)
#
# Creates layers necessary to simulate the posterior occupancy probability for 
# each species across Colombia at 2 km resolution
#
# Note: rasterising Ayerbe maps and saving them also seems redundant.. could 
# just run one loop and store prediction info..

# housekeeping
library(sf); library(reticulate); library(raster)

`%ni%` <- Negate(`%in%`)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 
    +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

##### Get elevation raster across Colombia #####
# Set up GEE session
use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()             # Trigger the authentication
np <- import("numpy")       # Import Numpy        needed for converting gee raster to R raster object
pd <- import("pandas")      # Import Pandas       ditto the above

# Get elevations for Colombia
countries <- ee$FeatureCollection('USDOS/LSIB_SIMPLE/2017')
roi <- countries$filterMetadata('country_na', 'equals', 'Colombia')  
DEM_full <- ee$Image("JAXA/ALOS/AW3D30/V2_2")
DEM <- DEM_full$select(list("AVE_DSM"),list("elevation"))
latlng <- ee$Image$pixelLonLat()$addBands(DEM)
latlng <- latlng$reduceRegion(reducer = ee$Reducer$toList(),
                              geometry = roi,
                              maxPixels = 10^9,
                              scale=2000)
# Convert to arrays
lats <- np$array((ee$Array(latlng$get("latitude"))$getInfo()))
lngs <- np$array((ee$Array(latlng$get("longitude"))$getInfo()))
elevs <- np$array((ee$Array(latlng$get("elevation"))$getInfo()))
# Convert to elevation raster
elevation <- data.frame(x = lngs, y = lats, elevation = elevs)
raster_elev <- raster::rasterFromXYZ(elevation,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# Reproject
raster_elev_AEA <- raster::projectRaster(raster_elev, crs = sp::CRS(AEAstring))

# Crop
mask <- readRDS("outputs/WWF_terr_ecoregions_mask.rds")
source("code/GIS_processing/get_mainland.R")
mainland <- st_transform(mainland, AEAstring)
raster_elev_AEA <- crop(raster_elev_AEA, extent(mainland))
# Save raster
raster::writeRaster(raster_elev, "outputs/elev_raster/raster_elev.grd", overwrite = T)
raster::writeRaster(raster_elev_AEA, "outputs/elev_raster/raster_elev_AEA.grd", overwrite = T)

raster_elev_AEA_masked <- st_crop(raster_elev_AEA, mask)
write_stars(raster_elev_AEA_masked, 
            "outputs/elev_raster/raster_elev_AEA_masked.grd", 
            overwrite = T)

##### Rasterize the Ayerbe maps #####
ayerbe_list_updated <- readRDS('outputs/ayerbe_list_updated.RDS')
ayerbe_buffered_ranges_updated <- readRDS("outputs/ayerbe_buffered_ranges_updated.RDS")

n_files <- length(list.files('outputs/Ayerbe_rasters/buffered/', '.grd'))
if(length(ayerbe_list_updated) != n_files) { 
    
    pb <- txtProgressBar(min = 0, max = length(ayerbe_buffered_ranges_updated), initial = 0, style = 3) 
    for(i in 1:length(ayerbe_buffered_ranges_updated)){
        setTxtProgressBar(pb,i)
        ayerbe_raster <- fasterize::fasterize(st_sf(st_union(ayerbe_list_updated[[i]])), raster_elev_AEA)
        
        raster::writeRaster(ayerbe_raster, 
                            paste0('outputs/Ayerbe_rasters/regular/', names(ayerbe_list_updated)[i], '.grd'), 
                            overwrite = TRUE)
        ayerbe_buffer_raster <- fasterize::fasterize(st_sf(ayerbe_buffered_ranges_updated[[i]]), raster_elev_AEA)
        raster::writeRaster(ayerbe_buffer_raster, 
                            paste0('outputs/Ayerbe_rasters/buffered/', names(ayerbe_list_updated)[i], '_buffered.grd'), 
                            overwrite = TRUE)
    }
}

# extract indexing from stars object ----
col_index_stars <- raster_elev_AEA_masked %>%
    setNames("elevation") %>%
    mutate(id_cell = 1:n(), id_x=NA, id_y = NA, elevation=NULL)

xy_dim <- dim(col_index_stars)
col_index_stars[["id_x"]] <- matrix(rep(1:xy_dim[2], each=xy_dim[1]), 
                                    nrow=xy_dim[2], ncol=xy_dim[1])
col_index_stars[["id_y"]] <- matrix(rep(1:xy_dim[1], xy_dim[2]), 
                                    nrow=xy_dim[2], ncol=xy_dim[1])

col_index_sf <- st_as_sf(col_index_stars, as_points = TRUE)
col_index_dt <- as.data.table(col_index_sf)
col_index_dt[, geometry := NULL]

## save indexing stars object and data.table
saveRDS(col_index_dt, "outputs/xy_lookup.rds")
saveRDS(col_index_stars, "outputs/xy_lookup_stars.rds")

# calculate distance-to-range ----
# get CO mainland 
source("code/GIS_processing/get_mainland.R")
mainland <- st_transform(mainland, AEAstring)
mainland_inward <- st_buffer(mainland, -7000)

# get elevation scaling info etc.
ayerbe_list_updated <- readRDS('outputs/ayerbe_list_updated.RDS')
bird_data <- readRDS("outputs/bird_stan_data6_package.RDS")
birds <- readRDS("outputs/birds.RDS")
sp_list <- unique(birds$species)

dtr_offset <-bird_data$means_and_sds$distance_to_range_offset
dtr_sd <- bird_data$means_and_sds$distance_to_range_sd
dtr_rescale <- bird_data$means_and_sds$distance_to_range_logit_rescale

# convert elevation raster to points
points <- st_as_sf(raster_elev_AEA_masked, as_points = TRUE)

## create elevation sf ----
elev_sf <- st_as_sf(raster_elev_AEA_masked, as_points=TRUE, na.rm = FALSE)

# add elevation to this for use below
dt_template <- copy(col_index_dt)
dt_template[,`:=`(id_x = NULL, id_y = NULL)]
dt_template[,elev := elev_sf$raster_elev_AEA.grd]

## blank grid for templating ----
blank_grid <- raster_elev_AEA_masked %>%
    setNames("elevation") %>%
    mutate(elevation = NA)

## make write directory ----
dir.create("outputs/pred_dt_species/", showWarnings = FALSE)

# do calculation ----
pb <- txtProgressBar(min = 0, max = length(sp_list), initial = 0, style = 3) 
for(i in 1:length(sp_list)){
  setTxtProgressBar(pb,i)
  sp <- sp_list[i]
  ayerbe_buffer_polygon <- st_union(ayerbe_buffered_ranges_updated[[gsub("_", " ", sp)]])
  indices_sp <- st_intersects(points, ayerbe_buffer_polygon, sparse = F)
  indices_m <- which(indices_sp)
  points_m <- points[indices_m,]

  ayerbe_polygon <- st_union(ayerbe_list_updated[[gsub("_", " ", sp)]])
  range_linestring <- st_cast(ayerbe_polygon, "MULTILINESTRING")
  
  range_linestring_cropped <- st_intersection(mainland_inward, range_linestring) 
  
  points_m <- points_m %>%
      mutate(distance = st_distance(., range_linestring_cropped),
             inside = st_intersects(., ayerbe_polygon, sparse = F), 
             distance2 = distance * c(1, -1)[inside + 1])

  ayerbe_dist_raster_i <- cbind(st_coordinates(points_m), points_m$distance2) %>%
      as.data.frame %>%
      raster::rasterFromXYZ(.,crs=AEAstring)

  
  ayerbe_dist_raster_i <- st_rasterize(points_m["distance2"], template = blank_grid)

  ayerbe_dist_sf_i <- st_as_sf(ayerbe_dist_raster_i, as_points = TRUE, na.rm = FALSE)
    
  sp_lower <- unique(birds$lower[birds$species == sp])
  sp_upper <- unique(birds$upper[birds$species == sp])
  sp_breadth <- sp_upper - sp_lower
  
  dt_template[, distance := round(ayerbe_dist_sf_i$distance2, 0)]  
  
  dt_i <- dt_template[!is.na(distance)]
  dt_i[,elev := ifelse(elev > (sp_lower - sp_breadth) & 
                           elev < (sp_upper + sp_breadth), 
                       elev, NA)]
  dt_i[, relev := ((elev - sp_lower)/sp_breadth - relev_offset)/relev_sd]
  dt_i[, species := sp]  
  
  saveRDS(dt_i[!is.na(relev)], paste0("outputs/pred_dt_species/", sp, ".rds"))
}

# source("code/analysis_and_plotting/helper_functions.R")
# ggplot(dt_i[!is.na(relev)],
#        aes(cell_to_x_pos(id_cell, c(712, 517)), 
#            cell_to_y_pos(id_cell, c(712, 517)), fill=relev)) +
#     geom_tile() +
#     coord_equal()

## bind all dts together ----
