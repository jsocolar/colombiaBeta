library(sf)
library(reticulate)

`%ni%` <- Negate(`%in%`)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"


##### Set up GEE session #####
use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()             # Trigger the authentication
np <- import("numpy")       # Import Numpy        needed for converting gee raster to R raster object
pd <- import("pandas")      # Import Pandas       ditto the above
countries <- ee$FeatureCollection('USDOS/LSIB_SIMPLE/2017')
roi <- countries$filterMetadata('country_na', 'equals', 'Colombia')  
DEM_full <- ee$Image("JAXA/ALOS/AW3D30/V2_2")
DEM <- DEM_full$select(list("AVE_DSM"),list("elevation"))
latlng <- ee$Image$pixelLonLat()$addBands(DEM)
latlng <- latlng$reduceRegion(reducer = ee$Reducer$toList(),
                                geometry = roi,
                                maxPixels = 10^9,
                                scale=500)
  
# Convert to arrays
lats <- np$array((ee$Array(latlng$get("latitude"))$getInfo()))
lngs <- np$array((ee$Array(latlng$get("longitude"))$getInfo()))
elevs <- np$array((ee$Array(latlng$get("elevation"))$getInfo()))
  
# Convert to elevation raster
elevation <- data.frame(x = lngs, y = lats, elevation = elevs)
raster_elev <- raster::rasterFromXYZ(elevation,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

raster_elev_AEA <- raster::projectRaster(raster_elev, crs = sp::CRS(AEAstring))

plot(raster_elev_AEA)

###### Rasterize the Ayerbe maps
ayerbe_list_updated <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/ayerbe_maps/ayerbe_list_updated.RDS')
ayerbe_buffered_ranges_updated <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/ayerbe_maps/ayerbe_buffered_ranges_updated.RDS")

ayerbe_raster <- ayerbe_buffer_raster <- list()
for(i in 1:length(ayerbe_buffered_ranges_updated)){
  print(i)
  ayerbe_raster[[i]] <- fasterize::fasterize(st_transform(ayerbe_list_updated[[i]], 4326), raster_elev)
  ayerbe_buffer_raster[[i]] <- fasterize::fasterize(st_transform(st_sf(ayerbe_buffered_ranges_updated[[i]]), 4326), raster_elev)
}

ayerbe_rasters <- list(species = names(ayerbe_list_updated), ayerbe_raster = ayerbe_raster, ayerbe_buffer_raster = ayerbe_buffer_raster)

saveRDS(ayerbe_rasters, "/Users/JacobSocolar/Dropbox/Work/Colombia/Data/Analysis/ayerbe_rasters.RDS")









birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_data_trimmed.RDS")
all_pts <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Points/all_pts.RDS")


birds_q <- birds[birds$Q == 1, ]
quantile(birds_q$elev_sp_standard, .0155)
quantile(birds_q$elev_sp_standard, .95562)

# 94% of species-point combinations with detections are within Ayerbe elevational ranges