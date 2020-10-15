library(sf)
library(reticulate)

setwd("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS")


`%ni%` <- Negate(`%in%`)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# Load point locations
all_pts <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Points/all_pts.RDS")

##### Set up GEE session #####
use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()             # Trigger the authentication
np <- import("numpy")       # Import Numpy        needed for converting gee raster to R raster object
pd <- import("pandas")      # Import Pandas       ditto the above


##### Get ALOS elevation raster #####
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
rasterElev <- raster::rasterFromXYZ(elevation,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


##### Get continent outlines #####
#Download the continents shapefile
download.file("http://faculty.baruch.cuny.edu/geoportal/data/esri/world/continent.zip",
              destfile = "continents/cont.zip")
#Unzip it
unzip("continents/cont.zip", exdir="continents")
#Load it
cont <- st_read("continents/continent.shp")

n_am <- st_cast(cont[cont$CONTINENT == "North America", ]$geometry, "POLYGON")
n_am_main <- n_am[which(st_area(n_am) == max(st_area(n_am)))]

s_am <- st_cast(cont[cont$CONTINENT == "South America", ]$geometry, "POLYGON")
s_am_main <- s_am[which(st_area(s_am) == max(st_area(s_am)))]

americas <- st_union(n_am_main, s_am_main)


##### Plotting #####
col <- colorspace::lighten(viridis::viridis(501), amount = .4)

plot(americas, col = 'gray90', border = 'gray90', xlim = c(-80, -66), ylim = c(-4.6, 13.1))
plot(rasterElev, col = col, legend = T, add = T)

points(all_pts$lon, all_pts$lat, pch = 16, col = rgb(0, 0, 0, max = 255, alpha = 15))
points(all_pts$lon, all_pts$lat, pch = 1, col = rgb(0, 0, 0, max = 255, alpha = 255))

