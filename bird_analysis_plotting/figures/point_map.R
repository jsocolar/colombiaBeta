library(sf)
library(reticulate)

source('/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/GIS_processing/get_mainland.R')
mainland2 <- st_simplify(mainland, dTolerance = .01)

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


source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/GIS_processing/hydrosheds_extraction.R")

##### Plotting #####
col <- colorspace::lighten(viridis::viridis(501), amount = .4)

plot(americas, col = 'gray90', border = 'gray90', xlim = c(-80, -66), ylim = c(-4.6, 13.1))
plot(rasterElev, col = col, legend = T, add = T)

points(all_pts$lon, all_pts$lat, pch = 16, col = rgb(0, 0, 0, max = 255, alpha = 15))
points(all_pts$lon, all_pts$lat, pch = 1, col = rgb(0, 0, 0, max = 255, alpha = 255))

regions <- list(pacific = pacific, cauca_west = cauca_west, cauca_east = cauca_east, magdalena_west = magdalena_west,
                magdalena_east = magdalena_east, amazon_orinoco = amazon_orinoco, snsm = snsm, guajira_perija = guajira_perija,
                catatumbo = catatumbo, pasto = pasto, tacarcuna = tacarcuna)

colors <- viridis::viridis(11)[c(5,11,3,10,1,4,9,6,8,7,2)]
plot(americas, col = 'gray90', border = 'gray90', xlim = c(-80, -66), ylim = c(-4.6, 13.1))
for(i in 1:length(regions)){
  plot(st_intersection(regions[[i]], mainland2), col = colors[i], border = colors[i], add = T)
}
# west Andes west; Cauca Valley west; Cauca Valley east; Magdalena Valley west; Magdalena Valley east & north;
# east Andes east; Tacarcuna foothills; Sierra Nevada de Santa Marta; Guajira/Perijá; Catatumbo; Pasto

points(all_pts$lon, all_pts$lat, pch = 16, col = rgb(0, 0, 0, max = 255, alpha = 15))
points(all_pts$lon, all_pts$lat, pch = 1, col = rgb(0, 0, 0, max = 255, alpha = 255))

