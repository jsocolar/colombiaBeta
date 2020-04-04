# Experimenting with getting data analyzed via GEE from the R console and bringing the result directly back into R.

# IMPORTANT: In order for this to work at all, it is necessary to first create a GEE account, download Anaconda, and 
# run GEE_setup.sh from a bash shell. I make no guarantees about whether GEE_setup.sh works properly on Windows.
# I experimented with running the entire bash script from R via system2(), but ultimately we'll have to manually input 
# our GEE authentication keys either way, so it's quite cumbersome to do entirely via R, and we don't get the benefit of
# complete automation anyway, since we'll need to manually input the GEE key no matter what.

# Some parts of code adapted from https://philippgaertner.github.io/2019/12/earth-engine-rstudio-reticulate/

library(reticulate)

# Use machine identifier to automatically set directory with points file
jorgen.desktop <- file.exists('C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD\\machine_identifier_lu847jp1o.txt')
socolar.desktop <- file.exists('/Users/jacobsocolar/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
socolar.laptop <- file.exists('/Users/jacob/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')

if(jorgen.desktop){
  dir.path <- "C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD\\Data_raw\\elevations"
}else if(socolar.desktop){
  dir.path <- '/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Points'
}else if(socolar.laptop){
  dir.path <- '/Users/jacob/Dropbox/Work/Colombia/Data/GIS/Points'
}


##### Set up GEE session #####
#Sys.setenv(PATH = paste(c("/Applications/anaconda3/bin", Sys.getenv("PATH")), collapse = .Platform$path.sep))
use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()             # Trigger the authentication
np <- import("numpy")       # Import Numpy        needed for converting gee raster to R raster object
pd <- import("pandas")      # Import Pandas       ditto the above

# Read coordinate file
if(jorgen.desktop){
  pts <- read.csv(paste(dir.path,"\\CO_sampling_points_metafile.csv",sep=""),stringsAsFactors = F)
}else if(socolar.desktop | socolar.laptop){
  pts <- read.csv(paste0(dir.path,"/CO_sampling_points_metafile.csv"),stringsAsFactors = F)
}
pts <- pts[-which(pts$long > 0),]   # Needed due to a cluster with error in the current metafile

# Load raster files
ALOS <- ee$Image('JAXA/ALOS/AW3D30/V2_2')
SRTM30 <- ee$Image("USGS/SRTMGL1_003")
RESOLVE <- ee$FeatureCollection("RESOLVE/ECOREGIONS/2017")
reg_ras <- RESOLVE$reduceToImage(reducer = ee$Reducer$first(),properties = list("ECO_ID"))$unmask(0)$reproject('epsg:4326',NULL,ee$Number(30.922080775909325))
biom_ras <- RESOLVE$reduceToImage(reducer = ee$Reducer$first(),properties = list("BIOME_NUM"))$unmask(0)$reproject('epsg:4326',NULL,ee$Number(30.922080775909325))

##### Extract ALOS elevation (image) from a given point #####
point_name <- "RBP1"
point_of_interest <- ee$Geometry$Point(pts[which(pts$point_id==point_name),]$long, pts[which(pts$point_id==point_name),]$lat)
poi_elev <- ALOS$reduceRegions(point_of_interest, ee$Reducer$mean())$getInfo()
poi_elev$features[[1]]$properties$AVE_DSM

##### Extract RESOLVE ecoregion (featurecollection) from a given point #####
point_name <- "MOP1"
point_of_interest <- ee$Geometry$Point(pts[which(pts$point_id==point_name),]$long, pts[which(pts$point_id==point_name),]$lat)
poi_ecoreg <- RESOLVE$filterBounds(point_of_interest)$getInfo()$features[[1]]$properties$ECO_NAME
poi_biome <- RESOLVE$filterBounds(point_of_interest)$getInfo()$features[[1]]$properties$BIOME_NAME

##### Extract RESOLVE ecoregion/biomes for multiple points (without raster conversion)#####
pts_ecoreg <- pts_biome <- rep(NA, nrow(pts))
for(i in 1:nrow(pts)){
  poi_info <- RESOLVE$filterBounds(ee$Geometry$Point(pts$long[i], pts$lat[i]))$getInfo()$features[[1]]$properties
  pts_ecoreg[i] <- poi_info$ECO_NAME
  pts_biome[i] <- poi_info$BIOME_NAME
}

##### Extract raster values from all points #####
# Featurecollection of point geometries
geompts <- sapply(1:nrow(pts),function(x)ee$Geometry$Point(c(pts$long[x],pts$lat[x])))
geompts <- ee$FeatureCollection(c(unlist(geompts)))

# Extract raster values for all points
pts_ALOS <- ALOS$select('AVE_DSM')$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
ALOSelev <- sapply(c(1:length(pts_ALOS$features)),function(x)pts_ALOS$features[[x]]$properties$mean)

pts_SRTM <- SRTM30$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
SRTMelev <- sapply(c(1:length(pts_SRTM$features)),function(x)pts_SRTM$features[[x]]$properties$mean)

pts_REG <- reg_ras$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
ecoreg <- sapply(c(1:length(pts_REG$features)),function(x)pts_REG$features[[x]]$properties$mean)

pts_BIOM <- biom_ras$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
biome <- sapply(c(1:length(pts_BIOM$features)),function(x)pts_BIOM$features[[x]]$properties$mean)

# Combined dataframe
cbind.data.frame(point_id=pts$point_id,ALOSelev,SRTMelev,ecoreg,biome)

##### Extract mean elevation from 100-m buffer around point (or many points) #####
buffer.width <- 100      # radius, in meters
max.error <- 1             # maximum error (controls number of vertices), in meters
geomcircs <- sapply(1:nrow(pts),function(x)ee$Geometry$Point(c(pts$long[x],pts$lat[x]))$buffer(buffer.width,max.error))
geomcircs <- ee$FeatureCollection(c(unlist(geomcircs)))
circs_elev <- ALOS$reduceRegions(geomcircs, ee$Reducer$mean())$getInfo()
ALOSelev_circ <- sapply(c(1:length(circs_elev$features)),function(x)circs_elev$features[[x]]$properties$AVE_DSM)
ALOSelev_circ <- cbind.data.frame(point_id=as.character(pts$point_id),ALOSelev_circ)

##### Do some fancier calculation on raster at a point (or many points) in GEE; import result to R #####
# For example, get slope and aspect data from ALOS

# Extract slope for all points - combine into dataframe - NB: slope algorithm has NOT been checked
ALOS_slope <- ee$Terrain$slope(ALOS)
pts_slope <- ALOS_slope$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
ALOSslope <- sapply(c(1:length(pts_slope$features)),function(x)pts_slope$features[[x]]$properties$mean)
ALOSslope <- cbind.data.frame(point_id=as.character(pts$point_id),ALOSslope)

# Extract aspect for all points - combine into dataframe - NB: aspect algorithm has NOT been checked
ALOS_aspect <- ee$Terrain$aspect(ALOS)
pts_aspect <- ALOS_aspect$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
ALOSaspect <- sapply(c(1:length(pts_aspect$features)),function(x)pts_aspect$features[[x]]$properties$mean)
ALOSaspect <- cbind.data.frame(point_id=as.character(pts$point_id),ALOSaspect)


##### Extract actual elevation raster over some defined polygon (e.g. a 10 km buffer around a point), import to R as raster object #####
# (This might be more relevant if we want to automatically pull in GFC rasters at some future date to derive complicated fragstat-style measures that we can't figure out how to code directly in GEE)
buffer.width <- 5000
max.error <- 1
point_name <- "SAF1"
poi <- ee$Geometry$Point(pts[which(pts$point_id==point_name),]$long, pts[which(pts$point_id==point_name),]$lat)
latlng <- ee$Image$pixelLonLat()$addBands(ALOS)

latlng = latlng$reduceRegion(reducer = ee$Reducer$toList(),
                             geometry = poi$buffer(buffer.width,max.error),
                             maxPixels = 1e10,
                             scale=30.922080775909325)

# Convert to numpy arrays
lats <- np$array((ee$Array(latlng$get("latitude"))$getInfo()))
lngs <- np$array((ee$Array(latlng$get("longitude"))$getInfo()))
ras_vals = np$array((ee$Array(latlng$get("AVE_DSM"))$getInfo()))

# Convert to raster
ras <- data.frame(x = lngs, y = lats, ras = ras_vals)
RasterR <- raster::rasterFromXYZ(ras)
raster::plot(RasterR)
raster::crs(RasterR) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "

# Write raster to disk
#raster::writeRaster(RasterR,"RasterR.tiff",overwrite=T)


##### Visualise elevation raster with rayshader #####
library(rayshader)
rasval_m = raster_to_matrix(RasterR$ras)
pixscale <- 30.922080775909325

rasval_m %>%
  sphere_shade(texture = "desert") %>%
  add_shadow(ray_shade(rasval_m, zscale = 3), 0.5) %>%
  add_shadow(ambient_shade(rasval_m), 0.5) %>%
  plot_3d(rasval_m, zscale = pixscale, fov = 30, theta = 45, phi = 35, windowsize = c(1000,800), zoom = 0.6)

render_depth(focus = 0.6, focallength = 200, clear = TRUE)
render_snapshot(clear = TRUE)
