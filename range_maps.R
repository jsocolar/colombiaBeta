library(rgdal)
library(sp)
library(magrittr)

setwd('/Users/jacobsocolar/Dropbox/Work/Colombia')

# read in GADM colombia shapefile
colombia <- rgdal::readOGR('Data/GIS/colombia_maps/gadm36_COL_shp/gadm36_COL_0.shp')

# file consists of many disjoint polygons, representing the mainland and numerous islands. Figure out which is the mainland
# and extract it.
npoly <- length(colombia@polygons[[1]]@Polygons)
size <- rep(0,npoly)
for(i in 1:npoly){
  size[i] <- dim(colombia@polygons[[1]]@Polygons[[i]]@coords)[1]
}

mainland <- colombia@polygons[[1]]@Polygons[[which(size == max(size))]] %>% 
  list() %>% Polygons(1) %>% list() %>% SpatialPolygons()
  
proj4string(mainland) <- proj4string(colombia)
plot(mainland)

