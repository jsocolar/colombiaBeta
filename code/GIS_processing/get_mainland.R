# get shapefile for mainland Colombia

# housekeeping ----
library(sf)

# read in GADM colombia shapefile
colombia <- st_read('inputs/GIS/colombia_maps/gadm36_COL_shp/gadm36_COL_0.shp')

# File consists of many disjoint polygons, representing the mainland and numerous islands. Figure out which is the mainland
# and extract it.
npoly <- length(colombia$geometry[[1]])
size <- rep(0,npoly)
for(i in 1:npoly){
  size[i] <- dim(colombia$geometry[[1]][[i]][[1]])[1]
}

mainland <- colombia$geometry[[1]][[which(size == max(size))]] %>%
  st_polygon() %>% st_sfc() %>% st_sf()
st_crs(mainland) <- st_crs(colombia)

# no outputs, where does this get called from?