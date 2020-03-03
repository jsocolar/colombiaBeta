library(sf)
library(magrittr)

setwd('/Users/jacobsocolar/Dropbox/Work/Colombia')

##### Basic Colombia map #####
# read in GADM colombia shapefile
colombia <- st_read('Data/GIS/colombia_maps/gadm36_COL_shp/gadm36_COL_0.shp')

# file consists of many disjoint polygons, representing the mainland and numerous islands. Figure out which is the mainland
# and extract it.
npoly <- length(colombia$geometry[[1]])
size <- rep(0,npoly)
for(i in 1:npoly){
  size[i] <- dim(colombia$geometry[[1]][[i]][[1]])[1]
}

mainland <- colombia$geometry[[1]][[which(size == max(size))]] %>%
                st_polygon() %>% st_sfc() %>% st_sf()
st_crs(mainland) <- st_crs(colombia)

# Transform to AEA conic centered on Colombia
mainland <- st_transform(mainland, "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
plot(mainland)

# create new polygon with 100 km buffer
b_mainland <- st_buffer(mainland, dist = 100000)
plot(b_mainland)

##### Birdlife maps #####
botw <- "Data/GIS/birdlife_maps/BOTW/BOTW.gdb"
fc_list = ogrListLayers(botw)
range_maps <- sf::st_read(dsn=botw,layer="All_Species")
birdlife.species <- unique(range_maps$SCINAME)
nsp <- length(birdlife.species)

map1 <- range_maps[1,]
plot(map1)
#dim(range_maps)
#str(range_maps)


##### Hydrosheds #####
#This code is nothing at the moment--just a place to save some bits that will potentially be useful to have in the future

basins <- rgdal::readOGR('/Users/jacobsocolar/Downloads/hybas_sa_lev01-06_v1c/hybas_sa_lev03_v1c.shp')
plot(basins[1, ], col='red', xlim=c(-79.2, -66.7), ylim=c(-4.3, 12.5))
plot(basins[2, ], col='blue', add = T) # Magdalena sensu lato
plot(basins[3, ], col='green', add = T) # Maracaibo sensu lato
plot(basins[4, ], col='purple', add = T) # Orinoco
plot(basins[7, ], col='orange', add = T) # Amazon
plot(basins[25, ], col='brown', add = T) # Pacific

basins <- rgdal::readOGR('/Users/jacobsocolar/Downloads/hybas_sa_lev01-06_v1c/hybas_sa_lev05_v1c.shp')
plot(basins[3, ], col='gray', add = T)
plot(basins[4, ], col='gray', add = T)
plot(basins[5, ], col='gray', add = T) # Cauca

basins <- rgdal::readOGR('/Users/jacobsocolar/Downloads/hybas_sa_lev01-06_v1c/hybas_sa_lev06_v1c.shp')
#plot(basins)
plot(basins[6, ], add =T) # Sinu
