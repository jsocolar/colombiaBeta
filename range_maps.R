library(sf)
library(magrittr)

`%ni%` <- Negate(`%in%`)

##### For collaborative projects--figure out what machine we're on and automatically set the working directory ####
socolar.desktop <- file.exists('/Users/jacobsocolar/Dropbox/Work/Code/machine_identifier_n5L8paM.txt')
socolar.laptop <- file.exists('/Users/jacob/Dropbox/Work/Code/machine_identifier_n5L8paM.txt')
if(socolar.desktop){
  dir.path <- "/Users/JacobSocolar/Dropbox/Work/Colombia"
}else if(socolar.laptop){
  dir.path <- "/Users/jacob/Dropbox/Work/Colombia"
}# else if(){dir.path <- }
# Edit the above for whatever computer(s) you use.  Just make absolutely sure that the if condition is something that definitely
# wouldn't possibly evaluate as true on anybody else's system, and that none of the preceding conditions could possibly evaluate
# to TRUE on your system!  (This isn't just about making sure that we get the right working directories; in some cases we might
# conceivably invoke system commands for file management that depend on dir.path.)
setwd(dir.path)
############################

AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
'%ni%' <- Negate('%in%')

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
mainland <- st_transform(mainland, AEAstring)
plot(mainland)

# create new polygon with 100 km buffer
b_mainland <- st_buffer(mainland, dist = 100000)
plot(b_mainland)

##### Birdlife maps #####
botw <- "Data/GIS/birdlife_maps/BOTW/BOTW.gdb"
fc_list = ogrListLayers(botw)
orig_range_maps <- st_read(dsn=botw,layer="All_Species")
recast_range_maps <- st_cast(orig_range_maps, to = "MULTIPOLYGON")
range_maps <- st_transform(recast_range_maps, AEAstring)

birdlife.species <- unique(range_maps$SCINAME)
nsp <- length(birdlife.species)

##### Getting a unified list of species of interest #####
b_mainland_latlon <- st_transform(b_mainland, st_crs(recast_range_maps))

colombia_overlaps <- st_intersects(recast_range_maps, b_mainland_latlon)
# although coordinates are longitude/latitude, st_intersects assumes that they are planar
# I think the above is equivalent to checking for intersections on a Mercator projection.  It won't be a problem here.
colombia_species <- unique(range_maps$SCINAME[unlist(lapply(colombia_overlaps, FUN = function(i){return(length(i) != 0)}))])

taxonomy <- read.csv("Data/Birds/HBW_eBird_taxonomy.csv", stringsAsFactors = F)
t2 <- taxonomy[taxonomy$HBW_CAT == "sp" | taxonomy$CLEM_CAT_2019 == "sp", ] # remove taxa that are not treated as species by either eBird or HBW
t2$HBW_LATIN[t2$HBW_LATIN == "Nyctibius bracteatus"] <- "Phyllaemulor bracteatus"
t2$HBW_LATIN[t2$HBW_LATIN == "Gallinula melanops"] <- "Porphyriops melanops"
t2$HBW_LATIN[t2$HBW_LATIN == "Claravis mondetoura"] <- "Paraclaravis mondetoura"
t2$HBW_LATIN[t2$CLEM_SCI_2019 == "Scytalopus alvarezlopezi"] <- "Scytalopus alvarezlopezi"
t2$HBW_LATIN[t2$CLEM_SCI_2019 == "Vireo chivi"] <- "Vireo chivi"
t2$HBW_LATIN[t2$CLEM_SCI_2019 == "Megascops gilesi"] <- "Megascops gilesi"
t2$HBW_LATIN[t2$CLEM_SCI_2019 == "Anthocephala berlepschi"] <- "Anthocephala berlepschi"

t2 <- rbind(t2, data.frame(orig_sort = NA, concept_8 = NA, latin_name = NA, TAXON_nid_MATCHED = NA, HBW_CAT = "sp",
                           HBW_LATIN = c("Catharus maculatus", "Pyrrhura chapmani"), HBW_LATIN_rev = NA, HBW_EN_EBIRD = NA,
                           HBW.Clem = "sp-ssp", cat_match_2019 = "HBW split", sci_match = NA, CLEM_SORT = NA,
                           CLEM_SCI_2019 = c("Catharus dryas", "Pyrrhura melanura"), CLEM_ENG_2019 = NA, notes = NA,
                           CLEM_CAT_2019 = NA, SPECIES_CODE_2019 = NA, CLEM_SPECIES_CODE_2018 = NA, range = NA, 
                           SORT_INTEGRATED = NA))

csp2 <- t2[t2$HBW_LATIN %in% colombia_species, ]

elevation_dataset <- read.csv("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/Bird_elevations_modified.csv", stringsAsFactors = F)

species_to_check <- csp2$HBW_LATIN[csp2$HBW_LATIN %ni% elevation_dataset$Scientific & 
                                     !((csp2$CLEM_SCI_2019 %in% elevation_dataset$Scientific) & csp2$HBW.Clem == "sp-sp")]
# These are the species from the range-maps that are not included in the compiled elevation dataset.
# Many of these do not occur in Colombia. Those that do are mostly splits or naming issues.

write.csv(species_to_check, file = '/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/species_to_check.csv')

elevation_species <- elevation_dataset$Scientific[is.na(elevation_dataset$coastal_pelagic) & is.na(elevation_dataset$outside_birdlife)]
s2c2 <- elevation_species[(elevation_species %ni% csp2$HBW_LATIN) & (elevation_species %ni% csp2$CLEM_SCI_2019[csp2$HBW.Clem == 'sp-sp'])]

changes <- read.csv('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/species_to_check_checked.csv')
s2c2[(s2c2 %ni% changes$synonym) & (s2c2 %ni% changes$clean.split.lump) & (s2c2 %ni% changes$messy.split.lump)]

colombia_species_extra <- t2$CLEM_SCI_2019[(t2$HBW_LATIN %in% colombia_species) & (t2$HBW.Clem == "sp-sp")]

colombia_species[colombia_species %ni% t2$CLEM_SCI_2019]

elevation_dataset$Nombre.cientifico...Scientific.name[elevation_dataset$Nombre.cientifico...Scientific.name %ni% t2$CLEM_SCI_2019]
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
