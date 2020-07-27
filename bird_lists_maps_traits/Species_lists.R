library(sf)
library(magrittr)

`%ni%` <- Negate(`%in%`)

##### For collaborative projects--figure out what machine we're on and automatically set the working directory ####
socolar.desktop <- file.exists('/Users/jacobsocolar/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
socolar.laptop <- file.exists('/Users/jacob/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
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

# This script is the workhorse that I used to get a unified species list based on Birdlife + Donegan, using HBW Taxonomy
# Parts of the process were conducted manually using saved files.
# Annotations in this script explain the full process, including the manual parts

##### Basic Colombia map #####
# read in GADM colombia shapefile
colombia <- st_read('Data/GIS/colombia_maps/gadm36_COL_shp/gadm36_COL_0.shp')

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

# Transform to AEA conic centered on Colombia, for accurate 100 km buffering
mainland <- st_transform(mainland, AEAstring)
#plot(mainland)

# create new polygon with 100 km buffer
b_mainland <- st_buffer(mainland, dist = 100000)
plot(b_mainland)

##### Birdlife maps #####
botw <- "Data/GIS/birdlife_maps/BOTW/BOTW.gdb"
st_layers(botw)
orig_range_maps <- st_read(dsn=botw,layer="All_Species")
recast_range_maps <- st_cast(orig_range_maps, to = "MULTIPOLYGON")
save(recast_range_maps, file = "Data/GIS/birdlife_maps/recast_range_maps.Rdata")
load("Data/GIS/birdlife_maps/recast_range_maps.Rdata")

birdlife.species <- unique(recast_range_maps$SCINAME)
nsp <- length(birdlife.species)

##### Getting a unified list of species of interest #####
# Transform the buffered Colombia shapefile back to the lat-lon that the range maps use
b_mainland_latlon <- st_transform(b_mainland, st_crs(recast_range_maps))

# Extract a list of all species that ovelap the buffered polygon
colombia_overlaps <- st_intersects(recast_range_maps, b_mainland_latlon)
# "although coordinates are longitude/latitude, st_intersects assumes that they are planar"
# I think the above is equivalent to checking for intersections on a Mercator projection.  It won't be a problem here.
save(colombia_overlaps, file = "Data/Birds/species_list_creation/HBW_colombia_overlaps.Rdata")
load("Data/Birds/species_list_creation/HBW_colombia_overlaps.Rdata")
# Get the list of species
colombia_species <- unique(recast_range_maps$SCINAME[unlist(lapply(colombia_overlaps, FUN = function(i){return(length(i) != 0)}))])
save(colombia_species, file = "Data/Birds/species_list_creation/colombia_species.Rdata")



load("Data/Birds/species_list_creation/colombia_species.Rdata")

# Read in the HBW/eBird taxonomic interconversion that Marshall Iliff prepared
taxonomy <- read.csv("Data/Birds/species_list_creation/HBW_eBird_taxonomy.csv", stringsAsFactors = F)
t2 <- taxonomy[taxonomy$HBW_CAT == "sp" | taxonomy$CLEM_CAT_2019 == "sp", ] # remove taxa that are not treated as species by either eBird or HBW

# Make a few changes in the HBW taxonomy that were not yet current in the file version provided by Marshall
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

colombia_species[colombia_species %ni% t2$HBW_LATIN]
# The ebird/HBW interconversion now contains all species in with range-maps overlapping the buffered Colombia polygon

csp2 <- t2[t2$HBW_LATIN %in% colombia_species, ]
# Extract just the part of the interconversion file that is relevant to Colombia


elevation_dataset <- read.csv("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/species_list_creation/Bird_elevations_initial.csv", stringsAsFactors = F)
# pull in the elevation data based on the Donegan checklist. These data are described in 
# /Users/jacobsocolar/Dropbox/Work/Colombia/Writing/Bird_elevation methods.docx, except that they do not yet
# have the various HBW splits implemented.  (Those are in Bird_elevations_final.csv)

species_to_check <- csp2$HBW_LATIN[csp2$HBW_LATIN %ni% elevation_dataset$Scientific & 
                                     !((csp2$CLEM_SCI_2019 %in% elevation_dataset$Scientific) & csp2$HBW.Clem == "sp-sp")]
# These are the species from the range-maps that are not included in the compiled elevation dataset.
# These species either do not occur in the Donegan checklist (despite being mapped within 100 km of Colombia
# by BirdLife, or are species with taxonomic or nomenclatural variation between HBW and Donegan)

write.csv(species_to_check, file = '/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/species_list_creation/species_to_check.csv')
# This is an initial save of the species to check, which is the basis for a manually edited file, saved as 
# "species_to_check_checked.csv".  
# Each species was hand-checked against BirdLife maps and various taxonomic sources to assign it to one of four categories:
# 1) Not Colombia: species that do not occur in Colombia, following the Donegan checklist
# 2) Synonym: simple nomenclatural difference (e.g. different authorties place in different genera)
# 3) Clean split/lump: HBW splits or lumps species in Donegan, but in such a way that all colombian populations belong to 
#    a single species in both HBW and Donegan
# 4) Messy split/lump: HBW splits or lumps species in Donegan, in such a way that one source treats as heterospecific two 
#    populations, both occurring in Colombia, that the other source treats as conspecific
# Additionally, there is a column "really.messy" which contains a 1 if multiple species (following one authority or the 
# other) from a messy split/lump are both represented in our field-collected data. And a column "mapped.2.border"
# contains a 1 if a species not on the Donegan list is nevertheless mapped to the immediate vicinity of the Colombian border,
# excluding species that are confined to the south bank of the Amazon or to Venezuelan highlands.

# This next block of code was used to identify names that are present in Donegan (and relevant to our analysis)
# but absent from HBW/BirdLife. Those issues were added manually at lines 283-end of species_to_check_checked.csv
elevation_species <- elevation_dataset$Scientific[is.na(elevation_dataset$coastal_pelagic) & 
                                                    is.na(elevation_dataset$outside_birdlife) &
                                                    is.na(elevation_dataset$vagrant)]
# Extract species from the Donegan data that are relevant to analysis (not vagrant or strictly coastal)

s2c2 <- elevation_species[(elevation_species %ni% csp2$HBW_LATIN) & (elevation_species %ni% csp2$CLEM_SCI_2019[csp2$HBW.Clem == 'sp-sp'])]
# Get species absent from HBW that don't correspond to species-species synonymy in Marshall's file

changes <- read.csv('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/species_list_creation/species_to_check_checked.csv', stringsAsFactors = F)
s2c2[(s2c2 %ni% changes$synonym) & (s2c2 %ni% changes$clean.split.lump) & (s2c2 %ni% changes$messy.split.lump)]
# Check for extra Donegan species that aren't accounted for in species_to_check_checked.csv
# All changes accounted for! (This is because changes not accounted for were added manually to species_to_check_checked.csv)
# The species that still show up here are parts of complicated splits/lumps that are contained in species_to_check_checked.csv

# Because species_to_check_checked.csv covers everything, we can be confident that there's no serious error of Colombia-wide
# omission in the BirdLife maps, since I had checked every single one by hand.

# Get a list of species of interest for analysis, but starting with colombia_species and:
# 1) removing species marked in species_to_check_checked.csv as not occurring in Colombia
HBW <- as.character(colombia_species[colombia_species %ni% changes$hbw.name[changes$not.colombia == 1]])

# 2) removing species that are flagged as vagrant or coastal in bird_elevations_initial.csv
for(i in length(HBW):1){
  if(HBW[i] %in% elevation_dataset$Scientific){
    if(!is.na(elevation_dataset$outside_birdlife[elevation_dataset$Scientific == HBW[i]])){
      print(paste("outside", HBW[i]))
    }
   if(!is.na(elevation_dataset$coastal_pelagic[elevation_dataset$Scientific == HBW[i]])){
      print(HBW[i])
      HBW <- HBW[-i]
   }
    else if(!is.na(elevation_dataset$vagrant[elevation_dataset$Scientific == HBW[i]])){
      print(HBW[i])
      HBW <- HBW[-i]
    }
  }
  else if((t2$HBW.Clem[t2$HBW_LATIN == HBW[i]] == 'sp-sp') & 
            (t2$CLEM_SCI_2019[t2$HBW_LATIN == HBW[i]] %in% elevation_dataset$Scientific)){
    if(!is.na(elevation_dataset$coastal_pelagic[elevation_dataset$Scientific == t2$CLEM_SCI_2019[t2$HBW_LATIN == HBW[i]]])){
        print(HBW[i])
        HBW <- HBW[-i]
    }
  }
}

HBW[HBW %in% changes$not.colombia]
HBW[HBW %in% changes$synonym]
HBW[HBW %in% changes$clean.split.lump]

# Create a dataframe to hold our analyzeable species set
initial_species_list <- data.frame(HBW = HBW, Donegan = NA, Donegan2 = NA, eBird = NA, eBird2 = NA, stringsAsFactors = F)
for(i in 1:nrow(initial_species_list)){
  hbw <- initial_species_list$HBW[i]
  if(hbw %in% elevation_species){initial_species_list$Donegan[i] <- hbw}
  else if(hbw %in% changes$hbw.name){
    if(changes$synonym[changes$hbw.name == hbw] %in% elevation_species){
      initial_species_list$Donegan[i] <- changes$synonym[changes$hbw.name == hbw]
    }
    else if(changes$clean.split.lump[changes$hbw.name == hbw] %in% elevation_species){
      initial_species_list$Donegan[i] <- changes$clean.split.lump[changes$hbw.name == hbw]
    }
  }
  else if((t2$HBW.Clem[t2$HBW_LATIN == hbw] == "sp-sp") & (t2$CLEM_SCI_2019[t2$HBW_LATIN == hbw] %in% elevation_species)){
    initial_species_list$Donegan[i] <- t2$CLEM_SCI_2019[(t2$HBW_LATIN == hbw) & (t2$HBW.Clem == 'sp-sp')]
  }
  else{print(hbw)}
}

for(i in 1:288){
  if(changes$messy.split.lump[i] %in% elevation_dataset$Scientific){
    initial_species_list$Donegan[which(initial_species_list$HBW == changes$hbw.name[i])] <-
      changes$messy.split.lump[i]
  }
}

for(i in 1:nrow(initial_species_list)){
  donegan <- initial_species_list$Donegan[i]
  if(donegan %in% t2$CLEM_SCI_2019){
    initial_species_list$eBird[i] <- donegan
  }else{print(donegan)}
}

initial_species_list$Donegan[initial_species_list$HBW == "Megascops vermiculatus"] <- "Megascops centralis"
initial_species_list$Donegan2[initial_species_list$HBW == "Megascops vermiculatus"] <- "Megascops roraimae"
initial_species_list$Donegan[initial_species_list$HBW == "Trogon violaceus"] <- "Trogon ramonianus"
initial_species_list$Donegan2[initial_species_list$HBW == "Trogon violaceus"] <- "Trogon caligatus"
initial_species_list$Donegan2[initial_species_list$HBW == "Corapipo leucorrhoa"] <- "Corapipo altera"
initial_species_list$Donegan2[initial_species_list$HBW == "Zimmerius chrysops"] <- "Zimmerius minimus"
initial_species_list$Donegan2[initial_species_list$HBW == "Butorides striatus"] <- "Butorides virescens"
initial_species_list$Donegan2[initial_species_list$HBW == "Chaetura chapmani"] <- "Chaetura viridipennis"
initial_species_list$Donegan2[initial_species_list$HBW == "Chlorostilbon mellisugus"] <- "Chlorostilbon melanorhynchus"
initial_species_list$Donegan2[initial_species_list$HBW == "Turdus albicollis"] <- "Turdus assimilis"
initial_species_list$eBird[initial_species_list$HBW == "Megascops vermiculatus"] <- "Megascops centralis"
initial_species_list$eBird2[initial_species_list$HBW == "Megascops vermiculatus"] <- "Megascops roraimae"
initial_species_list$eBird[initial_species_list$HBW == "Trogon violaceus"] <- "Trogon ramonianus"
initial_species_list$eBird2[initial_species_list$HBW == "Trogon violaceus"] <- "Trogon caligatus"
initial_species_list$eBird2[initial_species_list$HBW == "Corapipo leucorrhoa"] <- "Corapipo altera"
initial_species_list$eBird2[initial_species_list$HBW == "Zimmerius chrysops"] <- "Zimmerius minimus"
initial_species_list$eBird2[initial_species_list$HBW == "Butorides striatus"] <- "Butorides virescens"
initial_species_list$eBird2[initial_species_list$HBW == "Chaetura chapmani"] <- "Chaetura viridipennis"
initial_species_list$eBird2[initial_species_list$HBW == "Chlorostilbon mellisugus"] <- "Chlorostilbon melanorhynchus"
initial_species_list$eBird2[initial_species_list$HBW == "Turdus albicollis"] <- "Turdus assimilis"
initial_species_list$eBird[initial_species_list$HBW == "Paraclaravis mondetoura"] <- "Paraclaravis mondetoura"
initial_species_list$eBird[initial_species_list$HBW == "Juliamyia julie"] <- "Juliamyia julie"
initial_species_list$eBird[initial_species_list$HBW == "Hapalocrex flaviventer"] <- "Hapalocrex flaviventer"
initial_species_list$eBird[initial_species_list$HBW == "Leuconotopicus fumigatus"] <- "Dryobates fumigatus"
initial_species_list$eBird[initial_species_list$HBW == "Leistes militaris"] <- "Leistes militaris"
initial_species_list$eBird[initial_species_list$HBW == "Leistes bellicosus"] <- "Leistes bellicosus"
initial_species_list$eBird[initial_species_list$HBW == "Anisognathus somptuosus"] <- "Anisognathus somptuosus"
initial_species_list$eBird[initial_species_list$HBW == "Anisognathus notabilis"] <- "Anisognathus notabilis"
initial_species_list$eBird[initial_species_list$HBW == "Geospizopsis unicolor"] <- "Geospizopsis unicolor"
initial_species_list$eBird[initial_species_list$HBW == "Spodiornis rusticus"] <- "Spodiornis rusticus"
initial_species_list$eBird[initial_species_list$HBW == "Epinecrophylla haematonota"] <- "Epinecrophylla haematonota"
initial_species_list$eBird[initial_species_list$HBW == "Sclerurus mexicanus"] <- "Sclerurus mexicanus"
initial_species_list$eBird[initial_species_list$HBW == "Manacus vitellinus"] <- "Manacus vitellinus"
initial_species_list$eBird[initial_species_list$HBW == "Pyrrhura chapmani"] <- "Pyrrhura melanura"
initial_species_list$eBird[initial_species_list$HBW == "Pyrrhura pacifica"] <- "Pyrrhura melanura"
initial_species_list$eBird[initial_species_list$HBW == "Coeligena conradii"] <- "Coeligena torquata"
initial_species_list$eBird[initial_species_list$HBW == "Forpus spengeli"] <- "Forpus xanthopterygius"
initial_species_list$eBird[initial_species_list$HBW == "Heliangelus clarisse"] <- "Heliangelus amethysticollis"
initial_species_list$eBird[initial_species_list$HBW == "Hypnelus bicinctus"] <- "Hypnelus ruficollis"
initial_species_list$eBird[initial_species_list$HBW == "Veniliornis affinis"] <- "Dryobates affinis"
initial_species_list$eBird[initial_species_list$HBW == "Veniliornis callonotus"] <- "Dryobates callonotus"
initial_species_list$eBird[initial_species_list$HBW == "Veniliornis chocoensis"] <- "Dryobates chocoensis"
initial_species_list$eBird[initial_species_list$HBW == "Veniliornis dignus"] <- "Dryobates dignus"
initial_species_list$eBird[initial_species_list$HBW == "Atlapetes nigrifrons"] <- "Atlapetes latinuchus"
initial_species_list$eBird[initial_species_list$HBW == "Veniliornis kirkii"] <- "Dryobates kirkii"
initial_species_list$eBird[initial_species_list$HBW == "Veniliornis nigriceps"] <- "Dryobates nigriceps"
initial_species_list$eBird[initial_species_list$HBW == "Veniliornis passerinus"] <- "Dryobates passerinus"
initial_species_list$eBird[initial_species_list$HBW == "Atlapetes crassus"] <- "Atlapetes tricolor"
initial_species_list$eBird[initial_species_list$HBW == "Automolus virgatus"] <- "Automolus subulatus"
initial_species_list$eBird[initial_species_list$HBW == "Pachyramphus xanthogenys"] <- "Pachyramphus viridis"
initial_species_list$eBird[initial_species_list$HBW == "Amazona bodini"] <- "Amazona festiva"
initial_species_list$eBird[initial_species_list$HBW == "Calliphlox mitchellii"] <- "Philodice mitchellii"
initial_species_list$eBird[initial_species_list$HBW == "Campephilus splendens"] <- "Campephilus haematogaster"
initial_species_list$eBird[initial_species_list$HBW == "Tephrophilus wetmorei"] <- "Buthraupis wetmorei"
initial_species_list$eBird[initial_species_list$HBW == "Islerothraupis cristata"] <- "Tachyphonus cristatus"
initial_species_list$eBird[initial_species_list$HBW == "Islerothraupis luctuosa"] <- "Tachyphonus luctuosus"
initial_species_list$eBird[initial_species_list$HBW == "Sporathraupis cyanocephala"] <- "Thraupis cyanocephala"
initial_species_list$eBird[initial_species_list$HBW == "Grallaria fenwickorum"] <- "Grallaria urraoensis"
initial_species_list$HBW[is.na(initial_species_list$Donegan)]
initial_species_list$HBW[is.na(initial_species_list$eBird)]



# Now match to Pulido tree names
pulido <- ape::read.tree(file = "Data/Birds/phylogeny/Pulido_phylogeny/JETZ TREES/All_birds_MaxCladeCredTree.txt")
pulido.spp <- gsub("_", " ", pulido$tip.label)

initial_species_list$HBW[(initial_species_list$HBW %ni% pulido.spp) &
                           (initial_species_list$Donegan %ni% pulido.spp) &
                           (initial_species_list$eBird %ni% pulido.spp)]

initial_species_list$HBW[(initial_species_list$HBW %ni% pulido.spp) &
                           (initial_species_list$Donegan %ni% pulido.spp) &
                           (initial_species_list$eBird %ni% pulido.spp)]

initial_species_list$Pulido <- NA
for(i in 1:nrow(initial_species_list)){
  if(initial_species_list$HBW[i] %in% pulido.spp){
    initial_species_list$Pulido[i] <- initial_species_list$HBW[i]
  }else if((initial_species_list$Donegan[i] %in% pulido.spp) & (initial_species_list$Donegan2[i] %ni% pulido.spp)){
    initial_species_list$Pulido[i] <- initial_species_list$Donegan[i]
  }else if((initial_species_list$eBird[i] %in% pulido.spp) & (initial_species_list$eBird2[i] %ni% pulido.spp)){
    initial_species_list$Pulido[i] <- initial_species_list$eBird[i]
  }
}

initial_species_list$Pulido[initial_species_list$HBW == "Aramides cajaneus"] <- "Aramides cajanea"
initial_species_list$Pulido[initial_species_list$HBW == "Ardea alba"] <- "Casmerodius albus"
initial_species_list$Pulido[initial_species_list$HBW == "Spatula clypeata"] <- "Anas clypeata"
initial_species_list$Pulido[initial_species_list$HBW == "Machaeropterus striolatus"] <- "Machaeropterus regulus"
initial_species_list$Pulido[initial_species_list$HBW == "Porphyriops melanops"] <- "Gallinula melanops"
initial_species_list$Pulido[initial_species_list$HBW == "Chordeiles nacunda"] <- "Podager nacunda"
initial_species_list$Pulido[initial_species_list$HBW == "Ciccaba albitarsis"] <- "Strix albitarsis"
initial_species_list$Pulido[initial_species_list$HBW == "Ciccaba huhula"] <- "Strix huhula"
initial_species_list$Pulido[initial_species_list$HBW == "Ciccaba nigrolineata"] <- "Strix nigrolineata"
initial_species_list$Pulido[initial_species_list$HBW == "Ciccaba virgata"] <- "Strix virgata"
initial_species_list$Pulido[initial_species_list$HBW == "Circus hudsonius"] <- "Circus cyaneus"
initial_species_list$Pulido[initial_species_list$HBW == "Cryptoleucopteryx plumbea"] <- "Leucopternis plumbeus"
initial_species_list$Pulido[initial_species_list$HBW == "Eriocnemis aline"] <- "Eriocnemis alinae"
initial_species_list$Pulido[initial_species_list$HBW == "Oxypogon stuebelii"] <- "Oxypogon guerinii"
initial_species_list$Pulido[initial_species_list$HBW == "Gallinago delicata"] <- "Gallinago gallinago"
initial_species_list$Pulido[initial_species_list$HBW == "Gallinula galeata"] <- "Gallinula chloropus"
initial_species_list$Pulido[initial_species_list$HBW == "Geranoaetus polyosoma"] <- "Buteo polyosoma"
initial_species_list$Pulido[initial_species_list$HBW == "Geotrygon purpurata"] <- "Geotrygon saphirina"
initial_species_list$Pulido[initial_species_list$HBW == "Geranoaetus albicaudatus"] <- "Buteo albicaudatus"
initial_species_list$Pulido[initial_species_list$HBW == "Spinus magellanicus"] <- "Carduelis magellanica"
initial_species_list$Pulido[initial_species_list$HBW == "Hydropsalis cayennensis"] <- "Caprimulgus cayennensis"
initial_species_list$Pulido[initial_species_list$HBW == "Hydropsalis maculicaudus"] <- "Caprimulgus maculicaudus"
initial_species_list$Pulido[initial_species_list$HBW == "Uromyias agilis"] <- "Anairetes agilis"
initial_species_list$Pulido[initial_species_list$HBW == "Myiothlypis basilica"] <- "Basileuterus basilicus"
initial_species_list$Pulido[initial_species_list$HBW == "Myiothlypis cinereicollis"] <- "Basileuterus cinereicollis"
initial_species_list$Pulido[initial_species_list$HBW == "Myiothlypis conspicillata"] <- "Basileuterus conspicillatus"
initial_species_list$Pulido[initial_species_list$HBW == "Myiothlypis coronata"] <- "Basileuterus coronatus"
initial_species_list$Pulido[initial_species_list$HBW == "Myiothlypis nigrocristata"] <- "Basileuterus nigrocristatus"
initial_species_list$Pulido[initial_species_list$HBW == "Mareca americana"] <- "Anas americana"
initial_species_list$Pulido[initial_species_list$HBW == "Oxypogon cyanolaemus"] <- "Oxypogon guerinii"
initial_species_list$Pulido[initial_species_list$HBW == "Conirostrum binghami"] <- "Oreomanes fraseri"
initial_species_list$Pulido[initial_species_list$HBW == "Cercomacroides tyrannina"] <- "Cercomacra tyrannina"
initial_species_list$Pulido[initial_species_list$HBW == "Cyanoloxia brissonii"] <- "Cyanocompsa brissonii"
initial_species_list$Pulido[initial_species_list$HBW == "Setophaga fusca"] <- "Dendroica fusca"
initial_species_list$Pulido[initial_species_list$HBW == "Setophaga striata"] <- "Dendroica striata"
initial_species_list$Pulido[initial_species_list$HBW == "Setophaga pensylvanica"] <- "Dendroica pensylvanica"
initial_species_list$Pulido[initial_species_list$HBW == "Cyanoloxia cyanoides"] <- "Cyanocompsa cyanoides"
initial_species_list$Pulido[initial_species_list$HBW == "Cyanoloxia rothschildii"] <- "Cyanocompsa cyanoides"
initial_species_list$Pulido[initial_species_list$HBW == "Sphenopsis frontalis"] <- "Hemispingus frontalis"
initial_species_list$Pulido[initial_species_list$HBW == "Thlypopsis superciliaris"] <- "Hemispingus superciliaris"
initial_species_list$Pulido[initial_species_list$HBW == "Pachysylvia hypoxantha"] <- "Hylophilus hypoxanthus"
initial_species_list$Pulido[initial_species_list$HBW == "Leptotila cassinii"] <- "Leptotila cassini"
initial_species_list$Pulido[initial_species_list$HBW == "Arremon assimilis"] <- "Arremon torquatus"
initial_species_list$Pulido[initial_species_list$HBW == "Leptotrygon veraguensis"] <- "Geotrygon veraguensis"
initial_species_list$Pulido[initial_species_list$HBW == "Orochelidon flavipes"] <- "Notiochelidon flavipes"
initial_species_list$Pulido[initial_species_list$HBW == "Myrmophylax atrothorax"] <- "Myrmeciza atrothorax"
initial_species_list$Pulido[initial_species_list$HBW == "Sipia berlepschi"] <- "Myrmeciza berlepschi"
initial_species_list$Pulido[initial_species_list$HBW == "Akletos melanoceps"] <- "Myrmeciza melanoceps"
initial_species_list$Pulido[initial_species_list$HBW == "Hafferia fortis"] <- "Myrmeciza fortis"
initial_species_list$Pulido[initial_species_list$HBW == "Sipia palliata"] <- "Myrmeciza laemosticta"
initial_species_list$Pulido[initial_species_list$HBW == "Myrmelastes hyperythrus"] <- "Myrmeciza hyperythra"
initial_species_list$Pulido[initial_species_list$HBW == "Isleria hauxwelli"] <- "Myrmotherula hauxwelli"
initial_species_list$Pulido[initial_species_list$HBW == "Cacicus latirostris"] <- "Ocyalus latirostris"
initial_species_list$Pulido[initial_species_list$HBW == "Orochelidon murina"] <- "Notiochelidon murina"
initial_species_list$Pulido[initial_species_list$HBW == "Geothlypis formosa"] <- "Oporornis formosus"
initial_species_list$Pulido[initial_species_list$HBW == "Sporophila crassirostris"] <- "Oryzoborus crassirostris"
initial_species_list$Pulido[initial_species_list$HBW == "Myiothlypis fulvicauda"] <- "Phaeothlypis fulvicauda"
initial_species_list$Pulido[initial_species_list$HBW == "Phelpsia inornata"] <- "Phelpsia inornatus"
initial_species_list$Pulido[initial_species_list$HBW == "Ceratopipra erythrocephala"] <- "Pipra erythrocephala"
initial_species_list$Pulido[initial_species_list$HBW == "Ceratopipra mentalis"] <- "Pipra mentalis"
initial_species_list$Pulido[initial_species_list$HBW == "Pseudopipra pipra"] <- "Pipra pipra"
initial_species_list$Pulido[initial_species_list$HBW == "Myrmelastes schistaceus"] <- "Schistocichla schistacea"
initial_species_list$Pulido[initial_species_list$HBW == "Scytalopus opacus"] <- "Scytalopus canus"
initial_species_list$Pulido[initial_species_list$HBW == "Mazaria propinqua"] <- "Synallaxis propinqua"
initial_species_list$Pulido[initial_species_list$HBW == "Euchrepomis spodioptila"] <- "Terenura spodioptila"
initial_species_list$Pulido[initial_species_list$HBW == "Cantorchilus leucopogon"] <- "Thryothorus leucopogon"
initial_species_list$Pulido[initial_species_list$HBW == "Asthenes perijana"] <- "Schizoeaca perijana"
initial_species_list$Pulido[initial_species_list$HBW == "Cantorchilus nigricapillus"] <- "Thryothorus nigricapillus"
initial_species_list$Pulido[initial_species_list$HBW == "Pheugopedius rutilus"] <- "Thryothorus rutilus"
initial_species_list$Pulido[initial_species_list$HBW == "Leiothlypis peregrina"] <- "Vermivora peregrina"
initial_species_list$Pulido[initial_species_list$HBW == "Chloropipo flavicapilla"] <- "Xenopipo flavicapilla"
initial_species_list$Pulido[initial_species_list$HBW == "Megascops vermiculatus"] <- "Megascops guatemalae"
initial_species_list$Pulido[initial_species_list$HBW == "Systellura longirostris"] <- "Caprimulgus longirostris"
initial_species_list$Pulido[initial_species_list$HBW == "Momotus subrufescens"] <- "Momotus momota"
initial_species_list$Pulido[initial_species_list$HBW == "Asemospiza obscura"] <- "Tiaris obscurus"
initial_species_list$Pulido[initial_species_list$HBW == "Morphnarchus princeps"] <- "Leucopternis princeps"
initial_species_list$Pulido[initial_species_list$HBW == "Nyctipolus nigrescens"] <- "Caprimulgus nigrescens"
initial_species_list$Pulido[initial_species_list$HBW == "Ortalis columbiana"] <- "Ortalis guttata"
initial_species_list$Pulido[initial_species_list$HBW == "Arremon perijanus"] <- "Arremon torquatus"
initial_species_list$Pulido[initial_species_list$HBW == "Pseudastur albicollis"] <- "Leucopternis albicollis"
initial_species_list$Pulido[initial_species_list$HBW == "Psittacara leucophthalmus"] <- "Aratinga leucophthalma"
initial_species_list$Pulido[initial_species_list$HBW == "Psittacara wagleri"] <- "Aratinga wagleri"
initial_species_list$Pulido[initial_species_list$HBW == "Schistes albogularis"] <- "Schistes geoffroyi"
initial_species_list$Pulido[initial_species_list$HBW == "Setopagis heterura"] <- "Caprimulgus heterurus"
initial_species_list$Pulido[initial_species_list$HBW == "Spatula cyanoptera"] <- "Anas cyanoptera"
initial_species_list$Pulido[initial_species_list$HBW == "Spatula discors"] <- "Anas discors"
initial_species_list$Pulido[initial_species_list$HBW == "Sternula superciliaris"] <- "Sterna superciliaris"
initial_species_list$Pulido[initial_species_list$HBW == "Urochroa leucura"] <- "Urochroa bougueri"
initial_species_list$Pulido[initial_species_list$HBW == "Zentrygon frenata"] <- "Geotrygon frenata"
initial_species_list$Pulido[initial_species_list$HBW == "Zentrygon goldmani"] <- "Geotrygon goldmani"
initial_species_list$Pulido[initial_species_list$HBW == "Zentrygon linearis"] <- "Geotrygon linearis"
initial_species_list$Pulido[initial_species_list$HBW == "Arremon basilicus"] <- "Arremon torquatus"
initial_species_list$Pulido[initial_species_list$HBW == "Cercomacroides fuscicauda"] <- "Cercomacra nigrescens"
initial_species_list$Pulido[initial_species_list$HBW == "Drymophila striaticeps"] <- "Drymophila caudata"
initial_species_list$Pulido[initial_species_list$HBW == "Elaenia brachyptera"] <- "Elaenia chiriquensis"
initial_species_list$Pulido[initial_species_list$HBW == "Gymnopithys bicolor"] <- "Gymnopithys leucaspis"
initial_species_list$Pulido[initial_species_list$HBW == "Henicorhina anachoreta"] <- "Henicorhina leucophrys"
initial_species_list$Pulido[initial_species_list$HBW == "Drymophila hellmayri"] <- "Drymophila caudata"
initial_species_list$Pulido[initial_species_list$HBW == "Automolus virgatus"] <- "Hyloctistes subulatus"
initial_species_list$Pulido[initial_species_list$HBW == "Automolus subulatus"] <- "Hyloctistes subulatus"
initial_species_list$Pulido[initial_species_list$HBW == "Tunchiornis ochraceiceps"] <- "Hylophilus ochraceiceps"
initial_species_list$Pulido[initial_species_list$HBW == "Lepidocolaptes duidae"] <- "Lepidocolaptes albolineatus"
initial_species_list$Pulido[initial_species_list$HBW == "Aprositornis disjuncta"] <- "Myrmeciza disjuncta"
initial_species_list$Pulido[initial_species_list$HBW == "Hafferia zeledoni"] <- "Myrmeciza immaculata"
initial_species_list$Pulido[initial_species_list$HBW == "Sciaphylax castanea"] <- "Myrmeciza castanea"
initial_species_list$Pulido[initial_species_list$HBW == "Hafferia immaculata"] <- "Myrmeciza immaculata"
initial_species_list$Pulido[initial_species_list$HBW == "Ammonastes pelzelni"] <- "Myrmeciza pelzelni"
initial_species_list$Pulido[initial_species_list$HBW == "Nephelomyias pulcher"] <- "Myiophobus pulcher"
initial_species_list$Pulido[initial_species_list$HBW == "Premnornis guttuliger"] <- "Premnornis guttuligera"
initial_species_list$Pulido[initial_species_list$HBW == "Pseudocolaptes johnsoni"] <- "Pseudocolaptes lawrencii johnsoni"
initial_species_list$Pulido[initial_species_list$HBW == "Pheugopedius columbianus"] <- "Thryothorus sclateri"
initial_species_list$Pulido[initial_species_list$HBW == "Thamnophilus melanonotus"] <- "Sakesphorus melanonotus"
initial_species_list$Pulido[initial_species_list$HBW == "Schiffornis aenea"] <- "Schiffornis turdina"
initial_species_list$Pulido[initial_species_list$HBW == "Schiffornis stenorhyncha"] <- "Schiffornis turdina"
initial_species_list$Pulido[initial_species_list$HBW == "Schiffornis veraepacis"] <- "Schiffornis turdina"
initial_species_list$Pulido[initial_species_list$HBW == "Sirystes albocinereus"] <- "Sirystes sibilator"
initial_species_list$Pulido[initial_species_list$HBW == "Sporophila funerea"] <- "Oryzoborus funereus"
initial_species_list$Pulido[initial_species_list$HBW == "Pheugopedius euophrys"] <- "Thryothorus euophrys"
initial_species_list$Pulido[initial_species_list$HBW == "Thryophilus nicefori"] <- "Thryothorus nicefori"
initial_species_list$Pulido[initial_species_list$HBW == "Thryophilus sernai"] <- "Thryothorus rufalbus"
initial_species_list$Pulido[initial_species_list$HBW == "Turdus sanchezorum"] <- "Turdus haplochrous"
initial_species_list$Pulido[initial_species_list$HBW == "Cryptopipo holochlora"] <- "Xenopipo holochlora"
initial_species_list$Pulido[initial_species_list$HBW == "Cryptopipo litae"] <- "Xenopipo holochlora"
initial_species_list$Pulido[initial_species_list$HBW == "Xenops genibarbis"] <- "Xenops minutus remoratus"
initial_species_list$Pulido[initial_species_list$HBW == "Myiothlypis chlorophrys"] <- "Basileuterus chlorophrys"
initial_species_list$Pulido[initial_species_list$HBW == "Myiothlypis luteoviridis"] <- "Basileuterus luteoviridis"
initial_species_list$Pulido[initial_species_list$HBW == "Kleinothraupis atropileus"] <- "Hemispingus atropileus"
initial_species_list$Pulido[initial_species_list$HBW == "Tangara argentea"] <- "Tangara cyanoptera"
initial_species_list$Pulido[initial_species_list$HBW == "Sphenopsis melanotis"] <- "Hemispingus melanotis"
initial_species_list$Pulido[initial_species_list$HBW == "Setophaga tigrina"] <- "Dendroica tigrina"
initial_species_list$Pulido[initial_species_list$HBW == "Pachysylvia aurantiifrons"] <- "Hylophilus aurantiifrons"
initial_species_list$Pulido[initial_species_list$HBW == "Pachysylvia decurtata"] <- "Hylophilus decurtatus"
initial_species_list$Pulido[initial_species_list$HBW == "Sipia nigricauda"] <- "Myrmeciza nigricauda"
initial_species_list$Pulido[initial_species_list$HBW == "Geothlypis philadelphia"] <- "Oporornis philadelphia"
initial_species_list$Pulido[initial_species_list$HBW == "Cardellina canadensis"] <- "Wilsonia canadensis"
initial_species_list$Pulido[initial_species_list$HBW == "Orthopsittaca manilatus"] <- "Orthopsittaca manilata"
initial_species_list$Pulido[initial_species_list$HBW == "Parabuteo leucorrhous"] <- "Buteo leucorrhous"
initial_species_list$Pulido[initial_species_list$HBW == "Rupornis magnirostris"] <- "Buteo magnirostris"
initial_species_list$Pulido[initial_species_list$HBW == "Poliocrania maculifer"] <- "Myrmeciza exsul"
initial_species_list$Pulido[initial_species_list$HBW == "Anabacerthia ruficaudata"] <- "Philydor ruficaudatum"
initial_species_list$Pulido[initial_species_list$HBW == "Euchrepomis callinota"] <- "Terenura callinota"
initial_species_list$Pulido[initial_species_list$HBW == "Setophaga petechia"] <- "Dendroica petechia"
initial_species_list$Pulido[initial_species_list$HBW == "Setophaga coronata"] <- "Dendroica coronata"
initial_species_list$Pulido[initial_species_list$HBW == "Dendrocolaptes certhia"] <- "Dendrocolaptes certhia certhia"
initial_species_list$Pulido[initial_species_list$HBW == "Frederickena unduliger"] <- "Frederickena unduligera"
initial_species_list$Pulido[initial_species_list$HBW == "Parkesia motacilla"] <- "Seiurus motacilla"
initial_species_list$Pulido[initial_species_list$HBW == "Parkesia noveboracensis"] <- "Seiurus noveboracensis"
initial_species_list$Pulido[initial_species_list$HBW == "Calidris subruficollis"] <- "Tryngites subruficollis"
initial_species_list$Pulido[initial_species_list$HBW == "Drymophila klagesi"] <- "Drymophila caudata"
initial_species_list$Pulido[initial_species_list$HBW == "Aglaiocercus kingii"] <- "Aglaiocercus kingi"
initial_species_list$Pulido[initial_species_list$HBW == "Amazilia saucerottei"] <- "Amazilia saucerrottei"
initial_species_list$Pulido[initial_species_list$HBW == "Amazona mercenarius"] <- "Amazona mercenaria"
initial_species_list$Pulido[initial_species_list$HBW == "Antrostomus carolinensis"] <- "Caprimulgus carolinensis"
initial_species_list$Pulido[initial_species_list$HBW == "Asio clamator"] <- "Pseudoscops clamator"
initial_species_list$Pulido[initial_species_list$HBW == "Buteogallus schistaceus"] <- "Leucopternis schistaceus"
initial_species_list$Pulido[initial_species_list$HBW == "Buteogallus solitarius"] <- "Harpyhaliaetus solitarius"
initial_species_list$Pulido[initial_species_list$HBW == "Asthenes fuliginosa"] <- "Schizoeaca fuliginosa"
initial_species_list$Pulido[initial_species_list$HBW == "Clibanornis rubiginosus"] <- "Automolus rubiginosus nigricauda"
initial_species_list$Pulido[initial_species_list$HBW == "Cercomacroides nigrescens"] <- "Cercomacra nigrescens"
initial_species_list$Pulido[initial_species_list$HBW == "Cercomacroides parkeri"] <- "Cercomacra parkeri"
initial_species_list$Pulido[initial_species_list$HBW == "Cercomacroides serva"] <- "Cercomacra serva"
initial_species_list$Pulido[initial_species_list$HBW == "Drymotoxeres pucheranii"] <- "Drymotoxeres pucherani"
initial_species_list$Pulido[initial_species_list$HBW == "Pachysylvia semibrunnea"] <- "Hylophilus semibrunneus"
initial_species_list$Pulido[initial_species_list$HBW == "Scytalopus perijanus"] <- "Scytalopus meridanus"
initial_species_list$Pulido[initial_species_list$HBW == "Pygochelidon melanoleuca"] <- "Atticora melanoleuca"
initial_species_list$Pulido[initial_species_list$HBW == "Atticora tibialis"] <- "Neochelidon tibialis"
initial_species_list$Pulido[initial_species_list$HBW == "Myrmelastes leucostigma"] <- "Schistocichla leucostigma"
initial_species_list$Pulido[initial_species_list$HBW == "Antrostomus rufus"] <- "Caprimulgus rufus"
initial_species_list$Pulido[initial_species_list$HBW == "Myiothlypis flaveola"] <- "Basileuterus flaveolus"
initial_species_list$Pulido[initial_species_list$HBW == "Cnemathraupis eximia"] <- "Buthraupis eximia"
initial_species_list$Pulido[initial_species_list$HBW == "Spinus spinescens"] <- "Carduelis spinescens"
initial_species_list$Pulido[initial_species_list$HBW == "Spinus xanthogastrus"] <- "Carduelis xanthogastra"
initial_species_list$Pulido[initial_species_list$HBW == "Setophaga virens"] <- "Dendroica virens"
initial_species_list$Pulido[initial_species_list$HBW == "Chlorospingus flavopectus"] <- "Chlorospingus ophthalmicus"
initial_species_list$Pulido[initial_species_list$HBW == "Cacicus oseryi"] <- "Clypicterus oseryi"
initial_species_list$Pulido[initial_species_list$HBW == "Setophaga castanea"] <- "Dendroica castanea"
initial_species_list$Pulido[initial_species_list$HBW == "Sporophila fringilloides"] <- "Dolospingus fringilloides"
initial_species_list$Pulido[initial_species_list$HBW == "Pseudospingus verticalis"] <- "Hemispingus verticalis"
initial_species_list$Pulido[initial_species_list$HBW == "Pheugopedius spadix"] <- "Thryothorus spadix"
initial_species_list$Pulido[initial_species_list$HBW == "Melanospiza bicolor"] <- "Tiaris bicolor"
initial_species_list$Pulido[initial_species_list$HBW == "Pheugopedius fasciatoventris"] <- "Thryothorus fasciatoventris"
initial_species_list$Pulido[initial_species_list$HBW == "Cantorchilus leucotis"] <- "Thryothorus leucotis"
initial_species_list$Pulido[initial_species_list$HBW == "Pheugopedius mystacalis"] <- "Thryothorus mystacalis"
initial_species_list$Pulido[initial_species_list$HBW == "Thryophilus rufalbus"] <- "Thryothorus rufalbus"
initial_species_list$Pulido[initial_species_list$HBW == "Asemospiza fuliginosa"] <- "Tiaris fuliginosus"
initial_species_list$Pulido[initial_species_list$HBW == "Setophaga pitiayumi"] <- "Parula pitiayumi"
initial_species_list$Pulido[initial_species_list$HBW == "Spinus psaltria"] <- "Carduelis psaltria"
initial_species_list$Pulido[initial_species_list$HBW == "Dives warczewiczi"] <- "Dives warszewiczi"
initial_species_list$Pulido[initial_species_list$HBW == "Pheugopedius coraya"] <- "Thryothorus coraya"
initial_species_list$Pulido[initial_species_list$HBW == "Paroaria nigrogenis"] <- "Paroaria gularis"
initial_species_list$Pulido[initial_species_list$HBW == "Clibanornis rufipectus"] <- "Automolus rufipectus"
initial_species_list$Pulido[initial_species_list$HBW == "Eupsittula pertinax"] <- "Aratinga pertinax"
initial_species_list$Pulido[initial_species_list$HBW == "Arremon atricapillus"] <- "Arremon torquatus"
initial_species_list$Pulido[initial_species_list$HBW == "Psittacara acuticaudatus"] <- "Aratinga acuticaudata"
initial_species_list$Pulido[initial_species_list$HBW == "Setophaga cerulea"] <- "Dendroica cerulea"
initial_species_list$Pulido[initial_species_list$HBW == "Spinus cucullatus"] <- "Carduelis cucullata"
initial_species_list$Pulido[initial_species_list$HBW == "Sporophila angolensis"] <- "Oryzoborus angolensis"
initial_species_list$Pulido[initial_species_list$HBW == "Atlapetes blancae"] <- "Atlapetes schistaceus"
initial_species_list$Pulido[initial_species_list$HBW == "Anthocephala berlepschi"] <- "Anthocephala floriceps"
initial_species_list$Pulido[initial_species_list$HBW == "Vireo chivi"] <- "Vireo olivaceus"
initial_species_list$Pulido[initial_species_list$HBW == "Grallaria fenwickorum"] <- "Grallaria milleri"

# The following species are here lumped with an extralimital near-relative, so that they will not appear as duplicates
# in the Pulido column.
initial_species_list$Pulido[initial_species_list$HBW == "Podiceps andinus"] <- "Podiceps nigricollis"
initial_species_list$Pulido[initial_species_list$HBW == "Megascops gilesi"] <- "Megascops roboratus"
initial_species_list$Pulido[initial_species_list$HBW == "Scytalopus alvarezlopezi"] <- "Scytalopus robbinsi"
initial_species_list$Pulido[initial_species_list$HBW == "Synallaxis beverlyae"] <- "Synallaxis hypospodia"
initial_species_list$Pulido[initial_species_list$HBW == "Synallaxis fuscorufa"] <- "Synallaxis castanea"

initial_species_list$Pulido[initial_species_list$Pulido %ni% pulido.spp]

initial_species_list[duplicated(initial_species_list$Pulido), ]



##### Import Eltontraits data
#download.file("https://ndownloader.figshare.com/files/5631081",
#              destfile = 'Data/Birds/traits/elton.txt')
traits <- read.delim('Data/Birds/traits/elton.txt', header=T, stringsAsFactors = F)

initial_species_list$Pulido[initial_species_list$Pulido %ni% traits$Scientific]


initial_species_list$eltontraits <- NA
for(i in 1:nrow(initial_species_list)){
  if(initial_species_list$Pulido[i] %in% traits$Scientific){
    initial_species_list$eltontraits[i] <- initial_species_list$Pulido[i]
  }
}

initial_species_list$eltontraits[initial_species_list$Pulido == "Pseudocolaptes lawrencii johnsoni"] <- "Pseudocolaptes johnsoni"
initial_species_list$eltontraits[initial_species_list$Pulido == "Xiphorhynchus chunchotambo"] <- "Xiphorhynchus ocellatus"
initial_species_list$eltontraits[initial_species_list$Pulido == "Xenops minutus remoratus"] <- "Xenops minutus"
initial_species_list$eltontraits[initial_species_list$Pulido == "Dendrocolaptes certhia certhia"] <- "Dendrocolaptes certhia"
initial_species_list$eltontraits[initial_species_list$Pulido == "Microxenops milleri"] <- "Xenops milleri"
initial_species_list$eltontraits[initial_species_list$Pulido == "Automolus rubiginosus nigricauda"] <- "Automolus rubiginosus"
initial_species_list$eltontraits[initial_species_list$Pulido == "Certhiasomus stictolaemus"] <- "Deconychura stictolaema"
initial_species_list$eltontraits[initial_species_list$Pulido == "Drymotoxeres pucherani"] <- "Campylorhamphus pucherani"
initial_species_list$eltontraits[initial_species_list$Pulido == "Automolus rufipectus"] <- "Automolus rubiginosus"

initial_species_list[nrow(initial_species_list) + 1, ] <- c('Columba livia', 'Columba livia', NA, 'Columba livia', NA, 'Columba livia', 'Columba livia')
initial_species_list[nrow(initial_species_list) + 1, ] <- c('Grallaria satuata', 'Grallaria rufula', NA, 'Grallaria rufula', NA, 'Grallaria rufula', 'Grallaria rufula')
initial_species_list[nrow(initial_species_list) + 1, ] <- c('Grallaria spatiator', 'Grallaria rufula', NA, 'Grallaria rufula', NA, 'Grallaria rufula', 'Grallaria rufula')
initial_species_list[nrow(initial_species_list) + 1, ] <- c('Grallaria alvarezi', 'Grallaria rufula', NA, 'Grallaria rufula', NA, 'Grallaria rufula', 'Grallaria rufula')

# Remove Pteroglossus aracari, which is included in McMullan in error.
initial_species_list <- initial_species_list[-which(initial_species_list$HBW == "Pteroglossus aracari"), ]

write.csv(initial_species_list, file = "Data/Birds/species_list_creation/initial_species_list.csv")

