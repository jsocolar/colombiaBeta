`%ni%` <- Negate(`%in%`)

##### For collaborative projects--figure out what machine we're on and automatically set the working directory ####
socolar.desktop <- file.exists('/Users/jacobsocolar/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
socolar.laptop <- file.exists('/Users/jacob/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
if(socolar.desktop){
  dir.path <- "/Users/JacobSocolar/Dropbox/Work/Colombia/Data"
  simon.file.path <- "/Users/jacobsocolar/Google Drive/Simon_data/data/bird data_Jan&Jun2019/data_Jan&Jun2019_currentVersion.xlsx"
  simon.RDS.path <- "/Users/jacobsocolar/Google Drive/Simon_data/data/bird data_Jan&Jun2019/formattedDataset_JanJun2019.rds"
}else if(socolar.laptop){
  dir.path <- "/Users/jacob/Dropbox/Work/Colombia/Data"
}# else if(){dir.path <- }
# Edit the above for whatever computer(s) you use.  Just make absolutely sure that the if condition is something that definitely
# wouldn't possibly evaluate as true on anybody else's system, and that none of the preceding conditions could possibly evaluate
# to TRUE on your system!  (This isn't just about making sure that we get the right working directories; in some cases we might
# conceivably invoke system commands for file management that depend on dir.path.)
setwd(dir.path)
############################

# read in united eBird/HBW taxonomy
taxonomy <- read.csv("Birds/species_list_creation/HBW_eBird_taxonomy.csv")
t2 <- taxonomy[taxonomy$HBW_CAT == "sp" | taxonomy$CLEM_CAT_2019 == "sp", ] # remove taxa that are not treated as species by either eBird or HBW





########
##### Jacob #####
jacob1 <- read.csv("Birds/Jacob_data_v1.1.csv")

# some of the below consists of data checks that won't be necessary for final analysis but are retained for now
# (until the data file is finalized) to guard against typos and entry errors.
all.equal(jacob1$Point == "", is.na(jacob1$Take)) #  make sure that NAs in jacob1$Take universally match ""'s in jacob1$Point
unique(jacob1$Dis)
length(unique(jacob1$Dis))
unique(jacob1$Dist)
length(unique(jacob1$Dist))
unique(jacob1$FO)
length(unique(jacob1$FO))

# Fill in point and visit identifiers
for(i in 2:nrow(jacob1)){
  if(is.na(jacob1$Take[i])){
    jacob1$Take[i] <- jacob1$Take[i - 1]
    jacob1$Point[i] <- jacob1$Point[i - 1]
  }
}

load("Birds/species_list_creation/colombia_species.Rdata")
unique(jacob1$Species[gsub("_", " ", jacob1$Species) %ni% colombia_species])
length(unique(jacob1$Species[gsub("_", " ", jacob1$Species) %ni% colombia_species]))

# Extract the analyzeable records, but do so in several steps to implement typo checks and to keep track of 
# how many species/entries are unanalyzeable.
jacob2 <- droplevels(jacob1[jacob1$Species != "FLOCK" & is.na(jacob1$no_birds) & (is.na(jacob1$Dis) | jacob1$Point_time == "LVG"), ])
unique(jacob2$Species[gsub("_", " ", jacob2$Species) %ni% colombia_species])
length(unique(jacob2$Species[gsub("_", " ", jacob2$Species) %ni% colombia_species]))
unique(jacob1$Species[which(jacob1$Species %ni% jacob2$Species)])
length(unique(jacob1$Species[which(jacob1$Species %ni% jacob2$Species)]))
unique(jacob2$Dist)
length(unique(jacob2$Dist))

jacob_all <- droplevels(jacob2[jacob2$Species %ni% c("Sono", "Visu"), ])

jacob_all$Species[gsub("_", " ", jacob_all$Species) %ni% t2$HBW_LATIN]

jacob_all$ebird <- NA

for(i in 1:nrow(jacob_all)){
  if(gsub("_", " ", jacob_all$Species[i]) %in% t2$HBW_LATIN){
    jacob_all$ebird[i] <- t2$CLEM_SCI_2019[min(which(t2$HBW_LATIN == gsub("_", " ", jacob_all$Species[i])))]
  }else if(jacob_all$Species[i] %in% c("Grallaria_spatiator", "Grallaria_saturata")){
    jacob_all$ebird[i] <- "Grallaria rufula"
  }else if(jacob_all$Species[i] == "Vireo_chivi"){
    jacob_all$ebird[i] <- "Vireo chivi"
  }else if(jacob_all$Species[i] == "Catharus_maculatus"){
    jacob_all$ebird[i] <- "Catharus dryas"
  }
}

unique(jacob_all$ebird)

jacob_all$ebird[jacob_all$ebird == "Xiphorhynchus ocellatus beauperthuysii/lineatocapilla"] <- "Xiphorhynchus ocellatus"
jacob_all$ebird[jacob_all$ebird == "Catharus ustulatus [swainsoni Group]"] <- "Catharus ustulatus"
jacob_all$ebird[jacob_all$ebird == "Butorides virescens/striata"] <- "Butorides striata"
jacob_all$ebird[jacob_all$ebird == "Pionus menstruus menstruus/rubrigularis"] <- "Pionus menstruus"
jacob_all$ebird[jacob_all$ebird == "Amazona farinosa farinosa"] <- "Amazona farinosa"
jacob_all$ebird[jacob_all$ebird == "Amazona festiva festiva"] <- "Amazona festiva"
jacob_all$ebird[jacob_all$ebird == "Aramides albiventris/cajaneus"] <- "Aramides cajaneus"
jacob_all$ebird[jacob_all$ebird == "Coeligena bonapartei bonapartei"] <- "Coeligena bonapartei"

jacob_all$ebird[jacob_all$Species == "Myiodynastes_hemichrysus"] <- "Myiodynastes chrysocephalus"
jacob_all$ebird[jacob_all$Species == "Piaya_cayana"] <- "Piaya cayana"
jacob_all$ebird[jacob_all$Species == "Troglodytes_aedon"] <- "Troglodytes aedon"
jacob_all$ebird[jacob_all$Species == "Colibri_thalassinus"] <- "Colibri cyanotus"
jacob_all$ebird[jacob_all$Species == "Basileuterus_tristriatus"] <- "Basileuterus tristriatus"
jacob_all$ebird[jacob_all$Species == "Chondrohierax_uncinatus"] <- "Chondrohierax uncinatus"
jacob_all$ebird[jacob_all$Species == "Aulacorhynchus_albivitta"] <- "Aulacorhynchus albivitta"
jacob_all$ebird[jacob_all$Species == "Cistothorus_platensis"] <- "Cistothorus platensis"
jacob_all$ebird[jacob_all$Species == "Tangara_aurulenta"] <- "Tangara arthus"
jacob_all$ebird[jacob_all$Species == "Setophaga_pitiayumi"] <- "Setophaga pitiayumi"
jacob_all$ebird[jacob_all$Species == "Leptopogon_superciliaris"] <- "Leptopogon superciliaris"
jacob_all$ebird[jacob_all$Species == "Picumnus_squamulatus"] <- "Picumnus squamulatus"
jacob_all$ebird[jacob_all$Species == "Turdus_ignobilis"] <- "Turdus ignobilis"
jacob_all$ebird[jacob_all$Species == "Grallaria_rufula"] <- "Grallaria rufula"
jacob_all$ebird[jacob_all$Species == "Trogon_violaceus" & !grepl("PLF", jacob_all$Point)
                & !grepl("PLP", jacob_all$Point)
                & !grepl("PSF", jacob_all$Point)
                & !grepl("PSP", jacob_all$Point)
                & !grepl("SGP", jacob_all$Point)] <- "Trogon violaceus"
jacob_all$ebird[jacob_all$Species == "Trogon_violaceus" & jacob_all$ebird == ""] <- "Trogon ramonianus"
jacob_all$ebird[jacob_all$Species == "Phaeomyias_murina"] <- "Phaeomyias murina"
jacob_all$ebird[jacob_all$Species == "Amazona_autumnalis"] <- "Amazona autumnalis"
jacob_all$ebird[jacob_all$Species == "Automolus_ochrolaemus"] <- "Automolus ochrolaemus"
jacob_all$ebird[jacob_all$Species == "Dendrocincla_fuliginosa"] <- "Dendrocincla fuliginosa"
jacob_all$ebird[jacob_all$Species == "Anthracothorax_nigricollis"] <- "Anthracothorax nigricollis"
jacob_all$ebird[jacob_all$Species == "Chalybura_buffonii"] <- "Chalybura buffonii"
jacob_all$ebird[jacob_all$Species == "Sittasomus_griseus"] <- "Sittasomus griseicapillus"
jacob_all$ebird[jacob_all$Species == "Amazilia_saucerottei"] <- "Amazilia saucerottei"
jacob_all$ebird[jacob_all$Species == "Pyrocephalus_rubinus"] <- "Pyrocephalus rubinus"
jacob_all$ebird[jacob_all$Species == "Contopus_bogotensis"] <- "Contopus cinereus"
jacob_all$ebird[jacob_all$Species == "Polioptila_plumbea"] <- "Polioptila plumbea"
jacob_all$ebird[jacob_all$Species == "Formicarius_analis"] <- "Formicarius analis"
jacob_all$ebird[jacob_all$Species == "Corapipo_leucorrhoa"] <- "Corapipo leucorrhoa"
jacob_all$ebird[jacob_all$Species == "Xenops_genibarbis"] <- "Xenops minutus"
jacob_all$ebird[jacob_all$Species == "Turdus_albicollis"] <- "Turdus albicollis"
jacob_all$ebird[jacob_all$Species == "Colaptes_rubiginosus"] <- "Colaptes rubiginosus"
jacob_all$ebird[jacob_all$Species == "Chaetura_chapmani"] <- "Chaetura chapmani"
jacob_all$ebird[jacob_all$Species == "Xiphorhynchus_guttatoides"] <- "Xiphorhynchus guttatus"
jacob_all$ebird[jacob_all$Species == "Tolmomyias_assimilis"] <- "Tolmomyias assimilis"
jacob_all$ebird[jacob_all$Species == "Sittasomus_griseicapillus"] <- "Sittasomus griseicapillus"
jacob_all$ebird[jacob_all$Species == "Myrmothera_campanisona"] <- "Myrmothera campanisona"
jacob_all$ebird[jacob_all$Species == "Tunchiornis_ochraceiceps"] <- "Tunchiornis ochraceiceps"
jacob_all$ebird[jacob_all$Species == "Polioptila_guianensis"] <- "Polioptila facilis"
jacob_all$ebird[jacob_all$Species == "Pyrrhura_melanura"] <- "Pyrrhura melanura"
jacob_all$ebird[jacob_all$Species == "Threnetes_leucurus"] <- "Threnetes leucurus"
jacob_all$ebird[jacob_all$Species == "Glaucidium_brasilianum"] <- "Glaucidium brasilianum"
jacob_all$ebird[jacob_all$Species == "Psophia_crepitans"] <- "Psophia crepitans"


jacob_all_new <- jacob_all
jacob_all_new$pt_take <- paste(jacob_all_new$Point, jacob_all_new$Take, sep = "__")
jacob_all_new$pt_take_sp <- paste(jacob_all_new$pt_take, jacob_all_new$ebird, sep = "__")
upts <- unique(jacob_all_new$pt_take_sp)


for(i in 1:length(upts)){
  jacob_all_new$Count[jacob_all_new$pt_take_sp == upts[i]] <- jacob_all_new$Count[jacob_all_new$pt_take_sp == upts[i]]
}

jacob_all_new <- jacob_all_new[!duplicated(jacob_all_new$pt_take_sp),]


all_pts <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Points/all_pts.RDS")


ebird <- data.frame(Common_name = rep("", nrow(jacob_all_new)))
ebird <- cbind(ebird, reshape2::colsplit(jacob_all_new$ebird ," ",c("Genus","species")))
ebird$number <- jacob_all_new$Count
ebird$species_comment <- ""
ebird$location_name <- jacob_all_new$Point
ebird$Latitude <- NA
ebird$Longitude <- NA
for(i in 1:nrow(ebird)){
  ebird$Latitude[i] <- all_pts$lat[all_pts$point == ebird$location_name[i]]
  ebird$Longitude[i] <- all_pts$lon[all_pts$point == ebird$location_name[i]]
}

ebird$Date <- NA
ebird$Start_time <- NA
for(i in 1:nrow(ebird)){
  ttt <- jacob_all_new$Take[i]
  if(ttt == 1){d <- all_pts$posix1[all_pts$point == ebird$location_name[i]]}
  if(ttt == 2){d <- all_pts$posix2[all_pts$point == ebird$location_name[i]]}
  if(ttt == 3){d <- all_pts$posix3[all_pts$point == ebird$location_name[i]]}
  if(ttt == 4){d <- all_pts$posix4[all_pts$point == ebird$location_name[i]]}
  ebird$Date[i] <- format(as.Date(d), "%m/%d/%Y")
  ebird$Start_time[i] <- strftime(d, format="%H:%M")
}


ebird$state <- ""
ebird$country <- "CO"
ebird$protocol <- "stationary"
ebird$n_obs <- 1
ebird$duration <- NA
for(i in 1:nrow(ebird)){
  times <- jacob_all_new$Point_time[jacob_all_new$pt_take == jacob_all_new$pt_take[i]]
  times[times == "LVG"] <- 1100
  times <- as.numeric(times)
  ebird$duration[i] <- max(10, max(ceiling(times/100)), na.rm = T)
}
ebird$all_obs_rep <- "Y"
ebird$distance <- ""
ebird$acres <- ""
ebird$sub_comments <- "part of bulk upload S1"

ebird <- ebird[-which(((jacob_all_new$Point %in% c("IGF10", "IGF11", "IGF12")) & jacob_all_new$Take == 4) |
                 ((jacob_all_new$Point %in% c("IGF10", "IGF11")) & jacob_all_new$Take == 3)), ]

ebird$Common_name[1] <- "White-lored Warbler"

names(ebird) <- NULL



write.csv(ebird, "/Users/Jacobsocolar/Desktop/ebird_bulk_S1.csv", row.names = F)

