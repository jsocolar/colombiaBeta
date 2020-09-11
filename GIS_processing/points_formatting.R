##### Script dependencies: coord_processing.R, GEE_setup.sh, bird_import_and_cleaning.R #####

library(reticulate)

`%ni%` <- Negate(`%in%`)

##### For collaborative projects--figure out what machine we're on and automatically set the working directory ####
socolar.desktop <- file.exists('/Users/jacobsocolar/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
socolar.laptop <- file.exists('/Users/jacob/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
if(socolar.desktop){
  dir.path <- "/Users/JacobSocolar/Dropbox/Work/Colombia/Data"
  simon.file.path <- "/Users/jacobsocolar/Google Drive/Simon_data/data/bird data_Jan&Jun2019/data_Jan&Jun2019_currentVersion.xlsx"
  simon.RDS.path <- "/Users/jacobsocolar/Google Drive/Simon_data/data/bird data_Jan&Jun2019/formattedDataset_JanJun2019.rds"
  simon.points.path <- "/Users/jacobsocolar/Google Drive/Simon_data/data/points/points_EasternCordillera.rds"
}else if(socolar.laptop){
  dir.path <- "/Users/jacob/Dropbox/Work/Colombia/Data"
}# else if(){dir.path <- }
# Edit the above for whatever computer(s) you use.  Just make absolutely sure that the if condition is something that definitely
# wouldn't possibly evaluate as true on anybody else's system, and that none of the preceding conditions could possibly evaluate
# to TRUE on your system!  (This isn't just about making sure that we get the right working directories; in some cases we might
# conceivably invoke system commands for file management that depend on dir.path.)
setwd(dir.path)
############################

points_meta <- read.csv('GIS/Points/CO_sampling_points_metafile.csv')
wandes_pts <- read.csv("GIS/Points/James_andes_points.csv")
wandes2_pts <- points_meta[points_meta$region == "occidental", ]
llanos_pts <- readxl::read_excel("GIS/Points/Llanos master habitat data.xlsx")
eandes_pts <- readRDS(simon.points.path)
socolar_pts <- read.csv("GIS/Points/socolar_points_v1.csv")
socolar_pts <- socolar_pts[socolar_pts$name %ni% c("MOF1", "MOF2", "MOF3", "MOF4", "MOF5", "MOF6", "MOP1", "MOP2", "MOP3", "MOP4", "MOP5", "MOP6"), ]

eandes_pts$geometry <- NULL
eandes_pts <- eandes_pts[eandes_pts$group == "SM", ]


socolar_part <- data.frame(point = socolar_pts$name, lat = socolar_pts$lat, lon = socolar_pts$lon, site = substr(socolar_pts$name, start = 1, stop = 2), cluster = socolar_pts$cluster,
                              birds = 1, beetles = 1, habitat = NA)
socolar_part$birds[socolar_part$cluster %in% c("cluster_RCF_2", "cluster_RCP_2")] <- 0
socolar_part$habitat[substr(socolar_pts$name, start = 3, stop = 3) == "F"] <- "Forest"
socolar_part$habitat[substr(socolar_pts$name, start = 3, stop = 3) == "P"] <- "Pasture"
socolar_part$habitat[socolar_part$point %in% c("IGF10", "IGF11", "IGF12", "PNF1", "PNF2", "PNF3", "PUF4", "PUF5", "PUF6")] <- "Paramo"

simon_part <- data.frame(point = eandes_pts$point, lat = eandes_pts$lat, lon = eandes_pts$lon, site = eandes_pts$site, cluster = eandes_pts$cluster, 
                         birds = 1, beetles = 1, habitat = eandes_pts$habitat)
simon_part <- simon_part[simon_part$point != 'TAP1', ]

llanos_part <- data.frame(point = llanos_pts$Point, lat = llanos_pts$Latitude, lon = llanos_pts$Longitude, site = llanos_pts$Site, cluster = llanos_pts$Square,
                          birds = as.numeric(llanos_pts$`bird_point?` == "Y"), beetles = 1, habitat = llanos_pts$Habitat)

wandes_part <- data.frame(point = wandes_pts$Point, lat = wandes_pts$lat, lon = wandes_pts$lon, site = wandes_pts$Site, cluster = stringr::str_extract(wandes_pts$Point, ".*\\D(?=\\d)"),
                          birds = 1, beetles = 1, habitat = wandes_pts$Habcode)

wandes2_part <- data.frame(point = wandes2_pts$point_id, lat = wandes2_pts$lat, lon = wandes2_pts$long, site = wandes2_pts$site, cluster = stringr::str_extract(wandes2_pts$point_id, ".*\\D(?=\\d)"),
                           birds = as.numeric(wandes2_pts$point_id %in% wandes_pts$Point), beetles = 1, habitat = NA)

for(i in 1:nrow(wandes2_part)){
  if(wandes2_part$point[i] %in% wandes_part$point){
    wandes2_part$habitat[i] <- wandes_part$habitat[wandes_part$point == wandes2_part$point[i]]
  }
}

# Check against wandes from james against wandes2 from metafile
check_wandes <- merge(wandes_part, wandes2_part, by.x = "point", by.y = "point")
discrepancies <- check_wandes[check_wandes$lat.x != check_wandes$lat.y | check_wandes$lon.x != check_wandes$lon.y, ]
discrepancies$distance <- 111000 * sqrt((discrepancies$lat.x - discrepancies$lat.y)^2 + (discrepancies$lon.x - discrepancies$lon.y)^2)
discrepancies <- discrepancies[order(discrepancies$distance), ]
View(discrepancies)


# Assemble all points
all_pts <- rbind(socolar_part, simon_part, llanos_part, wandes2_part)
length(unique(all_pts$point)) == nrow(all_pts)

# Check other points against metafile
points_meta$point_id <- gsub("d", "D", points_meta$point_id)
sum(all_pts$point %ni% points_meta$point_id)
View(points_meta[points_meta$point_id %ni% all_pts$point,])
check_all <- merge(points_meta, all_pts, by.x = "point_id", by.y = "point")
discrepancies <- check_all[check_all$lat.x != check_all$lat.y | check_all$lon != check_all$long, ]
discrepancies$distance <- 111000 * sqrt((discrepancies$lat.x - discrepancies$lat.y)^2 + (discrepancies$lon - discrepancies$long)^2)
discrepancies <- discrepancies[order(discrepancies$distance), ]
View(discrepancies)



unique(all_pts$habitat)

all_pts$natural <- as.numeric(all_pts$habitat %in% c("Forest", "Paramo", "FOREST", "P", "Sm"))
all_pts$paramo <- as.numeric(all_pts$habitat == "Paramo")
all_pts$pasture <- as.numeric(all_pts$habitat %in% c("Pasture", "PASTURE", "G"))
all_pts$other <- as.numeric(all_pts$habitat %in% c("PALM", "Sy"))

all_pts$mixed_cluster <- 0
for(i in 1:nrow(all_pts)){
  cluster <- all_pts$cluster[i]
  if(length(unique(all_pts$habitat[all_pts$cluster == cluster & !is.na(all_pts$habitat)])) != 1){all_pts$mixed_cluster[i] <- 1}
}

# get ALOS elevations
use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()             # Trigger the authentication

ALOS <- ee$Image('JAXA/ALOS/AW3D30/V2_2')
ALOS_elev <- ALOS$select('AVE_DSM')
# Featurecollection of point geometries
geompts <- sapply(1:nrow(all_pts),function(x)ee$Geometry$Point(c(all_pts$lon[x],all_pts$lat[x])))
geompts <- ee$FeatureCollection(c(unlist(geompts)))
# Extract ALOS elevations for all points - combine into dataframe
pts_elev <- ALOS$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
ALOSelev <- sapply(c(1:length(pts_elev$features)),function(x)pts_elev$features[[x]]$properties$AVE_DSM)

all_pts$elev_ALOS <- ALOSelev

wandes1 <- read.csv("Birds/James_WAndes_all_birds.csv")
llanos <- read.csv("Birds/James_llanos_all_birds.csv")
simon1 <- data.frame(readxl::read_excel(simon.file.path))
simon_visit_data <- simon1[!is.na(simon1$Time), names(simon1) %in% c("Point", "Visit", "Date", "Wind", "Vis", "Rain", "Sun", "Time", "Observer")]

# Correct minor errors in data file
simon_visit_data$Date[simon_visit_data$Date == 42237] <- "22/08/2019"
simon_visit_data$Date[simon_visit_data$Date == "31/06/2019"] <- "31/07/2019"
simon_visit_data$Visit[simon_visit_data$Point == "CHA3" & simon_visit_data$Date == "23/07/2019"] <- 2
simon_visit_data$Visit[simon_visit_data$Point == "CHA2" & simon_visit_data$Date == "23/07/2019"] <- 2

# format times/dates to POSIXct times, with correct timezone, etc
simon_visit_data$t2 <- lubridate::hms(sapply(strsplit(as.character(simon_visit_data$Time), " "), "[[", 2))
simon_visit_data$posix <- lubridate::force_tz(simon_visit_data$t2 + as.Date(simon_visit_data$Date, format = "%d/%m/%Y"), "America/Bogota")

simon_visit_data$Observer[is.na(simon_visit_data$Observer)] <- "SCM"

jacob1 <- read.csv("Birds/Jacob_data_v1.1.csv")
jacob_visit_data <- jacob1[!is.na(jacob1$Time),]
jacob_visit_data$Date <- gsub("/18$", "/2018", jacob_visit_data$Date)
jacob_visit_data$Date <- gsub("/19$", "/2019", jacob_visit_data$Date)
jacob_visit_data$d2 <- as.Date(jacob_visit_data$Date, format = "%d/%m/%Y")
jacob_visit_data$t2 <- paste(stringr::str_sub(jacob_visit_data$Time, 1, nchar(jacob_visit_data$Time) - 2),
                             stringr::str_sub(jacob_visit_data$Time, nchar(jacob_visit_data$Time) - 1, nchar(jacob_visit_data$Time)),
                             sep = ":")
jacob_visit_data$posix <- lubridate::force_tz(jacob_visit_data$d2 + lubridate::hm(jacob_visit_data$t2), "America/Bogota")

llanos$pt_ct <- paste(llanos$Point, llanos$Count_event, sep = "__")
llanos_visit_data <- llanos[!duplicated(llanos$pt_ct), ]
llanos_visit_data$posix <- lubridate::force_tz(as.Date(llanos_visit_data$Date, format = "%m/%d/%y") + lubridate::hm(llanos_visit_data$Time), "America/Bogota")

wandes1$pt_ct <- paste(wandes1$Point, wandes1$Count, sep = "__")
wandes_visit_data <- wandes1[!duplicated(wandes1$pt_ct), ]
wandes_visit_data$Date[grepl("CORE1_", wandes_visit_data$Point) & wandes_visit_data$Count == 1] <- "07/03/2012"
wandes_visit_data$Date[grepl("CORE1_", wandes_visit_data$Point) & wandes_visit_data$Count == 2] <- "08/03/2012"
wandes_visit_data$Date[grepl("CORE1_", wandes_visit_data$Point) & wandes_visit_data$Count == 3] <- "09/03/2012"
wandes_visit_data$Date[grepl("CORE1_", wandes_visit_data$Point) & wandes_visit_data$Count == 4] <- "10/03/2012"

wandes_visit_data$Date[grepl("CORE2_", wandes_visit_data$Point) & wandes_visit_data$Count == 1] <- "07/03/2012"
wandes_visit_data$Date[grepl("CORE2_", wandes_visit_data$Point) & wandes_visit_data$Count == 2] <- "08/03/2012"
wandes_visit_data$Date[grepl("CORE2_", wandes_visit_data$Point) & wandes_visit_data$Count == 3] <- "09/03/2012"
wandes_visit_data$Date[grepl("CORE2_", wandes_visit_data$Point) & wandes_visit_data$Count == 4] <- "10/03/2012"

wandes_visit_data$Date[grepl("CORE3_", wandes_visit_data$Point) & wandes_visit_data$Count == 1] <- "07/03/2012"
wandes_visit_data$Date[grepl("CORE3_", wandes_visit_data$Point) & wandes_visit_data$Count == 2] <- "08/03/2012"
wandes_visit_data$Date[grepl("CORE3_", wandes_visit_data$Point) & wandes_visit_data$Count == 3] <- "09/03/2012"
wandes_visit_data$Date[grepl("CORE3_", wandes_visit_data$Point) & wandes_visit_data$Count == 4] <- "10/03/2012"

wandes_visit_data$Date[grepl("EXTRA9_", wandes_visit_data$Point) & wandes_visit_data$Count == 1] <- "10/03/2012"
wandes_visit_data$Date[grepl("EXTRA9_", wandes_visit_data$Point) & wandes_visit_data$Count == 2] <- "11/03/2012"
wandes_visit_data$Date[grepl("EXTRA9_", wandes_visit_data$Point) & wandes_visit_data$Count == 3] <- "12/03/2012"
wandes_visit_data$Date[grepl("EXTRA9_", wandes_visit_data$Point) & wandes_visit_data$Count == 4] <- "12/03/2012"

wandes_visit_data$Date[wandes_visit_data$Point == "EXTRA9_3" & wandes_visit_data$Count == 4] <- "11/03/2012"

wandes_visit_data$posix <- lubridate::force_tz(as.Date(wandes_visit_data$Date, format = "%d/%m/%Y") + lubridate::hm(wandes_visit_data$Time), "America/Bogota")

for(i in 1:nrow(all_pts)){
  counter <- 0
  if(all_pts$birds[i] == 1){
    if(all_pts$point[i] %in% jacob_visit_data$Point){counter <- counter + 1}
    if(all_pts$point[i] %in% simon_visit_data$Point){counter <- counter + 1}
    if(all_pts$point[i] %in% wandes1$Point){counter <- counter + 1}
    if(all_pts$point[i] %in% llanos$Point){counter <- counter + 1}
    if(counter != 1 & length(grep('MOF', all_pts$point[i])) == 0){stop()}
  }
}

posix_nas <- rep(NA, nrow(all_pts))
class(posix_nas) <- c("POSIXct", "POSIXt")
posix_nas <- lubridate::force_tz(posix_nas, "America/Bogota")
all_pts$posix4 <- all_pts$posix3 <- all_pts$posix2 <- all_pts$posix1 <- posix_nas

all_pts$obs4 <- all_pts$obs3 <- all_pts$obs2 <- all_pts$obs1 <- NA

for(i in 1:nrow(all_pts)){
  if(all_pts$birds[i] == 1){
    if(all_pts$point[i] %in% jacob_visit_data$Point){
      vd <- jacob_visit_data[jacob_visit_data$Point == all_pts$point[i], ]
      if(1 %in% vd$Take){
        all_pts$posix1[i] <- vd$posix[vd$Take == 1]
        all_pts$obs1[i] <- "JBS"
      }
      if(2 %in% vd$Take){
        all_pts$posix2[i] <- vd$posix[vd$Take == 2]
        all_pts$obs2[i] <- "JBS"
      }
      if(3 %in% vd$Take){
        all_pts$posix3[i] <- vd$posix[vd$Take == 3]
        all_pts$obs3[i] <- "JBS"
      }
      if(4 %in% vd$Take){
        all_pts$posix4[i] <- vd$posix[vd$Take == 4]
        all_pts$obs4[i] <- "JBS"
      }
    }
    
    
    if(all_pts$point[i] %in% simon_visit_data$Point){
      vd <- simon_visit_data[simon_visit_data$Point == all_pts$point[i], ]
      if(1 %in% vd$Visit){
        all_pts$posix1[i] <- vd$posix[vd$Visit == 1]
        all_pts$obs1[i] <- vd$Observer[vd$Visit == 1]
      }
      if(2 %in% vd$Visit){
        all_pts$posix2[i] <- vd$posix[vd$Visit == 2]
        all_pts$obs2[i] <- vd$Observer[vd$Visit == 2]
      }
      if(3 %in% vd$Visit){
        all_pts$posix3[i] <- vd$posix[vd$Visit == 3]
        all_pts$obs3[i] <- vd$Observer[vd$Visit == 3]
      }
      if(4 %in% vd$Visit){
        all_pts$posix4[i] <- vd$posix[vd$Visit == 4]
        all_pts$obs4[i] <- vd$Observer[vd$Visit == 4]
      }
    }

    
    
    if(all_pts$point[i] %in% wandes_visit_data$Point){
      vd <- wandes_visit_data[wandes_visit_data$Point == all_pts$point[i], ]
      if(1 %in% vd$Count){
        all_pts$posix1[i] <- vd$posix[vd$Count == 1]
        all_pts$obs1[i] <- vd$Observer[vd$Count == 1]
      }
      if(2 %in% vd$Count){
        all_pts$posix2[i] <- vd$posix[vd$Count == 2]
        all_pts$obs2[i] <- vd$Observer[vd$Count == 2]
      }
      if(3 %in% vd$Count){
        all_pts$posix3[i] <- vd$posix[vd$Count == 3]
        all_pts$obs3[i] <- vd$Observer[vd$Count == 3]
      }
      if(4 %in% vd$Count){
        all_pts$posix4[i] <- vd$posix[vd$Count == 4]
        all_pts$obs4[i] <- vd$Observer[vd$Count == 4]
      }
    }
    
    if(all_pts$point[i] %in% llanos_visit_data$Point){
      vd <- llanos_visit_data[llanos_visit_data$Point == all_pts$point[i], ]
      if(1 %in% vd$Count_event){
        all_pts$posix1[i] <- vd$posix[vd$Count_event == 1]
        all_pts$obs1[i] <- "JJG"
      }
      if(2 %in% vd$Count_event){
        all_pts$posix2[i] <- vd$posix[vd$Count_event == 2]
        all_pts$obs2[i] <- "JJG"
      }
      if(3 %in% vd$Count_event){
        all_pts$posix3[i] <- vd$posix[vd$Count_event == 3]
        all_pts$obs3[i] <- "JJG"
      }
      if(4 %in% vd$Count_event){
        all_pts$posix4[i] <- vd$posix[vd$Count_event == 4]
        all_pts$obs4[i] <- "JJG"
      }
    }
  }
}


all_pts$obs3[all_pts$point %in% c("IGF10", "IGF11")] <- "DPE"
all_pts$obs3[all_pts$point %in% c("IGF10", "IGF11", "IGF12")] <- "DPE"

all_pts$oday1 <- lubridate::yday(all_pts$posix1)
all_pts$oday2 <- lubridate::yday(all_pts$posix2)
all_pts$oday3 <- lubridate::yday(all_pts$posix3)
all_pts$oday4 <- lubridate::yday(all_pts$posix4)

# function to calculate hours-post-sunrise for a given lat, lon, date, and time
hps_calc <- function(lat, lon, posix){
  sc_data <- data.frame(date = as.Date(posix),
                     lat = lat,
                     lon = lon)
  sunrise_t <- suncalc::getSunlightTimes(data = sc_data, keep = "sunrise")
  return(difftime(posix, sunrise_t$sunrise, units = "hours"))
}

all_pts$hps1 <- as.numeric(hps_calc(all_pts$lat, all_pts$lon, all_pts$posix1))
all_pts$hps2 <- as.numeric(hps_calc(all_pts$lat, all_pts$lon, all_pts$posix2))
all_pts$hps3 <- as.numeric(hps_calc(all_pts$lat, all_pts$lon, all_pts$posix3))
all_pts$hps4 <- as.numeric(hps_calc(all_pts$lat, all_pts$lon, all_pts$posix4))


saveRDS(all_pts, file = "GIS/Points/all_pts.RDS")




