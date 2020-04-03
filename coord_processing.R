`%ni%` <- Negate(`%in%`)

##### For collaborative projects--figure out what machine we're on and automatically set the working directory ####
socolar.desktop <- file.exists('/Users/jacobsocolar/Dropbox/Work/Code/machine_identifier_n5L8paM.txt')
socolar.laptop <- file.exists('/Users/jacob/Dropbox/Work/Code/machine_identifier_n5L8paM.txt')
if(socolar.desktop){
  dir.path <- "/Users/JacobSocolar/Dropbox/Work/Colombia/Data"
  googleAPIkey <- readLines('/Users/jacobsocolar/Dropbox/Work/Code/code_keychain/GoogleAPIkey.txt')
}else if(socolar.laptop){
  dir.path <- "/Users/jacob/Dropbox/Work/Colombia/Data"
  googleAPIkey <- readLines('/Users/jacob/Dropbox/Work/Code/code_keychain/GoogleAPIkey.txt')
}# else if(){dir.path <- }
# Edit the above for whatever computer(s) you use.  Just make absolutely sure that the if condition is something that definitely
# wouldn't possibly evaluate as true on anybody else's system, and that none of the preceding conditions could possibly evaluate
# to TRUE on your system!  (This isn't just about making sure that we get the right working directories; in some cases we might
# conceivably invoke system commands for file management that depend on dir.path.)
setwd(dir.path)
############################

jacob1 <- read.csv("Birds/Jacob_data_v1.1.csv")
points_list <- droplevels(unique(jacob1$Point[jacob1$Point != ""]))

files1 <- list.files('GIS/GPS/GPS1')            # GPX folder of my GPS
files2 <- list.files('GIS/GPS/GPS2')       # GPX folder of GPS used by Diego & Marcela

# Older versions of my GPX folder: confirm that they don't have files that have been deleted
gXX <- list.files('GIS/GPS/GPXDump/GPX_old')
gX <- list.files('GIS/GPS/GPXDump/1_3_2018/GPX')
gX2 <- list.files('GIS/GPS/GPXDump/6_3_2019')
gX3 <- list.files('GIS/GPS/GPXDump/10_5_2019')
gXX[gXX %ni% files1]
gX[gX %ni% files1]
gX2[gX2 %ni% files1]
gX3[gX3 %ni% files1]
# files1 contains all of gX3 plus additional points from August and September 2019

# Extract files that potentially contain relevant waypoints:
gpx1 <- paste0('GIS/GPS/GPS1/', files1[grep('\\.gpx', files1)])
gpx2 <- paste0('GIS/GPS/GPS2/', files2[grep('\\.gpx', files2)])
gpx1 <- c(gpx1[grep('Waypoints_', gpx1)], gpx1[grep('/Waypoints.gpx', gpx1)], gpx1[grep('TAB', gpx1)], 
          gpx1[grep('TAP', gpx1)], gpx1[grep('ABF', gpx1)], gpx1[grep('PSF', gpx1)])
gpx2 <- gpx2[grep('Waypoints_', gpx2)]

# Read in waypoints
points1 <- plotKML::readGPX(gpx1[1])$waypoints[ , c('lat','lon','ele','time','name')]
for(i in 2:length(gpx1)){
  points1 <- rbind(points1, plotKML::readGPX(gpx1[i])$waypoints[ , c('lat','lon','ele','time','name')])
}
#points1 <- points1[-c(294, 293), ]
points1$name <- gsub('CCFX', 'CCF', points1$name)
points1$name <- gsub('CFX', 'CCF', points1$name)
points1$name[points1$name == '455'] <- 'PUF12'
points1$name[points1$name == '454'] <- 'PUF11'
points1$name[points1$name == '459'] <- 'PUF9'
points1$name[points1$name == '458'] <- 'PUF8'
points1$name[points1$name == '457'] <- 'PUF7'
points1$name[points1$name == '342'] <- 'BRF2'
points1$name[points1$name == 'BF1'] <- 'BRF1'
points1$name[points1$name == 'BF3'] <- 'BRF3'
points1$name[points1$name == 'BF4'] <- 'BRF4'
points1$name[points1$name == 'BF5'] <- 'BRF5'
points1$name[points1$name == 'BF6'] <- 'BRF6'
points1$name[points1$name == '334'] <- 'BRP4'
points1$name[points1$name == '335'] <- 'BRP5'
points1$name[points1$name == '337'] <- 'BRP6'
points1$name[points1$name == 'BP1'] <- 'BRP1'
points1$name[points1$name == 'BP2'] <- 'BRP2'
points1$name[points1$name == 'BP3'] <- 'BRP3'
points1$name[points1$name == 'BP3'] <- 'BRP3'
points1$name[points1$name == '522'] <- 'RAP4'
points1$name[points1$name == '523'] <- 'RAP6'
points1$name[points1$name == 'SOF 1'] <- 'SOF1'
points1$name[points1$name == 'EPP3N'] <- 'EPP3'
points1$name[points1$name == 'EPP9N'] <- 'EPP9'


points2 <- plotKML::readGPX(gpx2[1])$waypoints[ , c('lat','lon','ele','time','name')]
for(i in 2:length(gpx2)){
  points2 <- rbind(points2, plotKML::readGPX(gpx2[i])$waypoints[ , c('lat','lon','ele','time','name')])
}
points2$name <- gsub('ALB', 'ALF', points2$name)
points2 <- points2[-grep('ALPO', points2$name), ]
points2$name[points2$name == 'PU1'] <- 'PUF1'
points2$name[points2$name == 'PU2'] <- 'PUF2'
points2$name[points2$name == 'PU3'] <- 'PUF3'
points2$name[points2$name == 'PU4'] <- 'PUF4'
points2$name[points2$name == 'PU5'] <- 'PUF5'
points2$name[points2$name == 'PU6'] <- 'PUF6'
points2 <- points2[-which(points2$name == 'ENP1'), ]


DiegoExtras <- data.frame(lat = c(4.99644, 4.99488, 4.9954, 4.99558, 4.99741), 
                          lon = c(-75.40913, -75.41056, -75.42219, -75.424, -75.42269), 
                          ele = NA,
                          time = NA, 
                          name = c('RCF3', 'RCF2', 'RCP1', 'RCP2', 'RCP3'))

points <- rbind(points1, points2, DiegoExtras)
points <- points[(points$lat > -1) & (points$lat < 13), ]

pointsF <- points[grep('^[[:alnum:]][[:alnum:]]F[[:digit:]]', points$name), ]
pointsF$forest <- 1
pointsP <- points[grep('^[[:alnum:]][[:alnum:]]P[[:digit:]]', points$name), ]
pointsP$forest <- 0
pts <- rbind(pointsF, pointsP)
pts <- pts[order(pts$name), ]
pts$ele <- as.numeric(pts$ele)

#############################################################3
points <- pts[1:465,-4]   # remove María José's transect points (rows) and the "time" column

# Handle nonstandard names
points <- points[points$name %ni% c('MOF3', 'SLF2'), ]
points$name[points$name == 'MOF3NEW'] <- 'MOF3'
points$name[points$name == 'SLF2BIEN'] <- 'SLF2'
points$name[points$name == 'IGP1 1'] <- 'IGP1'
points$name[points$name == 'SLP3BUENO'] <- 'SLP3'
points$name[points$name == 'SLF4MIO'] <- 'SLF4'
points <- points[points$name %ni% c('ALF2FIN', 'IGF8A'),]
points <- points[-grep('ENT', points$name), ]

# Remove exact duplicate rows
points <- points[!duplicated(points), ]
row.names(points) <- seq(nrow(points))
View(points)

# Check all elevations against google maps API
gelevs1 <- googleway::google_elevation(df_locations = points[1,], location_type = "individual", key = googleAPIkey)
gelevs <- gelevs1$results[ , c(1,3)]
for(i in 2:nrow(points)){
  gelevs1 <- googleway::google_elevation(df_locations = points[i,], location_type = "individual", key = googleAPIkey)
  if(gelevs1$status != "OK"){c(print("NOT OK"), i)}
  gelevs <- rbind(gelevs, gelevs1$results[ , c(1,3)])
}

points$g_elev_init <- gelevs$elevation

# And against ALOS (JAXA)
library(reticulate)

use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()             # Trigger the authentication

ALOS <- ee$Image('JAXA/ALOS/AW3D30/V2_2')
ALOS_elev <- ALOS$select('AVE_DSM')
# Featurecollection of point geometries
geompts <- sapply(1:nrow(points),function(x)ee$Geometry$Point(c(points$lon[x],points$lat[x])))
geompts <- ee$FeatureCollection(c(unlist(geompts)))
# Extract ALOS elevations for all points - combine into dataframe
pts_elev <- ALOS$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
ALOSelev <- sapply(c(1:length(pts_elev$features)),function(x)pts_elev$features[[x]]$properties$AVE_DSM)

points$ALOS_init <- ALOSelev

# Where are the large discrepancies?
points[which(abs(points$g_elev_init - points$ALOS_init) > 30),]
# ALOS is a WAY better match to our field measurements at the really bad SAF points, and is substantially better
# at the moderately bad CCF points

points[which(abs(points$g_elev_init - points$ALOS_init) <= 30 & abs(points$g_elev_init - points$ALOS_init) > 20),]
# ALOS and google earth API are similar in quality at lower discrepancies, which is probably to be expected since 
# at this level the GPS errors are likely to be large enough to matter. Interestingly, at the PLP points ALOS and 
# field-measurements are in good agreement and google earth tends to run > 20 meters higher.

points$ALOSdiff_init <- abs(points$ele - points$ALOS_init)

# Remove elevation when one point has bad elevation but Felecity's file has the average of the coordinates
points$ele[c(429, 431, 435, 437, 439)] <- NA

# Remove points with bad coordinates and points with identical coordinates (in latter case, retain point with
# good elevation).  Also, if points are duplicated, equally acceptable-seeming, and Felicity's metapoints file
# has one or the other (rather than the average), retain only the one from the metapoints file.
points <- points[-c(15, 17, 19, 20, 23, 107, 109, 111, 113, 115, 117, 125, 127, 130, 237, 239, 241, 292, 293,
                    298, 300, 302, 304, 306, 317, 387, 390, 394, 398, 400, 402), ]

# Potential elevation problems at BRF4-6; BRP1-6, CCF3, all CCP points (problems in extreme heat?), EPP7-8, HTF3,
# IGF1, MOF4-6, MOP1-6, PPP7-9, PSF23, PUF11-12, PUP10-12, RAP2, RBP1-3, RCF1-6, RCP1-6



getCount <- function(x) {
  u <- unique(x);
  data.frame(
    value=u,
    count=sapply(u, function(v) { length(which(x==v)) } )
  )
}
unique(points$name)[getCount(points$name)$count == 3]


p2rm <- vector()
for(i in 2:nrow(points)){
  if(points$name[i] == points$name[i-1]){
    points$lat[i-1] <- mean(c(points$lat[i-1], points$lat[i]))
    points$lon[i-1] <- mean(c(points$lon[i-1], points$lon[i]))
    points$ele[i-1] <- mean(c(points$ele[i-1], points$ele[i]), na.rm = T)
    p2rm <- c(p2rm, i)
  }
}

points <- points[-p2rm, ]

points$cluster <- paste('cluster', gsub('[[:digit:]]', '', points$name), ceiling(as.numeric(gsub('[[:alpha:]]', '', points$name))/3), sep = "_")
points$cluster[235:241] <- c(rep('cluster_PSP_1', 2), rep('cluster_PSP_2', 2), rep('cluster_PSP_3', 3))

# within-cluster distances
for(i in 1:length(unique(cluster)))



dists <- raster::pointDistance(p1 = points[,c(2,1)], p2 = points[,c(2,1)], lonlat = T, allpairs = T)
diag(dists) <- NA

hist(dists[dists < 200])

inds <- which(dists < 194)
bads <- cbind((inds %/% nrow(points)) + 1, inds %% nrow(points))
bads <- bads[which(bads[,1] >= bads[,2]),]

cbind(as.character(points$name[bads[,1]]), as.character(points$name[bads[,2]]))


plot(points$ele ~ points$g_elev_init)
plot(points$ele ~ points$ALOS_init)
plot(points$g_elev_init ~ points$ALOS_init)
hist(points$g_elev_init - points$ALOS_init)
points[points$g_elev_init - points$ALOS_init < -20,]
points[points$g_elev_init - points$ALOS_init > 20,]

write.csv(points, file = "GIS/Points/socolar_points_v1.csv")

