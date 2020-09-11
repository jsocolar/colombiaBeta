`%ni%` <- Negate(`%in%`)
initial_species_list <- read.csv("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/species_list_creation/initial_species_list.csv")
migDates <- read.csv("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/traits/migratory.csv")

initial_species_list$HBW[initial_species_list$HBW %ni% migDates$latin]
migDates2 <- migDates[migDates$latin %in% initial_species_list$HBW,]

mig_codes <- c("B", "A", "P", "B2", "A2", "P2")

for(i in 1:nrow(migDates2)){
  if(migDates2$status[i] %in% c("B", "B2")){
    if(any(migDates2[i, 3:6] != "")){stop()}
    migDates2$start1[i] <- "Sep_1"
    migDates2$end1[i] <- "May_15"
  }
  if(migDates2$status[i] %in% c("A", "A2")){
    if(any(migDates2[i, 3:6] != "")){stop()}
    migDates2$start1[i] <- "Mar_1"
    migDates2$end1[i] <- "Nov_1"
  }
  if(migDates2$status[i] %in% mig_codes){
    if(migDates2$start1[i] == "" | migDates2$end1[i] == ""){stop()}
    if((migDates2$start2[i] == "" | migDates2$end2[i] == "") & migDates2$start2[i] != migDates2$end2[i]){stop()}
  }
  if(migDates2$start2[i] == ""){
    migDates2$start2[i] <- migDates2$start1[i]
    migDates2$end2[i] <- migDates2$end1[i]
  }
}

alldates <- unique(c(migDates$start1, migDates$start2, migDates$end1, migDates$end2))[2:length(unique(c(migDates$start1, migDates$start2, migDates$end1, migDates$end2)))]
alldates

ad2 <- strsplit(alldates, "_")
month_lookup <- day_lookup <- oday_lookup <- vector()
date_lookup <- list()
for(i in 1:length(ad2)){
  month_lookup[i] <- which(month.abb == ad2[[i]][1])
  day_lookup[i] <- ad2[[i]][2]
  date_lookup[[i]] <- lubridate::mdy(paste0(month_lookup[i], "-", day_lookup[i], "-", "2019"))
  oday_lookup[i] <- lubridate::yday(date_lookup[[i]])
}

alldates_lookup <- data.frame(abbr = alldates, oday = oday_lookup)
