# Read in table of migratory dates compiled by JBS.  
mig_dates <- read.csv("inputs/migratory.csv")
initial_species_list <- read.csv("outputs/initial_species_list.csv")

mig_dates <- mig_dates[mig_dates$latin %in% initial_species_list$HBW, ]

mig_dates$start1[mig_dates$status == "B"] <- mig_dates$start2[mig_dates$status == "B"] <- "Sep_1"
mig_dates$end1[mig_dates$status == "B"] <- mig_dates$end2[mig_dates$status == "B"] <- "May_15"
mig_dates$start1[mig_dates$status == "A"] <- mig_dates$start2[mig_dates$status == "A"] <- "Mar_15"
mig_dates$end1[mig_dates$status == "A"] <- mig_dates$end2[mig_dates$status == "A"] <- "Nov_1"

# get a lookup table to translate between abbreviations in table and ordinal days
unique_dates1 <- unique(c(mig_dates$start1, mig_dates$start2, mig_dates$end1, mig_dates$end2))
unique_dates <- unique_dates1[unique_dates1 != ""]
ud_month <- match(sapply(strsplit(unique_dates, "_"), "[[", 1), month.abb)
ud_day <- sapply(strsplit(unique_dates, "_"), "[[", 2)
ud_year <- 2018
dmy <- paste(ud_day, ud_month, ud_year, sep = "/")
date_lookup <- data.frame(abbrev = unique_dates, oday = lubridate::yday(as.Date(dmy, format = "%d/%m/%Y")))

# Replace abbreviations with ordinal days
for(i in 1:nrow(mig_dates)){
  if(mig_dates$start1[i] != ""){
    mig_dates$start1[i] <- date_lookup$oday[date_lookup$abbrev == mig_dates$start1[i]]
    mig_dates$end1[i] <- date_lookup$oday[date_lookup$abbrev == mig_dates$end1[i]]
  }
  if(mig_dates$start2[i] != ""){
    mig_dates$start2[i] <- date_lookup$oday[date_lookup$abbrev == mig_dates$start2[i]]
    mig_dates$end2[i] <- date_lookup$oday[date_lookup$abbrev == mig_dates$end2[i]]
  }
}

mig_dates$start1 <- as.numeric(mig_dates$start1)
mig_dates$start2 <- as.numeric(mig_dates$start2)
mig_dates$end1 <- as.numeric(mig_dates$end1)
mig_dates$end2 <- as.numeric(mig_dates$end2)

saveRDS(mig_dates, "outputs/mig_dates.RDS")
