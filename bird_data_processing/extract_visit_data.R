# This code was written and run on 2 March 2020 to extract visit data from the version of Jacob_data_v1.1 that was
# current on that date. Future changes to Jacob_data_v1.1 are not guaranteed to remain compatible with this script.

# This script is saved for purposes of internal documentation of the workflow only.

# At the time of execution, Socolar performed detailed manual checks to verify that all/only the correct rows of 
# the data file are being extracted and saved.

setwd("/Users/JacobSocolar/Dropbox/Work/Colombia")
jacob1 <- read.csv("Data/Birds/Jacob_data_v1.1.csv", stringsAsFactors = F)
jj <- jacob1[!is.na(jacob1$Time), 1:12]
jj$observer <- "Socolar"
jj$observer[paste(jj$Point, jj$Take, sep="_") %in% c('IGF10_3', 'IGF11_3', 'IGF10_4', 'IGF11_4', 'IGF12_4')] <- "Edwards"
write.csv(jj, file = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/visit_data.csv")

pts <- unique(jj$Point)
for(i in pts){
  r <- jj[jj$Point == i, ]
  if(nrow(r) != 4){
    print(r)
  }
}

# Points with 1 rep: MOF7-18 (remaining reps are Simon's)
# Points with 2 reps: SAP4-6 (guerrillas)
# Points with 3 reps: ATP3 (bulls), 
#                     PSP (plane crash), 
#                     ABP, ABF (illness)
#                     SAP8-12 (weather + guerrillas)



# still need to add cluster designations to all points