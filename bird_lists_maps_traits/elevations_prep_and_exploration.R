##### Script dependencies: Species_lists.R #####
`%ni%` <- Negate(`%in%`)

elevations <- read.csv("/Users/jacobsocolar/Google Drive/Simon_data/data/elevational_ranges_Quinones.csv")
elevations$latin <- paste(elevations$genus, elevations$species, sep = ' ')
elevations <- elevations[elevations$latin != "Frederickena unduliger", ]
initial_species_list <- read.csv("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/species_list_creation/initial_species_list.csv")

# One-to-one (from a Colombian perspective) synonymy that is resolved by HBW/eBird/EltonTraits lookup
for(i in 1:nrow(elevations)){
  if((elevations$latin[i] %ni% initial_species_list$HBW) &             #  if the Ayerbe name is missing from HBW    
     (sum(initial_species_list$eBird == elevations$latin[i]) == 1)){  # but is present exactly once in the eBird synonymy
    elevations$latin[i] <- initial_species_list$HBW[initial_species_list$eBird == elevations$latin[i]]   # Replace Ayerbe name with corresponding HBW name
  }
}
for(i in 1:nrow(elevations)){
  if((elevations$latin[i] %ni% c(initial_species_list$HBW, initial_species_list$eBird)) &   #  if the Ayerbe name is missing from HBW & eBird
     (sum(initial_species_list$eltontraits == elevations$latin[i]) == 1)){  # but is present exactly once in the EltonTraits synonymy
    elevations$latin[i] <- initial_species_list$HBW[initial_species_list$eltontraits == elevations$latin[i]]   # Replace Ayerbe name with corresponding HBW name
  }
}

# A weird case: Zimmerius vilissimus gets split to parvus and improbus, but because Elton, eBird, and HBW all split improbus but
# keep parvus lumped with vilissimus, the autoconversion changes vilissimus directly to parvus.  Here we change it back so that
# we can handle the splitting later on.
elevations$latin[elevations$latin == "Zimmerius parvus"] <- "Zimmerius vilissimus"

elevations$latin[elevations$latin == "Dolichonyx orzyivorus"] <- "Dolichonyx oryzivorus"  # Fix spelling error in elevations table
elevations$latin[elevations$latin == "Asemospiza fuliginosus"] <- "Asemospiza fuliginosa" # Fix spelling error in elevations table
elevations$latin[elevations$latin == "Passerina cuanea"] <- "Passerina cyanea" # Fix spelling error in elevations table


##### What HBW species are now missing elevations? #####
initial_species_list$HBW[initial_species_list$HBW %ni% elevations$latin]

missing_species <- data.frame(english = NA, genus = c("Troglodytes", "Grallaria", "Hemitriccus"), species = c("ochraceus", "rufocinerea", "inornatus"),
                              lower = c(1000, 2000, 0), upper = c(2200, 3100, 200), is_na = "", notes = "from McMullan")
missing_species$latin <- paste(missing_species$genus, missing_species$species, sep = " ")
elevations2 <- rbind(elevations, missing_species)
elevations2$lower[is.na(elevations2$lower & !is.na(elevations2$is_na))] <- 0
write.csv(elevations2, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/traits/elevations/elevations2.csv")


# elevations2.csv is modified by hand to elevations_final.csv to account for specific elevational limits for splits,
# inclusive elevational limits for lumps, and insertion of McMullan data for species with presences in our data at
# species-standardized elevations outside of [-1,2].  Affected species are Scytalopus vicinior, Syrigma sibilatrix,
# Mitrephanes, Dendroplex picus, Myrmoborus leucophrys, Creurgops verticalis.
elevations_final <- read.csv("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/traits/elevations/elevations_final.csv")
elevations_final <- elevations_final[!is.na(elevations_final$is_na),]


elevations_final$latin[elevations_final$latin %ni% initial_species_list$HBW]
initial_species_list$HBW[initial_species_list$HBW %ni% elevations_final$latin]



