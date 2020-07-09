load('/Users/jacobsocolar/Dropbox/Work/Useful_data/BirdlifeTraits/nf_species.Rdata')
birdlife_list <- readxl::read_xlsx('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/species_list_creation/HBW-BirdLife_Checklist_v4_Dec19/HBW-BirdLife_List_of_Birds_v4.xlsx')

neotrop_pages <- list()
for(i in 1:length(nf_species)){
  print(i)
  sn <- nf_species[i]
  cn <- birdlife_list$`Common name`[birdlife_list$`Scientific name` == sn]
  cn <- gsub("'", "", cn)
  cn <- gsub('á', 'a', cn)
  cn <- gsub('ç', 'c', cn)
  urlname <- paste(c(strsplit(cn, ' ')[[1]], strsplit(sn, ' ')[[1]]), sep = '-', collapse = '-')
  if(sn == 'Vireo chivi'){urlname <- 'vireo-chivi'}
  theurl <- paste0('http://datazone.birdlife.org/species/factsheet/', urlname, '/details')
  neotrop_pages[[i]] = readLines(theurl)
}

save(neotrop_pages, file = "/Users/jacobsocolar/Dropbox/Work/Useful_data/BirdlifeTraits/neotrop_pages.Rdata")

load('/Users/jacobsocolar/Dropbox/Work/Useful_data/BirdlifeTraits/neotrop_pages.Rdata')

# Confirm that all species pages include the following two lines
for(i in 1:length(neotrop_pages)){
  begin_line <- grep('Habitat \\(level 1\\)', neotrop_pages[[i]])
  end_line <- grep('Occasional altitudinal limits', neotrop_pages[[i]])
  if(length(end_line) == 0 | length(begin_line) == 0){stop()}
}


habitats <- as.list(rep(NA, length(nf_species)))
names(habitats) <- nf_species
for(i in 1:length(nf_species)){
  begin_line <- grep('Habitat \\(level 1\\)', neotrop_pages[[i]])
  end_line <- grep('Occasional altitudinal limits', neotrop_pages[[i]])
  lines_suitable1 <- grep('suitable', neotrop_pages[[i]])
  lines_suitable <- lines_suitable1[lines_suitable1 > begin_line & lines_suitable1 < end_line]
  lines_major1 <- grep('major', neotrop_pages[[i]])
  lines_major <- lines_major1[lines_major1 > begin_line & lines_major1 < end_line]
  
  a <- list()
  for(k in 1:4){
    thestring1 <- trimws(neotrop_pages[[i]][lines_suitable + (k-3)])
    thestring2 <- stringr::str_remove(thestring1, "<td>")
    thestring3 <- stringr::str_remove(thestring2, "</td>")
    
    thestring4 <- trimws(neotrop_pages[[i]][lines_major + (k-3)])
    thestring5 <- stringr::str_remove(thestring4, "<td>")
    thestring6 <- stringr::str_remove(thestring5, "</td>")
    
    a[[k]] <- c(thestring3, thestring6)
  }
  
  habitats[[i]] <- paste(a[[1]], a[[2]], a[[3]], a[[4]], sep = "; ")
}

altitude <- as.list(rep(NA, length(nf_species)))
names(altitude) <- nf_species
for(i in 1:length(nf_species)){
  end_line <- grep('Occasional altitudinal limits', neotrop_pages[[i]])
  lines_altitude1 <- grep('Altitude', neotrop_pages[[i]])
  line_altitude <- max(lines_altitude1[lines_altitude1 < end_line])
  thestring1 <- trimws(neotrop_pages[[i]][line_altitude + 2])
  thestring2 <- stringr::str_remove(thestring1, "<td>")
  thestring3 <- stringr::str_remove(thestring2, '</td>')
  a <- strsplit(thestring3, "-")
  a1 <- as.numeric(gsub("\\D", "", a[[1]][1]))
  a2 <- as.numeric(gsub("\\D", "", a[[1]][2]))
  altitude[[i]] <- c(a1, a2)
}


generation <- as.list(rep(NA, length(nf_species)))
names(generation) <- nf_species
for(i in 1:length(nf_species)){
  line_gen <- grep('Generation length \\(yrs\\)', neotrop_pages[[i]])
  if(length(line_gen) > 1){stop()}
  thestring1 <- trimws(neotrop_pages[[i]][line_gen + 1])
  thestring2 <- stringr::str_remove(thestring1, "<td>")
  thestring3 <- stringr::str_remove(thestring2, '</td>')
  generation[[i]] <- as.numeric(thestring3)
}

forest_dep <- as.list(rep(NA, length(nf_species)))
names(forest_dep) <- nf_species
for(i in 1:length(nf_species)){
  line_forest <- grep('Forest dependency', neotrop_pages[[i]])
  if(length(line_forest) != 1){stop()}
  thestring1 <- trimws(neotrop_pages[[i]][line_forest + 1])
  thestring2 <- stringr::str_remove(thestring1, "<td>")
  thestring3 <- stringr::str_remove(thestring2, '</td>')
  forest_dep[[i]] <- thestring3
}

migratory_status <- as.list(rep(NA, length(nf_species)))
names(migratory_status) <- nf_species
for(i in 1:length(nf_species)){
  line_migstat <- grep('Migratory status', neotrop_pages[[i]])
  if(length(line_migstat) != 1){stop()}
  thestring1 <- trimws(neotrop_pages[[i]][line_migstat + 1])
  thestring2 <- stringr::str_remove(thestring1, "<td>")
  thestring3 <- stringr::str_remove(thestring2, '</td>')
  migratory_status[[i]] <- thestring3
}


body_mass <- as.list(rep(NA, length(nf_species)))
names(body_mass) <- nf_species
for(i in 1:length(nf_species)){
  line_avgmass <- grep('Average mass', neotrop_pages[[i]])
  if(length(line_avgmass) != 1){stop()}
  thestring1 <- trimws(neotrop_pages[[i]][line_avgmass + 1])
  thestring2 <- stringr::str_remove(thestring1, "<td>")
  thestring3 <- stringr::str_remove(thestring2, '</td>')
  thestring4 <- stringr::str_remove(thestring3, ' g')
  if(thestring3 != "-"){body_mass[[i]] <- as.numeric(thestring4)}
}


birdlife_traits <- list(names = nf_species, habitats = habitats, altitude = altitude,
                        generation = generation, forest_dep = forest_dep, 
                        migratory_status = migratory_status, body_mass = body_mass)
save(birdlife_traits, file = '/Users/jacobsocolar/Dropbox/Work/Useful_data/BirdlifeTraits/birdlife_traits.Rdata')
