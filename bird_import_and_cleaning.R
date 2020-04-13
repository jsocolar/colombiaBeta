`%ni%` <- Negate(`%in%`)

##### For collaborative projects--figure out what machine we're on and automatically set the working directory ####
socolar.desktop <- file.exists('/Users/jacobsocolar/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
socolar.laptop <- file.exists('/Users/jacob/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
if(socolar.desktop){
  dir.path <- "/Users/JacobSocolar/Dropbox/Work/Colombia/Data"
  simon.file.path <- "/Users/jacobsocolar/Google_drive/Simon_data/data/bird data_Jan&Jun2019/data_Jan&Jun2019_currentVersion.xlsx"
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

##### Western Andes #####
wandes1 <- read.csv("Birds/James_WAndes_all_birds.csv")
wandes <- wandes1[wandes1$Distance %in% c('A', 'B', 'C'), ]
wandes_pts <- read.csv("GIS/Points/James_andes_points.csv")
unique(wandes_pts$Habcode)
goodpts <- wandes_pts$Point[wandes_pts$Habcode != "Sy"]

wandes$Species <- as.character(wandes$Species)
wandes$Species <- gsub("_", " ", wandes$Species)

unique(wandes$Species[wandes$Species %ni% t2$HBW_LATIN & wandes$Species %ni% t2$CLEM_SCI_2019])
p1 <- unique(wandes$Species[wandes$Species %ni% t2$HBW_LATIN])
p2 <- unique(wandes$Species[wandes$Species %ni% t2$CLEM_SCI_2019])
p1[p1 %ni% p2]
p2[p2 %ni% p1]

# Align to HBW taxonomy.
# Most differences correspond to different names for the same taxon, but a minority of changes involve splits/lumps
# Only one HBW split potentially has both daughters represented (Ramphcelus flammigerus/icteronotus)
wandes$Species[wandes$Species == "Myadestes ralliodes"] <- "Myadestes ralloides"
wandes$Species[wandes$Species == "Basileuterus coronatus"] <- "Myiothlypis coronata"
wandes$Species[wandes$Species == "Carduelis xanthogastra"] <- "Spinus xanthogastrus"
wandes$Species[wandes$Species == "Phyllomyias nigricapilla"] <- "Phyllomyias nigrocapillus"
wandes$Species[wandes$Species == "Tyrranus melancholicus"] <- "Tyrannus melancholicus"
wandes$Species[wandes$Species == "Lepidicolaptes lacrymiger"] <- "Lepidocolaptes lacrymiger"
wandes$Species[wandes$Species == "Momotus aequitorialis"] <- "Momotus aequatorialis"
wandes$Species[wandes$Species == "Tangara labradoides"] <- "Tangara labradorides"
wandes$Species[wandes$Species == "Phyllomyias nigricapilla"] <- "Phyllomyias nigrocapillus"
wandes$Species[wandes$Species == "Pseudocollaptes boissonneautii"] <- "Pseudocolaptes boissonneautii"
wandes$Species[wandes$Species == "Calliphlox mulsant"] <- "Chaetocercus mulsant"
wandes$Species[wandes$Species == "Ocreatus underwoodi"] <- "Ocreatus underwoodii"
wandes$Species[wandes$Species == "Pipreola riefferi"] <- "Pipreola riefferii"
wandes$Species[wandes$Species == "Calliphlox mitchelli"] <- "Calliphlox mitchellii"
wandes$Species[wandes$Species == "Picoides fumigatus"] <- "Leuconotopicus fumigatus"
wandes$Species[wandes$Species == "Carduelis spinescens"] <- "Spinus spinescens"
wandes$Species[wandes$Species == "Icterus chrysaster"] <- "Icterus chrysater"
wandes$Species[wandes$Species == "Aglaiocercus kingi"] <- "Aglaiocercus kingii"
wandes$Species[wandes$Species == "Geotrygon frenata"] <- "Zentrygon frenata"
wandes$Species[wandes$Species == "Trogoldytes solstitalis"] <- "Troglodytes solstitialis"
wandes$Species[wandes$Species == "Anisognathus lacrimosus"] <- "Anisognathus lacrymosus"
wandes$Species[wandes$Species == "Cinnycerthia sharpei"] <- "Cinnycerthia olivascens"
wandes$Species[wandes$Species == "Oreothraupis arremenops"] <- "Oreothraupis arremonops"
wandes$Species[wandes$Species == "Henicorhina negretti"] <- "Henicorhina negreti"
wandes$Species[wandes$Species == "Cyclaris nigrirostris"] <- "Cyclarhis nigrirostris"
wandes$Species[wandes$Species == "Margorornis stellatus"] <- "Margarornis stellatus"
wandes$Species[wandes$Species == "Chlorornis riefferi"] <- "Chlorornis riefferii"
wandes$Species[wandes$Species == "Amazona mercenaria"] <- "Amazona mercenarius"
wandes$Species[wandes$Species == "Rupicola peruviana"] <- "Rupicola peruvianus"
wandes$Species[wandes$Species == "Ciphorhinus thoracicus"] <- "Cyphorhinus thoracicus"
wandes$Species[wandes$Species == "Basileuterus luteoviridis"] <- "Myiothlypis luteoviridis"
wandes$Species[wandes$Species == "Phylloscartes opthalmicus"] <- "Pogonotriccus ophthalmicus"
wandes$Species[wandes$Species == "Calliphlox mulsant"] <- "Chaetocercus mulsant"
wandes$Species[wandes$Species == "Tangara ruficervix"] <- "Chalcothraupis ruficervix"
wandes$Species[wandes$Species == "Tangara artus"] <- "Tangara arthus"
wandes$Species[wandes$Species == "Chlorospingus opthalmicus"] <- "Chlorospingus flavopectus"
wandes$Species[wandes$Species == "Phyllosacrtes poecilotis"] <- "Pogonotriccus poecilotis"
wandes$Species[wandes$Species == "Hemispingus frontalis"] <- "Sphenopsis frontalis"
wandes$Species[wandes$Species == "Carduelis psaltria"] <- "Spinus psaltria"
wandes$Species[wandes$Species == "Premnornis guttuligera"] <- "Premnornis guttuliger"
wandes$Species[wandes$Species == "Myiarchus tubericulifer"] <- "Myiarchus tuberculifer"
wandes$Species[wandes$Species == "Myiozetetes simillis"] <- "Myiozetetes similis"
wandes$Species[wandes$Species == "Oryzoborus funereus"] <- "Sporophila funerea"
wandes$Species[wandes$Species == "Parula pitiayumi"] <- "Setophaga pitiayumi"
wandes$Species[wandes$Species == "Leiothlypis peregrinea"] <- "Leiothlypis peregrina"
wandes$Species[wandes$Species == "Sturnella militaris"] <- "Leistes militaris"
wandes$Species[wandes$Species == "Atilla spadiceus"] <- "Attila spadiceus"
wandes$Species[wandes$Species == "Cercomacra parkeri"] <- "Cercomacroides parkeri"
wandes$Species[wandes$Species == "Contopus cireneus"] <- "Contopus cinereus"
wandes$Species[wandes$Species == "Thalurannia fannyi"] <- "Thalurania colombica"
wandes$Species[wandes$Species == "Scytalopus altopisones"] <- "Scytalopus alvarezlopezi"
wandes$Species[wandes$Species == "Myrmeciza immaculata"] <- "Hafferia zeledoni"
wandes$Species[wandes$Species == "Leptopogon supercilliaris"] <- "Leptopogon superciliaris"
wandes$Species[wandes$Species == "Basileuterus chrysogaster"] <- "Myiothlypis chrysogaster"
wandes$Species[wandes$Species == "Dendrocincla tyrannia"] <- "Dendrocincla tyrannina"
wandes$Species[wandes$Species == "Hylophilus semibrunneus"] <- "Pachysylvia semibrunnea"
wandes$Species[wandes$Species == "Hemispingus supercilliaris"] <- "Thlypopsis superciliaris"
wandes$Species[wandes$Species == "Terenura callinota"] <- "Euchrepomis callinota"
wandes$Species[wandes$Species == "Pharomarchus antisianus"] <- "Pharomachrus antisianus"
wandes$Species[wandes$Species == "Malacoptila mysticalis"] <- "Malacoptila mystacalis"
wandes$Species[wandes$Species == "Aratinga wagleri"] <- "Psittacara wagleri"
wandes$Species[wandes$Species == "Gliphorynchus spirurus"] <- "Glyphorynchus spirurus"
wandes$Species[wandes$Species == "Urochroa boucheri"] <- "Urochroa bougueri"
wandes$Species[wandes$Species == "Automolus rubiginosus"] <- "Clibanornis rubiginosus"
wandes$Species[wandes$Species == "Melanerpes formicarivorus"] <- "Melanerpes formicivorus"
wandes$Species[wandes$Species == "Cyanocompsa cyanoides"] <- "Cyanoloxia cyanoides"
wandes$Species[wandes$Species == "Aulacorhynchus haematopygas"] <- "Aulacorhynchus haematopygus"
wandes$Species[wandes$Species == "Phaeothlypis fulvicauda"] <- "Myiothlypis fulvicauda"
wandes$Species[wandes$Species == "Hemispingus atropileus"] <- "Kleinothraupis atropileus"
wandes$Species[wandes$Species == "Ortalis colombianus"] <- "Ortalis columbiana"
wandes$Species[wandes$Species == "Leptoptila verreauxi"] <- "Leptotila verreauxi"
wandes$Species[wandes$Species == "Chamaeza mollisima"] <- "Chamaeza mollissima"
wandes$Species[wandes$Species == "Hemispingus melanotis"] <- "Sphenopsis melanotis"
wandes$Species[wandes$Species == "Atticora tobialis"] <- "Atticora tibialis"
wandes$Species[wandes$Species == "Nephilomyias fasciatus"] <- "Myiophobus fasciatus"
wandes$Species[wandes$Species == "Colibri delphinnae"] <- "Colibri delphinae"
wandes$Species[wandes$Species == "Dendroica fusca"] <- "Setophaga fusca"
wandes$Species[wandes$Species == "Cathatres aura"] <- "Cathartes aura"
wandes$Species[wandes$Species == "Dendroica striata"] <- "Setophaga striata"
wandes$Species[wandes$Species == "Wilsonia canadensis"] <- "Cardellina canadensis"
wandes$Species[wandes$Species == "Contopus virenus"] <- "Contopus virens"
wandes$Species[wandes$Species == "Leucopternis princeps"] <- "Morphnarchus princeps"
wandes$Species[wandes$Species == "Buteo albicaudatus"] <- "Geranoaetus albicaudatus"
wandes$Species[wandes$Species == "Buteo magnirostris"] <- "Rupornis magnirostris"
wandes$Species[wandes$Species == "Buteo polysoma"] <- "Geranoaetus polyosoma"
wandes$Species[wandes$Species == "Schistes geoffroyi"] <- "Schistes albogularis"
wandes$Species[wandes$Species == "Tangara parzudakii"] <- "Tangara lunigera"
wandes$Species[wandes$Species == "Patagioenas fasciata"] <- "Patagioenas albilinea"
wandes$Species[wandes$Species == "Aulacorhynchus prasinus"] <- "Aulacorhynchus albivitta"
wandes$Species[wandes$Species == "Pseudocolaptes lawrencii"] <- "Pseudocolaptes johnsoni"
wandes$Species[wandes$Species == "Platyrinchus mystaceus"] <- "Platyrinchus albogularis"
wandes$Species[wandes$Species == "Mionectes olivaceus"] <- "Mionectes galbinus"
wandes$Species[wandes$Species == "Contopus cinereus"] <- "Contopus bogotensis"
wandes$Species[wandes$Species == "Myiodynastes chrysocephalus"] <- "Myiodynastes hemichrysus"
wandes$Species[wandes$Species == "Cyanolyca armillata"] <- "Cyanolyca quindiuna"
wandes$Species[wandes$Species == "Cyphorhinus thoracicus"] <- "Cyphorhinus dichrous"
wandes$Species[wandes$Species == "Cacicus chrysonotus"] <- "Cacicus leucoramphus"
wandes$Species[wandes$Species == "Myiothlypis chrysogaster"] <- "Myiothlypis chlorophrys"
wandes$Species[wandes$Species == "Myioborus ornatus"] <- "Myioborus chrysops"
wandes$Species[wandes$Species == "Piranga flava"] <- "Piranga hepatica"
wandes$Species[wandes$Species == "Chalcothraupis ruficervix"] <- "Tangara ruficervix"
wandes$Species[wandes$Species == "Tangara arthus"] <- "Tangara aurulenta"
wandes$Species[wandes$Species == "Sporophila corvina"] <- "Sporophila ophthalmica"
wandes$Species[wandes$Species == "Chlorostilbon melanorhynchus"] <- "Chlorostilbon mellisugus"
wandes$Species[wandes$Species == "Ochthoeca diadema"] <- "Silvicultrix diadema"
wandes$Species[wandes$Species == "Pseudocolaptes boissonneautii"] <- "Pseudocolaptes boissonneauii"
wandes$Species[wandes$Species == "Thraupis cyanocephala"] <- "Sporathraupis cyanocephala"
wandes$Species[wandes$Species == "Thraupis episcopus"] <- "Tangara episcopus"
wandes$Species[wandes$Species == "Thraupis palmarum"] <- "Tangara palmarum"
wandes$Species[wandes$Species == "Xenops rutilans"] <- "Xenops rutilus"
wandes$Species[wandes$Species == "Atlapetes tricolor"] <- "Atlapetes crassus"
wandes$Species[wandes$Species == "Chlorothraupis stolzmanni"] <- "Habia stolzmanni"

wandes <- wandes[-which(wandes$Species %in% c("Stelgidopteryx ruficollis", "Atticora tibialis", "Hirundo rustica",
                                              "Pygochelidon cyanoleuca", "Orochelidon murina")), ]
# The above species (all Swallows) are not noted as flyovers but it is certain that most records were flyovers.
# JGG confirms that there is no good way to separate the flyovers from the perched birds, so we are discarding
# all swallows from the West Andes data

wandes <- wandes[-which(wandes$Species %in% c("Hummingbird sp", "Chaetura sp", "Unkn.", "Chaetura_cinereiventris", "Vireo_olivaceus")), ]
# The above either are not species level taxa, or were seen exclusively in flight (Chaetura cinereiventris, pers. comm. JGG),
# or in the case of Vireo olivaceus cannot be reliably assigned to olivaceus vs. chivi (pers. comm. JGG).

unique(wandes$Species[wandes$Migrant.or.transient. == "Y"])
wandes <- wandes[-which(wandes$Migrant.or.transient. == "Y" & wandes$Species %in% c('Bubulcus ibis', 'Streptoprocne zonaris',
                                                                                'Coragyps atratus', 'Falco sparverius',
                                                                                'Geranoaetus polyosoma', 'Cathartes aura',
                                                                                'Caracara cheriway', 'Rupornis magnirostris',
                                                                                'Streptoprocne rutila', 'Milvago chimachima',
                                                                                'Elanoides forficatus', 'Morphnarchus princeps',
                                                                                'Geranoaetus albicaudatus', 'Accipiter striatus')), ]
# Remove all species that have this designation because they are flyovers, but not species that have this designation
# because they are migrants


unique(wandes$Species[wandes$Species %ni% t2$HBW_LATIN & wandes$Species %ni% t2$CLEM_SCI_2019])
p1 <- unique(wandes$Species[wandes$Species %ni% t2$HBW_LATIN])
p2 <- unique(wandes$Species[wandes$Species %ni% t2$CLEM_SCI_2019])
p1
p2


dim(wandes)
wandes2 <- wandes
wandes <- wandes[wandes$Point %in% goodpts,]

j1sp <- unique(wandes$Species)
j2sp <- unique(wandes2$Species)

WAndesSpp <- gsub(" ", "_", unique(wandes$Species))


##### Llanos #####
llanos <- read.csv("Birds/James_llanos_all_birds.csv")
llanos$Species <- as.character(llanos$Species)

# Align to HBW taxonomy.
# Most differences correspond to different names for the same taxon, but a minority of changes involve splits/lumps
# It appears that no HBW splits have both daughters represented, but borderline cases include Myiodynastes maculatus/solitarius, 
# Turdus ignobilis/debilis, and Tolmomyias flaviventris/viridiceps.  The Myiodynastes are assumed to all be maculatus
# based on date.  The Turdus are assumed to all be ignobilis based on photos from the region in eBird.  The Tolmomyias
# are assumed to all be flaviventris based on recordings from the region on Xeno-Canto.

llanos$Species[llanos$Species == "Buteo magnirostris"] <- "Rupornis magnirostris"
llanos$Species[llanos$Species == "Penelope jacuae"] <- "Penelope jacquacu"
llanos$Species[llanos$Species == "Trogolodytes aedon"] <- "Troglodytes aedon"
llanos$Species[llanos$Species == "Dendroica striata"] <- "Setophaga striata"
llanos$Species[llanos$Species == "Aratinga pertinax"] <- "Eupsittula pertinax"
llanos$Species[llanos$Species == "Orthopsittaca manilata"] <- "Orthopsittaca manilatus"
llanos$Species[llanos$Species == "Tyrannus melancolchicus"] <- "Tyrannus melancholicus"
llanos$Species[llanos$Species == "Forpus conspiciliatus"] <- "Forpus conspicillatus"
llanos$Species[llanos$Species == "Pitangus sulpuratus"] <- "Pitangus sulphuratus"
llanos$Species[llanos$Species == "Patagioenas cayannensis"] <- "Patagioenas cayennensis"
llanos$Species[llanos$Species == "Myiozetetes cayannensis"] <- "Myiozetetes cayanensis"
llanos$Species[llanos$Species == "Bulbuculus ibis"] <- "Bubulcus ibis"
llanos$Species[llanos$Species == "Phimosus infusciatus"] <- "Phimosus infuscatus"
llanos$Species[llanos$Species == "Setophaga ruticula"] <- "Setophaga ruticilla"
llanos$Species[llanos$Species == "Eurypygia helias"] <- "Eurypyga helias"
llanos$Species[llanos$Species == "Parkesia noveboravenis"] <- "Parkesia noveboracensis"
llanos$Species[llanos$Species == "Megarhynchus pitangua"] <- "Megarynchus pitangua"
llanos$Species[llanos$Species == "Mesembrinibis cayannensis"] <- "Mesembrinibis cayennensis"
llanos$Species[llanos$Species == "Syrigma silbilatrix"] <- "Syrigma sibilatrix"
llanos$Species[llanos$Species == "Leptoptilla rufaxilla"] <- "Leptotila rufaxilla"
llanos$Species[llanos$Species == "Procne tapera"] <- "Progne tapera"
llanos$Species[llanos$Species == "Egretta alba"] <- "Ardea alba"
llanos$Species[llanos$Species == "Sturnella militaris"] <- "Leistes militaris"
llanos$Species[llanos$Species == "Molothrus bonairensis"] <- "Molothrus bonariensis"
llanos$Species[llanos$Species == "Volatina jacinaria"] <- "Volatinia jacarina"
llanos$Species[llanos$Species == "Myrmeciza atrothorax"] <- "Myrmophylax atrothorax"
llanos$Species[llanos$Species == "Caprimulgus maculicaudus"] <- "Hydropsalis maculicaudus"
llanos$Species[llanos$Species == "Phacellodromus rufifrons"] <- "Phacellodomus rufifrons"
llanos$Species[llanos$Species == "Myiarchus tubericulifer"] <- "Myiarchus tuberculifer"
llanos$Species[llanos$Species == "Dendroica aestiva"] <- "Setophaga petechia"
llanos$Species[llanos$Species == "Aramides cayanea"] <- "Aramides cajaneus"
llanos$Species[llanos$Species == "Tringa solitarius"] <- "Tringa solitaria"
llanos$Species[llanos$Species == "Eucometis penicilliata"] <- "Eucometis penicillata"
llanos$Species[llanos$Species == "Helicolestes halmatus"] <- "Helicolestes hamatus"
llanos$Species[llanos$Species == "Pteroglossus puricinctus"] <- "Pteroglossus pluricinctus"
llanos$Species[llanos$Species == "Aratinga leucopthalma"] <- "Psittacara leucophthalmus"
llanos$Species[llanos$Species == "Thamnophilus punctuata"] <- "Thamnophilus punctatus"
llanos$Species[llanos$Species == "Schistocichla leucostigma"] <- "Myrmelastes leucostigma"
llanos$Species[llanos$Species == "Procne chalybea"] <- "Progne chalybea"
llanos$Species[llanos$Species == "Pipra erythrocephala"] <- "Ceratopipra erythrocephala"
llanos$Species[llanos$Species == "Dacnis cyanea"] <- "Dacnis cayana"
llanos$Species[llanos$Species == "Arremonops taciturnus"] <- "Arremon taciturnus"
llanos$Species[llanos$Species == "Campsiempis flaveola"] <- "Capsiempis flaveola"
llanos$Species[llanos$Species == "Pionus menstrurus"] <- "Pionus menstruus"
llanos$Species[llanos$Species == "Theristictus caudatus"] <- "Theristicus caudatus"
llanos$Species[llanos$Species == "Sporophila cinereum"] <- "Sporophila plumbea"
llanos$Species[llanos$Species == "Dendroica petechia"] <- "Setophaga petechia"
llanos$Species[llanos$Species == "Neochen jubata"] <- "Oressochen jubatus"
llanos$Species[llanos$Species == "Oryzoborus angolensis"] <- "Sporophila angolensis"
llanos$Species[llanos$Species == "Eudocimus alba"] <- "Eudocimus albus"
llanos$Species[llanos$Species == "Clavaris pretiosa"] <- "Claravis pretiosa"
llanos$Species[llanos$Species == "Philherodias pileatus"] <- "Pilherodius pileatus"
llanos$Species[llanos$Species == "Tringites subruficollis"] <- "Calidris subruficollis"
llanos$Species[llanos$Species == "Tityra squamiger"] <- "Tityra cayana"
llanos$Species[llanos$Species == "Momotus subrufescens"] <- "Momotus momota"
llanos$Species[llanos$Species == "Pteroglossus inscriptus"] <- "Pteroglossus humboldti"
llanos$Species[llanos$Species == "Ramphastos tucanus"] <- "Ramphastos cuvieri"
llanos$Species[llanos$Species == "Dryocopus lineatus"] <- "Hylatomus lineatus"
llanos$Species[llanos$Species == "Phacellodomus rufifrons"] <- "Phacellodomus inornatus"
llanos$Species[llanos$Species == "Arremon taciturnus"] <- "Arremon axillaris"
llanos$Species[llanos$Species == "Icterus cayanensis"] <- "Icterus chrysocephalus"
llanos$Species[llanos$Species == "Pitangus lictor"] <- "Philohydor lictor"
llanos$Species[llanos$Species == "Tachyphonus luctuosus"] <- "Islerothraupis luctuosa"
llanos$Species[llanos$Species == "Thraupis episcopus"] <- "Tangara episcopus"
llanos$Species[llanos$Species == "Thraupis palmarum"] <- "Tangara palmarum"
llanos$Species[llanos$Species == "Hypnelus ruficollis"] <- "Hypnelus bicinctus"

# We assume that all of the llanos Myiodynastes maculatus are nominate maculatus, given the Jan-March timeframe.
llanos <- llanos[llanos$Species != "Vireo olivaceus", ] # Records are in February and March, when there's no way to separate true olivaceus from resident chivi after the fact


llanos1 <- llanos
llanos2 <- llanos[llanos$HAB != "PALM",]
llanos <- llanos2[llanos2$Flying. == "N", ]

L1sp <- unique(llanos1$Species)
L2sp <- unique(llanos2$Species)
L3sp <- unique(llanos$Species)

L1sp[L1sp %ni% L2sp]
L2sp[L2sp %ni% L3sp]


##### East Andes #####
simon1 <- data.frame(readxl::read_excel(simon.file.path))

# some of the below consists of data checks that won't be necessary for final analysis but are retained for now
# (until the data file is finalized) to guard against typos and entry errors.
all.equal(is.na(simon1$Point), is.na(simon1$Visit)) #  make sure that NAs in jacob1$Take universally match ""'s in jacob1$Point
which(is.na(simon1$Point) & !is.na(simon1$Visit))
which(!is.na(simon1$Point) & is.na(simon1$Visit))
#View(simon1[which(!is.na(simon1$Point) & is.na(simon1$Visit)), ])

simon1$Visit[which(!is.na(simon1$Point) & is.na(simon1$Visit))] <- c(4, 4, 4, 3, 3, 3, 3, 3, 4, 4, 4, 4)
all.equal(is.na(simon1$Point), is.na(simon1$Visit))

unique(simon1$Dist)
length(unique(simon1$Dist)) # make sure that unique doesn't include "" as a 12th item

unique(simon1$FO)
length(unique(simon1$FO)) # make sure that unique doesn't include "" as a 3rd item

unique(simon1$Visit) # Just need to make sure later that 5's are sites genuinely visited 5 times and not typos for 4.

# Fill in point and visit identifiers
for(i in 2:nrow(simon1)){
  if(is.na(simon1$Visit[i])){
    simon1$Take[i] <- simon1$Take[i - 1]
    simon1$Point[i] <- simon1$Point[i - 1]
  }
}

# Read in my lookup table for Simon's dataset
simon_list <- read.csv("Birds/Simon_list_28-02-2019.csv", stringsAsFactors = F)
simon_list$English.Ebird[229] <- "White-banded tyrannulet"
simon_list$English.Simon[simon_list$English.Simon == "Black-collared Jay"] <- "Black-collared jay"

simon_species <- unique(simon1$Species)
simon_species[simon_species %ni% simon_list$English.Simon]

#View(simon1[simon1$Species %ni% c(simon_list$English.Simon, 'Sono', 'Vis', 'Note', '-', '"') & !grepl(" sp$", simon1$Species) & !(is.na(simon1$Species)), ])
# Most of these hits are species that aren't really in the dataset (unknown IDs, flyovers, D-band)
# However, we need to deal with one:
simon1$Species[simon1$Species == 'Smoky bush tyrant'] <- 'Smoky bush-tyrant'

# Extract the analyzeable records, but do so in several steps to implement typo checks and to keep track of 
# how many species/entries are unanalyzeable.
simon2 <- simon1[simon1$Species %in% simon_list$English.Simon, ]
simon2$Dist[is.na(simon2$Dist)] <- 'C'
simon <- droplevels(simon2[is.na(simon2$FO) & simon2$Dist %in% c('A', 'B', 'C', 'C/B', 'B/C', 'A/B'), ])

# The below is ugly ugly ugly and slow slow slow, but it does the job
simon$Latin <- NA
for(i in 1:nrow(simon)){
  simon$Latin[i] <- as.character(t2$CLEM_SCI_2019[grep(paste0("^",simon_list$English.Ebird[which(simon_list$English.Simon == simon$Species[i])],"$"), t2$CLEM_ENG_2019, ignore.case = T)])
}

simon$Latin <- gsub(" ", "_", simon$Latin)

simon$Latin[simon$Latin == "Grallaria_quitensis"] <- "Grallaria_alticola"
simon$Latin[simon$Latin == "Patagioenas_fasciata"] <- "Patagioenas_albilinea"
simon$Latin[simon$Latin == "Dryobates_fumigatus"] <- "Leuconotopicus_fumigatus"
simon$Latin[simon$Latin == "Dryocopus_lineatus"] <- "Hylatomus_lineatus"
simon$Latin[simon$Latin == "Mionectes_olivaceus"] <- "Mionectes_galbinus"
simon$Latin[simon$Latin == "Myiodynastes_chrysocephalus"] <- "Myiodynastes_hemichrysus"
simon$Latin[simon$Latin == "Catharus_ustulatus"] <- "Catharus_swainsoni"
simon$Latin[simon$Latin == "Cacicus_chrysonotus"] <- "Cacicus_leucoramphus"
simon$Latin[simon$Latin == "Anisognathus_igniventris"] <- "Anisognathus_lunulatus"
simon$Latin[simon$Latin == "Tangara_arthus"] <- "Tangara_aurulenta"
simon$Latin[simon$Latin == "Colibri_cyanotus"] <- "Colibri_thalassinus"
simon$Latin[simon$Latin == "Ochthoeca_diadema"] <- "Silvicultrix_diadema"
simon$Latin[simon$Latin == "Ochthoeca_frontalis"] <- "Silvicultrix_frontalis"
simon$Latin[simon$Latin == "Phylloscartes_poecilotis"] <- "Pogonotriccus_poecilotis"
simon$Latin[simon$Latin == "Pseudocolaptes_boissonneautii"] <- "Pseudocolaptes_boissonneauii"
simon$Latin[simon$Latin == "Stilpnia_heinei"] <- "Tangara_heinei"
simon$Latin[simon$Latin == "Stilpnia_cyanicollis"] <- "Tangara_cyanicollis"
simon$Latin[simon$Latin == "Stilpnia_vitriolina"] <- "Tangara_vitriolina"
simon$Latin[simon$Latin == "Thraupis_cyanocephala"] <- "Sporathraupis_cyanocephala"
simon$Latin[simon$Latin == "Thraupis_episcopus"] <- "Tangara_episcopus"
simon$Latin[simon$Latin == "Thraupis_palmarum"] <- "Tangara_palmarum"
simon$Latin[simon$Latin == "Xenops_rutilans"] <- "Xenops_rutilus"


simonSpp <- unique(simon$Latin)

load("Birds/species_list_creation/colombia_species.Rdata")
unique(simonSpp[gsub("_", " ", simonSpp) %ni% colombia_species])


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
jacob2 <- droplevels(jacob1[jacob1$Species != "FLOCK" & is.na(jacob1$no_birds) & is.na(jacob1$Dis), ])
unique(jacob2$Species[gsub("_", " ", jacob2$Species) %ni% colombia_species])
length(unique(jacob2$Species[gsub("_", " ", jacob2$Species) %ni% colombia_species]))
unique(jacob1$Species[which(jacob1$Species %ni% jacob2$Species)])
length(unique(jacob1$Species[which(jacob1$Species %ni% jacob2$Species)]))
unique(jacob2$Dist)
length(unique(jacob2$Dist))

jacob3 <- droplevels(jacob2[is.na(jacob2$FO) & (jacob2$Dist %in% c("A", "B", "C")), ])
unique(jacob2$Species[which(jacob2$Species %ni% jacob3$Species)])
length(unique(jacob2$Species[which(jacob2$Species %ni% jacob3$Species)]))

jacob <- droplevels(jacob3[jacob3$Species %ni% c("Sono", "Visu"), ])


# Get data file of visit-specific information
jacob_visit_data <- jacob1[!is.na(jacob1$Time), names(jacob1) %in% c("Point", "Take", "Date", "Wind", "Vis", "Sky", "Sun", "Drip", "Time")]
jacob_visit_data$Obs <- "JBS"
jacob_visit_data$Obs[jacob_visit_data$Point == IGF]

table(jacob_visit_data$Point)[table(jacob_visit_data$Point) < 4]

##### jacob dataset: change from entry-convenient format to analysis-convenient format #####
jacobACF1 <- jacob[, names(jacob) %in% c("Point", "Take", "Species", "Count")]

jacobACF <- doBy::summaryBy(Count ~ Point + Take + Species, data = jacobACF1, FUN = sum)

jacobACF[1:20,]

jacobPOINTSUMMARY <- doBy::summaryBy(Count ~ Point + Species, data = jacobACF1, FUN = sum)



##### Tools to monitor progress on dataset as it's entered--not part of final analysis #####
# 237 Llanos spp
# 319 WAndes spp
# 247 EAndes spp (Simon+David only)
# 757 Jacob spp
## 281 Amazon
## 557 Jacob outside Amazon
# 607 other spp (all data except Jacob)
# 959 total spp
## 802 total spp outside of Amazon

dim(jacob)
dim(jacob3)[1] - dim(jacob)[1]

jacobSpp <- as.character(unique(jacob$Species)[order(unique(jacob$Species))])
jacobSpp[gsub("_", " ", jacobSpp) %ni% colombia_species]
jacobSpp

##### Combining Species Lists #####
LlanosSpp <- gsub(" ", "_", unique(llanos$Species))
WAndesSpp <- gsub(" ", "_", unique(wandes$Species))

allSpp <- unique(c(jacobSpp, simonSpp, LlanosSpp, WAndesSpp))
allSpp
#961

otherSpp <- unique(c(simonSpp, LlanosSpp, WAndesSpp))
otherSpp
#607

allSpp[gsub("_", " ", allSpp) %ni% colombia_species]

grepphrase <- "Pogono"

allSpp[grep(grepphrase, allSpp)]
jacobSpp[grep(grepphrase, jacobSpp)]
simonSpp[grep(grepphrase, simonSpp)]
WAndesSpp[grep(grepphrase, WAndesSpp)]
LlanosSpp[grep(grepphrase, LlanosSpp)]


##### Additional tools to monitor my portion of the dataset
for(i in 1:length(jacobSpp)){
  ti <- agrep(jacobSpp[i], jacobSpp)
  if(length(ti) > 1){
    print(jacobSpp[ti])
  }
}

jacobSpp[grep("Mecocerculus_s", jacobSpp)]

sum(jacob$Species == "Mecocerculus_stictopterus")

# counting number of point-visits per species
v <- vector()
for(i in 1:length(jacobSpp)){
  v[i] <- sum(jacobACF$Species == jacobSpp[i])
}

hist(v)
jacobSpp[which(v > 50)]
hist(v[v <= 50])
hist(v[v <= 10])
sum(v==1)
sum(v==2)
sum(v==3)
sum(v==4)
sum(v==5)

jacobSpp[which(v<5)]
# 385
jacobSpp[which(v>=5)]


# counting number of points per species
v <- vector()
for(i in 1:length(jacobSpp)){
  v[i] <- sum(jacobPOINTSUMMARY$Species == jacobSpp[i])
}

hist(v)
jacobSpp[which(v > 50)]
hist(v[v <= 50])
hist(v[v <= 10])
sum(v==1)
sum(v==2)
sum(v==3)
sum(v==4)
sum(v==5)

jacobSpp[which(v<5)]
# 430
jacobSpp[which(v>=5)]

sum(jacob$Species %in% jacobSpp[which(v<5)])
sum(jacob$Species %in% jacobSpp[which(v>=5)])





NonAmazonSpp <- as.character(unique(jacob$Species[1:(min(which(jacob$Point == 'SGP12'))-1)]))
NonAmazonSpp

AllNonAmazon <- unique(c(NonAmazonSpp, simonSpp, LlanosSpp, WAndesSpp))
AllNonAmazon



##### For Felicity ######
allSpp # all species currently in the dataset (HBW/BirdLife taxonomy)

AmazonSpp <- unique(jacob$Species[min(which(jacob$Point == 'SGP12')):nrow(jacob)])
AmazonSpp # all species recorded from Amazonian points



