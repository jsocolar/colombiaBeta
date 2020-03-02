`%ni%` <- Negate(`%in%`)
setwd("/Users/JacobSocolar/Dropbox/Work/Colombia")

# read in united eBird/HBW taxonomy
taxonomy <- read.csv("Data/Birds/HBW_eBird_taxonomy.csv")
t2 <- taxonomy[taxonomy$HBW_CAT == "sp" | taxonomy$CLEM_CAT_2019 == "sp", ] # remove taxa that are not treated as species by either eBird or HBW

##### Western Andes #####
wandes1 <- read.csv("Data/Birds/James_WAndes_all_birds.csv")
wandes <- wandes1[wandes1$Distance %in% c('A', 'B', 'C'), ]
wandes_pts <- read.csv("Data/Points/James_andes_points.csv")
unique(wandes_pts$Habcode)
goodpts <- wandes_pts$Point[wandes_pts$Habcode != "Sy"]

wandes$Species <- as.character(wandes$Species)
wandes$Species <- gsub("_", " ", wandes$Species)

unique(wandes$Species[wandes$Species %ni% t2$HBW_LATIN & wandes$Species %ni% t2$CLEM_SCI_2019])
p1 <- unique(wandes$Species[wandes$Species %ni% t2$HBW_LATIN])
p2 <- unique(wandes$Species[wandes$Species %ni% t2$CLEM_SCI_2019])
p1[p1 %ni% p2]
p2[p2 %ni% p1]

# Most differences correspond to different names for the same taxon.  However, Chlorostilbon melanorhynchus is lumped 
# with C. mellisugus in HBW/BirdLife. Here we will align with eBird:

wandes$Species[wandes$Species == "Myadestes ralliodes"] <- "Myadestes ralloides"
wandes$Species[wandes$Species == "Basileuterus coronatus"] <- "Myiothlypis coronata"
wandes$Species[wandes$Species == "Carduelis xanthogastra"] <- "Spinus xanthogastrus"
wandes$Species[wandes$Species == "Tangara heinei"] <- "Stilpnia heinei"
wandes$Species[wandes$Species == "Phyllomyias nigricapilla"] <- "Phyllomyias nigrocapillus"
wandes$Species[wandes$Species == "Tyrranus melancholicus"] <- "Tyrannus melancholicus"
wandes$Species[wandes$Species == "Lepidicolaptes lacrymiger"] <- "Lepidocolaptes lacrymiger"
wandes$Species[wandes$Species == "Momotus aequitorialis"] <- "Momotus aequatorialis"
wandes$Species[wandes$Species == "Tangara labradoides"] <- "Tangara labradorides"
wandes$Species[wandes$Species == "Phyllomyias nigricapilla"] <- "Phyllomyias nigrocapillus"
wandes$Species[wandes$Species == "Pseudocollaptes boissonneautii"] <- "Pseudocolaptes boissonneautii"
wandes$Species[wandes$Species == "Calliphlox mulsant"] <- "Chaetocercus mulsant"
wandes$Species[wandes$Species == "Tangara vitriolina"] <- "Stilpnia vitriolina"
wandes$Species[wandes$Species == "Ocreatus underwoodi"] <- "Ocreatus underwoodii"
wandes$Species[wandes$Species == "Pipreola riefferi"] <- "Pipreola riefferii"
wandes$Species[wandes$Species == "Calliphlox mitchelli"] <- "Philodice mitchellii"
wandes$Species[wandes$Species == "Picoides fumigatus"] <- "Dryobates fumigatus"
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
wandes$Species[wandes$Species == "Veniliornis dignus"] <- "Dryobates dignus"
wandes$Species[wandes$Species == "Rupicola peruviana"] <- "Rupicola peruvianus"
wandes$Species[wandes$Species == "Ciphorhinus thoracicus"] <- "Cyphorhinus thoracicus"
wandes$Species[wandes$Species == "Basileuterus luteoviridis"] <- "Myiothlypis luteoviridis"
wandes$Species[wandes$Species == "Phylloscartes opthalmicus"] <- "Phylloscartes ophthalmicus"
wandes$Species[wandes$Species == "Calliphlox mulsant"] <- "Chaetocercus mulsant"
wandes$Species[wandes$Species == "Tangara ruficervix"] <- "Chalcothraupis ruficervix"
wandes$Species[wandes$Species == "Tangara artus"] <- "Tangara arthus"
wandes$Species[wandes$Species == "Chlorospingus opthalmicus"] <- "Chlorospingus flavopectus"
wandes$Species[wandes$Species == "Phyllosacrtes poecilotis"] <- "Phylloscartes poecilotis"
wandes$Species[wandes$Species == "Hemispingus frontalis"] <- "Sphenopsis frontalis"
wandes$Species[wandes$Species == "Carduelis psaltria"] <- "Spinus psaltria"
wandes$Species[wandes$Species == "Premnornis guttuligera"] <- "Premnornis guttuliger"
wandes$Species[wandes$Species == "Myiarchus tubericulifer"] <- "Myiarchus tuberculifer"
wandes$Species[wandes$Species == "Myiozetetes simillis"] <- "Myiozetetes similis"
wandes$Species[wandes$Species == "Oryzoborus funereus"] <- "Sporophila funerea"
wandes$Species[wandes$Species == "Parula pitiayumi"] <- "Setophaga pitiayumi"
wandes$Species[wandes$Species == "Leiothlypis peregrinea"] <- "Leiothlypis peregrina"
wandes$Species[wandes$Species == "Tangara cyanicollis"] <- "Stilpnia cyanicollis"
wandes$Species[wandes$Species == "Tangara larvata"] <- "Stilpnia larvata"
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
wandes$Species[wandes$Species == "Tangara rufigula"] <- "Ixothraupis rufigula"
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
wandes$Species[wandes$Species == "Colibri thalassinus"] <- "Colibri cyanotus"
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
llanos <- read.csv("Data/Birds/James_llanos_all_birds.csv")
llanos$Species <- as.character(llanos$Species)

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
llanos$Species[llanos$Species == "Tangara cayana"] <- "Stilpnia cayana"
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
llanos$Species[llanos$Species == "Veniliornis passerinus"] <- "Dryobates passerinus"
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
llanos <- llanos[llanos$Species != "Vireo olivaceus", ] # Records are in February and March, when there's no way to separate true olivaceus from resident chivi after the fact

llanos1 <- llanos
llanos2 <- llanos[llanos$HAB != "PALM",]
llanos <- llanos2[llanos2$Flying. == "N", ]

L1sp <- unique(llanos1$Species)
L2sp <- unique(llanos2$Species)
L3sp <- unique(llanos$Species)

L1sp[L1sp %ni% L2sp]
L2sp[L2sp %ni% L3sp]

##### Jacob #####
jacob1 <- read.csv("Data/Birds/Jacob_data_v1.1.csv")
jacob2 <- jacob1[jacob1$Dist %in% c('A', 'B', 'C') & is.na(jacob1$Dis) & is.na(jacob1$FO), ]
jacob <- jacob2[jacob2$Species %ni% c("Sono", "Visu", "", "Henicorhina_sp"),]

##### Tools to monitor progress on dataset as it's entered--not part of final analysis #####
# 231 Llanos spp
# 319 WAndes spp
# 746 Jacob spp
# 243 EAndes spp (Simon+David only)
# 952 total spp

dim(jacob)
dim(jacob2)[1] - dim(jacob)[1]

jacobSpp <- as.character(unique(jacob$Species)[order(unique(jacob$Species))])
jacobSpp[gsub("_", " ", jacobSpp) %ni% t2$CLEM_SCI_2019]
jacobSpp

simon <- read.csv("Data/Birds/Simon_list_28-02-2019.csv", stringsAsFactors = F)
simon_species <- simon$English.Ebird
for(i in 1:length(simon_species)){
  if(length(grep(simon_species[i], t2$CLEM_ENG_2019, ignore.case = T)) == 0){
    print(c(i, simon_species[i]))
  }
}
simon_species[229] <- simon$English.Ebird[229] <- "White-banded tyrannulet"
simon$scientific <- NA

for(i in 1:nrow(simon)){
  simon$scientific[i] <- as.character(t2$CLEM_SCI_2019[grep(paste0("^",simon_species[i],"$"), t2$CLEM_ENG_2019, ignore.case = T)])
}

simonSpp <- gsub(" ", "_", simon$scientific)
simonSpp[simonSpp == "Grallaria_quitensis"] <- "Grallaria_alticola"
LlanosSpp <- gsub(" ", "_", unique(llanos$Species))
WAndesSpp <- gsub(" ", "_", unique(wandes$Species))

allSpp <- unique(c(jacobSpp, simonSpp, LlanosSpp, WAndesSpp))
allSpp <- allSpp[-which(allSpp %in% c('Henicorhina_bangsi', 'Grallaria_rufula_spatiator'))]
allSpp
#952

allSpp[grep("Ammod", allSpp)]


jacobSpp[grep("wagleri", jacobSpp)]
simonSpp[grep("Schistes_", simonSpp)]
WAndesSpp[grep("Schistes_", WAndesSpp)]


##### Additional tools to monitor my portion of the dataset
for(i in 1:length(jacobSpp)){
  ti <- agrep(jacobSpp[i], jacobSpp)
  if(length(ti) > 1){
    print(jacobSpp[ti])
  }
}

jacobSpp[grep("Schistes_", jacobSpp)]

sum(jacob$Species == "Pachysylvia_semibrunnea")

v <- vector()
for(i in 1:length(jacobSpp)){
  v[i] <- sum(jacob$Species == jacobSpp[i])
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

sum(v<5)
jacobSpp[which(v>=5)]


