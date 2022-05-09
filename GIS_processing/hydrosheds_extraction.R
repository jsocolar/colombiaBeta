socolar.desktop <- file.exists('/Users/jacobsocolar/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
socolar.laptop <- file.exists('/Users/jacob/Dropbox/Work/Code/code_keychain/machine_identifier_n5L8paM.txt')
if(socolar.desktop){
  dir.path <- "/Users/JacobSocolar/Dropbox"
}else if(socolar.laptop){
  dir.path <- "/Users/jacob/Dropbox"
}


library(sf)

hs2 <- st_read(paste0(dir.path, '/Work/Colombia/Data/GIS/hydrosheds/hybas_sa_lev01-06_v1c/hybas_sa_lev02_v1c.shp'))
hs3 <- st_read(paste0(dir.path, '/Work/Colombia/Data/GIS/hydrosheds/hybas_sa_lev01-06_v1c/hybas_sa_lev03_v1c.shp'))
hs4 <- st_read(paste0(dir.path, '/Work/Colombia/Data/GIS/hydrosheds/hybas_sa_lev01-06_v1c/hybas_sa_lev04_v1c.shp'))
hs5 <- st_read(paste0(dir.path, '/Work/Colombia/Data/GIS/hydrosheds/hybas_sa_lev01-06_v1c/hybas_sa_lev05_v1c.shp'))

amazon_orinoco <- st_make_valid(st_union(st_make_valid(hs2[2, ]), st_make_valid(hs3[4, ])))
pacific_prelim <- st_make_valid(st_union(st_make_valid(hs3[1,]), st_make_valid(hs3[25,])))
gen_magdalena <- st_make_valid(hs3[2, ])
cauca <- st_make_valid(hs5[5, ])
valledupar <- st_make_valid(hs5[8, ])
magdalena <- st_make_valid(st_difference(gen_magdalena, st_make_valid(st_union(cauca, valledupar))))[,c('HYBAS_ID', 'geometry')]
catatumbo <- st_make_valid(hs5[12, ])
maracaibo <- st_make_valid(hs3[3, ])
guajira_valledupar <- st_make_valid(st_union(valledupar, st_difference(maracaibo, catatumbo))[c('HYBAS_ID', 'geometry')])

central <- st_make_valid(st_zm(st_read(paste0(dir.path, '/Work/Colombia/Data/GIS/biogeographic_clips/mountains_clips/centralAndes__cauca__magdalena.kml'))))
st_is_valid(central) # see https://github.com/r-spatial/sf/issues/1771
sf_use_s2(FALSE)
central <- st_make_valid(central)
sf_use_s2(TRUE)
st_is_valid(central)

sm <- st_make_valid(st_zm(st_read(paste0(dir.path, '/Work/Colombia/Data/GIS/biogeographic_clips/mountains_clips/SNSM__guajira_valledupar.kml'))))
pasto <- st_make_valid(st_zm(st_read(paste0(dir.path, '/Work/Colombia/Data/GIS/biogeographic_clips/mountains_clips/pasto__pacific.kml'))))
tacarcuna <- st_make_valid(st_zm(st_read(paste0(dir.path, '/Work/Colombia/Data/GIS/biogeographic_clips/mountains_clips/tacarcuna__pacific.kml'))))


pacific_prelim2 <- st_make_valid(st_difference(pacific_prelim, pasto))
pacific <- st_make_valid(st_difference(pacific_prelim2, tacarcuna))
pasto <- st_make_valid(st_intersection(pacific_prelim, pasto))
tacarcuna <- st_make_valid(st_intersection(pacific_prelim, tacarcuna))

cauca_west <- st_make_valid(st_difference(cauca, central))
cauca_east <- st_make_valid(st_intersection(cauca, central))

magdalena_west <- st_make_valid(st_intersection(magdalena, central))
magdalena_east <- st_make_valid(st_difference(magdalena, central))

guajira_perija <- st_make_valid(st_difference(guajira_valledupar, sm))
snsm <- st_make_valid(st_intersection(guajira_valledupar, sm))

