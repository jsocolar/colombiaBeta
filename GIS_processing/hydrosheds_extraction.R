library(sf)

hs2 <- st_read('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/hydrosheds/hybas_sa_lev01-06_v1c/hybas_sa_lev02_v1c.shp')
hs3 <- st_read('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/hydrosheds/hybas_sa_lev01-06_v1c/hybas_sa_lev03_v1c.shp')
hs4 <- st_read('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/hydrosheds/hybas_sa_lev01-06_v1c/hybas_sa_lev04_v1c.shp')
hs5 <- st_read('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/hydrosheds/hybas_sa_lev01-06_v1c/hybas_sa_lev05_v1c.shp')

amazon_orinoco <- st_make_valid(st_union(hs2[2, ], hs3[4, ]))
pacific_prelim <- st_make_valid(st_union(hs3[1,], hs3[25,]))
gen_magdalena <- st_make_valid(hs3[2, ])
cauca <- st_make_valid(hs5[5, ])
valledupar <- st_make_valid(hs5[8, ])
magdalena <- st_make_valid(st_difference(gen_magdalena, st_union(cauca, valledupar)))[,c('HYBAS_ID', 'geometry')]
catatumbo <- st_make_valid(hs5[12, ])
maracaibo <- st_make_valid(hs3[3, ])
guajira_valledupar <- st_make_valid(st_union(valledupar, st_difference(maracaibo, catatumbo))[c('HYBAS_ID', 'geometry')])

central <- st_make_valid(st_zm(st_read('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/biogeographic_clips/mountains_clips/centralAndes__cauca__magdalena.kml')))
sm <- st_make_valid(st_zm(st_read('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/biogeographic_clips/mountains_clips/SNSM__guajira_valledupar.kml')))
pasto <- st_make_valid(st_zm(st_read('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/biogeographic_clips/mountains_clips/pasto__pacific.kml')))

pacific <- st_make_valid(st_difference(pacific_prelim, pasto))
pasto <- st_make_valid(st_intersection(pacific_prelim, pasto))

cauca_west <- st_make_valid(st_difference(cauca, central))
cauca_east <- st_make_valid(st_intersection(cauca, central))

magdalena_west <- st_make_valid(st_intersection(magdalena, central))
magdalena_east <- st_make_valid(st_difference(magdalena, central))

guajira_perija <- st_make_valid(st_difference(guajira_valledupar, sm))
snsm <- st_make_valid(st_intersection(guajira_valledupar, sm))


# setwd('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/biogeographic_clips/clipping_polygons')
# st_write(amazon_orinoco, "amazon_orinoco.shp")
# st_write(pacific, "pacific.shp")
# st_write(pasto, 'pasto.shp')
# st_write(cauca_west, 'cauca_west.shp')
# st_write(cauca_east, 'cauca_east.shp')
# st_write(magdalena_west, 'magdalena_west.shp')
# st_write(magdalena_east, "magdalena_east.shp")
# st_write(guajira_perija, "guajira_perija.shp")
# st_write(snsm, 'snsm.shp')
# st_write(catatumbo, 'catatumbo.shp')

