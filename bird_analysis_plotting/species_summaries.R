fit1_data <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/fit1_data.RDS")
fit1 <- as.data.frame(fit1_data)

birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_data_trimmed.RDS")
birds$sp_cl <- paste(birds$species, birds$cluster, sep = "__")
birds$elev_median <- rowMeans(cbind(birds$lower, birds$upper))

det_data <- as.matrix(birds[,c("v1", "v2", "v3", "v4")])
det_data[is.na(det_data)] <- -1

time <- matrix((scale(c(birds$hps1, birds$hps2, birds$hps3, birds$hps4))), ncol = 4)
time[is.na(time)] <- 0

obsSM <- matrix(c(birds$obs1 == "SM", birds$obs2 == "SM", birds$obs3 == "SM", birds$obs4 == "SM"), ncol = 4)
obsSM[is.na(obsSM)] <- 0

obsDE <- matrix(c(birds$obs1 == "DE", birds$obs2 == "DE", birds$obs3 == "DE", birds$obs4 == "DE"), ncol = 4)
obsDE[is.na(obsDE)] <- 0

obsJG <- matrix(c(birds$obs1 == "JG", birds$obs2 == "JG", birds$obs3 == "JG", birds$obs4 == "JG"), ncol = 4)
obsJG[is.na(obsJG)] <- 0

# hack for now, to change later
birds$forest_present[is.na(birds$forest_present)] <- 0
birds$forest_specialist[is.na(birds$forest_specialist)] <- 0
birds$tf_specialist[is.na(birds$tf_specialist)] <- 0
birds$dry_forest_present[is.na(birds$dry_forest_present)] <- 0
birds$flood_dry_specialist[is.na(birds$flood_dry_specialist)] <- 0
birds$floodplain_specialist[is.na(birds$floodplain_specialist)] <- 0
birds$arid_present[is.na(birds$arid_present)] <- 0
birds$Family[is.na(birds$Family)] <- "Grallariidae"
birds$migrat_birdlife[is.na(birds$migrat_birdlife)] <- "not a migrant"

birds$id_spCl <- as.integer(as.factor(birds$sp_cl))
birds$id_sp <- as.integer(as.factor(birds$species))
birds$id_fam <- as.integer(as.factor(birds$Family))

birds$relev <- (birds$elev_sp_standard -.5)/sd(birds$elev_sp_standard)
birds$relev2 <- ((birds$elev_sp_standard -.5)/sd(birds$elev_sp_standard))^2
birds$lowland <- as.numeric(birds$lower == 0)

birds$elevMedian <- vscale(birds$elev_median)
birds$elevBreadth <- vscale(birds$elev_breadth)

birds$migratory <- as.numeric(birds$migrat_birdlife == "full migrant")           ## This needs to change for final version
birds$mass <- vscale(birds$BodyMass.Value)
birds$dietInvert <- as.numeric(birds$Diet.5Cat == "Invertebrate")
birds$dietCarn <- as.numeric(birds$Diet.5Cat == "VertFishScav")
birds$dietFruitNect <- as.numeric(birds$Diet.5Cat == "FruiNect")
birds$dietGran <- as.numeric(birds$Diet.5Cat == "PlantSeed")

traitdata <- birds[!duplicated(birds$species),]


forest_int <- pasture_offset <- matrix(data = NA, nrow = 1614, ncol = 500)

for(sp in 1:1614){
  forest_int[sp,] <- t(fit1[paste0('1.b0_sp[',sp,']')] + fit1[paste0('1.b0_fam[',traitdata$id_fam[sp],']')] +
                                 fit1['1.b3_eastOnly']*traitdata$east_only[sp] + fit1['1.b3_westOnly']*traitdata$west_only[sp] + fit1['1.b3_snsmOnly']*traitdata$snsm_only[sp] +
                                 fit1['1.b3_notWandes']*traitdata$wandes_absent[sp] + fit1['1.b3_notEandes']*traitdata$eandes_absent[sp] + 
                                 
                                 fit1['1.b3_elevMedian']*traitdata$elevMedian[sp] + fit1['1.b3_elevBreadth']*traitdata$elevBreadth[sp] +
                                 
                                 fit1['1.b3_forestPresent']*traitdata$forest_present[sp] + fit1['1.b3_forestSpecialist']*traitdata$forest_specialist[sp] + 
                                 fit1['1.b3_tfSpecialist']*traitdata$tf_specialist[sp] + fit1['1.b3_dryForestPresent']*traitdata$dry_forest_present[sp] +
                                 fit1['1.b3_floodDrySpecialist']*traitdata$flood_dry_specialist[sp] + fit1['1.b3_floodSpecialist']*traitdata$floodplain_specialist[sp] +
                                 fit1['1.b3_aridPresent']*traitdata$arid_present[sp] +
                                 
                                 fit1['1.b3_migratory']*traitdata$migratory[sp] + fit1['1.b3_mass']*traitdata$mass[sp] +
                                 
                                 fit1['1.b3_dietInvert']*traitdata$dietInvert[sp] + fit1['1.b3_dietCarn']*traitdata$dietCarn[sp] + 
                                 fit1['1.b3_dietFruitNect']*traitdata$dietFruitNect[sp] + fit1['1.b3_dietGran']*traitdata$dietGran[sp] +
                                 
                                 fit1['1.b3_x_elevMedian_forestPresent']*traitdata$elevMedian[sp]*traitdata$forest_present[sp]+
                                 fit1['1.b3_x_elevMedian_forestSpecialist']*traitdata$elevMedian[sp]*traitdata$forest_specialist[sp])
  
  
  pasture_offset[sp,] <- t(fit1[paste0('1.b2_pasture_sp[',sp,']')] + fit1[paste0('1.b2_pasture_fam[',traitdata$id_fam[sp],']')] +
                                     fit1['1.b4_eastOnly']*traitdata$east_only[sp] + fit1['1.b4_westOnly']*traitdata$west_only[sp] + fit1['1.b4_snsmOnly']*traitdata$snsm_only[sp] +
                                     fit1['1.b4_notWandes']*traitdata$wandes_absent[sp] + fit1['1.b4_notEandes']*traitdata$eandes_absent[sp] + 
                                     
                                     fit1['1.b4_elevMedian']*traitdata$elevMedian[sp] + fit1['1.b4_elevBreadth']*traitdata$elevBreadth[sp] +
                                     
                                     fit1['1.b4_forestPresent']*traitdata$forest_present[sp] + fit1['1.b4_forestSpecialist']*traitdata$forest_specialist[sp] + 
                                     fit1['1.b4_tfSpecialist']*traitdata$tf_specialist[sp] + fit1['1.b4_dryForestPresent']*traitdata$dry_forest_present[sp] +
                                     fit1['1.b4_floodDrySpecialist']*traitdata$flood_dry_specialist[sp] + fit1['1.b4_floodSpecialist']*traitdata$floodplain_specialist[sp] +
                                     fit1['1.b4_aridPresent']*traitdata$arid_present[sp] +
                                     
                                     fit1['1.b4_migratory']*traitdata$migratory[sp] + fit1['1.b4_mass']*traitdata$mass[sp] +
                                     
                                     fit1['1.b4_dietInvert']*traitdata$dietInvert[sp] + fit1['1.b4_dietCarn']*traitdata$dietCarn[sp] + 
                                     fit1['1.b4_dietFruitNect']*traitdata$dietFruitNect[sp] + fit1['1.b4_dietGran']*traitdata$dietGran[sp] +
                                     
                                     fit1['1.b4_x_elevMedian_forestPresent']*traitdata$elevMedian[sp]*traitdata$forest_present[sp]+
                                     fit1['1.b4_x_elevMedian_forestSpecialist']*traitdata$elevMedian[sp]*traitdata$forest_specialist[sp])
}

elev1_coef <- elev2_coef <- matrix(data = NA, nrow = 1614, ncol = 500)

for(sp in 1:1614){
  elev1_coef[sp, ] <- t(fit1[paste0('1.b1_relev_sp[',sp,']')])
  elev2_coef[sp, ] <- t(fit1[paste0('1.b1_relev2_sp[',sp,']')])
}

species_results_prelim <- list(forest_int = forest_int, pasture_offset = pasture_offset, elev1_coef = elev1_coef, elev2_coef = elev2_coef)

saveRDS(species_results_prelim, '/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/species_results_prelim.RDS')
  
pasture_CIs <- as.data.frame(cbind(traitdata$species, 
                     apply(pasture_offset, 1, function(x) quantile(x, .05)),
                     apply(pasture_offset, 1, function(x) quantile(x, .95)),
                     apply(pasture_offset, 1, function(x) quantile(x, .1)),
                     apply(pasture_offset, 1, function(x) quantile(x, .9)),
                     apply(pasture_offset, 1, function(x) quantile(x, .025)),
                     apply(pasture_offset, 1, function(x) quantile(x, .975))))
names(pasture_CIs) <- c('species', 'q05', 'q95', 'q10', 'q90', 'q025', 'q975')
pasture_CIs$species[pasture_CIs$q05 > 0]
pasture_CIs$species[pasture_CIs$q95 < 0]

pasture_CIs$species[pasture_CIs$q10 > 0]
pasture_CIs$species[pasture_CIs$q90 < 0]

pasture_CIs$species[pasture_CIs$q025 > 0]
pasture_CIs$species[pasture_CIs$q975 < 0]


pasture_CIs$species[pasture_CIs$q10 < 0 & pasture_CIs$q90 > 0]

#  b1_relev_sp[id_sp[r0+r]]*relev[r0+r] + b1_relev2_sp[id_sp[r0+r]]*relev2[r0+r] +





family_lookup <- traitdata[!duplicated(traitdata$Family),c("Family", "id_fam")]
family_lookup <- family_lookup[order(family_lookup$id_fam),]
fit1_summary <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/fit1_summary.RDS")

View(cbind(fit1_summary[grep('b2_pasture_fam\\[', fit1_summary$variable),], family_lookup$Family))