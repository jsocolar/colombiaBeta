fit1511_data <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/fit1511_data.RDS")
fit1511 <- as.data.frame(fit1511_data)

birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/birds.RDS")

birds_stan <- birds

birds_stan$id_spCl <- as.numeric(as.factor(birds$sp_cl))
birds_stan$id_sp <- as.numeric(as.factor(birds$species))
birds_stan$id_fam <- as.numeric(as.factor(birds$Family))
birds_stan$migratory <- as.numeric(!is.na(birds$start1))
birds_stan$dietInvert <- as.numeric(birds$Diet.5Cat == "Invertebrate")
birds_stan$dietCarn <- as.numeric(birds$Diet.5Cat == "VertFishScav")
birds_stan$dietFruitNect <- as.numeric(birds$Diet.5Cat == "FruiNect")
birds_stan$dietGran <- as.numeric(birds$Diet.5Cat == "PlantSeed")

bird_stan_data1_package <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data1_package.RDS")

traitdata <- birds_stan[!duplicated(birds_stan$species),]

forest_int <- pasture_offset <- matrix(data = NA, nrow = 1614, ncol = 500)

for(sp in 1:1614){
  forest_int[sp,] <- t(fit1511[paste0('1.b0_sp[',sp,']')] + fit1511[paste0('1.b0_fam[',traitdata$id_fam[sp],']')] +
                         
                                 fit1511['1.b1_lowland']*traitdata$lowland[sp] +
                         
                                 fit1511['1.b3_eastOnly']*traitdata$east_only[sp] + fit1511['1.b3_westOnly']*traitdata$west_only[sp] + fit1511['1.b3_snsmOnly']*traitdata$snsm_only[sp] +
                                 fit1511['1.b3_notWandes']*traitdata$wandes_absent[sp] + fit1511['1.b3_notEandes']*traitdata$eandes_absent[sp] + 
                                 
                                 fit1511['1.b3_elevMedian']*traitdata$elev_median_scaled[sp] + fit1511['1.b3_elevBreadth']*traitdata$elev_breadth_scaled[sp] +
                                 
                                 fit1511['1.b3_forestPresent']*traitdata$forest_present[sp] + fit1511['1.b3_forestSpecialist']*traitdata$forest_specialist[sp] + 
                                 fit1511['1.b3_tfSpecialist']*traitdata$tf_specialist[sp] + fit1511['1.b3_dryForestPresent']*traitdata$dry_forest_present[sp] +
                                 fit1511['1.b3_floodDrySpecialist']*traitdata$flood_dry_specialist[sp] + fit1511['1.b3_floodSpecialist']*traitdata$floodplain_specialist[sp] +
                                 fit1511['1.b3_aridPresent']*traitdata$arid_present[sp] +
                                 
                                 fit1511['1.b3_migratory']*traitdata$migratory[sp] + fit1511['1.b3_mass']*traitdata$log_mass_scaled[sp] +
                                 
                                 fit1511['1.b3_dietInvert']*traitdata$dietInvert[sp] + fit1511['1.b3_dietCarn']*traitdata$dietCarn[sp] + 
                                 fit1511['1.b3_dietFruitNect']*traitdata$dietFruitNect[sp] + fit1511['1.b3_dietGran']*traitdata$dietGran[sp] +
                                 
                                 fit1511['1.b3_x_elevMedian_forestPresent']*traitdata$elev_median_scaled[sp]*traitdata$forest_present[sp]+
                                 fit1511['1.b3_x_elevMedian_forestSpecialist']*traitdata$elev_median_scaled[sp]*traitdata$forest_specialist[sp])
  
  
  pasture_offset[sp,] <- t(fit1511[paste0('1.b2_pasture_sp[',sp,']')] + fit1511[paste0('1.b2_pasture_fam[',traitdata$id_fam[sp],']')] +
                                     fit1511['1.b4_eastOnly']*traitdata$east_only[sp] + fit1511['1.b4_westOnly']*traitdata$west_only[sp] + fit1511['1.b4_snsmOnly']*traitdata$snsm_only[sp] +
                                     fit1511['1.b4_notWandes']*traitdata$wandes_absent[sp] + fit1511['1.b4_notEandes']*traitdata$eandes_absent[sp] + 
                                     
                                     fit1511['1.b4_elevMedian']*traitdata$elev_median_scaled[sp] + fit1511['1.b4_elevBreadth']*traitdata$elev_breadth_scaled[sp] +
                                     
                                     fit1511['1.b4_forestPresent']*traitdata$forest_present[sp] + fit1511['1.b4_forestSpecialist']*traitdata$forest_specialist[sp] + 
                                     fit1511['1.b4_tfSpecialist']*traitdata$tf_specialist[sp] + fit1511['1.b4_dryForestPresent']*traitdata$dry_forest_present[sp] +
                                     fit1511['1.b4_floodDrySpecialist']*traitdata$flood_dry_specialist[sp] + fit1511['1.b4_floodSpecialist']*traitdata$floodplain_specialist[sp] +
                                     fit1511['1.b4_aridPresent']*traitdata$arid_present[sp] +
                                     
                                     fit1511['1.b4_migratory']*traitdata$migratory[sp] + fit1511['1.b4_mass']*traitdata$log_mass_scaled[sp] +
                                     
                                     fit1511['1.b4_dietInvert']*traitdata$dietInvert[sp] + fit1511['1.b4_dietCarn']*traitdata$dietCarn[sp] + 
                                     fit1511['1.b4_dietFruitNect']*traitdata$dietFruitNect[sp] + fit1511['1.b4_dietGran']*traitdata$dietGran[sp] +
                                     
                                     fit1511['1.b4_x_elevMedian_forestPresent']*traitdata$elev_median_scaled[sp]*traitdata$forest_present[sp]+
                                     fit1511['1.b4_x_elevMedian_forestSpecialist']*traitdata$elev_median_scaled[sp]*traitdata$forest_specialist[sp])
}

elev1_coef <- elev2_coef <- matrix(data = NA, nrow = 1614, ncol = 500)

for(sp in 1:1614){
  elev1_coef[sp, ] <- t(fit1511[paste0('1.b1_relev_sp[',sp,']')] + fit1511['1.b1_x_lowland_relev']*traitdata$lowland[sp])
  elev2_coef[sp, ] <- t(fit1511[paste0('1.b1_relev2_sp[',sp,']')] + fit1511['1.b1_x_lowland_relev2']*traitdata$lowland[sp])
}

species_results_prelim <- list(species = traitdata$species, forest_int = forest_int, pasture_offset = pasture_offset, elev1_coef = elev1_coef, elev2_coef = elev2_coef)

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
fit1511_summary <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/fit1511_summary.RDS")

View(cbind(fit1511_summary[grep('b2_pasture_fam\\[', fit1511_summary$variable),], family_lookup$Family))