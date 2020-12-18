model_fit <- cmdstanr::read_cmdstan_csv('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/CSVs/occupancyMod_ragged_parallel_v4_threads-202011251544-1-8e9bb7.csv')
model_fit2 <- cmdstanr::read_cmdstan_csv('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/CSVs/occupancyMod_ragged_parallel_v4_threads-202011251544-1-8e9bb7.csv', 
                                            variables = model_fit$metadata$stan_variables[c(1:3, 5:87)])
samples <- model_fit$post_warmup_draws
samples2 <- model_fit2$post_warmup_draws

source("/Users/jacobsocolar/Dropbox/Work/Code/stanstuff/parsummarise_draws.R") # parsummarise_draws function

posterior_summary <- parsummarise_draws(samples, 4, 1000)
posterior_summary2 <- parsummarise_draws(samples2, 4, 50)
print(posterior_summary2[posterior_summary2$rhat > 1.05,], n = 30)
print(posterior_summary2[posterior_summary2$ess_bulk < 40,], n = 30)

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

bird_stan_data3_package <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data3_package.RDS")

traitdata <- birds_stan[!duplicated(birds_stan$species),]

forest_int <- pasture_offset <- matrix(data = NA, nrow = 1614, ncol = 500)

for(sp in 1:1614){
  forest_int[sp,] <- t(as.numeric(samples[,,paste0('b0_sp[',sp,']')] + samples[,,paste0('b0_fam[',traitdata$id_fam[sp],']')] +
                         
                                 samples[,,'b1_lowland']*traitdata$lowland[sp] +
                         
                                 samples[,,'b3_eastOnly']*traitdata$east_only[sp] + samples[,,'b3_westOnly']*traitdata$west_only[sp] + samples[,,'b3_snsmOnly']*traitdata$snsm_only[sp] +
                                 samples[,,'b3_notWandes']*traitdata$wandes_absent[sp] + samples[,,'b3_notEandes']*traitdata$eandes_absent[sp] + 
                                 
                                 samples[,,'b3_elevMedian']*traitdata$elev_median_scaled[sp] + samples[,,'b3_elevBreadth']*traitdata$elev_breadth_scaled[sp] +
                                 
                                 samples[,,'b3_forestPresent']*traitdata$forest_present[sp] + samples[,,'b3_forestSpecialist']*traitdata$forest_specialist[sp] + 
                                 samples[,,'b3_tfSpecialist']*traitdata$tf_specialist[sp] + samples[,,'b3_dryForestPresent']*traitdata$dry_forest_present[sp] +
                                 samples[,,'b3_floodDrySpecialist']*traitdata$flood_dry_specialist[sp] + samples[,,'b3_floodSpecialist']*traitdata$floodplain_specialist[sp] +
                                 samples[,,'b3_aridPresent']*traitdata$arid_present[sp] +
                                 
                                 samples[,,'b3_migratory']*traitdata$migratory[sp] + samples[,,'b3_mass']*traitdata$log_mass_scaled[sp] +
                                 
                                 samples[,,'b3_dietInvert']*traitdata$dietInvert[sp] + samples[,,'b3_dietCarn']*traitdata$dietCarn[sp] + 
                                 samples[,,'b3_dietFruitNect']*traitdata$dietFruitNect[sp] + samples[,,'b3_dietGran']*traitdata$dietGran[sp] +
                                 
                                 samples[,,'b3_x_elevMedian_forestPresent']*traitdata$elev_median_scaled[sp]*traitdata$forest_present[sp]+
                                 samples[,,'b3_x_elevMedian_forestSpecialist']*traitdata$elev_median_scaled[sp]*traitdata$forest_specialist[sp]))
  
  
  pasture_offset[sp,] <- t(as.numeric(samples[,,paste0('b2_pasture_sp[',sp,']')] + samples[,,paste0('b2_pasture_fam[',traitdata$id_fam[sp],']')] +
                                     samples[,,'b4_eastOnly']*traitdata$east_only[sp] + samples[,,'b4_westOnly']*traitdata$west_only[sp] + samples[,,'b4_snsmOnly']*traitdata$snsm_only[sp] +
                                     samples[,,'b4_notWandes']*traitdata$wandes_absent[sp] + samples[,,'b4_notEandes']*traitdata$eandes_absent[sp] + 
                                     
                                     samples[,,'b4_elevMedian']*traitdata$elev_median_scaled[sp] + samples[,,'b4_elevBreadth']*traitdata$elev_breadth_scaled[sp] +
                                     
                                     samples[,,'b4_forestPresent']*traitdata$forest_present[sp] + samples[,,'b4_forestSpecialist']*traitdata$forest_specialist[sp] + 
                                     samples[,,'b4_tfSpecialist']*traitdata$tf_specialist[sp] + samples[,,'b4_dryForestPresent']*traitdata$dry_forest_present[sp] +
                                     samples[,,'b4_floodDrySpecialist']*traitdata$flood_dry_specialist[sp] + samples[,,'b4_floodSpecialist']*traitdata$floodplain_specialist[sp] +
                                     samples[,,'b4_aridPresent']*traitdata$arid_present[sp] +
                                     
                                     samples[,,'b4_migratory']*traitdata$migratory[sp] + samples[,,'b4_mass']*traitdata$log_mass_scaled[sp] +
                                     
                                     samples[,,'b4_dietInvert']*traitdata$dietInvert[sp] + samples[,,'b4_dietCarn']*traitdata$dietCarn[sp] + 
                                     samples[,,'b4_dietFruitNect']*traitdata$dietFruitNect[sp] + samples[,,'b4_dietGran']*traitdata$dietGran[sp] +
                                     
                                     samples[,,'b4_x_elevMedian_forestPresent']*traitdata$elev_median_scaled[sp]*traitdata$forest_present[sp]+
                                     samples[,,'b4_x_elevMedian_forestSpecialist']*traitdata$elev_median_scaled[sp]*traitdata$forest_specialist[sp]))
}

elev1_coef <- elev2_coef <- matrix(data = NA, nrow = 1614, ncol = 500)

for(sp in 1:1614){
  elev1_coef[sp, ] <- t(as.numeric(samples[,,paste0('b1_relev_sp[',sp,']')] + samples[,,'b1_x_lowland_relev']*traitdata$lowland[sp]))
  elev2_coef[sp, ] <- t(as.numeric(samples[,,paste0('b1_relev2_sp[',sp,']')] + samples[,,'b1_x_lowland_relev2']*traitdata$lowland[sp]))
}

cluster_sigma <- as.numeric(samples[,,"sigma_b0_spCl"])

species_results_prelim <- list(species = traitdata$species, forest_int = forest_int, pasture_offset = pasture_offset, elev1_coef = elev1_coef, elev2_coef = elev2_coef, cluster_sigma = cluster_sigma)

saveRDS(species_results_prelim, '/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/species_results_prelim.RDS')

species_results_prelim <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/species_results_prelim.RDS')

pasture_offset <- species_results_prelim$pasture_offset

pasture_effects <- as.data.frame(cbind(
                     apply(pasture_offset, 1, function(x) quantile(x, .05)),
                     apply(pasture_offset, 1, function(x) quantile(x, .95)),
                     apply(pasture_offset, 1, function(x) quantile(x, .1)),
                     apply(pasture_offset, 1, function(x) quantile(x, .9)),
                     apply(pasture_offset, 1, function(x) quantile(x, .025)),
                     apply(pasture_offset, 1, function(x) quantile(x, .975)),
                     apply(pasture_offset, 1, function(x) mean(x))))
pasture_effects <- cbind(traitdata$species, pasture_effects)
names(pasture_effects) <- c('species', 'q05', 'q95', 'q10', 'q90', 'q025', 'q975', 'mean')

pasture_effects$species[pasture_effects$q05 > 0]
pasture_effects$species[pasture_effects$q95 < 0]

pasture_effects$species[pasture_effects$q10 > 0]
pasture_effects$species[pasture_effects$q90 < 0]

pasture_effects$species[pasture_effects$q025 > 0]
pasture_effects$species[pasture_effects$q975 < 0]

pasture_effects$species[pasture_effects$q10 < 0 & pasture_effects$q90 > 0]

hist(pasture_effects$mean[pasture_effects$q05 > 0])
hist(pasture_effects$mean[pasture_effects$q95 < 0])

pasture_effects$species[pasture_effects$mean < -10]

##############




family_lookup <- traitdata[!duplicated(traitdata$Family),c("Family", "id_fam")]
family_lookup <- family_lookup[order(family_lookup$id_fam),]
samples_summary <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/samples_summary.RDS")

View(cbind(samples_summary[grep('b2_pasture_fam\\[', samples_summary$variable),], family_lookup$Family))
