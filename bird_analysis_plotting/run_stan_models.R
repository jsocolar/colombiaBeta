library("cmdstanr"); library("dplyr"); library("posterior")

bird_stan_data1_package <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data1_package.RDS")

# Run ragged model ----
mod_R_1211 <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/stan_files/full_colombia_model/occupancyMod_ragged_parallel_v1.stan",
                       cpp_options = list(stan_threads = TRUE))
samps_R_1211 <- mod_R_1211$sample(data = bird_stan_data1_package$data, 
                          chains = 1,
                          threads_per_chain = 4,
                          refresh = 1,
                          iter_sampling = 500,
                          iter_warmup = 500,
                          save_warmup = 1,
                          output_dir = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/CSVs")

saveRDS(samps_R_1211, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/samps_R_1211.RDS")

#############

fit1_data <- samps_R_2$draws()
fit1_summary <- samps_R_2$summary()

saveRDS(fit1_data, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/fit1_data.RDS")
saveRDS(fit1_summary, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/fit1_summary.RDS")  

fit1_data <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/fit1_data.RDS")
fit1_summary <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/fit1_summary.RDS")

stan_data_2 <- list(
  # Grainsize for reduce_sum
  grainsize = 1,
  
  # Dimensions
  n_spCl = length(unique(birds$sp_cl)),
  n_sp = length(unique(birds$species)),
  n_fam = length(unique(birds$Family)),
  n_tot = nrow(birds),
  n_visit_max = max(birds$nv),
  
  # Detection matrix
  det_data = det_data,
  
  # Q and nv
  Q = birds$Q,
  nv = birds$nv,
  
  # Random effect IDs
  id_spCl = as.numeric(as.factor(birds$sp_cl)),
  id_sp = as.numeric(as.factor(birds$species)),
  id_fam = as.numeric(as.factor(birds$Family)),
  
  # Covariates
  relev = (birds$elev_sp_standard -.5)/sd(birds$elev_sp_standard),
  relev2 = ((birds$elev_sp_standard -.5)/sd(birds$elev_sp_standard))^2,
  lowland = as.numeric(birds$lower == 0),
  pasture = birds$pasture,
  eastOnly = birds$east_only,
  westOnly = birds$west_only,
  snsmOnly = birds$snsm_only,
  notWandes = birds$wandes_absent,
  notEandes = birds$eandes_absent,
  elevMedian = vscale(birds$elev_median),
  elevBreadth = vscale(birds$elev_breadth),
  forestPresent = birds$forest_present,
  forestSpecialist = birds$forest_specialist,
  tfSpecialist = birds$tf_specialist,
  dryForestPresent = birds$dry_forest_present,
  floodDrySpecialist = birds$flood_dry_specialist,
  floodSpecialist = birds$floodplain_specialist,
  aridPresent = birds$arid_present,
  migratory = as.numeric(birds$migrat_birdlife == "full migrant"),           ## This needs to change for final version
  mass = vscale(birds$BodyMass.Value),
  dietInvert = as.numeric(birds$Diet.5Cat == "Invertebrate"),
  dietCarn = as.numeric(birds$Diet.5Cat == "VertFishScav"),
  dietFruitNect = as.numeric(birds$Diet.5Cat == "FruiNect"),
  dietGran = as.numeric(birds$Diet.5Cat == "PlantSeed"),
  time = time,
  obsSM = obsSM,
  obsJG = obsJG,
  obsDE = obsDE,
  
  # offsets and multipliers
  mu_b0_off = fit1_summary$median[fit1_summary$variable == "mu_b0"],
  mu_b0_mult = fit1_summary$mad[fit1_summary$variable == "mu_b0"],
  
  log_sigma_b0_sp_off = median(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_b0_sp")])),
  log_sigma_b0_sp_mult = mad(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_b0_sp")])),
  
  b0_sp_off = fit1_summary$median[grep("b0_sp\\[", fit1_summary$variable)],
  b0_sp_mult = fit1_summary$mad[grep("b0_sp\\[", fit1_summary$variable)],
  
  log_sigma_b0_fam_off = median(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_b0_fam")])),
  log_sigma_b0_fam_mult = mad(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_b0_fam")])),
  
  b0_fam_off = fit1_summary$median[grep("b0_fam\\[", fit1_summary$variable)],
  b0_fam_mult = fit1_summary$mad[grep("b0_fam\\[", fit1_summary$variable)],
  
  mu_b1_relev_off = fit1_summary$median[fit1_summary$variable == "mu_b1_relev"],
  mu_b1_relev_mult = fit1_summary$mad[fit1_summary$variable == "mu_b1_relev"],
  
  log_sigma_b1_relev_sp_off = median(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_b1_relev_sp")])),
  log_sigma_b1_relev_sp_mult = mad(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_b1_relev_sp")])),
  
  b1_relev_sp_off = fit1_summary$median[grep("b1_relev_sp\\[", fit1_summary$variable)],
  b1_relev_sp_mult = fit1_summary$mad[grep("b1_relev_sp\\[", fit1_summary$variable)],
  
  mu_b1_relev2_off = fit1_summary$median[fit1_summary$variable == "mu_b1_relev2"],
  mu_b1_relev2_mult = fit1_summary$mad[fit1_summary$variable == "mu_b1_relev2"],
  
  log_sigma_b1_relev2_sp_off = median(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_b1_relev2_sp")])),
  log_sigma_b1_relev2_sp_mult = mad(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_b1_relev2_sp")])),
  
  b1_relev2_sp_off = fit1_summary$median[grep("b1_relev2_sp\\[", fit1_summary$variable)],
  b1_relev2_sp_mult = fit1_summary$mad[grep("b1_relev2_sp\\[", fit1_summary$variable)],
  
  mu_b2_pasture_off = fit1_summary$median[fit1_summary$variable == "mu_b2_pasture"],
  mu_b2_pasture_mult = fit1_summary$mad[fit1_summary$variable == "mu_b2_pasture"],
  
  log_sigma_b2_pasture_sp_off = median(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_b2_pasture_sp")])),
  log_sigma_b2_pasture_sp_mult = mad(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_b2_pasture_sp")])),
  
  b2_pasture_sp_off = fit1_summary$median[grep("b2_pasture_sp\\[", fit1_summary$variable)],
  b2_pasture_sp_mult = fit1_summary$mad[grep("b2_pasture_sp\\[", fit1_summary$variable)],
  
  log_sigma_b2_pasture_fam_off = median(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_b2_pasture_fam")])),
  log_sigma_b2_pasture_fam_mult = mad(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_b2_pasture_fam")])),
  
  b2_pasture_fam_off = fit1_summary$median[grep("b2_pasture_fam\\[", fit1_summary$variable)],
  b2_pasture_fam_mult = fit1_summary$mad[grep("b2_pasture_fam\\[", fit1_summary$variable)],
  
  b3_eastOnly_off = fit1_summary$median[fit1_summary$variable == "b3_eastOnly"],
  b3_eastOnly_mult = fit1_summary$mad[fit1_summary$variable == "b3_eastOnly"],
  
  b3_westOnly_off = fit1_summary$median[fit1_summary$variable == "b3_westOnly"],
  b3_westOnly_mult = fit1_summary$mad[fit1_summary$variable == "b3_westOnly"],
  
  b3_snsmOnly_off = fit1_summary$median[fit1_summary$variable == "b3_snsmOnly"],
  b3_snsmOnly_mult = fit1_summary$mad[fit1_summary$variable == "b3_snsmOnly"],
  
  b3_notWandes_off = fit1_summary$median[fit1_summary$variable == "b3_notWandes"],
  b3_notWandes_mult = fit1_summary$mad[fit1_summary$variable == "b3_notWandes"],
  
  b3_notEandes_off = fit1_summary$median[fit1_summary$variable == "b3_notEandes"],
  b3_notEandes_mult = fit1_summary$mad[fit1_summary$variable == "b3_notEandes"],
  
  b3_elevMedian_off = fit1_summary$median[fit1_summary$variable == "b3_elevMedian"],
  b3_elevMedian_mult = fit1_summary$mad[fit1_summary$variable == "b3_elevMedian"],
  
  b3_elevBreadth_off = fit1_summary$median[fit1_summary$variable == "b3_elevBreadth"],
  b3_elevBreadth_mult = fit1_summary$mad[fit1_summary$variable == "b3_elevBreadth"],
  
  b3_forestPresent_off = fit1_summary$median[fit1_summary$variable == "b3_forestPresent"],
  b3_forestPresent_mult = fit1_summary$mad[fit1_summary$variable == "b3_forestPresent"],
  
  b3_forestSpecialist_off = fit1_summary$median[fit1_summary$variable == "b3_forestSpecialist"],
  b3_forestSpecialist_mult = fit1_summary$mad[fit1_summary$variable == "b3_forestSpecialist"],
  
  b3_tfSpecialist_off = fit1_summary$median[fit1_summary$variable == "b3_tfSpecialist"],
  b3_tfSpecialist_mult = fit1_summary$mad[fit1_summary$variable == "b3_tfSpecialist"],
  
  b3_dryForestPresent_off = fit1_summary$median[fit1_summary$variable == "b3_dryForestPresent"],
  b3_dryForestPresent_mult = fit1_summary$mad[fit1_summary$variable == "b3_dryForestPresent"],
  
  b3_floodDrySpecialist_off = fit1_summary$median[fit1_summary$variable == "b3_floodDrySpecialist"],
  b3_floodDrySpecialist_mult = fit1_summary$mad[fit1_summary$variable == "b3_floodDrySpecialist"],
  
  b3_floodSpecialist_off = fit1_summary$median[fit1_summary$variable == "b3_floodSpecialist"],
  b3_floodSpecialist_mult = fit1_summary$mad[fit1_summary$variable == "b3_floodSpecialist"],
  
  b3_aridPresent_off = fit1_summary$median[fit1_summary$variable == "b3_aridPresent"],
  b3_aridPresent_mult = fit1_summary$mad[fit1_summary$variable == "b3_aridPresent"],
  
  b3_migratory_off = fit1_summary$median[fit1_summary$variable == "b3_migratory"],
  b3_migratory_mult = fit1_summary$mad[fit1_summary$variable == "b3_migratory"],
  
  b3_mass_off = fit1_summary$median[fit1_summary$variable == "b3_mass"],
  b3_mass_mult = fit1_summary$mad[fit1_summary$variable == "b3_mass"],
  
  b3_dietInvert_off = fit1_summary$median[fit1_summary$variable == "b3_dietInvert"],
  b3_dietInvert_mult = fit1_summary$mad[fit1_summary$variable == "b3_dietInvert"],
  
  b3_dietCarn_off = fit1_summary$median[fit1_summary$variable == "b3_dietCarn"],
  b3_dietCarn_mult = fit1_summary$mad[fit1_summary$variable == "b3_dietCarn"],
  
  b3_dietFruitNect_off = fit1_summary$median[fit1_summary$variable == "b3_dietFruitNect"],
  b3_dietFruitNect_mult = fit1_summary$mad[fit1_summary$variable == "b3_dietFruitNect"],
  
  b3_dietGran_off = fit1_summary$median[fit1_summary$variable == "b3_dietGran"],
  b3_dietGran_mult = fit1_summary$mad[fit1_summary$variable == "b3_dietGran"],
  
  b3_x_elevMedian_forestPresent_off = fit1_summary$median[fit1_summary$variable == "b3_x_elevMedian_forestPresent"],
  b3_x_elevMedian_forestPresent_mult = fit1_summary$mad[fit1_summary$variable == "b3_x_elevMedian_forestPresent"],
  
  b3_x_elevMedian_forestSpecialist_off = fit1_summary$median[fit1_summary$variable == "b3_x_elevMedian_forestSpecialist"],
  b3_x_elevMedian_forestSpecialist_mult = fit1_summary$mad[fit1_summary$variable == "b3_x_elevMedian_forestSpecialist"],
  
  b4_eastOnly_off = fit1_summary$median[fit1_summary$variable == "b4_eastOnly"],
  b4_eastOnly_mult = fit1_summary$mad[fit1_summary$variable == "b4_eastOnly"],
  
  b4_westOnly_off = fit1_summary$median[fit1_summary$variable == "b4_westOnly"],
  b4_westOnly_mult = fit1_summary$mad[fit1_summary$variable == "b4_westOnly"],
  
  b4_snsmOnly_off = fit1_summary$median[fit1_summary$variable == "b4_snsmOnly"],
  b4_snsmOnly_mult = fit1_summary$mad[fit1_summary$variable == "b4_snsmOnly"],
  
  b4_notWandes_off = fit1_summary$median[fit1_summary$variable == "b4_notWandes"],
  b4_notWandes_mult = fit1_summary$mad[fit1_summary$variable == "b4_notWandes"],
  
  b4_notEandes_off = fit1_summary$median[fit1_summary$variable == "b4_notEandes"],
  b4_notEandes_mult = fit1_summary$mad[fit1_summary$variable == "b4_notEandes"],
  
  b4_elevMedian_off = fit1_summary$median[fit1_summary$variable == "b4_elevMedian"],
  b4_elevMedian_mult = fit1_summary$mad[fit1_summary$variable == "b4_elevMedian"],
  
  b4_elevBreadth_off = fit1_summary$median[fit1_summary$variable == "b4_elevBreadth"],
  b4_elevBreadth_mult = fit1_summary$mad[fit1_summary$variable == "b4_elevBreadth"],
  
  b4_forestPresent_off = fit1_summary$median[fit1_summary$variable == "b4_forestPresent"],
  b4_forestPresent_mult = fit1_summary$mad[fit1_summary$variable == "b4_forestPresent"],
  
  b4_forestSpecialist_off = fit1_summary$median[fit1_summary$variable == "b4_forestSpecialist"],
  b4_forestSpecialist_mult = fit1_summary$mad[fit1_summary$variable == "b4_forestSpecialist"],
  
  b4_tfSpecialist_off = fit1_summary$median[fit1_summary$variable == "b4_tfSpecialist"],
  b4_tfSpecialist_mult = fit1_summary$mad[fit1_summary$variable == "b4_tfSpecialist"],
  
  b4_dryForestPresent_off = fit1_summary$median[fit1_summary$variable == "b4_dryForestPresent"],
  b4_dryForestPresent_mult = fit1_summary$mad[fit1_summary$variable == "b4_dryForestPresent"],
  
  b4_floodDrySpecialist_off = fit1_summary$median[fit1_summary$variable == "b4_floodDrySpecialist"],
  b4_floodDrySpecialist_mult = fit1_summary$mad[fit1_summary$variable == "b4_floodDrySpecialist"],
  
  b4_floodSpecialist_off = fit1_summary$median[fit1_summary$variable == "b4_floodSpecialist"],
  b4_floodSpecialist_mult = fit1_summary$mad[fit1_summary$variable == "b4_floodSpecialist"],
  
  b4_aridPresent_off = fit1_summary$median[fit1_summary$variable == "b4_aridPresent"],
  b4_aridPresent_mult = fit1_summary$mad[fit1_summary$variable == "b4_aridPresent"],
  
  b4_migratory_off = fit1_summary$median[fit1_summary$variable == "b4_migratory"],
  b4_migratory_mult = fit1_summary$mad[fit1_summary$variable == "b4_migratory"],
  
  b4_mass_off = fit1_summary$median[fit1_summary$variable == "b4_mass"],
  b4_mass_mult = fit1_summary$mad[fit1_summary$variable == "b4_mass"],
  
  b4_dietInvert_off = fit1_summary$median[fit1_summary$variable == "b4_dietInvert"],
  b4_dietInvert_mult = fit1_summary$mad[fit1_summary$variable == "b4_dietInvert"],
  
  b4_dietCarn_off = fit1_summary$median[fit1_summary$variable == "b4_dietCarn"],
  b4_dietCarn_mult = fit1_summary$mad[fit1_summary$variable == "b4_dietCarn"],
  
  b4_dietFruitNect_off = fit1_summary$median[fit1_summary$variable == "b4_dietFruitNect"],
  b4_dietFruitNect_mult = fit1_summary$mad[fit1_summary$variable == "b4_dietFruitNect"],
  
  b4_dietGran_off = fit1_summary$median[fit1_summary$variable == "b4_dietGran"],
  b4_dietGran_mult = fit1_summary$mad[fit1_summary$variable == "b4_dietGran"],
  
  b4_x_elevMedian_forestPresent_off = fit1_summary$median[fit1_summary$variable == "b4_x_elevMedian_forestPresent"],
  b4_x_elevMedian_forestPresent_mult = fit1_summary$mad[fit1_summary$variable == "b4_x_elevMedian_forestPresent"],
  
  b4_x_elevMedian_forestSpecialist_off = fit1_summary$median[fit1_summary$variable == "b4_x_elevMedian_forestSpecialist"],
  b4_x_elevMedian_forestSpecialist_mult = fit1_summary$mad[fit1_summary$variable == "b4_x_elevMedian_forestSpecialist"],
  
  mu_d0_off = fit1_summary$median[fit1_summary$variable == "mu_d0"],
  mu_d0_mult = fit1_summary$mad[fit1_summary$variable == "mu_d0"],
  
  log_sigma_d0_sp_off = median(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_d0_sp")])),
  log_sigma_d0_sp_mult = mad(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_d0_sp")])),
  
  d0_sp_off = fit1_summary$median[grep("d0_sp\\[", fit1_summary$variable)],
  d0_sp_mult = fit1_summary$mad[grep("d0_sp\\[", fit1_summary$variable)],
  
  log_sigma_d0_fam_off = median(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_d0_fam")])),
  log_sigma_d0_fam_mult = mad(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_d0_fam")])),
  
  d0_fam_off = fit1_summary$median[grep("d0_fam\\[", fit1_summary$variable)],
  d0_fam_mult = fit1_summary$mad[grep("d0_fam\\[", fit1_summary$variable)],
  
  mu_d1_pasture_off = fit1_summary$median[fit1_summary$variable == "mu_d1_pasture"],
  mu_d1_pasture_mult = fit1_summary$mad[fit1_summary$variable == "mu_d1_pasture"],
  
  log_sigma_d1_pasture_sp_off = median(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_d1_pasture_sp")])),
  log_sigma_d1_pasture_sp_mult = mad(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_d1_pasture_sp")])),
  
  d1_pasture_sp_off = fit1_summary$median[grep("d1_pasture_sp\\[", fit1_summary$variable)],
  d1_pasture_sp_mult = fit1_summary$mad[grep("d1_pasture_sp\\[", fit1_summary$variable)],
  
  log_sigma_d1_pasture_fam_off = median(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_d1_pasture_fam")]*2)), # multiplying by 2 because accidentally had the raw values drawn from normal(0,2) in the previous model
  log_sigma_d1_pasture_fam_mult = mad(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_d1_pasture_fam")]*2)),
  
  d1_pasture_fam_off = fit1_summary$median[grep("d1_pasture_fam\\[", fit1_summary$variable)],
  d1_pasture_fam_mult = fit1_summary$mad[grep("d1_pasture_fam\\[", fit1_summary$variable)],
  
  mu_d2_time_off = fit1_summary$median[fit1_summary$variable == "mu_d2_time"],
  mu_d2_time_mult = fit1_summary$mad[fit1_summary$variable == "mu_d2_time"],
  
  log_sigma_d2_time_sp_off = median(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_d2_time_sp")])),
  log_sigma_d2_time_sp_mult = mad(log(fit1_data[,,which(dimnames(fit1_data)$variable == "sigma_d2_time_sp")])),
  
  d2_time_sp_off = fit1_summary$median[grep("d2_time_sp\\[", fit1_summary$variable)],
  d2_time_sp_mult = fit1_summary$mad[grep("d2_time_sp\\[", fit1_summary$variable)],
  
  d2_obsSM_off = fit1_summary$median[fit1_summary$variable == "d2_obsSM"],
  d2_obsSM_mult = fit1_summary$mad[fit1_summary$variable == "d2_obsSM"],
  
  d2_obsDE_off = fit1_summary$median[fit1_summary$variable == "d2_obsDE"],
  d2_obsDE_mult = fit1_summary$mad[fit1_summary$variable == "d2_obsDE"],
  
  d2_obsJG_off = fit1_summary$median[fit1_summary$variable == "d2_obsJG"],
  d2_obsJG_mult = fit1_summary$mad[fit1_summary$variable == "d2_obsJG"],
  
  d3_mass_off = fit1_summary$median[fit1_summary$variable == "d3_mass"],
  d3_mass_mult = fit1_summary$mad[fit1_summary$variable == "d3_mass"],
  
  d3_elevMedian_off = fit1_summary$median[fit1_summary$variable == "d3_elevMedian"],
  d3_elevMedian_mult = fit1_summary$mad[fit1_summary$variable == "d3_elevMedian"],
  
  d3_migratory_off = fit1_summary$median[fit1_summary$variable == "d3_migratory"],
  d3_migratory_mult = fit1_summary$mad[fit1_summary$variable == "d3_migratory"]
)



mod_R_3 <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/stan_files/full_colombia_model/occupancyMod_ragged_parallel_v3.stan",
                       cpp_options = list(stan_threads = TRUE))
fullmod_test_samples <- mod_R_3$sample(data = stan_data_2, 
                                    chains = 1,
                                    threads_per_chain = 4,
                                    refresh = 1,
                                    iter_sampling = 500,
                                    iter_warmup = 500,
                                    save_warmup = 1,
                                    step_size = .02,
                                    output_dir = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/CSVs")



