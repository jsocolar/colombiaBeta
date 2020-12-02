library("cmdstanr"); library("dplyr"); library("posterior")

bird_stan_data1_package <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data1_package.RDS")

# Run ragged model ----
mod_R_1511 <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/stan_files/full_colombia_model/occupancyMod_ragged_parallel_v1.stan",
                       cpp_options = list(stan_threads = TRUE))
samps_R_1511 <- mod_R_1511$sample(data = bird_stan_data1_package$data, 
                          chains = 1,
                          threads_per_chain = 4,
                          refresh = 1,
                          iter_sampling = 500,
                          iter_warmup = 500,
                          save_warmup = 1,
                          output_dir = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/CSVs")

saveRDS(samps_R_1511, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/samps_R_1511.RDS")

#############

fit1511_data <- samps_R_1511$draws()
fit1511_summary <- samps_R_1511$summary()

saveRDS(fit1511_data, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/fit1511_data.RDS")
saveRDS(fit1511_summary, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/fit1511_summary.RDS")  

fit1511_data <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/fit1511_data.RDS")
fit1511_summary <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/fit1511_summary.RDS")

offsets_and_multipliers <- list(
  # offsets and multipliers
  mu_b0_off = fit1511_summary$median[fit1511_summary$variable == "mu_b0"],
  mu_b0_mult = fit1511_summary$mad[fit1511_summary$variable == "mu_b0"],
  
  log_sigma_b0_sp_off = median(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_b0_sp")])),
  log_sigma_b0_sp_mult = mad(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_b0_sp")])),
  
  b0_sp_off = fit1511_summary$median[grep("b0_sp\\[", fit1511_summary$variable)],
  b0_sp_mult = fit1511_summary$mad[grep("b0_sp\\[", fit1511_summary$variable)],
  
  log_sigma_b0_fam_off = median(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_b0_fam")])),
  log_sigma_b0_fam_mult = mad(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_b0_fam")])),
  
  b0_fam_off = fit1511_summary$median[grep("b0_fam\\[", fit1511_summary$variable)],
  b0_fam_mult = fit1511_summary$mad[grep("b0_fam\\[", fit1511_summary$variable)],
  
  mu_b1_relev_off = fit1511_summary$median[fit1511_summary$variable == "mu_b1_relev"],
  mu_b1_relev_mult = fit1511_summary$mad[fit1511_summary$variable == "mu_b1_relev"],
  
  log_sigma_b1_relev_sp_off = median(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_b1_relev_sp")])),
  log_sigma_b1_relev_sp_mult = mad(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_b1_relev_sp")])),
  
  b1_relev_sp_off = fit1511_summary$median[grep("b1_relev_sp\\[", fit1511_summary$variable)],
  b1_relev_sp_mult = fit1511_summary$mad[grep("b1_relev_sp\\[", fit1511_summary$variable)],
  
  mu_b1_relev2_off = fit1511_summary$median[fit1511_summary$variable == "mu_b1_relev2"],
  mu_b1_relev2_mult = fit1511_summary$mad[fit1511_summary$variable == "mu_b1_relev2"],
  
  log_sigma_b1_relev2_sp_off = median(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_b1_relev2_sp")])),
  log_sigma_b1_relev2_sp_mult = mad(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_b1_relev2_sp")])),
  
  b1_relev2_sp_off = fit1511_summary$median[grep("b1_relev2_sp\\[", fit1511_summary$variable)],
  b1_relev2_sp_mult = fit1511_summary$mad[grep("b1_relev2_sp\\[", fit1511_summary$variable)],
  
  b1_lowland_off = fit1511_summary$median[fit1511_summary$variable == "b1_lowland"],
  b1_lowland_mult = fit1511_summary$mad[fit1511_summary$variable == "b1_lowland"],
  
  b1_x_lowland_relev_off = fit1511_summary$median[fit1511_summary$variable == "b1_x_lowland_relev"],
  b1_x_lowland_relev_mult = fit1511_summary$mad[fit1511_summary$variable == "b1_x_lowland_relev"],
  
  b1_x_lowland_relev2_off = fit1511_summary$median[fit1511_summary$variable == "b1_x_lowland_relev2"],
  b1_x_lowland_relev2_mult = fit1511_summary$mad[fit1511_summary$variable == "b1_x_lowland_relev2"],
  
  mu_b2_pasture_off = fit1511_summary$median[fit1511_summary$variable == "mu_b2_pasture"],
  mu_b2_pasture_mult = fit1511_summary$mad[fit1511_summary$variable == "mu_b2_pasture"],
  
  log_sigma_b2_pasture_sp_off = median(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_b2_pasture_sp")])),
  log_sigma_b2_pasture_sp_mult = mad(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_b2_pasture_sp")])),
  
  b2_pasture_sp_off = fit1511_summary$median[grep("b2_pasture_sp\\[", fit1511_summary$variable)],
  b2_pasture_sp_mult = fit1511_summary$mad[grep("b2_pasture_sp\\[", fit1511_summary$variable)],
  
  log_sigma_b2_pasture_fam_off = median(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_b2_pasture_fam")])),
  log_sigma_b2_pasture_fam_mult = mad(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_b2_pasture_fam")])),
  
  b2_pasture_fam_off = fit1511_summary$median[grep("b2_pasture_fam\\[", fit1511_summary$variable)],
  b2_pasture_fam_mult = fit1511_summary$mad[grep("b2_pasture_fam\\[", fit1511_summary$variable)],
  
  b3_eastOnly_off = fit1511_summary$median[fit1511_summary$variable == "b3_eastOnly"],
  b3_eastOnly_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_eastOnly"],
  
  b3_westOnly_off = fit1511_summary$median[fit1511_summary$variable == "b3_westOnly"],
  b3_westOnly_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_westOnly"],
  
  b3_snsmOnly_off = fit1511_summary$median[fit1511_summary$variable == "b3_snsmOnly"],
  b3_snsmOnly_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_snsmOnly"],
  
  b3_notWandes_off = fit1511_summary$median[fit1511_summary$variable == "b3_notWandes"],
  b3_notWandes_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_notWandes"],
  
  b3_notEandes_off = fit1511_summary$median[fit1511_summary$variable == "b3_notEandes"],
  b3_notEandes_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_notEandes"],
  
  b3_elevMedian_off = fit1511_summary$median[fit1511_summary$variable == "b3_elevMedian"],
  b3_elevMedian_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_elevMedian"],
  
  b3_elevBreadth_off = fit1511_summary$median[fit1511_summary$variable == "b3_elevBreadth"],
  b3_elevBreadth_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_elevBreadth"],
  
  b3_forestPresent_off = fit1511_summary$median[fit1511_summary$variable == "b3_forestPresent"],
  b3_forestPresent_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_forestPresent"],
  
  b3_forestSpecialist_off = fit1511_summary$median[fit1511_summary$variable == "b3_forestSpecialist"],
  b3_forestSpecialist_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_forestSpecialist"],
  
  b3_tfSpecialist_off = fit1511_summary$median[fit1511_summary$variable == "b3_tfSpecialist"],
  b3_tfSpecialist_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_tfSpecialist"],
  
  b3_dryForestPresent_off = fit1511_summary$median[fit1511_summary$variable == "b3_dryForestPresent"],
  b3_dryForestPresent_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_dryForestPresent"],
  
  b3_floodDrySpecialist_off = fit1511_summary$median[fit1511_summary$variable == "b3_floodDrySpecialist"],
  b3_floodDrySpecialist_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_floodDrySpecialist"],
  
  b3_floodSpecialist_off = fit1511_summary$median[fit1511_summary$variable == "b3_floodSpecialist"],
  b3_floodSpecialist_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_floodSpecialist"],
  
  b3_aridPresent_off = fit1511_summary$median[fit1511_summary$variable == "b3_aridPresent"],
  b3_aridPresent_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_aridPresent"],
  
  b3_migratory_off = fit1511_summary$median[fit1511_summary$variable == "b3_migratory"],
  b3_migratory_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_migratory"],
  
  b3_mass_off = fit1511_summary$median[fit1511_summary$variable == "b3_mass"],
  b3_mass_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_mass"],
  
  b3_dietInvert_off = fit1511_summary$median[fit1511_summary$variable == "b3_dietInvert"],
  b3_dietInvert_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_dietInvert"],
  
  b3_dietCarn_off = fit1511_summary$median[fit1511_summary$variable == "b3_dietCarn"],
  b3_dietCarn_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_dietCarn"],
  
  b3_dietFruitNect_off = fit1511_summary$median[fit1511_summary$variable == "b3_dietFruitNect"],
  b3_dietFruitNect_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_dietFruitNect"],
  
  b3_dietGran_off = fit1511_summary$median[fit1511_summary$variable == "b3_dietGran"],
  b3_dietGran_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_dietGran"],
  
  b3_x_elevMedian_forestPresent_off = fit1511_summary$median[fit1511_summary$variable == "b3_x_elevMedian_forestPresent"],
  b3_x_elevMedian_forestPresent_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_x_elevMedian_forestPresent"],
  
  b3_x_elevMedian_forestSpecialist_off = fit1511_summary$median[fit1511_summary$variable == "b3_x_elevMedian_forestSpecialist"],
  b3_x_elevMedian_forestSpecialist_mult = fit1511_summary$mad[fit1511_summary$variable == "b3_x_elevMedian_forestSpecialist"],
  
  b4_eastOnly_off = fit1511_summary$median[fit1511_summary$variable == "b4_eastOnly"],
  b4_eastOnly_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_eastOnly"],
  
  b4_westOnly_off = fit1511_summary$median[fit1511_summary$variable == "b4_westOnly"],
  b4_westOnly_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_westOnly"],
  
  b4_snsmOnly_off = fit1511_summary$median[fit1511_summary$variable == "b4_snsmOnly"],
  b4_snsmOnly_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_snsmOnly"],
  
  b4_notWandes_off = fit1511_summary$median[fit1511_summary$variable == "b4_notWandes"],
  b4_notWandes_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_notWandes"],
  
  b4_notEandes_off = fit1511_summary$median[fit1511_summary$variable == "b4_notEandes"],
  b4_notEandes_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_notEandes"],
  
  b4_elevMedian_off = fit1511_summary$median[fit1511_summary$variable == "b4_elevMedian"],
  b4_elevMedian_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_elevMedian"],
  
  b4_elevBreadth_off = fit1511_summary$median[fit1511_summary$variable == "b4_elevBreadth"],
  b4_elevBreadth_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_elevBreadth"],
  
  b4_forestPresent_off = fit1511_summary$median[fit1511_summary$variable == "b4_forestPresent"],
  b4_forestPresent_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_forestPresent"],
  
  b4_forestSpecialist_off = fit1511_summary$median[fit1511_summary$variable == "b4_forestSpecialist"],
  b4_forestSpecialist_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_forestSpecialist"],
  
  b4_tfSpecialist_off = fit1511_summary$median[fit1511_summary$variable == "b4_tfSpecialist"],
  b4_tfSpecialist_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_tfSpecialist"],
  
  b4_dryForestPresent_off = fit1511_summary$median[fit1511_summary$variable == "b4_dryForestPresent"],
  b4_dryForestPresent_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_dryForestPresent"],
  
  b4_floodDrySpecialist_off = fit1511_summary$median[fit1511_summary$variable == "b4_floodDrySpecialist"],
  b4_floodDrySpecialist_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_floodDrySpecialist"],
  
  b4_floodSpecialist_off = fit1511_summary$median[fit1511_summary$variable == "b4_floodSpecialist"],
  b4_floodSpecialist_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_floodSpecialist"],
  
  b4_aridPresent_off = fit1511_summary$median[fit1511_summary$variable == "b4_aridPresent"],
  b4_aridPresent_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_aridPresent"],
  
  b4_migratory_off = fit1511_summary$median[fit1511_summary$variable == "b4_migratory"],
  b4_migratory_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_migratory"],
  
  b4_mass_off = fit1511_summary$median[fit1511_summary$variable == "b4_mass"],
  b4_mass_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_mass"],
  
  b4_dietInvert_off = fit1511_summary$median[fit1511_summary$variable == "b4_dietInvert"],
  b4_dietInvert_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_dietInvert"],
  
  b4_dietCarn_off = fit1511_summary$median[fit1511_summary$variable == "b4_dietCarn"],
  b4_dietCarn_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_dietCarn"],
  
  b4_dietFruitNect_off = fit1511_summary$median[fit1511_summary$variable == "b4_dietFruitNect"],
  b4_dietFruitNect_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_dietFruitNect"],
  
  b4_dietGran_off = fit1511_summary$median[fit1511_summary$variable == "b4_dietGran"],
  b4_dietGran_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_dietGran"],
  
  b4_x_elevMedian_forestPresent_off = fit1511_summary$median[fit1511_summary$variable == "b4_x_elevMedian_forestPresent"],
  b4_x_elevMedian_forestPresent_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_x_elevMedian_forestPresent"],
  
  b4_x_elevMedian_forestSpecialist_off = fit1511_summary$median[fit1511_summary$variable == "b4_x_elevMedian_forestSpecialist"],
  b4_x_elevMedian_forestSpecialist_mult = fit1511_summary$mad[fit1511_summary$variable == "b4_x_elevMedian_forestSpecialist"],
  
  mu_d0_off = fit1511_summary$median[fit1511_summary$variable == "mu_d0"],
  mu_d0_mult = fit1511_summary$mad[fit1511_summary$variable == "mu_d0"],
  
  log_sigma_d0_sp_off = median(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_d0_sp")])),
  log_sigma_d0_sp_mult = mad(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_d0_sp")])),
  
  d0_sp_off = fit1511_summary$median[grep("d0_sp\\[", fit1511_summary$variable)],
  d0_sp_mult = fit1511_summary$mad[grep("d0_sp\\[", fit1511_summary$variable)],
  
  log_sigma_d0_fam_off = median(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_d0_fam")])),
  log_sigma_d0_fam_mult = mad(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_d0_fam")])),
  
  d0_fam_off = fit1511_summary$median[grep("d0_fam\\[", fit1511_summary$variable)],
  d0_fam_mult = fit1511_summary$mad[grep("d0_fam\\[", fit1511_summary$variable)],
  
  mu_d1_pasture_off = fit1511_summary$median[fit1511_summary$variable == "mu_d1_pasture"],
  mu_d1_pasture_mult = fit1511_summary$mad[fit1511_summary$variable == "mu_d1_pasture"],
  
  log_sigma_d1_pasture_sp_off = median(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_d1_pasture_sp")])),
  log_sigma_d1_pasture_sp_mult = mad(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_d1_pasture_sp")])),
  
  d1_pasture_sp_off = fit1511_summary$median[grep("d1_pasture_sp\\[", fit1511_summary$variable)],
  d1_pasture_sp_mult = fit1511_summary$mad[grep("d1_pasture_sp\\[", fit1511_summary$variable)],
  
  log_sigma_d1_pasture_fam_off = median(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_d1_pasture_fam")])),
  log_sigma_d1_pasture_fam_mult = mad(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_d1_pasture_fam")])),
  
  d1_pasture_fam_off = fit1511_summary$median[grep("d1_pasture_fam\\[", fit1511_summary$variable)],
  d1_pasture_fam_mult = fit1511_summary$mad[grep("d1_pasture_fam\\[", fit1511_summary$variable)],
  
  mu_d2_time_off = fit1511_summary$median[fit1511_summary$variable == "mu_d2_time"],
  mu_d2_time_mult = fit1511_summary$mad[fit1511_summary$variable == "mu_d2_time"],
  
  log_sigma_d2_time_sp_off = median(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_d2_time_sp")])),
  log_sigma_d2_time_sp_mult = mad(log(fit1511_data[,,which(dimnames(fit1511_data)$variable == "sigma_d2_time_sp")])),
  
  d2_time_sp_off = fit1511_summary$median[grep("d2_time_sp\\[", fit1511_summary$variable)],
  d2_time_sp_mult = fit1511_summary$mad[grep("d2_time_sp\\[", fit1511_summary$variable)],
  
  d2_obsSM_off = fit1511_summary$median[fit1511_summary$variable == "d2_obsSM"],
  d2_obsSM_mult = fit1511_summary$mad[fit1511_summary$variable == "d2_obsSM"],
  
  d2_obsDE_off = fit1511_summary$median[fit1511_summary$variable == "d2_obsDE"],
  d2_obsDE_mult = fit1511_summary$mad[fit1511_summary$variable == "d2_obsDE"],
  
  d2_obsJG_off = fit1511_summary$median[fit1511_summary$variable == "d2_obsJG"],
  d2_obsJG_mult = fit1511_summary$mad[fit1511_summary$variable == "d2_obsJG"],
  
  d3_mass_off = fit1511_summary$median[fit1511_summary$variable == "d3_mass"],
  d3_mass_mult = fit1511_summary$mad[fit1511_summary$variable == "d3_mass"],
  
  d3_elevMedian_off = fit1511_summary$median[fit1511_summary$variable == "d3_elevMedian"],
  d3_elevMedian_mult = fit1511_summary$mad[fit1511_summary$variable == "d3_elevMedian"],
  
  d3_migratory_off = fit1511_summary$median[fit1511_summary$variable == "d3_migratory"],
  d3_migratory_mult = fit1511_summary$mad[fit1511_summary$variable == "d3_migratory"],
  
  d3_dietCarn_off = fit1511_summary$median[fit1511_summary$variable == "d3_dietCarn"],
  d3_dietCarn_mult = fit1511_summary$mad[fit1511_summary$variable == "d3_dietCarn"],
  
  d3_x_time_elevMedian_off = fit1511_summary$median[fit1511_summary$variable == "d3_x_time_elevMedian"],
  d3_x_time_elevMedian_mult = fit1511_summary$mad[fit1511_summary$variable == "d3_x_time_elevMedian"]
)

bird_stan_data2_package <- list(data = c(bird_stan_data1_package$data, offsets_and_multipliers), means_and_sds = bird_stan_data1_package$means_and_sds)
saveRDS(bird_stan_data2_package, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data2_package.RDS")

mod_R_3 <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/stan_files/full_colombia_model/occupancyMod_ragged_parallel_v3.stan",
                       cpp_options = list(stan_threads = TRUE))
fullmod_samples <- mod_R_3$sample(data = bird_stan_data2_package$data, 
                                    chains = 1,
                                    threads_per_chain = 4,
                                    refresh = 1,
                                    iter_sampling = 500,
                                    iter_warmup = 500,
                                    save_warmup = 1,
                                    step_size = .0015,
                                    max_treedepth = 9,
                                    output_dir = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/CSVs")



saveRDS(fullmod_samples, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/fullmod_samples_2511.RDS")

#############
samples_2511 <- cmdstanr::read_cmdstan_csv('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/CSVs/occupancyMod_ragged_parallel_v3_threads-202011191247-1-deb0d8.csv')
samples2_2511 <- cmdstanr::read_cmdstan_csv('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/CSVs/occupancyMod_ragged_parallel_v3_threads-202011191247-1-deb0d8.csv', 
                                                        variables = samples_2511$metadata$stan_variables[c(1:3, 5:87)])
muSamples_2511 <- cmdstanr::read_cmdstan_csv('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/CSVs/occupancyMod_ragged_parallel_v3_threads-202011191247-1-deb0d8.csv', 
                                             variables = "mu_b0")

bird_stan_data3_package <- readRDS('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data2_package.RDS')
bird_stan_data3_package$data$mu_b0_off <- posterior::summarise_draws(muSamples_2511$post_warmup_draws)$median
bird_stan_data3_package$data$mu_b0_mult <- posterior::summarise_draws(muSamples_2511$post_warmup_draws)$mad
saveRDS(bird_stan_data3_package, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data3_package.RDS")

##############

mod_R_4 <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/stan_files/full_colombia_model/occupancyMod_ragged_parallel_v4.stan",
                         cpp_options = list(stan_threads = TRUE))
fullmod_samples <- mod_R_4$sample(data = bird_stan_data3_package$data, 
                                  chains = 1,
                                  threads_per_chain = 4,
                                  refresh = 1,
                                  iter_sampling = 500,
                                  iter_warmup = 500,
                                  save_warmup = 1,
                                  step_size = .0015,
                                  max_treedepth = 9,
                                  output_dir = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/CSVs")


saveRDS(fullmod_samples, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/fullmod_samples_0212.RDS")


