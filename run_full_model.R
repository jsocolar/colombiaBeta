# run pan-Colombia model on bessemer

library("cmdstanr")

bird_data <- readRDS("bird_stan_data2_package.RDS")

mod_R_3 <- cmdstan_model("stan_files/full_colombia_model/occupancyMod_ragged_parallel_v3.stan",
                         cpp_options = list(stan_threads = TRUE))

fullmod_samples <- mod_R_3$sample(data = bird_data$data,
                                  chains = 1,
                                  threads_per_chain = 40,
                                  refresh = 50,
                                  iter_sampling = 2500,
                                  iter_warmup = 2500,
                                  save_warmup = 1,
                                  step_size = .0015,
                                  max_treedepth = 9,
                                  output_dir = "stan_output")

# also save as RDS (..should be faster to read)
saveRDS(fullmod_samples, "stan_output/full_colombia_model_run1.rds")