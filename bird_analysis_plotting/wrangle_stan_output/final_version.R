remotes::install_github("stan-dev/cmdstanr")
install.packages("/Users/JacobSocolar/Dropbox/Work/Code/Stan/posterior", 
                 repos = NULL, 
                 type = "source")
library(posterior)
library(cmdstanr)

stan_csvs <- list.files("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/CSVs/Final", full.names = TRUE)

fc_long <- as_cmdstan_fit(stan_csvs, format = "draws_list")
fmsd <- fc_long$sampler_diagnostics()
unique(fmsd[,,"divergent__"])
ebfmi <- apply(fmsd[,,"energy__"], 2, function(x) {
  (sum(diff(x)^2)/length(x))/var(x)
})
ebfmi


final_draws <- fc_long$draws()

saveRDS(final_draws, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/v9_final/draws.RDS")
final_draws <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/v9_final/draws.RDS")
fdsum <- summarise_draws(final_draws, .cores = 3)
saveRDS(fdsum, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/v9_final/draws_summary.RDS")
fdsum <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/v9_final/draws_summary.RDS")

# The below uses Jacob's personal fork of posterior:
# rollup <- posterior::rollup_summary(fdsum)
# 
# rollup$unrolled_vars[rollup$unrolled_vars$variable == "sigma_b0_spSr", ]
