remotes::install_github("stan-dev/cmdstanr")
install.packages("/Users/JacobSocolar/Dropbox/Work/Code/Stan/posterior", 
                 repos = NULL, 
                 type = "source")
library(posterior)
library(cmdstanr)

long_csvs <- c("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/CSVs/Final/occupancy_v9-202105070926-1-3fd25b.csv",
               "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/CSVs/Final/occupancy_v9-202105070926-2-3fd25b.csv",
               "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/CSVs/Final/occupancy_v9-202105061701-1-2a4265.csv")

fc_long <- as_cmdstan_fit(long_csvs, format = "draws_list")
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
rollup <- posterior::rollup_summary(fdsum)

rollup$unrolled_vars[rollup$unrolled_vars$variable == "sigma_b0_spSr", ]
