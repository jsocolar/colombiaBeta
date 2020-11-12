library("cmdstanr")

jacobianTest_mod <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/stan_sandbox/test_jacobianadj_log.stan",
                         cpp_options = list(stan_threads = TRUE))
jacobianTest_samples <- jacobianTest_mod$sample(
                                       chains = 1,
                                       threads_per_chain = 1,
                                       refresh = 100000,
                                       iter_sampling = 1000000,
                                       iter_warmup = 10000,
                                       save_warmup = 1)

#jacobianTest_samples$summary()
quantile(as.numeric(exp(jacobianTest_samples$draws()[,,dimnames(jacobianTest_samples$draws())$variable == "log_sigma1"])),  probs = c(.05, .5, .95)) -
  quantile(as.numeric(jacobianTest_samples$draws()[,,dimnames(jacobianTest_samples$draws())$variable == "sigma"]),  probs = c(.05, .5, .95))
