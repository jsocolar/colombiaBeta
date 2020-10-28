parameters{
  real<lower=0> sigma;
  real log_sigma1;
  // real log_sigma2;
}
model{
  sigma ~ normal(0, 1);
  
  real sigma_adj = exp(log_sigma1);
  target += log_sigma1;
  sigma_adj ~ normal(0, 1);
        
  // real sigma_unadj = exp(log_sigma2);
  // sigma_unadj ~ normal(0, 1);
}