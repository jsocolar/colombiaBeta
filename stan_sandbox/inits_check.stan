parameters{
  real mu;
  real<lower=0> sigma;
  vector[500] v;
}
model{
  // priors
  sigma ~ normal(0, 10);
  mu ~ normal(0, 10);
  v ~ normal(mu, sigma);
}