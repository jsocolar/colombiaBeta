# Create pathological case, following a bunch of problems with unexpected behaviour
# fitting model to simulated data. 
#
# Simulation code inherits from Jacob's occupancysimulation2.R (accessed 27Apr2020)
# Only edits to simulation are minor reordering for legibility, and simplifying 
# the simulation so that:
#
# logit.occ[i, k] <- b0[k] + b1[clusterID[i],k] + b4[k]*pt_cov1[i] 
# logit.det[i, j, k] <- d0[k] + d5[k]*vis_cov1[i,j]
# (note d2 is no longer REff, but it's not in the simulation!)
#
# Only edits to stan code are so that:
# logit_psi[i,k] = b0[k] + b1[k, clusterID[i]] + b4[k]*pt_cov1[i]; 
# logit_theta[i,j,k] = d0[k] + d2*sp_cov1[k] + d1[k,i] + d5[k]*vis_cov1[i,j];
#
# Conclusion is that having d1[k,i] in the model completely screws up the inference. 
# In hindsight it seems obvious that a species:point varying intercept in detection 
# would be problematic, as it becomes *very* difficult to separately identify 
# detection from variation in point-level occupancy (I'm not convinced a much larger
# quantity of data would fix this either). It's bothersome that it will happily
# generate nonsense results rather than blowing up though. Note, after removing
# this term, models seem to largely behave sensibly. Note that behaviour can 
# still be concerning when some of the data dimensions are small. E.g. running
# a 30 points * 50 species case isn't so pathological but will struggle to converge
# and I'm not completely satisfied that mu_b0 isn't biased upwards in this 
# situation (though 95% intervals seem to consistently include 0).
#
# i indexes points, j indexes visits, k indexes species.

# Housekeeping ----
library(rstan)
rstan_options(auto_write = T)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
options(mc.cores = parallel::detectCores()).

# Seed ----
set.seed(101)

# Define data size ----
n_cluster <- 30
ppc <- 3
n_point <- n_cluster*ppc
n_species <- 12
n_visit <- rep(4, n_point)

# Define covariates ----
clusterID <- rep(c(1:n_cluster), ppc)
sp_cov1 <- runif(n_species) - .5
sp_cov2 <- runif(n_species) - .5
pt_cov1 <- runif(n_point) - .5
cl.cov.1 <- runif(n_cluster) - .5
pt_cov2 <- rep(NA, n_point)
for(i in 1:n_point){
    pt_cov2[i] <- cl.cov.1[clusterID[i]]
}

vis_cov1 <- matrix(data = NA, nrow = n_point, ncol = max(n_visit))
for(i in 1:n_point){
    vis_cov1[i, ] <- runif(n_visit[i]) - .5
}

exists_visit <- matrix(data = 0, nrow = n_point, ncol = max(n_visit))
for(i in 1:n_point){
    for(j in 1:n_visit[i]){
        exists_visit[i,j] <- 1
    }
}

# Define hyperparameters ----
# (minor reorganisation, bringing det.hyper to top)
occ.hyper <- list(b0 = c(0, .5), 
                  b1 = c(0, 1), 
                  b2 = c(1, .5), 
                  b3 = c(-1, 1), 
                  b4 = c(0, .2),
                  b5 = c(1,2), 
                  b6 = c(2, 2))

det.hyper <- list(d0 = c(-2, .5), 
                  d1 = c(0, 1), 
                  d2 = c(0, .5), 
                  d3 = c(1, 1), 
                  d4 = c(0, .2),
                  d5 = c(2, .5))

b0 <- rnorm(n_species, occ.hyper$b0[1], occ.hyper$b0[2])
b1 <- matrix(data = rnorm(n_cluster*n_species, occ.hyper$b1[1], occ.hyper$b1[2]), nrow = n_cluster)
b2 <- rnorm(n_species, occ.hyper$b2[1], occ.hyper$b2[2])
b3 <- rnorm(n_species, occ.hyper$b3[1], occ.hyper$b3[2])
b4 <- rnorm(n_species, occ.hyper$b4[1], occ.hyper$b4[2])
b5 <- rnorm(n_species, occ.hyper$b5[1], occ.hyper$b5[2])
b6 <- rnorm(n_species, occ.hyper$b6[1], occ.hyper$b6[2])

d0 <- rnorm(n_species, det.hyper$d0[1], det.hyper$d0[2])
d1 <- matrix(data = rnorm(n_point*n_species, det.hyper$d1[1], det.hyper$d1[2]), nrow = n_point)
d2 <- rnorm(n_species, det.hyper$d2[1], det.hyper$d2[2])
d3 <- rnorm(n_species, det.hyper$d3[1], det.hyper$d3[2])
d4 <- rnorm(n_species, det.hyper$d4[1], det.hyper$d4[2])
d5 <- rnorm(n_species, det.hyper$d5[1], det.hyper$d5[2])

# Simulate parameters from hyperparameters
logit.occ <- psi <- matrix(NA, nrow = n_point, ncol = n_species)
for(i in 1:n_point){
    for(k in 1:n_species){
        logit.occ[i, k] <- b0[k] + b1[clusterID[i],k] + b4[k]*pt_cov1[i] 
        psi[i, k] <- boot::inv.logit(logit.occ[i, k])
    }
}

logit.det <- theta <- array(NA, dim = c(n_point, max(n_visit), n_species))
for(i in 1:n_point){
    for(j in 1:n_visit[i]){
        for(k in 1:n_species){
            logit.det[i, j, k] <- d0[k] + d5[k]*vis_cov1[i,j]  
            theta[i, j, k] <- boot::inv.logit(logit.det[i, j, k])
        }
    }
}


# simulate data from parameters
Z <- matrix(NA, nrow = n_point, ncol = n_species)
for(i in 1:n_point){
    for(k in 1:n_species){
        Z[i, k] <- rbinom(1, 1, psi[i, k])
    }
}

det_data <- array(NA, dim = c(n_point, max(n_visit), n_species))
for(i in 1:n_point){
    for(j in 1:n_visit[i]){
        for(k in 1:n_species){
            det_data[i,j,k] <- Z[i, k] * rbinom(1, 1, theta[i,j,k])
        }
    }
}


# Stan analysis
Q <- apply(det_data, c(1,3), function(x){return(as.numeric(sum(x) > 0))})



library(rstan)
stan.model <- '
data {
    int<lower=1> n_point; //number of sites
    int<lower=1> n_visit; //fixed number of visits
    int<lower=1> n_species; //number of species
    int<lower=1> n_cluster; //number of clusters

    int<lower=0, upper=1> det_data[n_point, n_visit, n_species]; //detection history
    int<lower=0, upper=1> Q[n_point, n_species]; //at least one detection

    int<lower=0, upper=n_cluster> clusterID[n_point]; //cluster identifier (for random effects)
    vector[n_species] sp_cov1; //species covariate 1
    vector[n_species] sp_cov2; //species covariate 2
    vector[n_point] pt_cov1; //point covariate 1
    vector[n_point] pt_cov2; //point covariate 2
    matrix[n_point, n_visit] vis_cov1; //visit covariate 1
}

parameters {
    real mu_b0;
    real<lower=0> sigma_b0;
    vector[n_species] b0_raw;

    real<lower=0> sigma_b1;
    matrix[n_species, n_cluster] b1_raw;

    real b2;

    real b3;

    real mu_b4;
    real<lower=0> sigma_b4;
    vector[n_species] b4_raw;

    real mu_b5;
    real<lower=0> sigma_b5;
    vector[n_species] b5_raw;

    real mu_b6;
    real<lower=0> sigma_b6;
    vector[n_species] b6_raw;

    real mu_d0;
    real<lower=0> sigma_d0;
    vector[n_species] d0_raw;

    real<lower=0> sigma_d1;
    matrix[n_species, n_point] d1_raw;

    real d2;
    real d3;
    real mu_d4;
    real<lower=0> sigma_d4;
    vector[n_species] d4_raw;

    real mu_d5;
    real<lower=0> sigma_d5;
    vector[n_species] d5_raw;
}

transformed parameters{
    vector[n_species] b0 = mu_b0 + b0_raw * sigma_b0;
    matrix[n_species, n_cluster] b1 = b1_raw * sigma_b1;

    vector[n_species] b4 = mu_b4 + b4_raw * sigma_b4;
    vector[n_species] b5 = mu_b5 + b5_raw * sigma_b5;
    vector[n_species] b6 = mu_b6 + b6_raw * sigma_b6;

    vector[n_species] d0 = mu_d0 + d0_raw * sigma_d0;
    matrix[n_species, n_point] d1 = d1_raw * sigma_d1;

    vector[n_species] d4 = mu_d4 + d4_raw * sigma_d4;
    vector[n_species] d5 = mu_d5 + d5_raw * sigma_d5;

    real logit_psi[n_point, n_species];
    real logit_theta[n_point, n_visit, n_species];
    matrix[n_point, n_species] log_prob_increment;
for(i in 1:n_point){
        for(k in 1:n_species){
            logit_psi[i,k] = b0[k] + b1[k, clusterID[i]] + b4[k]*pt_cov1[i]; 
        }
    }
    for(i in 1:n_point){
        for(j in 1:n_visit){
            for(k in 1:n_species){
                logit_theta[i,j,k] = d0[k] + d2*sp_cov1[k] + d1[k,i] +
                    d5[k]*vis_cov1[i,j];
            }
        }
    }

for(i in 1:n_point){
    for(k in 1:n_species){
        if(Q[i,k] == 1)
            log_prob_increment[i,k] = log_inv_logit(logit_psi[i,k]) + 
                bernoulli_logit_lpmf(det_data[i,1,k] | logit_theta[i,1,k]) + 
                bernoulli_logit_lpmf(det_data[i,2,k] | logit_theta[i,2,k]) +
                bernoulli_logit_lpmf(det_data[i,3,k] | logit_theta[i,3,k]) +
                bernoulli_logit_lpmf(det_data[i,4,k] | logit_theta[i,4,k]);
        else
            log_prob_increment[i,k] = log_sum_exp(log_inv_logit(logit_psi[i,k]) + 
                    log1m_inv_logit(logit_theta[i,1,k]) + 
                    log1m_inv_logit(logit_theta[i,2,k]) + 
                    log1m_inv_logit(logit_theta[i,3,k]) + 
                    log1m_inv_logit(logit_theta[i,4,k]), 
                log1m_inv_logit(logit_psi[i,k]));
        }
    }
}

model {
    //Hyper-priors:
    mu_b0 ~ normal(0,10);
    b2 ~ normal(0,10);
    b3 ~ normal(0,10);
    mu_b4 ~ normal(0,10);
    mu_b5 ~ normal(0,10);
    mu_b6 ~ normal(0,10);

    mu_d0 ~ normal(0,10);
    d2 ~ normal(0,10);
    d3 ~ normal(0,10);
    mu_d4 ~ normal(0,10);
    mu_d5 ~ normal(0,10);

    sigma_b0 ~ normal(0,10);
    sigma_b1 ~ normal(0,10);
    sigma_b4 ~ normal(0,10);
    sigma_b5 ~ normal(0,10);
    sigma_b6 ~ normal(0,10);

    sigma_d0 ~ normal(0,10);
    sigma_d1 ~ normal(0,10);
    sigma_d4 ~ normal(0,10);
    sigma_d5 ~ normal(0,10);

    //Random Effects
    b0_raw ~ normal(0, 1);
    to_vector(b1_raw) ~ normal(0, 1);
    b4_raw ~ normal(0, 1);
    b5_raw ~ normal(0, 1);
    b6_raw ~ normal(0, 1);

    d0_raw ~ normal(0, 1);
    to_vector(d1_raw) ~ normal(0, 1);
    d4_raw ~ normal(0, 1);
    d5_raw ~ normal(0, 1);

    //Likelihood (data level)
    target += sum(log_prob_increment);
}'

# Fit model ----
stan.data <- list(n_point = n_point, n_species = n_species, 
                  n_cluster = n_cluster, n_visit = 4,
                  det_data = det_data, 
                  Q = Q, 
                  clusterID = clusterID, 
                  sp_cov1 = sp_cov1, sp_cov2 = sp_cov2, 
                  pt_cov1 = pt_cov1, pt_cov2 = pt_cov2,
                  vis_cov1 = vis_cov1)

nc <- 4

start.time <- Sys.time()
stan.samples <- stan(model_code = stan.model, data = stan.data, iter = 2000, chains = nc, cores = nc,
                     pars = c('logit_psi', 'logit_theta', 'log_prob_increment', 'd1_raw', 'd1', 'b1_raw', 'b1'),
                     include = FALSE)
elapsed <- Sys.time() - start.time

draws_fit <- as_draws(stan.samples) 
draws_fit %>% 
    subset_draws(., variable = c("mu_b0", "sigma_b0", "mu_d0", "sigma_d0", 
                                 "mu_d5", "sigma_d5")) %>%
    summarise_draws() %>%
    mutate(ACTUAL = c(occ.hyper["b0"], det.hyper[c("d0", "d5")]) %>% unlist)

## Model without d1[k,i] ----
stan.model_2 <- '
data {
int<lower=1> n_point; //number of sites
int<lower=1> n_visit; //fixed number of visits
int<lower=1> n_species; //number of species
int<lower=1> n_cluster; //number of clusters

int<lower=0, upper=1> det_data[n_point, n_visit, n_species]; //detection history
int<lower=0, upper=1> Q[n_point, n_species]; //at least one detection

int<lower=0, upper=n_cluster> clusterID[n_point]; //cluster identifier (for random effects)
vector[n_species] sp_cov1; //species covariate 1
vector[n_species] sp_cov2; //species covariate 2
vector[n_point] pt_cov1; //point covariate 1
vector[n_point] pt_cov2; //point covariate 2
matrix[n_point, n_visit] vis_cov1; //visit covariate 1
}

parameters {
real mu_b0;
real<lower=0> sigma_b0;
vector[n_species] b0_raw;

real<lower=0> sigma_b1;
matrix[n_species, n_cluster] b1_raw;

real b2;

real b3;

real mu_b4;
real<lower=0> sigma_b4;
vector[n_species] b4_raw;

real mu_b5;
real<lower=0> sigma_b5;
vector[n_species] b5_raw;

real mu_b6;
real<lower=0> sigma_b6;
vector[n_species] b6_raw;

real mu_d0;
real<lower=0> sigma_d0;
vector[n_species] d0_raw;

real<lower=0> sigma_d1;
matrix[n_species, n_point] d1_raw;

real d2;
real d3;
real mu_d4;
real<lower=0> sigma_d4;
vector[n_species] d4_raw;

real mu_d5;
real<lower=0> sigma_d5;
vector[n_species] d5_raw;
}

transformed parameters{
vector[n_species] b0 = mu_b0 + b0_raw * sigma_b0;
matrix[n_species, n_cluster] b1 = b1_raw * sigma_b1;

vector[n_species] b4 = mu_b4 + b4_raw * sigma_b4;
vector[n_species] b5 = mu_b5 + b5_raw * sigma_b5;
vector[n_species] b6 = mu_b6 + b6_raw * sigma_b6;

vector[n_species] d0 = mu_d0 + d0_raw * sigma_d0;
matrix[n_species, n_point] d1 = d1_raw * sigma_d1;

vector[n_species] d4 = mu_d4 + d4_raw * sigma_d4;
vector[n_species] d5 = mu_d5 + d5_raw * sigma_d5;

real logit_psi[n_point, n_species];
real logit_theta[n_point, n_visit, n_species];
matrix[n_point, n_species] log_prob_increment;
for(i in 1:n_point){
for(k in 1:n_species){
logit_psi[i,k] = b0[k] + b1[k, clusterID[i]] + b4[k]*pt_cov1[i]; 
}
}
for(i in 1:n_point){
for(j in 1:n_visit){
for(k in 1:n_species){
logit_theta[i,j,k] = d0[k] + d2*sp_cov1[k] +
d5[k]*vis_cov1[i,j];
}
}
}

for(i in 1:n_point){
for(k in 1:n_species){
if(Q[i,k] == 1)
log_prob_increment[i,k] = log_inv_logit(logit_psi[i,k]) + 
bernoulli_logit_lpmf(det_data[i,1,k] | logit_theta[i,1,k]) + 
bernoulli_logit_lpmf(det_data[i,2,k] | logit_theta[i,2,k]) +
bernoulli_logit_lpmf(det_data[i,3,k] | logit_theta[i,3,k]) +
bernoulli_logit_lpmf(det_data[i,4,k] | logit_theta[i,4,k]);
else
log_prob_increment[i,k] = log_sum_exp(log_inv_logit(logit_psi[i,k]) + 
log1m_inv_logit(logit_theta[i,1,k]) + 
log1m_inv_logit(logit_theta[i,2,k]) + 
log1m_inv_logit(logit_theta[i,3,k]) + 
log1m_inv_logit(logit_theta[i,4,k]), 
log1m_inv_logit(logit_psi[i,k]));
}
}
}

model {
//Hyper-priors:
mu_b0 ~ normal(0,10);
b2 ~ normal(0,10);
b3 ~ normal(0,10);
mu_b4 ~ normal(0,10);
mu_b5 ~ normal(0,10);
mu_b6 ~ normal(0,10);

mu_d0 ~ normal(0,10);
d2 ~ normal(0,10);
d3 ~ normal(0,10);
mu_d4 ~ normal(0,10);
mu_d5 ~ normal(0,10);

sigma_b0 ~ normal(0,10);
sigma_b1 ~ normal(0,10);
sigma_b4 ~ normal(0,10);
sigma_b5 ~ normal(0,10);
sigma_b6 ~ normal(0,10);

sigma_d0 ~ normal(0,10);
sigma_d1 ~ normal(0,10);
sigma_d4 ~ normal(0,10);
sigma_d5 ~ normal(0,10);

//Random Effects
b0_raw ~ normal(0, 1);
to_vector(b1_raw) ~ normal(0, 1);
b4_raw ~ normal(0, 1);
b5_raw ~ normal(0, 1);
b6_raw ~ normal(0, 1);

d0_raw ~ normal(0, 1);
to_vector(d1_raw) ~ normal(0, 1);
d4_raw ~ normal(0, 1);
d5_raw ~ normal(0, 1);

//Likelihood (data level)
target += sum(log_prob_increment);
}'

start.time <- Sys.time()
stan.samples_2 <- stan(model_code = stan.model_2, data = stan.data, iter = 2000, chains = nc, cores = nc,
                     pars = c('logit_psi', 'logit_theta', 'log_prob_increment', 'd1_raw', 'd1', 'b1_raw', 'b1'),
                     include = FALSE)
elapsed <- Sys.time() - start.time

draws_fit_2 <- as_draws(stan.samples_2) 
draws_fit_2 %>% 
    subset_draws(., variable = c("mu_b0", "sigma_b0", "mu_d0", "sigma_d0", 
                                 "mu_d5", "sigma_d5")) %>%
    summarise_draws() %>%
    mutate(ACTUAL = c(occ.hyper["b0"], det.hyper[c("d0", "d5")]) %>% unlist)
# better..

## Run second model with another seed ----
set.seed(1001)

clusterID <- rep(c(1:n_cluster), ppc)
sp_cov1 <- runif(n_species) - .5
sp_cov2 <- runif(n_species) - .5
pt_cov1 <- runif(n_point) - .5
cl.cov.1 <- runif(n_cluster) - .5
pt_cov2 <- rep(NA, n_point)
for(i in 1:n_point){
    pt_cov2[i] <- cl.cov.1[clusterID[i]]
}

vis_cov1 <- matrix(data = NA, nrow = n_point, ncol = max(n_visit))
for(i in 1:n_point){
    vis_cov1[i, ] <- runif(n_visit[i]) - .5
}

exists_visit <- matrix(data = 0, nrow = n_point, ncol = max(n_visit))
for(i in 1:n_point){
    for(j in 1:n_visit[i]){
        exists_visit[i,j] <- 1
    }
}

# (minor reorganisation, bringing det.hyper to top)
occ.hyper <- list(b0 = c(0, .5), 
                  b1 = c(0, 1), 
                  b2 = c(1, .5), 
                  b3 = c(-1, 1), 
                  b4 = c(0, .2),
                  b5 = c(1,2), 
                  b6 = c(2, 2))

det.hyper <- list(d0 = c(-2, .5), 
                  d1 = c(0, 1), 
                  d2 = c(0, .5), 
                  d3 = c(1, 1), 
                  d4 = c(0, .2),
                  d5 = c(2, .5))

b0 <- rnorm(n_species, occ.hyper$b0[1], occ.hyper$b0[2])
b1 <- matrix(data = rnorm(n_cluster*n_species, occ.hyper$b1[1], occ.hyper$b1[2]), nrow = n_cluster)
b2 <- rnorm(n_species, occ.hyper$b2[1], occ.hyper$b2[2])
b3 <- rnorm(n_species, occ.hyper$b3[1], occ.hyper$b3[2])
b4 <- rnorm(n_species, occ.hyper$b4[1], occ.hyper$b4[2])
b5 <- rnorm(n_species, occ.hyper$b5[1], occ.hyper$b5[2])
b6 <- rnorm(n_species, occ.hyper$b6[1], occ.hyper$b6[2])

d0 <- rnorm(n_species, det.hyper$d0[1], det.hyper$d0[2])
d1 <- matrix(data = rnorm(n_point*n_species, det.hyper$d1[1], det.hyper$d1[2]), nrow = n_point)
d2 <- rnorm(n_species, det.hyper$d2[1], det.hyper$d2[2])
d3 <- rnorm(n_species, det.hyper$d3[1], det.hyper$d3[2])
d4 <- rnorm(n_species, det.hyper$d4[1], det.hyper$d4[2])
d5 <- rnorm(n_species, det.hyper$d5[1], det.hyper$d5[2])

# Simulate parameters from hyperparameters
logit.occ <- psi <- matrix(NA, nrow = n_point, ncol = n_species)
for(i in 1:n_point){
    for(k in 1:n_species){
        logit.occ[i, k] <- b0[k] + b1[clusterID[i],k] + b4[k]*pt_cov1[i] 
        psi[i, k] <- boot::inv.logit(logit.occ[i, k])
    }
}

logit.det <- theta <- array(NA, dim = c(n_point, max(n_visit), n_species))
for(i in 1:n_point){
    for(j in 1:n_visit[i]){
        for(k in 1:n_species){
            logit.det[i, j, k] <- d0[k] + d5[k]*vis_cov1[i,j]  
            theta[i, j, k] <- boot::inv.logit(logit.det[i, j, k])
        }
    }
}


# simulate data from parameters
Z <- matrix(NA, nrow = n_point, ncol = n_species)
for(i in 1:n_point){
    for(k in 1:n_species){
        Z[i, k] <- rbinom(1, 1, psi[i, k])
    }
}

det_data <- array(NA, dim = c(n_point, max(n_visit), n_species))
for(i in 1:n_point){
    for(j in 1:n_visit[i]){
        for(k in 1:n_species){
            det_data[i,j,k] <- Z[i, k] * rbinom(1, 1, theta[i,j,k])
        }
    }
}


# Stan analysis
Q <- apply(det_data, c(1,3), function(x){return(as.numeric(sum(x) > 0))})

stan.data_3 <- list(n_point = n_point, n_species = n_species, 
                  n_cluster = n_cluster, n_visit = 4,
                  det_data = det_data, 
                  Q = Q, 
                  clusterID = clusterID, 
                  sp_cov1 = sp_cov1, sp_cov2 = sp_cov2, 
                  pt_cov1 = pt_cov1, pt_cov2 = pt_cov2,
                  vis_cov1 = vis_cov1)

# Model on diff sim ----
start.time <- Sys.time()
stan.samples_3 <- stan(model_code = stan.model_2, data = stan.data_3, iter = 2000, chains = nc, cores = nc,
                       pars = c('logit_psi', 'logit_theta', 'log_prob_increment', 'd1_raw', 'd1', 'b1_raw', 'b1'),
                       include = FALSE)
elapsed <- Sys.time() - start.time

draws_fit_3 <- as_draws(stan.samples_3) 
draws_fit_3 %>% 
    subset_draws(., variable = c("mu_b0", "sigma_b0", "mu_d0", "sigma_d0", 
                                 "mu_d5", "sigma_d5")) %>%
    summarise_draws() %>%
    mutate(ACTUAL = c(occ.hyper["b0"], det.hyper[c("d0", "d5")]) %>% unlist)
