# Checking that reduceSum is doing the same sampling as the non-parallelised 
# version ("occupancyMod_ragged"). 
#
# Simulating under:
# logit.occ[i, k] <- b0[k] + b1[clusterID[i],k] + b4[k]*pt_cov1[i]
# logit.det[i, j, k] <- d0[k] + d5[k]*vis_cov1[i,j] + d2*sp_cov1[k] 
#
# Stan models should match this. 
#
# Runs in cmdstanr
# i indexes points, j indexes visits, k indexes species.
#
# partial_sum has been tested in isolation ("compare_ragged_partialSum.R") and 
# produces correct result

# Housekeeping ----
library("cmdstanr"); library("dplyr"); library("posterior")
num_chains <- 4

set.seed(101)
# Define data size ----
n_cluster <- 30
ppc <- 3
n_point <- n_cluster*ppc
n_species <- 30
n_visit <- rep(4, n_point)
print(paste0("nrow: ", n_species*n_point))

# Define covariates
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

# Define parameters ----
# Hyperparameters
occ.hyper <- list(b0 = c(0, .5), b1 = c(0, 1), b2 = c(1, .5), b3 = c(-1, 1), 
                  b4 = c(0, .2), b5 = c(1,2), b6 = c(2, 2))
b0 <- rnorm(n_species, occ.hyper$b0[1], occ.hyper$b0[2])
b1 <- matrix(data = rnorm(n_cluster*n_species, occ.hyper$b1[1], occ.hyper$b1[2]), 
             nrow = n_cluster)
b2 <- rnorm(n_species, occ.hyper$b2[1], occ.hyper$b2[2])
b3 <- rnorm(n_species, occ.hyper$b3[1], occ.hyper$b3[2])
b4 <- rnorm(n_species, occ.hyper$b4[1], occ.hyper$b4[2])
b5 <- rnorm(n_species, occ.hyper$b5[1], occ.hyper$b5[2])
b6 <- rnorm(n_species, occ.hyper$b6[1], occ.hyper$b6[2])

det.hyper <- list(d0 = c(-2, .5), d1 = c(0, 1), d2 = c(0, .5), d3 = c(1, 1), 
                  d4 = c(0, .2), d5 = c(2, .5))
d0 <- rnorm(n_species, det.hyper$d0[1], det.hyper$d0[2])
d1 <- matrix(data = rnorm(n_point*n_species, det.hyper$d1[1], det.hyper$d1[2]),
             nrow = n_point)
d2 <- .3 #rnorm(n_species, det.hyper$d2[1], det.hyper$d2[2])
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
            logit.det[i, j, k] <- d0[k] + d5[k]*vis_cov1[i,j] + d2*sp_cov1[k] 
            theta[i, j, k] <- boot::inv.logit(logit.det[i, j, k])
        }
    }
}

# Simulate data ----
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

Q <- apply(det_data, c(1,3), function(x){return(as.numeric(sum(x) > 0))})

stan.data_J <- list(n_point = n_point, n_species = n_species, 
                    n_cluster = n_cluster, n_visit = 4,
                    det_data = det_data, 
                    Q = Q, 
                    clusterID = clusterID, 
                    sp_cov1 = sp_cov1, sp_cov2 = sp_cov2, 
                    pt_cov1 = pt_cov1, pt_cov2 = pt_cov2,
                    vis_cov1 = vis_cov1)

# Format for Stan ----
det_df <- apply(det_data, 3, as_data_frame) %>%
    bind_rows(., .id="id_species") %>%
    rename(det_1 = V1, det_2 = V2, det_3 = V3, det_4 = V4) %>%
    mutate(id_species = as.numeric(id_species),
           id_point = rep(1:n_point, n_species), 
           id_cluster = clusterID[id_point], 
           id_sp_cl = as.integer(as.factor(paste0("s", id_species, "c", id_cluster))),
           id_sp_pt = as.integer(as.factor(paste0("s", id_species, "c", id_point))),
           vis_cov1_1 = vis_cov1[id_point, 1], 
           vis_cov1_2 = vis_cov1[id_point, 2], 
           vis_cov1_3 = vis_cov1[id_point, 3],
           vis_cov1_4 = vis_cov1[id_point, 4], 
           sp_cov1 = sp_cov1[id_species], 
           sp_cov2 = sp_cov2[id_species], 
           pt_cov1 = pt_cov1[id_point], 
           pt_cov2 = pt_cov2[id_point]) %>%
    mutate(Q = rowSums(select(.,det_1, det_2, det_3, det_4)) > 0) %>%
    arrange(desc(Q), id_species, id_point)

stan_data <- list(
    n_visit = 4, 
    n_species = length(unique(det_df$id_species)),
    n_sp_cl = length(unique(det_df$id_sp_cl)),
    n_sp_pt = length(unique(det_df$id_sp_pt)),
    n_pt = length(unique(det_df$id_point)),
    n_visit_max = 4, 
    n_tot = nrow(det_df),
    n_1 = sum(det_df$Q),
    id_pt = det_df$id_point,
    id_cl = det_df$id_cluster,
    id_sp = det_df$id_species,
    id_sp_cl = det_df$id_sp_cl,
    # note: currently just give it a vector of length nrows, with replicated 
    # entries for covariates. This could be changed, but don't think it will 
    # meaningfully alter memory requirements
    sp_cov1 = det_df$sp_cov1,
    sp_cov2 = det_df$sp_cov2,
    pt_cov1 = det_df$pt_cov1,
    pt_cov2 = det_df$pt_cov2,
    vis_cov1 = as.matrix(det_df[,paste0("vis_cov1_", 1:4)]),
    det_data = as.matrix(det_df[,paste0("det_", 1:4)]), 
    Q = det_df$Q, 
    grainsize = 1350)

# Run ragged model ----
mod_R <- cmdstan_model("stan_files/occupancyMod_ragged.stan")
time_R <- system.time(samps_R <- mod_R$sample(data = stan_data,
                                              num_chains = num_chains,
                                              num_cores = num_chains))

# Run 1 thread ----
set_num_threads(1)
mod_P <- cmdstan_model("stan_files/occupancyMod_ragged_parallel.stan", 
                       threads=TRUE)
time_P <- system.time(samps_P <- mod_P$sample(data = stan_data,
                                              num_chains = num_chains,
                                              num_cores = num_chains))

# Get sample subsets ----
sub_R <- samps_R$draws() %>%
    subset_draws(., variable = c("mu_b0", "sigma_b0", "mu_d0", "sigma_d0"))
sub_P <- samps_P$draws() %>%
    subset_draws(., variable = c("mu_b0", "sigma_b0", "mu_d0", "sigma_d0"))

# Save ----
# Reduce filesize by dropping full sample object and save workspace
rm(list = c("samps_R", "samps_P"))
save.image("workspace/check_reduceSum_1thread.RData")
