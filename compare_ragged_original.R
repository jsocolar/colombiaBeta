# Comparing "ragged" model formulation, and the original written by Jacob.
#
# Inherits from Jacob's occupancy_simulation2.R (version accessed 27Apr2020)- 
# simulation code largely follows the original, but model is different. 
#
# Simulating under:
# logit.occ[i, k] <- b0[k] + b1[clusterID[i],k] + b4[k]*pt_cov1[i]
# logit.det[i, j, k] <- d0[k] + d5[k]*vis_cov1[i,j] + d2*sp_cov1[k] 
#
# Stan models should match this. 
#
# Runs in rstan
# i indexes points, j indexes visits, k indexes species.
#
# I'm satisfied that they are producing the same result (and ragged code actually
# appears to run slightly faster). 90 points * 30 species has  8.07 mins vs 10.35
# 
# Note: final check that makes array ragged and checks correctly runs hasn't been 
# run recently (and code out of date). Not concerned about this, but for 
# completeness should be uncommented and updated. 

# Housekeeping ----
library("rstan"); library("dplyr"); library("reshape2");
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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
    det_data = as.matrix(det_df[,paste0("det_", 1:4)]))

# Run model ----
stan.model <- stan_model("stan_files/occupancyMod_ragged.stan")
start.time <- Sys.time()
stan.samples <- sampling(stan.model, 
                         data = stan_data,
                         iter = 2000,
                         chains = 4,
                         pars = c('logit_psi', 'logit_theta', 'd1_raw', 'd1', 'b1_raw', 'b1'),
                         include = FALSE)
(elapsed <- Sys.time() - start.time)

## Compare with original model ----
stan.data_J <- list(n_point = n_point, n_species = n_species, 
                    n_cluster = n_cluster, n_visit = 4,
                    det_data = det_data, 
                    Q = Q, 
                    clusterID = clusterID, 
                    sp_cov1 = sp_cov1, sp_cov2 = sp_cov2, 
                    pt_cov1 = pt_cov1, pt_cov2 = pt_cov2,
                    vis_cov1 = vis_cov1)

stan.model_J <- stan_model("stan_files/occupancyMod_original.stan")
start.time_J <- Sys.time()
stan.samples_J <- sampling(stan.model_J, 
                           data = stan.data_J,
                           iter = 2000,
                           chains = 4)
(elapsed_J <- Sys.time() - start.time_J)

library(posterior)
draws_R <- as_draws(stan.samples) 
draws_J <- as_draws(stan.samples_J) 
sub_R <- draws_R %>% 
    subset_draws(., variable = c("mu_b0", "sigma_b0", "mu_d0", "sigma_d0"))
sub_J <- draws_J %>% 
    subset_draws(., variable = c("mu_b0", "sigma_b0", "mu_d0", "sigma_d0"))
bind_draws(sub_R, sub_J, along="chain") %>% summarise_draws

summarise_draws(sub_R)
summarise_draws(sub_J)


# # Check model on ragged dataset ----
# det_df_2 <- det_df %>%
#     # reduce size a bit
#     filter(!id_cluster %in% c(15:30)) %>%
#     group_by(id_species) %>%
#     # remove 1 cluster from the first 10 species
#     filter(ifelse(id_species %in% 1:10, id_cluster!=14, id_cluster!=Inf)) %>%
#     ungroup %>%
#     arrange(desc(Q), id_species) %>%
#     # recalculate unique id:cluster combinations in dataset
#     mutate(id_sp_cl = as.integer(as.factor(paste0("s", id_species, "c", id_cluster))))
# 
# stan_data <- list(
#     n_visit = 4, 
#     n_species = length(unique(det_df_2$id_species)),
#     n_sp_cl = length(unique(det_df_2$id_sp_cl)),
#     n_sp_pt = length(unique(det_df_2$id_sp_pt)),
#     n_pt = length(unique(det_df_2$id_point)),
#     n_visit_max = 4, 
#     n_tot = nrow(det_df_2),
#     n_1 = sum(det_df_2$Q),
#     id_pt = det_df_2$id_point,
#     id_cl = det_df_2$id_cluster,
#     id_sp = det_df_2$id_species,
#     id_sp_cl = det_df_2$id_sp_cl,
#     sp_cov1 = det_df_2$sp_cov1,
#     sp_cov2 = det_df_2$sp_cov2,
#     pt_cov1 = det_df_2$pt_cov1,
#     pt_cov2 = det_df_2$pt_cov2,
#     vis_cov1 = as.matrix(det_df_2[,paste0("vis_cov1_", 1:4)]),
#     det_data = as.matrix(det_df_2[,paste0("det_", 1:4)]))
# 
# stan_data$id_sp_cl %>% unique %>% length
# stan.samples_3 <- sampling(stan.model, 
#                            data = stan_data,
#                            iter = 2000,
#                            chains = 1,
#                            pars = c('logit_psi', 'logit_theta', 'd1_raw', 'd1', 'b1_raw', 'b1'),
#                            include = FALSE)
# # runs fine! 
