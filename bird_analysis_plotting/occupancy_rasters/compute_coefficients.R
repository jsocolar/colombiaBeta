remotes::install_github("jsocolar/flocker")

library(flocker)
library("dplyr"); library("brms")
library(data.table)

# Assemble dataframe to predict the spatially constant parts of the response
uc1 <- readRDS("outputs/unit_covariates.RDS") %>%
  filter(!duplicated(species)) %>%
  select(
    c(
      "species", "Family", "lowland",
      "mountain_barrier", "valley_barrier",
      "forestPresent", "forestSpecialist",
      "tfSpecialist", "dryForestPresent",
      "floodDrySpecialist", "aridPresent",
      "migratory", "dietInvert", "dietCarn",
      "dietFruitNect", "dietGran",
      "elevMedian", "elevBreadth", "mass",
      "elevMedian_x_forestPresent", "elevMedian_x_forestSpecialist"
    )
  )

# Read in fitted flocker model
fm3 <- readRDS("outputs/ffit3_copy.RDS")
# Model was fit from `more-models` branch; we will switch so that we can
# use `fitted_flocker()` which is not yet updated to work with that branch
attributes(fm3)$lik_type <- "V"
# Munging because when we built fm3 we thinned manually
fm3$fit@sim$thin <- 100
fm3$fit@sim$n_save[1] <- 10

# can save 1.6GB by only retaining essential info in $data
fm3$data <- fm3$data %>%
    group_by(species) %>%
    slice(1)

##### prediction function #####
# helper for monotonic effects, translated from brms-generated stan
stan_mo <- function(scale, i) {
  out <- rep(NA, length(i))
  for(j in 1:length(i)){
    if (i[j] == 0) {
      out[j] <- 0
    } else {
      out[j] <- length(scale) * sum(scale[1 : i[j]])
    }
  }
  out
}

compute_coefs <- function(uc1, fm3, iter = NULL){
    ## TODO: update to run for varying number of iterations (at moment set to 10)
    # Linear predictor part that is constant within species 
    uc <- uc1 %>%
        mutate(
            distance_bin = 1,
            relev = 0,
            relev2 = 0,
            lowland_x_relev = 0,
            lowland_x_relev2 = 0, pasture = 1
        ) 
    uc <- bind_rows(uc %>% mutate(pasture = -1),
                    uc %>% mutate(pasture = 1)) %>%
        mutate(
            mountainBarrier_x_pasture = mountain_barrier * pasture,
            valleyBarrier_x_pasture = valley_barrier * pasture,
            forestPresent_x_pasture = forestPresent * pasture,
            forestSpecialist_x_pasture = forestSpecialist * pasture,
            tfSpecialist_x_pasture = tfSpecialist * pasture,
            dryForestPresent_x_pasture = dryForestPresent * pasture,
            floodDrySpecialist_x_pasture = floodDrySpecialist * pasture,
            aridPresent_x_pasture = aridPresent * pasture,
            migratory_x_pasture = migratory * pasture,
            dietInvert_x_pasture = dietInvert * pasture,
            dietCarn_x_pasture = dietCarn * pasture,
            dietFruitNect_x_pasture = dietFruitNect * pasture,
            dietGran_x_pasture = dietGran * pasture,
            elevMedian_x_pasture = elevMedian * pasture,
            elevBreadth_x_pasture = elevBreadth * pasture,
            mass_x_pasture = mass * pasture,
            elevMedian_x_forestPresent_x_pasture = elevMedian_x_forestPresent * pasture,
            elevMedian_x_forestSpecialist_x_pasture = elevMedian_x_forestSpecialist * pasture
        )
    
    # The detection components do not matter since we will not be conditioning on any
    # observed histories.
    ec <- list(
        time = matrix(0, nrow(uc), 4),
        time_x_elev = matrix(0, nrow(uc), 4),
        id_spObs = matrix(0, nrow(uc), 4),
        obsSM = matrix(0, nrow(uc), 4),
        obsDE = matrix(0, nrow(uc), 4),
        obsJG = matrix(0, nrow(uc), 4)
    )
    obs <- matrix(0, nrow(uc), 4)
    newdata <- flocker::make_flocker_data(obs, uc, ec)
    
    lpo <- t(
        brms::posterior_linpred(
            fm3, 
            dpar = "occ",
            draw_ids = 1:10,
            ##### TODO: make sure that when fm3 is multi-chain, the draw_ids is getting out the same
            ##### iterations as iter does elsewhere!
            
            newdata = newdata$data[1:nrow(uc), ], 
            re_formula = ~ (1 + pasture |g1| species) + (1 + pasture |g2| Family),
            allow_new_levels = FALSE
        )
    )
    
    # Component of linear predictor that is variable across space
    ## extract subregion and cluster effects
    sd_subregion <- rstan::extract(fm3$fit,"sd_species_subregion__occ_Intercept", permuted = FALSE)
    sd_cluster <- rstan::extract(fm3$fit,"sd_species_cluster__occ_Intercept", permuted = FALSE)
    
    ## extract relev terms
    uc_index <- 1:nrow(uc1)
    template_mat <- matrix(1, nrow = length(uc_index), ncol = 10)
   
    relev_term1 <- template_mat %*%
        diag(as.vector(rstan::extract(fm3$fit, "b_occ_relev", permuted = FALSE))) 
    
    relev_term2 <- matrix(uc$lowland[uc_index], ncol = 10, nrow = length(uc_index)) %*% 
        diag(as.vector(rstan::extract(fm3$fit, "b_occ_lowland_x_relev", permuted = FALSE)))
    
    relev_term3 <- t(rstan::extract(fm3$fit, 
                           paste0("r_species__occ[", uc$species[uc_index], ",relev]"), 
                           permuted = FALSE)[1:10, 1, 1:1614])

    relev_term <- relev_term1 + relev_term2 + relev_term3
    
    ## extract relev2 terms
    relev2_term1 <- template_mat %*%
        diag(as.vector(rstan::extract(fm3$fit, "b_occ_relev2", permuted = FALSE))) 
    
    relev2_term2 <- matrix(uc$lowland[uc_index], ncol = 10, nrow = length(uc_index)) %*% 
        diag(as.vector(rstan::extract(fm3$fit, "b_occ_lowland_x_relev2", permuted = FALSE)))
    
    relev2_term3 <- t(rstan::extract(fm3$fit, 
                                  paste0("r_species__occ[", uc$species[uc_index], ",relev2]"), 
                                  permuted = FALSE)[1:10, 1, 1:1614])
    
    relev2_term <- relev2_term1 + relev2_term2 + relev2_term3
    
    # extract monotonic distance effect
    mos <- matrix(NA, ncol = 10, nrow = 12)
    for(i in 1:10) {
        bsp <- rstan::extract(fm3$fit, "bsp_occ_modistance_bin", permuted = FALSE)[i]
        simo <- rstan::extract(fm3$fit, paste0("simo_occ_modistance_bin1[", c(1:11), "]"), permuted = FALSE)[i,1,]
        mos[,i] <- bsp * stan_mo(simo, c(0:11)) 
    }
    
    # format and return list
    rownames(relev_term) <- rownames(relev2_term) <- uc$species[uc_index]
    list(species = uc$species[uc_index],
         lpo_forest = lpo[1:(nrow(lpo)/2),],
         lpo_pasture = lpo[(nrow(lpo)/2 + 1):nrow(lpo),],
         relev_term = relev_term,
         relev2_term = relev2_term, 
         mos = mos,
         sd_subregion = sd_subregion,
         sd_cluster = sd_cluster)
}

out <- compute_coefs(uc1, fm3)
saveRDS(out, "outputs/lpo_and_coefs.rds")
