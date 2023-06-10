library(cmdstanr)

library(flocker)
library("dplyr"); library("brms")

bird_stan_data9_1_package <- readRDS("/home/ec2-user/inputs/bird_stan_data9_package.RDS")

# Get out the integer data (response and integer covariates) as a dataframe
integer_data <- do.call(cbind, bird_stan_data9_1_package$data$integer_data)
colnames(integer_data) <- names(bird_stan_data9_1_package$data$integer_data)
integer_data <- as.data.frame(integer_data)

length(unique(integer_data$id_sp[integer_data$Q == 1]))

# Get out the observations
obs <- integer_data[, paste0("det_data.v", 1:4)] %>%  as.matrix()
obs[obs == -1] <- NA

# Get out the real unit data as a dataframe
real_data_unit_cols <- c(
  "relev", "relev2",
  "lowland_x_relev", "lowland_x_relev2",
  "elevMedian", "elevBreadth",
  "mass",
  "elevMedian_x_forestPresent", "elevMedian_x_forestSpecialist",
  "elevMedian_x_pasture", "elevBreadth_x_pasture",
  "mass_x_pasture", "elevMedian_x_forestPresent_x_pasture", 
  "elevMedian_x_forestSpecialist_x_pasture"
)

uc_real <- do.call(cbind, bird_stan_data9_1_package$data[real_data_unit_cols])
names(uc_real) <- real_data_unit_cols
uc_real <- as.data.frame(uc_real)

# Get out the integer uc data 
uc_integer_cols <- c(
#  "id_spCl", "id_spSr", "id_sp", "id_fam", 
  "lowland",
  "pasture", 
  "mountain_barrier", "valley_barrier",
  "forestPresent", "forestSpecialist", "tfSpecialist", "dryForestPresent",
  "floodDrySpecialist","aridPresent",
  "migratory",
  "dietInvert", "dietCarn", "dietFruitNect", "dietGran", 
  "mountainBarrier_x_pasture", "valleyBarrier_x_pasture",
  "forestPresent_x_pasture", "forestSpecialist_x_pasture", 
  "tfSpecialist_x_pasture", "dryForestPresent_x_pasture", 
  "floodDrySpecialist_x_pasture", "aridPresent_x_pasture", 
  "migratory_x_pasture", 
  "dietInvert_x_pasture", "dietCarn_x_pasture", "dietFruitNect_x_pasture",
  "dietGran_x_pasture"
)
uc_integer <- integer_data[,uc_integer_cols]
uc <- cbind(uc_integer, uc_real)
for (i in 1:ncol(uc)) {
  if (identical(unique(uc[,i])[order(unique(uc[,i]))], c(1, 2))) {
    uc[,i] <- c(-1, 1)[uc[,i]]
  }
}

birds <- readRDS("/home/ec2-user/inputs/birds.RDS")
assertthat::assert_that(identical(birds$elev_breadth_scaled, uc$elevBreadth))

# a <- brm(
#   Q ~ s(distance_from_range) + (1|species), 
#   family = bernoulli(), 
#   data = birds,
#   chains = 2,
#   cores = 2,
#   threads = 4,
#   backend = "cmdstanr",
#   refresh = 1)
# saveRDS(a, "/Users/jacob/Dropbox/Work/Colombia/Data/Analysis/checking_distance_relationships/a.RDS")
# marginal_smooths(a)

# b <- brm(
#   Q ~ s(distance_from_range) + (1|species), 
#   family = bernoulli(), 
#   data = birds[birds$distance_from_range > -150000, ],
#   chains = 2,
#   cores = 2,
#   threads = 4,
#   backend = "cmdstanr",
#   refresh = 1)
# saveRDS(b, "/Users/jacob/Dropbox/Work/Colombia/Data/Analysis/checking_distance_relationships/b.RDS")
# conditional_smooths(b)

birds$distance_bin <- 
  lapply(birds$distance_from_range, function(x){max(which(c(-Inf, seq(-60000, 140000, 20000)) < x))}) %>%
  unlist() %>%
  factor(ordered = TRUE)


uc <- cbind(birds[, c("species", "Family", "subregion", "cluster", "distance_bin")], uc)
uc$species_subregion <- paste0(uc$species, "__", uc$subregion)
uc$species_cluster <- paste0(uc$species, "__", uc$cluster)



# Get out event covariates
obsSM <- integer_data[,paste0("obsSM.", 1:4)]
obsJG <- integer_data[,paste0("obsJG.", 1:4)]
obsDE <- integer_data[,paste0("obsDE.", 1:4)]

obsSM[obsSM == 1] <- -1
obsSM[obsSM == 2] <- 1
obsDE[obsDE == 1] <- -1
obsDE[obsDE == 2] <- 1
obsJG[obsJG == 1] <- -1
obsJG[obsJG == 2] <- 1

ec <- list(
  time = bird_stan_data9_1_package$data$time,
  time_x_elev = bird_stan_data9_1_package$data$time_x_elev,
  id_spObs = birds[,paste0("sp_obs", 1:4)],
  obsSM = obsSM,
  obsDE = obsDE,
  obsJG = obsJG
)

fd <- make_flocker_data(obs, uc, ec)

bprior1 <- 
  prior(normal(-3,1), class = "b", coef = "Intercept", dpar = "mu") +
  prior(normal(0, 0.75), class = "b", coef = "pasture", dpar = "mu") +
  prior(normal(0, 0.5), class = "b", coef = "time", dpar = "mu") +
  prior(normal(0, 0.5), class = "b", coef = "elevMedian", dpar = "mu") +
  prior(normal(0, 0.5), class = "b", coef = "time_x_elev", dpar = "mu") +
  prior(normal(0, 0.5), class = "b", coef = "dietCarn", dpar = "mu") +
  prior(normal(0, 0.5), class = "b", coef = "mass", dpar = "mu") +
  prior(normal(0, 1), class = "b", coef = "migratory", dpar = "mu") +
  prior(normal(0, 0.25), class = "b", coef = "obsDE", dpar = "mu") +
  prior(normal(0, 0.25), class = "b", coef = "obsJG", dpar = "mu") +
  prior(normal(0, 0.25), class = "b", coef = "obsSM", dpar = "mu") +
  prior(normal(0, 2), class = "sd", coef = "Intercept", group = "species", dpar = "mu") +
  prior(normal(0, 2), class = "sd", coef = "Intercept", group = "Family", dpar = "mu") +
  prior(std_normal(), class = "sd", coef = "pasture", group = "species", dpar = "mu") +
  prior(std_normal(), class = "sd", coef = "pasture", group = "Family", dpar = "mu") +
  prior(std_normal(), class = "sd", coef = "time", group = "species", dpar = "mu")
bprior2 <-   
  prior(normal(-7,2.5), class = "b", coef = "Intercept", dpar = "occ") +
  prior(normal(0,5), class = "b", coef = "relev", dpar = "occ") +
  prior(normal(0,5), class = "b", coef = "relev2", dpar = "occ") +
  prior(normal(0,5), class = "b", coef = "lowland_x_relev", dpar = "occ") +
  prior(normal(0,5), class = "b", coef = "lowland_x_relev2", dpar = "occ") +
  prior("", class = "b", coef = "modistance_bin", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "lowland", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "pasture", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "aridPresent", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "aridPresent_x_pasture", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "dryForestPresent", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "dryForestPresent_x_pasture", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "elevBreadth", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "elevBreadth_x_pasture", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "elevMedian", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "elevMedian_x_pasture", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "elevMedian_x_forestPresent", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "elevMedian_x_forestPresent_x_pasture", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "elevMedian_x_forestSpecialist", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "elevMedian_x_forestSpecialist_x_pasture", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "floodDrySpecialist", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "floodDrySpecialist_x_pasture", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "forestPresent", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "forestPresent_x_pasture", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "forestSpecialist", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "forestSpecialist_x_pasture", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "migratory", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "migratory_x_pasture", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "mountain_barrier", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "mountainBarrier_x_pasture", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "tfSpecialist", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "tfSpecialist_x_pasture", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "valley_barrier", dpar = "occ") +
  prior(std_normal(), class = "b", coef = "valleyBarrier_x_pasture", dpar = "occ")
bprior3 <- 
  prior(normal(0, 0.5), class = "b", coef = "dietCarn", dpar = "occ") +
  prior(normal(0, 0.5), class = "b", coef = "dietCarn_x_pasture", dpar = "occ") +
  prior(normal(0, 0.5), class = "b", coef = "dietFruitNect", dpar = "occ") +
  prior(normal(0, 0.5), class = "b", coef = "dietFruitNect_x_pasture", dpar = "occ") +
  prior(normal(0, 0.5), class = "b", coef = "dietGran", dpar = "occ") +
  prior(normal(0, 0.5), class = "b", coef = "dietGran_x_pasture", dpar = "occ") +
  prior(normal(0, 0.5), class = "b", coef = "dietInvert", dpar = "occ") +
  prior(normal(0, 0.5), class = "b", coef = "dietInvert_x_pasture", dpar = "occ") +
  prior(normal(0, 0.5), class = "b", coef = "mass", dpar = "occ") +
  prior(normal(0, 0.5), class = "b", coef = "mass_x_pasture", dpar = "occ") +  
  prior(normal(0, 2), class = "sd", coef = "Intercept", group = "species", dpar = "occ") +
  prior(normal(0, 2), class = "sd", coef = "Intercept", group = "Family", dpar = "occ") +
  prior(std_normal(), class = "sd", coef = "pasture", group = "species", dpar = "occ") +
  prior(std_normal(), class = "sd", coef = "pasture", group = "Family", dpar = "occ") +
  prior(normal(0, 3), class = "sd", coef = "Intercept", group = "species_cluster", dpar = "occ") +
  prior(normal(0, 3), class = "sd", coef = "Intercept", group = "species_subregion", dpar = "occ") +
  prior(normal(0, 2), class = "sd", coef = "relev", group = "species", dpar = "occ") +
  prior(normal(0, 2), class = "sd", coef = "relev2", group = "species", dpar = "occ")

bprior <- bprior1 + bprior2 + bprior3

fm <- flock(
  f_occ = ~ 0 + Intercept + 
    # Species range
    mo(distance_bin) +
    relev + relev2 + 
    lowland + lowland_x_relev + lowland_x_relev2 +
    # Pasture
    pasture + 
    # Biogeography
    #  Barriers
    mountain_barrier + valley_barrier + 
    #  Elevations
    elevMedian + elevBreadth +
    #  Habitats
    forestPresent + forestSpecialist +
    tfSpecialist + dryForestPresent + 
    floodDrySpecialist + aridPresent +
    #  Migratory
    migratory + 
    #  Interactions
    elevMedian_x_forestPresent + elevMedian_x_forestSpecialist +
    # Functional (mass & diet)
    mass + dietInvert + dietCarn + dietFruitNect + dietGran +
    # Pasture interactions
    mountainBarrier_x_pasture + valleyBarrier_x_pasture + 
    elevMedian_x_pasture + elevBreadth_x_pasture +
    forestPresent_x_pasture + forestSpecialist_x_pasture +
    tfSpecialist_x_pasture + dryForestPresent_x_pasture + 
    floodDrySpecialist_x_pasture + aridPresent_x_pasture +
    migratory_x_pasture + 
    elevMedian_x_forestPresent_x_pasture + 
    elevMedian_x_forestSpecialist_x_pasture +
    mass_x_pasture + 
    dietInvert_x_pasture + dietCarn_x_pasture + 
    dietFruitNect_x_pasture + dietGran_x_pasture +
    (1 + relev + relev2 + pasture |g1| species) +
    (1 + pasture |g2| Family) +
    (1 | species_subregion) +
    (1 | species_cluster),

  f_det = ~ 0 + Intercept + 
    pasture +
    mass + elevMedian + migratory + dietCarn +
    time + time_x_elev + 
    obsSM + obsDE + obsJG +
    (1 + pasture + time |g1| species) +
    (1 + pasture |g2| Family) +
    (1 | id_spObs),
  
  flocker_data = fd, 
  prior = bprior,
  cores = 4,
  chains = 4,
  refresh = 1,
  backend = "cmdstanr",
  stan_model_args=list(stanc_options = list("O1")),
  output_dir = "/home/ec2-user/outputs",
  save_warmup = TRUE
)

saveRDS(fm, "/home/ec2-user/outputs/fm.RDS")


# On a machine with more memory (>100 GB), run the above thru line 
#   195 (bprior <- ...) and then run the below after transfering 
#   over the csv outputs from the above

# transfer with something like
# scp -i  /home/ec2-user/inputs/morphnus-ec2-test.cer  ec2-user@3.140.186.105:/home/ec2-user/outputs/model_43804f8e7996969de02f9357488b6f0f-202305101859-4-7c6866.csv /home/ec2-user/outputs/model_43804f8e7996969de02f9357488b6f0f-202305101859-4-7c6866.csv

fm2 <- flock(
  f_occ = ~ 0 + Intercept + 
    # Species range
    mo(distance_bin) +
    relev + relev2 + 
    lowland + lowland_x_relev + lowland_x_relev2 +
    # Pasture
    pasture + 
    # Biogeography
    #  Barriers
    mountain_barrier + valley_barrier + 
    #  Elevations
    elevMedian + elevBreadth +
    #  Habitats
    forestPresent + forestSpecialist +
    tfSpecialist + dryForestPresent + 
    floodDrySpecialist + aridPresent +
    #  Migratory
    migratory + 
    #  Interactions
    elevMedian_x_forestPresent + elevMedian_x_forestSpecialist +
    # Functional (mass & diet)
    mass + dietInvert + dietCarn + dietFruitNect + dietGran +
    # Pasture interactions
    mountainBarrier_x_pasture + valleyBarrier_x_pasture + 
    elevMedian_x_pasture + elevBreadth_x_pasture +
    forestPresent_x_pasture + forestSpecialist_x_pasture +
    tfSpecialist_x_pasture + dryForestPresent_x_pasture + 
    floodDrySpecialist_x_pasture + aridPresent_x_pasture +
    migratory_x_pasture + 
    elevMedian_x_forestPresent_x_pasture + 
    elevMedian_x_forestSpecialist_x_pasture +
    mass_x_pasture + 
    dietInvert_x_pasture + dietCarn_x_pasture + 
    dietFruitNect_x_pasture + dietGran_x_pasture +
    (1 + relev + relev2 + pasture |g1| species) +
    (1 + pasture |g2| Family) +
    (1 | species_subregion) +
    (1 | species_cluster),
  
  f_det = ~ 0 + Intercept + 
    pasture +
    mass + elevMedian + migratory + dietCarn +
    time + time_x_elev + 
    obsSM + obsDE + obsJG +
    (1 + pasture + time |g1| species) +
    (1 + pasture |g2| Family) +
    (1 | id_spObs),
  
  flocker_data = fd, 
  prior = bprior,
  backend = "cmdstanr",
  empty = TRUE
)

new_fit <- brms:::read_csv_as_stanfit(
  list.files("/home/ec2-user/outputs", full.names = TRUE)
)
fm2$fit <- new_fit
fm3 <- brms::rename_pars(fm2)
saveRDS(fm3, "/Users/JacobSocolar/Desktop/fm3.RDS")
