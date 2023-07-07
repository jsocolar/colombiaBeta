library(flocker); library("dplyr"); library("brms")

bird_stan_data9_1_package <- readRDS("outputs/bird_stan_data9_package.RDS")

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

birds <- readRDS("outputs/birds.RDS")
assertthat::assert_that(identical(birds$elev_breadth_scaled, uc$elevBreadth))

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

saveRDS(fd, "outputs/flocker_data.RDS")
saveRDS(uc, "outputs/unit_covariates.RDS")
