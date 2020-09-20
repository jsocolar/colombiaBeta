library("cmdstanr"); library("dplyr"); library("posterior")


vscale <- function(x){return(as.vector(scale(x)))}


birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_data_trimmed.RDS")
birds$sp_cl <- paste(birds$species, birds$cluster, sep = "__")
birds$elev_median <- rowMeans(cbind(birds$lower, birds$upper))

det_data <- as.matrix(birds[,c("v1", "v2", "v3", "v4")])
det_data[is.na(det_data)] <- -1

time <- matrix((scale(c(birds$hps1, birds$hps2, birds$hps3, birds$hps4))), ncol = 4)
time[is.na(time)] <- 0

obsSM <- matrix(c(birds$obs1 == "SM", birds$obs2 == "SM", birds$obs3 == "SM", birds$obs4 == "SM"), ncol = 4)
obsSM[is.na(obsSM)] <- 0

obsDE <- matrix(c(birds$obs1 == "DE", birds$obs2 == "DE", birds$obs3 == "DE", birds$obs4 == "DE"), ncol = 4)
obsDE[is.na(obsDE)] <- 0

obsJG <- matrix(c(birds$obs1 == "JG", birds$obs2 == "JG", birds$obs3 == "JG", birds$obs4 == "JG"), ncol = 4)
obsJG[is.na(obsJG)] <- 0

# hack for now, to change later
birds$forest_present[is.na(birds$forest_present)] <- 0
birds$forest_specialist[is.na(birds$forest_specialist)] <- 0
birds$tf_specialist[is.na(birds$tf_specialist)] <- 0
birds$dry_forest_present[is.na(birds$dry_forest_present)] <- 0
birds$flood_dry_specialist[is.na(birds$flood_dry_specialist)] <- 0
birds$floodplain_specialist[is.na(birds$floodplain_specialist)] <- 0
birds$arid_present[is.na(birds$arid_present)] <- 0
birds$Family[is.na(birds$Family)] <- "Grallariidae"
birds$migrat_birdlife[is.na(birds$migrat_birdlife)] <- "not a migrant"


stan_data <- list(
# Grainsize for reduce_sum
  grainsize = 1,
  
# Dimensions
#  n_sp_cl = length(unique(birds$sp_cl)),
  n_sp = length(unique(birds$species)),
  n_fam = length(unique(birds$Family)),
  n_tot = nrow(birds),
  n_visit_max = max(birds$nv),

# Detection matrix
  det_data = det_data,

# Q and nv
  Q = birds$Q,
  nv = birds$nv,

# Random effect IDs
#  id_sp_cl = as.numeric(as.factor(birds$sp_cl))
  id_sp = as.numeric(as.factor(birds$species)),
  id_fam = as.numeric(as.factor(birds$Family)),

# Covariates
  relev = vscale(birds$elev_sp_standard),
  relev2 = vscale(birds$elev_sp_standard^2),
  pasture = birds$pasture,
  eastOnly = birds$east_only,
  westOnly = birds$west_only,
  snsmOnly = birds$snsm_only,
  notWandes = birds$wandes_absent,
  notEandes = birds$eandes_absent,
  elevMedian = vscale(birds$elev_median),
  elevBreadth = vscale(birds$elev_breadth),
  forestPresent = birds$forest_present,
  forestSpecialist = birds$forest_specialist,
  tfSpecialist = birds$tf_specialist,
  dryForestPresent = birds$dry_forest_present,
  floodDrySpecialist = birds$flood_dry_specialist,
  floodSpecialist = birds$floodplain_specialist,
  aridPresent = birds$arid_present,
  migratory = as.numeric(birds$migrat_birdlife == "full migrant"),           ## This needs to change for final version
  mass = vscale(birds$BodyMass.Value),
  dietInvert = as.numeric(birds$Diet.5Cat == "Invertebrate"),
  dietCarn = as.numeric(birds$Diet.5Cat == "VertFishScav"),
  dietFruitNect = as.numeric(birds$Diet.5Cat == "FruiNect"),
  dietGran = as.numeric(birds$Diet.5Cat == "PlantSeed"),
  time = time,
  obsSM = obsSM,
  obsJG = obsJG,
  obsDE = obsDE)



# Run ragged model ----
mod_R <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/stan_files/full_colombia_model/occupancyMod_ragged_parallel_v1.stan",
                       cpp_options = list(stan_threads = TRUE))
samps_R <- mod_R$sample(data = stan_data, 
                        chains = 1,
                        threads_per_chain = 4,
                        refresh = 1)


