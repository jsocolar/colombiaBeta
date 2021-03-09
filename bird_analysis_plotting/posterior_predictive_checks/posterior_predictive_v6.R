# This script implements posterior and mixed predictive checks on the model

##### Summary of results: #####
# sum of Q: passed
# Mackenzie-Bailey:
#   1111 underpredicts
#   1000 good
#   0100 good
#   0010 good
#   0001 underpredicts
#   1100 overpredicts
#   1010 overpredicts
#   1001 overpredicts
#   0110 overpredicts
#   0101 overpredicts
#   0011 overpredicts
#   1110 overpredicts
#   1101 overpredicts
#   1011 overpredicts
#   0111 underpredicts
#   0000 good
# Join count statistics of the data: failed
# gdm effect size estimates


# Phylogenetic signal in the species effects
# Normality in the random effects: passed



##### Functions for simulating posterior data replicates and computing discrepancy measures: #####
source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/bird_analysis_plotting/get_posterior/get_posterior_z_v6.R")
source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/bird_analysis_plotting/posterior_predictive_checks/discrepancy_functions.R")

##### Data import #####
# Run commented code once to load in the cmdstan CSV file, covert to draws_dataframe, and save
# v6 <- read_cmdstan_csv("/Users/jacobsocolar/Downloads/occupancy_v6_threads-202102151131-1-5f9d48.csv")
# draws <- posterior::as_draws_df(v6$post_warmup_draws[1:620,,])
# saveRDS(draws, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/v6_draws/draws.RDS")
# rm(v6)
# gc()

# Read in data
bird_data <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data6_package.RDS")
birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/birds.RDS")
draws <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/v6_draws/draws.RDS")

# create z_info object for computing posterior Z (see get_posterior_z.R)
z_info <- data.frame(bird_data$data[8:41])
z_info$point <- birds$point
z_info$species <- birds$species
z_info$cl_q_real <- cluster_q(z_info, z_info$Q)[z_info$id_spCl]

# additional data manipulation in preparation for computing join-count stats
species_list <- unique(birds$species)
bd2 <- data.frame(species = birds$species, id_sp = z_info$id_sp, id_spCl = z_info$id_spCl,
                  lon = birds$lon, lat = birds$lat, pasture = z_info$pasture, cl_q_real = z_info$cl_q_real)
bd3 <- bd2[!duplicated(bd2$id_spCl),]
bd3_forest <- bd3[bd3$pasture == 0, ]
bd3_pasture <- bd3[bd3$pasture == 1, ]

##### Perform posterior predictive checks #####
n_rep <- 100 # number of posterior iterations to use

# Create containers for the output
qsum_rep_include_post <- qsum_rep_resample_margin <- rep(0,n_rep)
mackenzie_bailey_counts <- matrix(nrow = n_rep, ncol = 16)
dist_threshold <- c(5e3, 1e4, 2e4, 5e4, 1e5)
jcs_rs <- jcs_forest_rs <- jcs_pasture_rs <-  jcs_p <- rep(list(matrix(nrow=n_rep, ncol=length(species_list))), length(dist_threshold))

for(i in 1:n_rep){
  print(i)
  iter = 6*i
  
  # Simulate data replicates
    # including spatial terms
  Z_probs_include <- get_Z_probs(draws, iter, z_info, spatial_effect = "include")
  data_rep_include <- get_data_rep(Z_probs_include)
    # resampling spatial terms
  Z_probs_resample <- get_Z_probs(draws, iter, z_info, spatial_effect = "resample")
  data_rep_resample <- get_data_rep(Z_probs_resample)
  
  # sum of Q
  qsum_rep_include_post[i] <- sum(data_rep_include$post)
  qsum_rep_resample_margin[i] <- sum(data_rep_resample$mixed)
  
  # Mackenzie-Bailey
  mackenzie_bailey_counts[i, ] <- mb_counts(data_rep_include$histories, z_info$nv)
  
  # Join count stats--posterior
  cq_include_post <- cluster_q(z_info, data_rep_include$post) # summarize simulated detection at the cluster level
  DRIP <- data_rep_include$post[!duplicated(bd2$id_spCl)] # DRIP = Data Rep Include Post
  DRIP_forest <- DRIP[bd3$pasture == 0]
  DRIP_pasture <- DRIP[bd3$pasture == 1]
  for(j in 1:length(species_list)){
    species <- species_list[j]
    # combined
    bd_sp <- bd3[bd3$species == species, ]
    DRIP_sp <- DRIP[bd3$species == species]
    if(sum(bd_sp$cl_q_real) > 1){
      for(k in 1:length(dist_threshold)){
        th <- dist_threshold[k]
        neighbor_weights <- get_neighbor_weights(bd_sp, function(x){1*(x < dist_threshold[k])})
        jcs_p[[k]][i,j] <- join_counts_stat(DRIP_sp, neighbor_weights)
      }
    }
  }
  
  # Join count stats--mixed
  cq_resample_mixed <- cluster_q(z_info, data_rep_resample$mixed) # summarize simulated detection at the cluster level
  DRRM <- data_rep_resample$mixed[!duplicated(bd2$id_spCl)] # DRRM = Data Rep Resample Mixed
  DRRM_forest <- DRRM[bd3$pasture == 0]
  DRRM_pasture <- DRRM[bd3$pasture == 1]
  for(j in 1:length(species_list)){
    species <- species_list[j]
    # combined
    bd_sp <- bd3[bd3$species == species, ]
    DRRM_sp <- DRRM[bd3$species == species]
    if(sum(bd_sp$cl_q_real) > 1){
      for(k in 1:length(dist_threshold)){
        th <- dist_threshold[k]
        neighbor_weights <- get_neighbor_weights(bd_sp, function(x){1*(x < dist_threshold[k])})
        jcs_rs[[k]][i,j] <- join_counts_stat(DRRM_sp, neighbor_weights)
      }
    }
    # # forest
    # bd_sp <- bd3_forest[bd3_forest$species == species,]
    # DRRM_forest_sp <- DRRM_forest[bd3_forest$species == species]
    # if(sum(bd_sp$cl_q_real) > 1){
    #   for(k in 1:length(dist_threshold)){
    #     th <- dist_threshold[k]
    #     neighbor_weights <- get_neighbor_weights(bd_sp, function(x){1*(x < dist_threshold[k])})
    #     jcs_forest_rs[[k]][i,j] <- join_counts_stat(DRRM_forest_sp, neighbor_weights)
    #   }
    # }
    # # pasture
    # bd_sp <- bd3_pasture[bd3_pasture$species == species,]
    # DRRM_pasture_sp <- DRRM_pasture[bd3_pasture$species == species]
    # if(sum(bd_sp$cl_q_real) > 1){
    #   for(k in 1:length(dist_threshold)){
    #     th <- dist_threshold[k]
    #     neighbor_weights <- get_neighbor_weights(bd_sp, function(x){1*(x < dist_threshold[k])})
    #     jcs_pasture_rs[[k]][i,j] <- join_counts_stat(DRRM_pasture_sp, neighbor_weights)
    #   }
    # }
  }
}

######



hist(qsum_rep_include_post)
hist(qsum_rep_resample_margin)
sum(birds$Q)
mean(qsum_rep_include_post > sum(birds$Q))
mean(qsum_rep_resample_margin > sum(birds$Q))

det_data <- birds[,c("v1","v2","v3","v4")]
real_histories <- apply(det_data, 1, function(x){paste(x, collapse = "")})
real_mb_counts <- mb_counts(real_histories, birds$nv)
dhs <- c("1111",
         "1000", "0100", "0010", "0001",
         "1100", "1010", "1001", "0110", "0101", "0011",
         "1110", "1101", "1011", "0111",
         "0000")
for(i in 1:16){
  hist(mackenzie_bailey_counts[1:16,i], main = paste(dhs[i], real_mb_counts[i]))
}



jcs_real <- jcs_real_forest <- jcs_real_pasture <- rep(list(vector()), length(threshold))
bd2 <- data.frame(species = birds$species, id_sp = z_info$id_sp, id_spCl = z_info$id_spCl,
                  lon = birds$lon, lat = birds$lat, pasture = z_info$pasture, cl_q_real = z_info$cl_q_real)

bd3 <- bd2[!duplicated(bd2$id_spCl),]
bd3_forest <- bd3[bd3$pasture == 0, ]
bd3_pasture <- bd3[bd3$pasture == 1, ]

for(j in 1:length(species_list)){
  species <- species_list[j]
  # combined
  bd_sp <- bd3[bd3$species == species, ]
  if(sum(bd_sp$cl_q_real) > 1){
    for(k in 1:length(threshold)){
      th <- threshold[k]
      neighbor_weights <- get_neighbor_weights(bd_sp, function(x){1*(x < threshold[k])})
      jcs_real[[k]][j] <- join_counts_stat(bd_sp$cl_q_real, neighbor_weights)
    }
  }
  
  # forest
  bd_sp <- bd3_forest[bd3_forest$species == species,]
  if(sum(bd_sp$cl_q_real) > 1){
    for(k in 1:length(threshold)){
      th <- threshold[k]
      neighbor_weights <- get_neighbor_weights(bd_sp, function(x){1*(x < threshold[k])})
      jcs_real_forest[[k]][j] <- join_counts_stat(bd_sp$cl_q_real, neighbor_weights)
    }
  }
  
  # pasture
  bd_sp <- bd3_pasture[bd3_pasture$species == species,]
  if(sum(bd_sp$cl_q_real) > 1){
    for(k in 1:length(threshold)){
      th <- threshold[k]
      neighbor_weights <- get_neighbor_weights(bd_sp, function(x){1*(x < threshold[k])})
      jcs_real_pasture[[k]][j] <- join_counts_stat(bd_sp$cl_q_real, neighbor_weights)
    }
  }
}

for(k in 1:length(threshold)){
  ilist <- ilist_forest <- ilist_pasture <- ng <- ng_forest <- ng_pasture <- vector()
  counter <- counter_f <- counter_p <- 0
  for(i in 1:length(jcs_real[[k]])){
    ref <- jcs_rs[[k]][1:n_rep,i]
    if(!is.na(jcs_real[[k]][i])){
      counter <- counter + 1
      ilist[counter] <- i
      # ng is the fraction of replicates that have lower spatial autocorrelation
      # than the real data.
      ng[counter] <- sum(ref > jcs_real[[k]][i], na.rm = T)/sum(!is.na(ref))
    }
    
    ref <- jcs_forest_rs[[k]][1:n_rep,i]
    if(!is.na(jcs_real_forest[[k]][i])){
      counter_f <- counter_f + 1
      ilist_forest[counter_f] <- i
      ng_forest[counter_f] <- sum(ref > jcs_real_forest[[k]][i], na.rm = T)/sum(!is.na(ref))
    }
    
    ref <- jcs_pasture_rs[[k]][1:n_rep,i]
    if(!is.na(jcs_real_pasture[i])){
      counter_p <- counter_p + 1
      ilist_pasture[counter_p] <- i
      ng_pasture[counter_p] <- sum(ref > jcs_real_pasture[[k]][i], na.rm = T)/sum(!is.na(ref))
    }
  }
  
  hist(ng, main = paste("all", threshold[k]))
  print(c(threshold[k], sum(ng>.95)))
  # hist(ng_forest, main = paste("forest", threshold[k]))
  # hist(ng_pasture, main = paste("pasture", threshold[k]))
  
  # for(i in 1:length(ng)){
  #   if(ng[i] > .95){
  #     print(i)
  #     hist(jcs_rs[[k]][,ilist[i]], main=paste(k, i, species_list[ilist[i]], jcs_real[ilist[i]], ng[i]))
  #   }
  # }
}




# sum of q in forest and pasture by species
forest_q <- pasture_q <- vector()
bd2 <- data.frame(species = z_info$species, id_sp = z_info$id_sp, pasture = birds$pasture,
                  cl_q_real = cq_real[z_info$id_spCl])
bd3 <- bd2[!duplicated(z_info$id_spCl),]
bd4_forest <- bd3[bd3$pasture == 0, ]
bd4_pasture <- bd3[bd3$pasture == 1, ]
for(i in 1:length(species_list)){
  species <- species_list[i]
  
  # forest
  bd_sp <- bd4_forest[bd4_forest$species == species,]
  forest_q[i] <- sum(bd_sp$cl_q_real)
  
  # pasture
  bd_sp <- bd4_pasture[bd4_pasture$species == species,]
  pasture_q[i] <- sum(bd_sp$cl_q_real)
}

plot()

plot(ng_forest ~ forest_q[ilist_forest], ylab = "rep_prop_greater", xlab = "n", main = "fitted clusters")
points(ng_pasture ~ pasture_q[ilist_pasture], col = "red")
abline(h = .95)
abline(h = .05)

plot(ng_forest_rs ~ forest_q[ilist_forest], ylab = "rep_prop_greater", xlab = "n", main = "resampled clusters")
points(ng_pasture_rs ~ pasture_q[ilist_pasture], col = "red")
abline(h = .95)
abline(h = .05)



mean(jcs_real, na.rm = T)
hist(apply(jcs_resample, 1, mean, na.rm = T))
rowSums(jcs, na.rm = T)
hist(rowSums(jcs, na.rm = T))
sum(jcs_real, na.rm = T)

cols_use <- which(!is.na(colSums(jcs_resample)) & !is.na(jcs_real))
mean(jcs_real[cols_use])
hist(apply(jcs_resample[,cols_use], 1, mean, na.rm = T))

cols_use <- which(!is.na(colSums(jcs)) & !is.na(jcs_real))
mean(jcs_real[cols_use])
hist(apply(jcs[,cols_use], 1, mean, na.rm = T))


##### Phylogenetic signal #####
species <- unique(birds$species)
n_pts_by_sp <- rep(0, length(species))
for(i in 1:length(species)){
  n_pts_by_sp[i] <- sum(z_info$Q[z_info$id_sp == i])
}
initial_species_list <- read.csv("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Birds/species_list_creation/initial_species_list.csv")
tree_files <- list.files("/Users/jacobsocolar/Downloads/doi_10.5061_dryad.2v462__v1/JETZ TREES/", full.names = T)[2:1001]
thresholds <- c(10,50,100)
params <- c("b0_sp_raw[", "b1_relev_sp_raw[", "b1_relev2_sp_raw[", "b2_pasture_sp_raw[", "b5_distance_to_range_sp_raw[",
            "d0_sp_raw[", "d1_pasture_sp_raw[", "d2_time_sp_raw[")

effect_sp_data <- data.frame(HBW = gsub("_", " ", species))
effect_sp_data$Pulido <- NA
for(i in 1:nrow(effect_sp_data)){
  if(length(initial_species_list$Pulido[initial_species_list$HBW == effect_sp_data$HBW[i]]) != 1){stop()}
  effect_sp_data$Pulido[i] <- gsub(" ", "_", initial_species_list$Pulido[initial_species_list$HBW == effect_sp_data$HBW[i]])
}
effect_sp_data$effect <- 0

n_iter <- 100

lambda_array <- lambda_null_array <- array(data = 0, dim = c(n_iter, length(params), length(thresholds)))


for(iter in 1:n_iter){
  print(iter)
  pulidoTree <- ape::read.tree(tree_files[iter*floor(1000/n_iter)])
  for(p in 1:length(params)){
    effect_sp_data$effect <- as.numeric(draws[floor(nrow(draws)/n_iter)*iter, paste0(params[p], 1:length(unique(birds$species)), "]")])
    for(th in 1:length(thresholds)){
      effect_sp_data_reduced <- effect_sp_data[n_pts_by_sp >= thresholds[th], ]
      effect_sdr_avg <- effect_sp_data_reduced[!duplicated(effect_sp_data_reduced$Pulido), ]
      for(i in 1:nrow(effect_sdr_avg)){
        effect_sdr_avg$effect[i] <- mean(effect_sp_data_reduced$effect[effect_sp_data_reduced$Pulido == effect_sdr_avg$Pulido[i]])
      }
      tree_pruned <- ape::drop.tip(pulidoTree, pulidoTree$tip.label[!(pulidoTree$tip.label %in% effect_sdr_avg$Pulido)])
      lambda_array[iter, p, th] <- phytools::phylosig(tree_pruned, effect_sdr_avg$effect, method = "lambda")$lambda
      lambda_null_array[iter, p, th] <- phytools::phylosig(tree_pruned, rnorm(nrow(effect_sdr_avg)), method = "lambda")$lambda
    }
  }
}

lambda_null_deviates <- lambda_array - lambda_null_array

summary(lambda_array[,1,1])
summary(lambda_array[,1,2])
summary(lambda_array[,1,3])

summary(lambda_array[,2,1])
summary(lambda_array[,2,2])
summary(lambda_array[,2,3])

summary(lambda_array[,3,1])
summary(lambda_array[,3,2])
summary(lambda_array[,3,3])

summary(lambda_array[,4,1])
summary(lambda_array[,4,2])
summary(lambda_array[,4,3])

summary(lambda_array[,5,1])
summary(lambda_array[,5,2])
summary(lambda_array[,5,3])

summary(lambda_array[,6,1])
summary(lambda_array[,6,2])
summary(lambda_array[,6,3])

summary(lambda_array[,7,1])
summary(lambda_array[,7,2])
summary(lambda_array[,7,3])
quantile(lambda_array[,7,2], .05)

summary(lambda_array[,8,1])
summary(lambda_array[,8,2])
summary(lambda_array[,8,3])


quantile(lambda_null_deviates[,1,1], c(.05, .95))
quantile(lambda_null_deviates[,1,2], c(.05, .95))
quantile(lambda_null_deviates[,1,3], c(.05, .95))

quantile(lambda_null_deviates[,2,1], c(.05, .95))
quantile(lambda_null_deviates[,2,2], c(.05, .95))
quantile(lambda_null_deviates[,2,3], c(.05, .95))

quantile(lambda_null_deviates[,3,1], c(.05, .95))
quantile(lambda_null_deviates[,3,2], c(.05, .95))
quantile(lambda_null_deviates[,3,3], c(.05, .95))

quantile(lambda_null_deviates[,4,1], c(.05, .95))
quantile(lambda_null_deviates[,4,2], c(.05, .95))
quantile(lambda_null_deviates[,4,3], c(.05, .95))

quantile(lambda_null_deviates[,5,1], c(.05, .95))
quantile(lambda_null_deviates[,5,2], c(.05, .95))
quantile(lambda_null_deviates[,5,3], c(.05, .95))

quantile(lambda_null_deviates[,6,1], c(.05, .95))
quantile(lambda_null_deviates[,6,2], c(.05, .95))
quantile(lambda_null_deviates[,6,3], c(.05, .95))

quantile(lambda_null_deviates[,7,1], c(.05, .95))
quantile(lambda_null_deviates[,7,2], c(.05, .95))
quantile(lambda_null_deviates[,7,3], c(.05, .95))

quantile(lambda_null_deviates[,8,1], c(.05, .95))
quantile(lambda_null_deviates[,8,2], c(.05, .95))
quantile(lambda_null_deviates[,8,3], c(.05, .95))

# Out of curiosity, how much phylogenetic signal is there in the species specific pasture terms (summing over the effects)?
for(i in 1:100){
  effect_sp_data$effect <- get_psi_components(draws, 6*i, z_info)$logit_psi_pasture_offset
  effect_sp_data_reduced <- effect_sp_data[n_pts_by_sp >= 10, ]
  effect_sdr_avg <- effect_sp_data_reduced[!duplicated(effect_sp_data_reduced$Pulido), ]
  for(i in 1:nrow(effect_sdr_avg)){
    effect_sdr_avg$effect[i] <- mean(effect_sp_data_reduced$effect[effect_sp_data_reduced$Pulido == effect_sdr_avg$Pulido[i]])
  }
  tree_pruned <- ape::drop.tip(pulidoTree, pulidoTree$tip.label[!(pulidoTree$tip.label %in% effect_sdr_avg$Pulido)])
  print(phytools::phylosig(tree_pruned, effect_sdr_avg$effect, method = "lambda")$lambda)
}

##### Normality of random effect terms #####
species <- unique(birds$species)
all.equal(species, birds$species[!duplicated(bird_data$data$id_sp)])
n_pts_by_sp <- rep(0, length(species))
for(i in 1:length(species)){
  n_pts_by_sp[i] <- sum(z_info$Q[z_info$id_sp == i])
}

shapiro_tests <- as.data.frame(matrix(0, nrow = n_rep, ncol = 17))
names(shapiro_tests) <- c("b0sp_min10", "b0sp_min50", "b0sp_min100", "b0_fam",
                          "b2pasture_sp_min10", "b2pasture_sp_min50", "b2pasture_sp_min100", "b2pasture_fam",
                          "d0sp_min10", "d0sp_min50", "d0sp_min100", "d0_fam",
                          "d1pasture_sp_min10", "d1pasture_sp_min50", "d1pasture_sp_min100", "d1pasture_fam",
                          "b0cluster")

for(i in 1:n_rep){
  iter = 6*i
  b0sp <- as.numeric(draws[iter, paste0("b0_sp_raw[", 1:length(unique(birds$species)), "]")])
  shapiro_tests$b0sp_min10[i] <- shapiro.test(b0sp[n_pts_by_sp > 9])$p.value
  shapiro_tests$b0sp_min50[i] <- shapiro.test(b0sp[n_pts_by_sp > 49])$p.value
  shapiro_tests$b0sp_min100[i] <- shapiro.test(b0sp[n_pts_by_sp > 99])$p.value
  
  b0fam <- as.numeric(draws[iter, paste0("b0_fam_raw[", 1:length(unique(birds$Family)), "]")])
  shapiro_tests$b0_fam[i] <- shapiro.test(b0fam)$p.value
  
  
  b2pasture_sp <- as.numeric(draws[iter, paste0("b2_pasture_sp_raw[", 1:length(unique(birds$species)), "]")])
  shapiro_tests$b2pasture_sp_min10[i] <- shapiro.test(b2pasture_sp[n_pts_by_sp > 9])$p.value
  shapiro_tests$b2pasture_sp_min50[i] <- shapiro.test(b2pasture_sp[n_pts_by_sp > 49])$p.value
  shapiro_tests$b2pasture_sp_min100[i] <- shapiro.test(b2pasture_sp[n_pts_by_sp > 99])$p.value
  
  b2pasture_fam <- as.numeric(draws[iter, paste0("b2_pasture_fam_raw[", 1:length(unique(birds$Family)), "]")])
  shapiro_tests$b2pasture_fam[i] <- shapiro.test(b2pasture_fam)$p.value
  
  
  d0sp <- as.numeric(draws[iter, paste0("d0_sp_raw[", 1:length(unique(birds$species)), "]")])
  shapiro_tests$d0sp_min10[i] <- shapiro.test(d0sp[n_pts_by_sp > 9])$p.value
  shapiro_tests$d0sp_min50[i] <- shapiro.test(d0sp[n_pts_by_sp > 49])$p.value
  shapiro_tests$d0sp_min100[i] <- shapiro.test(d0sp[n_pts_by_sp > 99])$p.value
  
  d0fam <- as.numeric(draws[iter, paste0("d0_fam_raw[", 1:length(unique(birds$Family)), "]")])
  shapiro_tests$d0_fam[i] <- shapiro.test(d0fam)$p.value
  
  
  d1pasture_sp <- as.numeric(draws[iter, paste0("d1_pasture_sp_raw[", 1:length(unique(birds$species)), "]")])
  shapiro_tests$d1pasture_sp_min10[i] <- shapiro.test(d1pasture_sp[n_pts_by_sp > 9])$p.value
  shapiro_tests$d1pasture_sp_min50[i] <- shapiro.test(d1pasture_sp[n_pts_by_sp > 49])$p.value
  shapiro_tests$d1pasture_sp_min100[i] <- shapiro.test(d1pasture_sp[n_pts_by_sp > 99])$p.value
  
  d1pasture_fam <- as.numeric(draws[iter, paste0("d1_pasture_fam_raw[", 1:length(unique(birds$Family)), "]")])
  shapiro_tests$d1pasture_fam[i] <- shapiro.test(d1pasture_fam)$p.value
  
  b0cluster <- as.numeric(draws[iter, paste0("b0_spCl_raw[", 1:length(unique(birds$sp_cl)), "]")])
  shapiro_tests$b0cluster[i] <- shapiro.test(sample(b0cluster, 5000))$p.value
}

for(i in 1:ncol(shapiro_tests)){
  hist(shapiro_tests[,i], main = names(shapiro_tests)[i])
  print(c(names(shapiro_tests)[i], sum(shapiro_tests[,i] > .1)))
}
