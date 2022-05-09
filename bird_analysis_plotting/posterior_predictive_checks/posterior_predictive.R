# This script implements posterior and mixed predictive checks on the model

v5 <- cmdstanr::read_cmdstan_csv("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/Stan_outputs/v5_first_run/occupancy_v5_threads-202012282018-1-261afe.csv")

# draws_summary <- parsummarise_draws(v5$post_warmup_draws[1:2000,,], n_cores = 2, n_chunks = 1000)
# print(draws_summary[draws_summary$rhat > 1.05,], n=400)
# print(draws_summary[draws_summary$ess_bulk < 20,], n=400)

draws <- posterior::as_draws_df(v5$post_warmup_draws[1:2000,,])
# Long format: 
bird_data <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/bird_stan_data4_package.RDS")
birds <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/birds.RDS")

z_info <- data.frame(bird_data$data[8:43])
z_info$point <- birds$point
z_info$species <- birds$species
z_info$cl_q_real <- cluster_q(z_info, z_info$Q)[z_info$id_spCl]

source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/bird_analysis_plotting/get_posterior_z.R")
source("/Users/jacobsocolar/Dropbox/Work/Code/colombiaBeta/bird_analysis_plotting/posterior_predictive_checks/discrepancy_functions.R")
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

species_list <- unique(birds$species)
bd2 <- data.frame(species = birds$species, id_sp = z_info$id_sp, id_spCl = z_info$id_spCl,
                  lon = birds$lon, lat = birds$lat, pasture = z_info$pasture, cl_q_real = z_info$cl_q_real)

qsum_rep_include_post <- qsum_rep_resample_margin <- rep(0,100)
mackenzie_bailey_counts <- matrix(nrow = 100, ncol = 16)
threshold <- c(5e3, 1e4, 2e4, 5e4, 1e5)
jcs_rs <- jcs_forest_rs <- jcs_pasture_rs <-  rep(list(matrix(nrow=100, ncol=length(species_list))), length(threshold))

for(i in 1:100){
  print(i)
  iter = 20*i
  Z_probs_include <- get_Z_probs(draws, iter, z_info, cluster_effect = "include")
  data_rep_include <- get_data_rep(Z_probs_include)
  Z_probs_resample <- get_Z_probs(draws, iter, z_info, cluster_effect = "resample")
  data_rep_resample <- get_data_rep(Z_probs_resample)
  
  # sum of Q
  qsum_rep_include_post[i] <- sum(data_rep_include$post)
  qsum_rep_resample_margin[i] <- sum(data_rep_resample$mixed)
  
  # Mackenzie-Bailey
  mackenzie_bailey_counts[i, ] <- mb_counts(data_rep_include$histories, z_info$nv)
  
  # Join count stats
  cq_resample_mixed <- cluster_q(z_info, data_rep_resample$mixed)
  bd3 <- bd2[!duplicated(bd2$id_spCl),]
  DRRM <- data_rep_resample$mixed[!duplicated(bd2$id_spCl)] # DRRM = Data Rep Resample Mixed
  bd3_forest <- bd3[bd3$pasture == 0, ]
  bd3_pasture <- bd3[bd3$pasture == 1, ]
  DRRM_forest <- DRRM[bd3$pasture == 0]
  DRRM_pasture <- DRRM[bd3$pasture == 1]

  for(j in 1:length(species_list)){
    species <- species_list[j]
    # combined
    bd_sp <- bd3[bd3$species == species, ]
    DRRM_sp <- DRRM[bd3$species == species]
    if(sum(bd_sp$cl_q_real) > 1){
      for(k in 1:length(threshold)){
        th <- threshold[k]
        neighbor_weights <- get_neighbor_weights(bd_sp, function(x){1*(x < threshold[k])})
        jcs_rs[[k]][i,j] <- join_counts_stat(DRRM_sp, neighbor_weights)
      }
    }

    # forest
    bd_sp <- bd3_forest[bd3_forest$species == species,]
    DRRM_forest_sp <- DRRM_forest[bd3_forest$species == species]
    if(sum(bd_sp$cl_q_real) > 1){
      for(k in 1:length(threshold)){
        th <- threshold[k]
        neighbor_weights <- get_neighbor_weights(bd_sp, function(x){1*(x < threshold[k])})
        jcs_forest_rs[[k]][i,j] <- join_counts_stat(DRRM_forest_sp, neighbor_weights)
      }
    }

    # pasture
    bd_sp <- bd3_pasture[bd3_pasture$species == species,]
    DRRM_pasture_sp <- DRRM_pasture[bd3_pasture$species == species]
    if(sum(bd_sp$cl_q_real) > 1){
      for(k in 1:length(threshold)){
        th <- threshold[k]
        neighbor_weights <- get_neighbor_weights(bd_sp, function(x){1*(x < threshold[k])})
        jcs_pasture_rs[[k]][i,j] <- join_counts_stat(DRRM_pasture_sp, neighbor_weights)
      }
    }
  }
}

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
  hist(mackenzie_bailey_counts[,i], main = paste(dhs[i], real_mb_counts[i]))
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
    ref <- jcs_rs[[k]][1:100,i]
    if(!is.na(jcs_real[[k]][i])){
      counter <- counter + 1
      ilist[counter] <- i
      ng[counter] <- sum(ref > jcs_real[[k]][i], na.rm = T)/sum(!is.na(ref))
    }
    
    ref <- jcs_forest_rs[[k]][1:100,i]
    if(!is.na(jcs_real_forest[[k]][i])){
      counter_f <- counter_f + 1
      ilist_forest[counter_f] <- i
      ng_forest[counter_f] <- sum(ref > jcs_real_forest[[k]][i], na.rm = T)/sum(!is.na(ref))
    }
    
    ref <- jcs_pasture_rs[[k]][1:100,i]
    if(!is.na(jcs_real_pasture[i])){
      counter_p <- counter_p + 1
      ilist_pasture[counter_p] <- i
      ng_pasture[counter_p] <- sum(ref > jcs_real_pasture[[k]][i], na.rm = T)/sum(!is.na(ref))
    }
  }
  
  hist(ng, main = paste("all", threshold[k]))
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




##### Normality of random effect terms #####
species <- unique(birds$species)
n_pts_by_sp <- rep(0, length(species))
for(i in 1:length(species)){
  n_pts_by_sp[i] <- sum(z_info$Q[z_info$id_sp == i])
}

shapiro_tests <- as.data.frame(matrix(0, nrow = 100, ncol = 17))
names(shapiro_tests) <- c("b0sp_min10", "b0sp_min50", "b0sp_min100", "b0_fam",
                          "b2pasture_sp_min10", "b2pasture_sp_min50", "b2pasture_sp_min100", "b2pasture_fam",
                          "d0sp_min10", "d0sp_min50", "d0sp_min100", "d0_fam",
                          "d1pasture_sp_min10", "d1pasture_sp_min50", "d1pasture_sp_min100", "d1pasture_fam",
                          "b0cluster")

for(i in 1:100){
  iter = 20*i
  b0sp <- as.numeric(draws[iter, paste0("b0_sp[", 1:length(unique(birds$species)), "]")])
  shapiro_tests$b0sp_min10[i] <- shapiro.test(b0sp[n_pts_by_sp > 9])$p.value
  shapiro_tests$b0sp_min50[i] <- shapiro.test(b0sp[n_pts_by_sp > 49])$p.value
  shapiro_tests$b0sp_min100[i] <- shapiro.test(b0sp[n_pts_by_sp > 99])$p.value
  
  b0fam <- as.numeric(draws[iter, paste0("b0_fam[", 1:length(unique(birds$Family)), "]")])
  shapiro_tests$b0_fam[i] <- shapiro.test(b0fam)$p.value
  
  
  b2pasture_sp <- as.numeric(draws[iter, paste0("b2_pasture_sp[", 1:length(unique(birds$species)), "]")])
  shapiro_tests$b2pasture_sp_min10[i] <- shapiro.test(b2pasture_sp[n_pts_by_sp > 9])$p.value
  shapiro_tests$b2pasture_sp_min50[i] <- shapiro.test(b2pasture_sp[n_pts_by_sp > 49])$p.value
  shapiro_tests$b2pasture_sp_min100[i] <- shapiro.test(b2pasture_sp[n_pts_by_sp > 99])$p.value
  
  b2pasture_fam <- as.numeric(draws[iter, paste0("b2_pasture_fam[", 1:length(unique(birds$Family)), "]")])
  shapiro_tests$b2pasture_fam[i] <- shapiro.test(b2pasture_fam)$p.value
  
  
  d0sp <- as.numeric(draws[iter, paste0("d0_sp[", 1:length(unique(birds$species)), "]")])
  shapiro_tests$d0sp_min10[i] <- shapiro.test(d0sp[n_pts_by_sp > 9])$p.value
  shapiro_tests$d0sp_min50[i] <- shapiro.test(d0sp[n_pts_by_sp > 49])$p.value
  shapiro_tests$d0sp_min100[i] <- shapiro.test(d0sp[n_pts_by_sp > 99])$p.value
  
  d0fam <- as.numeric(draws[iter, paste0("d0_fam[", 1:length(unique(birds$Family)), "]")])
  shapiro_tests$d0_fam[i] <- shapiro.test(d0fam)$p.value
  
  
  d1pasture_sp <- as.numeric(draws[iter, paste0("d1_pasture_sp[", 1:length(unique(birds$species)), "]")])
  shapiro_tests$d1pasture_sp_min10[i] <- shapiro.test(d1pasture_sp[n_pts_by_sp > 9])$p.value
  shapiro_tests$d1pasture_sp_min50[i] <- shapiro.test(d1pasture_sp[n_pts_by_sp > 49])$p.value
  shapiro_tests$d1pasture_sp_min100[i] <- shapiro.test(d1pasture_sp[n_pts_by_sp > 99])$p.value
  
  d1pasture_fam <- as.numeric(draws[iter, paste0("d1_pasture_fam[", 1:length(unique(birds$Family)), "]")])
  shapiro_tests$d1pasture_fam[i] <- shapiro.test(d1pasture_fam)$p.value
  
  b0cluster <- as.numeric(draws[iter, paste0("b0_spCl[", 1:length(unique(birds$sp_cl)), "]")])
  shapiro_tests$b0cluster[i] <- shapiro.test(sample(b0cluster, 5000))$p.value
}

for(i in 1:ncol(shapiro_tests)){
  hist(shapiro_tests[,i], main = names(shapiro_tests)[i])
  print(c(names(shapiro_tests)[i], sum(shapiro_tests[,i] > .1)))
}
