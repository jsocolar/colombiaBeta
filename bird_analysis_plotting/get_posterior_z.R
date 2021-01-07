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

get_z_components <- function(draws, iter, z_info){
  b <- z_info[!duplicated(z_info$id_sp),]
  psi_forest <- as.numeric(draws[iter, paste0("b0_sp[", b$id_sp, "]")]) + as.numeric(draws[iter, paste0("b0_fam[", b$id_fam, "]")]) +
      as.numeric(draws[iter, "b1_lowland"])*b$lowland + 
      as.numeric(draws[iter, "b3_mountain_barrier"])*b$mountain_barrier + as.numeric(draws[iter, "b3_valley_barrier"])*b$valley_barrier +
      as.numeric(draws[iter, "b3_eastOnly"])*b$eastOnly + as.numeric(draws[iter, "b3_westOnly"])*b$westOnly + as.numeric(draws[iter, "b3_snsmOnly"])*b$snsmOnly +
      as.numeric(draws[iter, "b3_notEandes"])*b$notEandes + as.numeric(draws[iter, "b3_notWandes"])*b$notWandes + 
      as.numeric(draws[iter, "b3_elevMedian"])*b$elevMedian + as.numeric(draws[iter, "b3_elevBreadth"])*b$elevBreadth +
      as.numeric(draws[iter, "b3_forestPresent"])*b$forestPresent + as.numeric(draws[iter, "b3_forestSpecialist"])*b$forestSpecialist +
      as.numeric(draws[iter, "b3_tfSpecialist"])*b$tfSpecialist + as.numeric(draws[iter, "b3_dryForestPresent"])*b$dryForestPresent +
      as.numeric(draws[iter, "b3_floodDrySpecialist"])*b$floodDrySpecialist + as.numeric(draws[iter, "b3_floodSpecialist"])*b$floodSpecialist +
      as.numeric(draws[iter, "b3_aridPresent"])*b$aridPresent + as.numeric(draws[iter, "b3_migratory"])*b$migratory + 
      as.numeric(draws[iter, "b3_mass"])*b$mass + as.numeric(draws[iter, "b3_dietInvert"])*b$dietInvert + as.numeric(draws[iter, "b3_dietCarn"])*b$dietCarn +
      as.numeric(draws[iter, "b3_dietFruitNect"])*b$dietFruitNect + as.numeric(draws[iter, "b3_dietGran"])*b$dietGran + 
      as.numeric(draws[iter, "b3_x_elevMedian_forestPresent"])*b$forestPresent*b$elevMedian +
      as.numeric(draws[iter, "b3_x_elevMedian_forestSpecialist"])*b$forestSpecialist*b$elevMedian
  
  psi_pasture_offset <- as.numeric(draws[iter, paste0("b2_pasture_sp[", b$id_sp, "]")]) + as.numeric(draws[iter, paste0("b2_pasture_fam[", b$id_fam, "]")]) +
      as.numeric(draws[iter, "b4_mountain_barrier"])*b$mountain_barrier + as.numeric(draws[iter, "b4_valley_barrier"])*b$valley_barrier +
      as.numeric(draws[iter, "b4_elevMedian"])*b$elevMedian + as.numeric(draws[iter, "b4_elevBreadth"])*b$elevBreadth +
      as.numeric(draws[iter, "b4_forestPresent"])*b$forestPresent + as.numeric(draws[iter, "b4_forestSpecialist"])*b$forestSpecialist +
      as.numeric(draws[iter, "b4_tfSpecialist"])*b$tfSpecialist + as.numeric(draws[iter, "b4_dryForestPresent"])*b$dryForestPresent +
      as.numeric(draws[iter, "b4_floodDrySpecialist"])*b$floodDrySpecialist + as.numeric(draws[iter, "b4_floodSpecialist"])*b$floodSpecialist +
      as.numeric(draws[iter, "b4_aridPresent"])*b$aridPresent + as.numeric(draws[iter, "b4_migratory"])*b$migratory + 
      as.numeric(draws[iter, "b4_mass"])*b$mass + as.numeric(draws[iter, "b4_dietInvert"])*b$dietInvert + as.numeric(draws[iter, "b4_dietCarn"])*b$dietCarn +
      as.numeric(draws[iter, "b4_dietFruitNect"])*b$dietFruitNect + as.numeric(draws[iter, "b4_dietGran"])*b$dietGran + 
      as.numeric(draws[iter, "b4_x_elevMedian_forestPresent"])*b$forestPresent*b$elevMedian +
      as.numeric(draws[iter, "b4_x_elevMedian_forestSpecialist"])*b$forestSpecialist*b$elevMedian

  theta_forest <- as.numeric(draws[iter, paste0("d0_sp[", b$id_sp, "]")]) + as.numeric(draws[iter, paste0("d0_fam[", b$id_fam, "]")]) +
      as.numeric(draws[iter, "d3_mass"])*b$mass + as.numeric(draws[iter, "d3_elevMedian"])*b$elevMedian +
      as.numeric(draws[iter, "d3_migratory"])*b$migratory + as.numeric(draws[iter, "d3_dietCarn"])*b$dietCarn
  
  theta_pasture_offset <- as.numeric(draws[iter, paste0("d1_pasture_sp[", b$id_sp, "]")]) + as.numeric(draws[iter, paste0("d1_pasture_fam[", b$id_fam, "]")])
  
  return(data.frame(psi_forest = psi_forest, psi_pasture_offset = psi_pasture_offset, theta_forest = theta_forest, theta_pasture_offset = theta_pasture_offset))
}

get_Z_probs <- function(draws, iter, z_info, cluster_effect = T){
  zc <- get_z_components(draws, iter, z_info)
  
  if(cluster_effect){
    logit_psi <- zc$psi_forest[z_info$id_sp] + zc$psi_pasture_offset[z_info$id_sp] * z_info$pasture +
      as.numeric(draws[iter, paste0("b0_spCl[", format(z_info$id_spCl, scientific = F, trim = "true", justify = "none"), "]")]) +
      as.numeric(draws[iter, paste0("b1_relev_sp[", z_info$id_sp, "]")]) * z_info$relev +
      as.numeric(draws[iter, paste0("b1_relev2_sp[", z_info$id_sp, "]")]) * z_info$relev2 +
      as.numeric(draws[iter, "b1_x_lowland_relev"]) * z_info$lowland * z_info$relev +
      as.numeric(draws[iter, "b1_x_lowland_relev2"]) * z_info$lowland * z_info$relev2 +
      as.numeric(draws[iter, paste0("b5_distance_to_range_sp[", z_info$id_sp, "]")]) * z_info$distance_to_range
  }else{
    logit_psi <- zc$psi_forest[z_info$id_sp] + zc$psi_pasture_offset[z_info$id_sp] * z_info$pasture +
      as.numeric(draws[iter, paste0("b1_relev_sp[", z_info$id_sp, "]")]) * z_info$relev +
      as.numeric(draws[iter, paste0("b1_relev2_sp[", z_info$id_sp, "]")]) * z_info$relev2 +
      as.numeric(draws[iter, "b1_x_lowland_relev"]) * z_info$lowland * z_info$relev +
      as.numeric(draws[iter, "b1_x_lowland_relev2"]) * z_info$lowland * z_info$relev2 +
      as.numeric(draws[iter, paste0("b5_distance_to_range_sp[", z_info$id_sp, "]")]) * z_info$distance_to_range
  }
  
  psi <- boot::inv.logit(logit_psi)
  
  logit_theta_matrix <- replicate(4, zc$theta_forest[z_info$id_sp] + zc$theta_pasture_offset[z_info$id_sp] * z_info$pasture) +
    sweep(as.matrix(z_info[,c("time.1", "time.2", "time.3", "time.4")]), MARGIN = 1, as.numeric(draws[iter, paste0("d2_time_sp[", z_info$id_sp, "]")]), `*`) +
    as.numeric(draws[iter, "d2_obsSM"]) * as.matrix(z_info[,c("obsSM.1", "obsSM.2", "obsSM.3", "obsSM.4")]) +
    as.numeric(draws[iter, "d2_obsDE"]) * as.matrix(z_info[,c("obsDE.1", "obsDE.2", "obsSM.3", "obsDE.4")]) +
    as.numeric(draws[iter, "d2_obsJG"]) * as.matrix(z_info[,c("obsJG.1", "obsJG.2", "obsJG.3", "obsJG.4")]) +
    sweep(as.matrix(z_info[,c("time.1", "time.2", "time.3", "time.4")]), MARGIN = 1, as.numeric(draws[iter, "d3_x_time_elevMedian"]) * z_info$elevMedian, `*`)
  
  theta_matrix <- boot::inv.logit(logit_theta_matrix)
  zero_det_lik <- apply(cbind(z_info$nv, theta_matrix), 1, FUN = function(x){prod(1 - x[2:(1+x[1])])})
  
  z1_lik <- psi*zero_det_lik
  z0_lik <- 1-psi
  z1_prob <- z1_lik/(z1_lik + z0_lik)
  
  df_out <- data.frame(point = z_info$point, species = z_info$species, Q = z_info$Q)
  df_out$Z_prob <- df_out$Q + (1 - df_out$Q)*z1_prob
  
  return(df_out)
}
