# generate posterior predicted occupancy across Colombia
# note: requires running compute_coefficients.R first, to generate 
# "lpo_and_coefs.rds" (the component of the linear predictor that doesn't vary 
# spatially)

# housekeeping ----
library(data.table); library(dplyr)
n_post_preds <- 100

xy_to_cell <- function(id_x, id_y, xy_dim = c(712, 517)) {
    (id_y - 1) * xy_dim[2] + id_x
}

# inputs ----
coefs <- readRDS("outputs/lpo_and_coefs_400draws.rds")
n_draws <- ncol(coefs$lpo_forest)

# Read in info for spatial parts of prediction
xy_info <- readRDS("outputs/xy_lookup.rds") %>%
    # calculate 20 km subregions
    dplyr::mutate(id_subregion_x = ceiling(id_x/10), 
                  id_subregion_y = ceiling(id_y/10), 
                  id_subregion = xy_to_cell(
                      id_subregion_x, id_subregion_y, 
                      c(max(id_subregion_x), max(id_subregion_y))
                  )
    ) %>%
    dplyr::select(c("id_cell", "id_subregion"))

# read in dataframe for prediction
pred_dt <- readRDS("outputs/prediction_info_dt.rds")

# calculate distance bins
pred_dt[, distance_bin := cut(distance/1000, 
                              c(-Inf, seq(-60, 140, 20), Inf), 
                              include.lowest = T, 
                              ordered_result = T, 
                              labels = F)]

# merge with xy_info for cell indexing
pred_dt <- xy_info[pred_dt, on = "id_cell"]
sr_lookup <- unique(pred_dt[,c("species", "id_subregion")])

# columns to sum N across posteriors (division at end)
pred_dt[,`:=`(p_forest_update = 0,
              p_pasture_update = 0)]

# calculate predicted occupancy (including spatial random effects)
# iterate across draws
for(i in 1:n_post_preds) {
    print(i)
    species_terms <- with(coefs,
                          data.table(
                              species = species,
                              lpo_pasture = lpo_pasture[,i],
                              lpo_forest = lpo_forest[,i],
                              relev_term = relev_term[,i],
                              relev2_term = relev2_term[,i]))
    
    mos <- coefs$mos[,i]      
    sr_lookup[,sr_effects := rnorm(.N) * coefs$sd_subregion[i]]

    pred_dt[species_terms, `:=`(logit_psi_pasture_partial = lpo_pasture + 
                                    mos[distance_bin] + 
                                    relev * relev_term + relev^2 * relev2_term, 
                                logit_psi_forest_partial = lpo_forest + 
                                    mos[distance_bin] + 
                                    relev * relev_term + relev^2 * relev2_term),
            on = "species"]
    
    pred_dt[sr_lookup, `:=`(logit_psi_pasture_partial = logit_psi_pasture_partial +
                                sr_effects,
                            logit_psi_forest_partial = logit_psi_forest_partial +
                                sr_effects), on = c("species", "id_subregion")]
    
    pred_dt[,`:=`(p_forest = 0,
                  p_pasture = 0)]
    
    # average across cluster effect
    for(cluster in seq_len(16)) {
        print(paste0("cluster_", cluster))
        pred_dt[, `:=`(
            p_forest = p_forest + 
                boot::inv.logit(logit_psi_forest_partial + rnorm(.N) * coefs$sd_cluster[i]),
            p_pasture = p_pasture + 
                boot::inv.logit(logit_psi_pasture_partial + rnorm(.N) * coefs$sd_cluster[i])
            )]   
    }
    pred_dt[,`:=`(p_forest = p_forest/16, 
                  p_pasture = p_pasture/16)]
    
    # save posterior
    print("saving")
    saveRDS(pred_dt[,.(id_cell, species, p_forest, p_pasture)], 
            paste0("outputs/predicted_occupancy_dts/posterior_", i, ".rds"), 
            compress = FALSE)
    
    # update sum
    pred_dt[,`:=`(p_forest_update = p_forest_update + p_forest, 
                  p_pasture_update = p_pasture_update + p_pasture)]
    gc()
}

pred_dt[,`:=`(p_forest_update = p_forest_update/n_post_preds, 
              p_pasture_update = p_pasture_update/n_post_preds)]

saveRDS(pred_dt[,.(id_cell, p_forest = p_forest_update, p_pasture = p_pasture_update)], 
        paste0("outputs/predicted_occupancy_dts/averaged_posterior_", n_post_preds, ".rds"), 
        compress = FALSE)
