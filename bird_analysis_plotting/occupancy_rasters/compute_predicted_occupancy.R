# generate posterior predicted occupancy across Colombia
library(data.table); library(dplyr)

# if(!file.exists())
coefs <- readRDS("outputs/lpo_and_coefs.rds")

# Read in info for spatial parts of prediction
# Spatial parts of the linear predictor: build covariate frame

xy_to_cell <- function(id_x, id_y, xy_dim = c(924, 679)) {
    (id_y - 1) * xy_dim[2] + id_x
}

xy_info <- readRDS("outputs/xy_info_lookup.rds") %>%
    # calculate 20 km subregions
    dplyr::mutate(id_subregion_x = ceiling(id_x/10), 
                  id_subregion_y = ceiling(id_y/10), 
                  id_subregion = xy_to_cell(
                      id_subregion_x, id_subregion_y, 
                      c(max(id_subregion_x), max(id_subregion_y))
                  )
    ) %>%
    dplyr::select(c("id_cell", "id_subregion"))

pred_dt <- readRDS("outputs/prediction_info_dt.rds")
pred_dt <- pred_dt[!is.na(tdist) & !is.na(relev), ]
pred_dt[, distance_bin := (boot::logit(tdist) * 14.9) %>%
            cut(c(-Inf, seq(-60, 140, 20), Inf), include.lowest = T, ordered_result = T, labels = F)]
pred_dt <- xy_info[pred_dt, on = "id_cell"]
pred_dt[, tdist := NULL]
sr_lookup <- unique(pred_dt[,c("species", "id_subregion")])
# cl_lookup now redundant
#cl_lookup <- unique(pred_dt[,c("species", "id_cell")])

# columns to sum N across posteriors (division at end)
pred_dt[,`:=`(N_forest_update = 0,
              N_pasture_update = 0)]

# extract
for(i in seq_len(10)) {
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

    pred_dt[species_terms, `:=`(logit_psi_pasture_partial = lpo_forest + 
                                    mos[distance_bin] + 
                                    relev * relev_term + relev^2 * relev2_term, 
                                logit_psi_forest_partial = lpo_pasture + 
                                    mos[distance_bin] + 
                                    relev * relev_term + relev^2 * relev2_term),
            on = "species"]
    
    pred_dt[sr_lookup, `:=`(logit_psi_pasture_partial = logit_psi_pasture_partial +
                                sr_effects,
                            logit_psi_forest_partial = logit_psi_forest_partial +
                                sr_effects), on = c("species", "id_subregion")]
    
    pred_dt[,`:=`(N_forest = 0,
                  N_pasture = 0)]
    
    for(cluster in seq_len(16)) {
        print(paste0("cluster_", cluster))
        pred_dt[, `:=`(
            N_forest = N_forest + 
                boot::inv.logit(logit_psi_forest_partial + rnorm(.N) * coefs$sd_cluster[i])*3,
            N_pasture = N_pasture + 
                boot::inv.logit(logit_psi_pasture_partial + rnorm(.N) * coefs$sd_cluster[i])*3
            )]   
    }
    
    # save posterior (only save 2 cols for space- change this later for 
    # safety/redundancy)
    print("saving")
    saveRDS(pred_dt[,"N_pasture"], 
            paste0("outputs/predicted_occupancy_dts/posterior_pasture_", i, ".rds"), 
            compress = FALSE)
    
    saveRDS(pred_dt[,"N_forest"], 
            paste0("outputs/predicted_occupancy_dts/posterior_forest_", i, ".rds"), 
            compress = FALSE)
    
    # update sum
    pred_dt[,`:=`(N_forest_update = N_forest_update + N_forest, 
                  N_pasture_update = N_pasture_update + N_pasture)]
    
    
    pred_dt[,`:=`(N_forest = NULL, 
                  N_pasture = NULL)]
    gc()
}

pred_dt[,`:=`(N_forest_update = N_forest_update/10, 
              N_pasture_update = N_pasture_update/10)]

saveRDS(pred_dt[,c("N_forest_update", "N_pasture_update")], 
        "outputs/predicted_occupancy_dts/averaged_posterior_10.rds", 
        compress = FALSE)
