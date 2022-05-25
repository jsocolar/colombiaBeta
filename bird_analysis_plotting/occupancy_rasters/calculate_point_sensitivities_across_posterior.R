# calculate point probabilities across Colombia, for each posterior draw
# note: requires less than 20GB RAM (approx.), and takes ~3-7 mins per draw

library(data.table)
source("bird_analysis_plotting/get_posterior/get_posterior_z_v6.R")

# read data
stan_output <- readRDS("data/occupancy_v9_202201080800_summary.rds")
draws <- posterior::as_draws_df(stan_output$draws_200)
birds<- readRDS("data/birds.RDS")
bird_data <- readRDS("data/bird_stan_data6_package.RDS")
species_names <- unique(birds$species)

# calculate z_info
z_info <- data.frame(bird_data$data[8:41])
z_info$point <- birds$point
z_info$species <- birds$species

# manage memory
rm(bird_data, birds, stan_output); gc()

# read in prediction info 
pred_info <- readRDS("outputs/prediction_info_dt.rds")
pred_info <- pred_info[complete.cases(pred_info),]

# function to calculate probs
# note: running via function requires just one join, so 2x speed
make_preds <- function(tdist, relev, logit_psi_0, logit_psi_pasture_offset,
                       b1_relev_sp, b1_relev2_sp, lowland, b1_x_lowland_relev,
                       b1_x_lowland_relev2, b5_distance_to_range_sp) {
    
    sp_forest_logit <-  logit_psi_0 - logit_psi_pasture_offset +
        relev * b1_relev_sp +
        relev^2 * b1_relev2_sp +
        relev * lowland * b1_x_lowland_relev +
        relev^2 * lowland * b1_x_lowland_relev2 +
        tdist * b5_distance_to_range_sp
    
    list(boot::inv.logit(sp_forest_logit),
         boot::inv.logit(sp_forest_logit + 2*logit_psi_pasture_offset))
}

# loop over draw ids, saving point preds and summary each time
draw_ids <- 10:50
for(i in draw_ids) {
    start.time <- Sys.time()
    
    # get prediction components
    pc <- as.data.table(get_prediction_components(draws, i, z_info))
    pc[, species := factor(species_names)]
    
    # calculate probs
    pred_info[pc, on="species", 
              c("p_forest", "p_pasture") := 
                  make_preds(tdist = tdist, relev = relev, 
                             logit_psi_0 = logit_psi_0, 
                             logit_psi_pasture_offset = logit_psi_pasture_offset, 
                             b1_relev_sp = b1_relev_sp, 
                             b1_relev2_sp = b1_relev2_sp, 
                             lowland = lowland, 
                             b1_x_lowland_relev = b1_x_lowland_relev, 
                             b1_x_lowland_relev2 = b1_x_lowland_relev2, 
                             b5_distance_to_range_sp = b5_distance_to_range_sp)]
    
    # save
    saveRDS(pred_info[,.(id_cell, species, p_forest, p_pasture)], 
            paste0("outputs/posterior_preds/draw_", i, ".rds"))
    
    # calculate summary
    p_summ <- pred_info[p_forest > .1 | p_pasture > .1, 
                        list(median_log_RR = median(log(p_forest/p_pasture)), 
                             median_log_OR = median(log((p_forest/(1-p_forest))/
                                                            (p_pasture/(1-p_pasture))))), 
                        by="id_cell"]
    
    # save
    saveRDS(p_summ, paste0("outputs/posterior_rel_diff/draw_", i, ".rds"))
    
    # remove prob cols to prevent column replication in next iteration
    pred_info[,`:=`(p_forest = NULL, p_pasture = NULL)]
    
    # manage memory
    rm(p_summ); gc()
    
    # print loop info
    print(paste0(i, ": "))
    print(Sys.time() - start.time)
}

