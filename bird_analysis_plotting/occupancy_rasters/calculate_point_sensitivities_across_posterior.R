# generate predicted occupancy for a single posterior iteration
# Note: there are several rm(list = ls()) calls: do not run with unsaved stuff 
# in memory!

start.time <- Sys.time()
library(dplyr); library(stars)
source("bird_analysis_plotting/get_posterior/get_posterior_z_v6.R")
source("bird_analysis_plotting/occupancy_rasters/helper_functions.R")

# read data
stan_output <- readRDS("data/occupancy_v9_202201080800_summary_sharc.rds")
draws <- posterior::as_draws_df(stan_output$draws_200)
birds <- readRDS("data/birds.RDS")
bird_data <- readRDS("data/bird_stan_data6_package.RDS")
species_names <- unique(birds$species)

# calculate z_info
z_info <- data.frame(bird_data$data[8:41])
z_info$point <- birds$point
z_info$species <- birds$species

# read and join relev and tdist, and then manage memory
# note: would be fastest to use pred_info directly, but doing subsetting 
# overflows memory 
pred_info <- readRDS("data/prediction_info.rds")

summ_pt <- vector("list", 200)
logodds_forest <- vector("list", 1614)
logodds_pasture <- vector("list", 1614)

for(i in 1:200) {
    # iterate across draws 
    print(i)
    # get prediction components
    pc <- get_prediction_components(draws, i, z_info)
    
    # calculate point-level occupancies
    for(j in 1:1614){
        sp_pc <- pc[j,]
        out <- pred_info %>%
            slice(j, along="species") %>%
            mutate(relev2 = relev * relev, 
                   sp_forest_logit = sp_pc$logit_psi_0 - sp_pc$logit_psi_pasture_offset +
                       relev * sp_pc$b1_relev_sp +
                       relev2 * sp_pc$b1_relev2_sp +
                       relev * sp_pc$lowland * sp_pc$b1_x_lowland_relev +
                       relev2 * sp_pc$lowland * sp_pc$b1_x_lowland_relev2 +
                       tdist * sp_pc$b5_distance_to_range_sp,
                   sp_pasture_logit = sp_forest_logit + 2*sp_pc$logit_psi_pasture_offset)
        
        logodds_forest[[j]] <- out %>%
            dplyr::select(sp_forest_logit)
        
        logodds_pasture[[j]] <- out %>% 
            dplyr::select(sp_pasture_logit)
    }
    
    forest_probs <- logodds_forest %>%
        lapply(., function(x) setNames(x, nm = paste0(sample(0:9, 3, T), collapse=""))) %>%
        do.call(c, .) %>%
        merge(.) %>%
        setNames("p") %>%
        mutate(p = boot::inv.logit(p)) %>%
        st_set_dimensions(3, values = species_names) %>%
        st_set_dimensions(names = c("x", "y", "species"))
    
    pasture_probs <- logodds_pasture %>%
        lapply(., function(x) setNames(x, nm = paste0(sample(0:9, 3, T), collapse=""))) %>%
        do.call(c, .) %>%
        merge(.) %>%
        setNames("p") %>%
        mutate(p = boot::inv.logit(p)) %>%
        st_set_dimensions(3, values = species_names) %>%
        st_set_dimensions(names = c("x", "y", "species"))
    
    summ_pt[[i]] <- calc_point_summary(forest_points, pasture_points, .1)
}

saveRDS(summ_pt, "outputs/summ_point_list_0.1_threshold.rds")

print(Sys.time() - start.time)