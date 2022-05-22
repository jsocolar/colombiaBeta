# generate predicted occupancy for a single posterior iteration and save
# subsequent read in times of ~ 1 minute per iteration. 

library(dplyr); library(stars)
source("bird_analysis_plotting/get_posterior/get_posterior_z_v6.R")

# helper fun
trim_stars <- function(x) {
    n_attr <- length(x)
    dim1_list <- vector("list", length(n_attr))
    dim2_list <- vector("list", length(n_attr))
    for(i in 1:n_attr) {
        dim1 = apply(x[[i]], 1, function(x) !all(is.na(x)))
        dim1 = which(dim1)
        dim1 = dim1[1]:dim1[length(dim1)]
        dim2 = apply(x[[i]], 2, function(x) !all(is.na(x)))
        dim2 = which(dim2)
        dim2 = dim2[1]:dim2[length(dim2)]
        dim1_list[[i]] <- dim1
        dim2_list[[i]] <- dim2
    }
    d1u <- unlist(dim1_list)
    d2u <- unlist(dim2_list)
    d1u <- sort(unique(d1u))
    d2u <- sort(unique(d2u))
    
    x = x[, d1u, d2u]
    x = st_normalize(x)
    return(x)
}

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

rm(birds, bird_data, stan_output)

# read pred_info
if(!file.exists("outputs/prediction_info_trimmed.rds")) {
    pred_info <- readRDS("data/prediction_info.rds")  
    pred_info <- trim_stars(pred_info)
    saveRDS(pred_info, "outputs/prediction_info_trimmed.rds")
} else {
    pred_info <- readRDS("outputs/prediction_info_trimmed.rds")
}

draw_ids <- 3:50
#mclapply(draw_ids, mc.cores = 4, function(i) {
for (i in draw_ids) {
    start.time <- Sys.time()
    summ_pt <- vector("list", max(draw_ids))
    logodds_forest <- vector("list", 1614)
    logodds_pasture <- vector("list", 1614)
    
    # get prediction components
    pc <- get_prediction_components(draws, i, z_info)
    
    # calculate point-level occupancies
    for(j in 1:1614){
        sp_pc <- pc[j,]
        out <- pred_info %>%
            slice(j, along="species") %>%
            mutate(sp_forest_logit = sp_pc$logit_psi_0 - sp_pc$logit_psi_pasture_offset +
                       relev * sp_pc$b1_relev_sp +
                       relev^2 * sp_pc$b1_relev2_sp +
                       relev * sp_pc$lowland * sp_pc$b1_x_lowland_relev +
                       relev^2 * sp_pc$lowland * sp_pc$b1_x_lowland_relev2 +
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
        setNames("p_forest") %>%
        mutate(p_forest = boot::inv.logit(p_forest)) %>%
        st_set_dimensions(3, values = species_names) %>%
        st_set_dimensions(names = c("x", "y", "species"))
    
    pasture_probs <- logodds_pasture %>%
        lapply(., function(x) setNames(x, nm = paste0(sample(0:9, 3, T), collapse=""))) %>%
        do.call(c, .) %>%
        merge(.) %>%
        setNames("p_pasture") %>%
        mutate(p_pasture = boot::inv.logit(p_pasture)) %>%
        st_set_dimensions(3, values = species_names) %>%
        st_set_dimensions(names = c("x", "y", "species"))
    
    print(paste0(i, " calc: "))
    print(Sys.time() - start.time)
    start.time <- Sys.time()
    
    write_stars(forest_probs, paste0("outputs/posterior_rasters/forest_occ_draw_", i, ".tif"))
    write_stars(pasture_probs, paste0("outputs/posterior_rasters/pasture_occ_draw_", i, ".tif"))
    
    print(paste0(i, " write: "))
    print(Sys.time() - start.time)
    
    rm(forest_probs, pasture_probs, logodds_forest, logodds_pasture)
    gc()    
}

#    saveRDS(c(forest_probs, pasture_probs), 
#            paste0("outputs/posterior_rasters/species_occupancies_draw_", i, ".rds"))
#})
