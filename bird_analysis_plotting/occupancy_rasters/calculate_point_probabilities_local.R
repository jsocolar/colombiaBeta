# generate predicted occupancy for a single posterior iteration
# Note: ultimately probably won't bother saving this output, but doing for now
# Note: this is a version that only works with a subset so can run locally 
# Can store about 2 species x point rasters locally (20GB RAM), but not 3
#
# Note: To run within 30GB RAM, there are several rm(list = ls()) calls: do not
# run with unsaved stuff in memory!

if(!exists("n_species_sample")) n_species_sample <- 200

start.time <- Sys.time()
library(dplyr); library(stars)
source("bird_analysis_plotting/get_posterior/get_posterior_z_v6.R")

# read data
stan_output <- readRDS("data/occupancy_v9_202201080800_summary.rds")
draws <- posterior::as_draws_df(stan_output$draws_200)
birds<- readRDS("data/birds.RDS")
bird_data <- readRDS("data/bird_stan_data6_package.RDS")
species_full <- unique(readRDS("data/birds.RDS")$species)

# calculate z_info
z_info <- data.frame(bird_data$data[8:41])
z_info$point <- birds$point
z_info$species <- birds$species

# create species subset ----
# to keep mem size down for working with locally
pasture_species <- c(Cattle_tyrant = "Machetornis_rixosa", 
                     Rufous_col_sparrow = "Zonotrichia_capensis", 
                     Tropical_kingbird = "Tyrannus_melancholicus", 
                     Sm_billed_ani = "Crotophaga_ani", 
                     Southern_lapwing = "Vanellus_chilensis", 
                     Black_faced_gquit = "Melanospiza_bicolor", 
                     Yellow_faced_gquit = "Tiaris_olivaceus", 
                     Many_striped_canastero = "Asthenes_flammulata", 
                     Ruddy_br_seedeater = "Sporophila_minuta", 
                     Bananaquit = "Coereba_flaveola",
                     Social_flycatcher = "Myiozetetes_similis", 
                     Carib_grackle = "Quiscalus_lugubris")

forest_species <- c("Nothocercus_bonapartei",
                    "Nothocercus_julius", 
                    "Grallaricula_cucullata",
                    "Grallaricula_ferrugineipectus",
                    "Grallaricula_flavirostris",
                    "Grallaricula_lineifrons",
                    "Grallaricula_nana", 
                    "Odontophorus_strophium", 
                    "Geotrygon_montana", 
                    "Crax_alberti", 
                    "Crax_alector", 
                    "Crax_daubentoni",
                    "Crax_globulosa",
                    "Crax_rubra")
sp_sub <- c(pasture_species, forest_species)  
species_names <- c(sp_sub, 
                   sample(species_full[!(species_full %in% sp_sub)], 
                          (n_species_sample - length(sp_sub))))

# calculate species index
species_index <- match(species_names, species_full)

# read and join relev and tdist, and then manage memory
# note: would be fastest to use pred_info directly, but doing subsetting 
# overflows memory 
tdist <- readRDS("outputs/tdist_stars.rds")[,,,species_names]
relev <- readRDS("outputs/relev_stars.rds")[,,,species_names]
pred_info <- c(tdist, relev)
rm(tdist, relev)

# get prediction components
pc <- get_prediction_components(draws, 1, z_info)[species_index,]

# calculate point-level occupancies
logodds_forest <- vector("list", n_species_sample)
logodds_pasture <- vector("list", n_species_sample)
for(i in 1:n_species_sample){
    print(i)
    sp_pc <- pc[i,]
    out <- pred_info %>%
        slice(i, along="species") %>%
        mutate(relev2 = relev * relev, 
               sp_forest_logit = sp_pc$logit_psi_0 - sp_pc$logit_psi_pasture_offset +
                   relev * sp_pc$b1_relev_sp +
                   relev2 * sp_pc$b1_relev2_sp +
                   relev * sp_pc$lowland * sp_pc$b1_x_lowland_relev +
                   relev2 * sp_pc$lowland * sp_pc$b1_x_lowland_relev2 +
                   tdist * sp_pc$b5_distance_to_range_sp,
               sp_pasture_logit = sp_forest_logit + 2*sp_pc$logit_psi_pasture_offset)
    
    logodds_forest[[i]] <- out %>%
        dplyr::select(sp_forest_logit)
    
    logodds_pasture[[i]] <- out %>% 
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

#saveRDS(forest_probs, "data/forest_probability_raster_iter1_subset.rds")
#saveRDS(pasture_probs, "data/pasture_probability_raster_iter1_subset.rds")

print(Sys.time() - start.time)