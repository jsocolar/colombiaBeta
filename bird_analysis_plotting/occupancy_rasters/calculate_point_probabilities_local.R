# generate predicted occupancy for a single posterior iteration
# Note: ultimately probably won't bother saving this output, but doing for now
# Note: this is a version that only works with a subset so can run locally 
# Can store about 2 species x point rasters locally (20GB RAM), but not 3

start.time <- Sys.time()
library(dplyr); library(stars)

# create subset for working with locally 
if(!exists(n_species_sample)) n_species_sample <- 200

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
sp_sub_full <- c(sp_sub, 
                 sample(species_names[!(species_names %in% sp_sub)], 
                        (n_species_sample - length(sp_sub))))

# calculate species index
species_full <- unique(readRDS("data/birds.RDS")$species)
species_index <- which(species_full %in% sp_sub_full)

# read and join relev and tdist, and then manage memory
relev_stars <- readRDS("outputs/relev_stars.rds")[,,,sp_sub_full]
tdist_stars <- readRDS("outputs/tdist_stars.rds")[,,,sp_sub_full]
pred_info <- c(relev_stars, tdist_stars)
rm(relev_stars, tdist_stars)

# get prediction components
pc <- get_prediction_components(draws, 1, z_info)[species_index,]

logodds_forest <- vector("list", n_species_sample)
logodds_pasture <- vector("list", n_species_sample)
for(i in 1:i){
    sp_pc <- pc[i,]
    pred_info[,,,i] %>%
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

species_names <- pred_info$species$values

forest_probs <- logodds_forest %>%
    do.call(c, .) %>%
    merge(.) %>%
    setNames("p") %>%
    mutate(p = boot::inv.logit(p)) %>%
    st_set_dimensions(3, values = species_names) %>%
    st_set_dimensions(names = c("x", "y", "species"))

pasture_probs <- logodds_pasture %>%
    do.call(c, .) %>%
    merge(.) %>%
    setNames("p") %>%
    mutate(p = boot::inv.logit(p)) %>%
    st_set_dimensions(3, values = species_names) %>%
    st_set_dimensions(names = c("x", "y", "species"))

saveRDS(forest_probs, "data/forest_probability_raster_iter1_subset.rds")
saveRDS(pasture_probs, "data/pasture_probability_raster_iter1_subset.rds")

print(Sys.time() - start.time)