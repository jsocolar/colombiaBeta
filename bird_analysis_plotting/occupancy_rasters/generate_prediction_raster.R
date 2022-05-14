# generate prediction stars object
# most of time spent per species is opening the .grd files. Looping across these, 
# extracting relev and tdist and then saving all species to a single rds file
# will be much faster (as pay this overhead once, rather than for each posterior)
#
# Note: To run within 30GB RAM, there are several rm(list = ls()) calls: do not
# run with unsaved stuff in memory!
# Would probably be able to compress saved stars files further by concatenating
# them, but no real advantage.. 

library(dplyr); library(stars)

# clear memory 
rm(list = ls())

# read files (think this aggregate step is redundant here.. remove)
raster_elev_AEA <- raster::raster("data/elev_raster/elev_raster/raster_elev_AEA.grd")
raster_agg <- raster::aggregate(raster_elev_AEA, 10)
raster::values(raster_agg) <- 1:(raster::ncell(raster_agg))

raster_disagg <- raster::disaggregate(raster_agg, 10)
elev_df <- raster::as.data.frame(raster_elev_AEA, xy = T) 
elev_df$cell_id <- 1:nrow(elev_df)
sr_df <- raster::as.data.frame(raster::crop(raster_disagg, raster_elev_AEA), xy = T)
all.equal(elev_df[,c("x","y")], sr_df[,c("x","y")])
elev_df$sr_id <- sr_df$elevation
names(elev_df)[3] <- "elevation"

source("bird_analysis_plotting/get_posterior/get_posterior_z_v6.R")
bird_data <- readRDS("data/bird_stan_data6_package.RDS")
birds <- readRDS("data/birds.RDS")
stan_output <- readRDS("data/occupancy_v9_202201080800_summary.rds")
class(stan_output$draws_200)
draws <- posterior::as_draws_df(stan_output$draws_200)

z_info <- data.frame(bird_data$data[8:41])
z_info$point <- birds$point
z_info$species <- birds$species

relev_sd <- bird_data$means_and_sds$relev_sd
relev_offset <- bird_data$means_and_sds$relev_offset
elev_stars <- st_as_stars(raster_elev_AEA)

tdist_list <- vector("list", 1614)
relev_list <- vector("list", 1614)
for(i in 1:1614){
    print(i)
    sp <- unique(birds$species[bird_data$data$id_sp == i])
    sp2 <- gsub("_", " ", sp)
    
    if(length(sp)!= 1){stop('sp does not have length 1')}
    
    tdist_stars <- read_stars(paste0("data/transformed_distance/transformed_distance/", sp, ".grd")) %>%
        setNames("tdist")
    
    buffer_stars <- read_stars(paste0('data/buffered/buffered/', sp2, '_buffered.grd')) %>%
        setNames("in_range") 
    
    # Get elevational min/max
    sp_lower <- unique(birds$lower[birds$species == sp])
    sp_upper <- unique(birds$upper[birds$species == sp])
    sp_breadth <- sp_upper - sp_lower
    
    out <- c(buffer_stars, elev_stars, tdist_stars) %>%
        mutate(tdist = ifelse(!is.na(in_range), tdist, NA), 
               elevation = ifelse(!is.na(in_range), elevation, NA), 
               # further remove elevations outside of buffered elevation
               elevation = ifelse(elevation > (sp_lower - sp_breadth) & 
                                      elevation < (sp_upper + sp_breadth), 
                                  elevation, NA),
               relev = ((elevation - sp_lower)/sp_breadth - relev_offset)/relev_sd)
    
    tdist_list[[i]] <- out %>% select(tdist)
    relev_list[[i]] <- out %>% select(relev)
}

# manage memory
saveRDS(tdist_list, "outputs/tdist_list.rds")
saveRDS(relev_list, "outputs/relev_list.rds")
rm(list = ls())

# species ids
species_names <- unique(readRDS("data/birds.RDS")$species)

# convert to stars object
# note: this bit is more memory hungry- RAM bumps up to ~30GB mark, hence all the 
# garbage collection. The lapply line is to avoid the error "variable names are 
# limited to 10000 bytes" at the merge stage. Not sure what this is about, but 
# this fixes it.  
tdist_list <- readRDS("outputs/tdist_list.rds")
tdist_stars <- tdist_list %>%
    lapply(., function(x) setNames(x, nm = paste0(sample(0:9, 3, T), collapse=""))) %>%
    do.call(c, .) %>%
    merge(.) %>%
    setNames("tdist") %>%
    st_set_dimensions("attributes", values = species_names) %>%
    st_set_dimensions(names = c("x", "y", "species"))

# save stars object & manage memory
saveRDS(tdist_stars, "outputs/tdist_stars.rds")
rm(tdist_stars, tdist_list)
relev_list <- readRDS("outputs/relev_list.rds")

# convert to stars object
relev_stars <- relev_list %>%
    lapply(., function(x) setNames(x, nm = paste0(sample(0:9, 3, T), collapse=""))) %>%
    do.call(c, .) %>%
    merge(.) %>%
    setNames("relev") %>%
    st_set_dimensions("attributes", values = species_names) %>%
    st_set_dimensions(names = c("x", "y", "species"))

# save stars object & manage memory
saveRDS(relev_stars, "outputs/relev_stars.rds")
rm(relev_list)

# concatenate to single stars object
tdist_stars <- readRDS("outputs/tdist_stars.rds")
saveRDS(c(tdist_stars, relev_stars), "prediction_info.rds")

rm(list=ls())
