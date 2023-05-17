# Generate prediction info data objects in varying formats (stars, data.table)
#
# Majority of time spent per species is opening the .grd files. Looping across 
# these, extracting relev and tdist and then saving all species to a single .rds 
# file is much faster for downstream processing (pay overhead once, rather than 
# for each posterior). The dt and xy lookup are mainly used downstream. 
# 
# Generates:
# outputs/prediction_info.rds: stars object with tdist and relev, species as bands
# outputs/prediction_info_dt.rds: data.table object, same info as above
# outputs/xy_info_lookup.rds: lookup table to relate data.table back to spatial
#   object
# 
# Note the rm(list = ls()) calls. Can't run this locally, though mileage may vary

# housekeeping ----
rm(list = ls())
library(dplyr); library(stars); library(data.table)
source("bird_analysis_plotting/occupancy_rasters/helper_functions.R")

# read files 
elev_stars <- read_stars("data/elev_raster/elev_raster/raster_elev_AEA.grd") %>%
    setNames("elevation")

bird_data <- readRDS("data/bird_stan_data6_package.RDS")
birds <- readRDS("data/birds.RDS")
species_names <- unique(birds$species)

relev_sd <- bird_data$means_and_sds$relev_sd
relev_offset <- bird_data$means_and_sds$relev_offset

# loop over species ----
pred_info_list <- vector("list", 1614)
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
     
     pred_info_list[[i]]  <- c(buffer_stars, elev_stars, tdist_stars) %>%
         mutate(tdist = ifelse(!is.na(in_range), tdist, NA), 
                elevation = ifelse(!is.na(in_range), elevation, NA), 
                # further remove elevations outside of buffered elevation
                elevation = ifelse(elevation > (sp_lower - sp_breadth) & 
                                       elevation < (sp_upper + sp_breadth), 
                                   elevation, NA),
                relev = ((elevation - sp_lower)/sp_breadth - relev_offset)/relev_sd) %>%
         select(tdist, relev)
}

# save list obj ----
saveRDS(pred_info_list, "outputs/pred_info_list.rds")

# manage memory
rm(bird_data, birds); gc()

# convert to stars ----
# note: this bit is the memory bottleneck 
pred_info <- do.call(c, c(pred_info_list, list(along="species"))) %>%
    st_set_dimensions("species", values = species_names)

pred_info <- trim_stars(pred_info)

# save
saveRDS(pred_info, "outputs/prediction_info.rds")

# convert to dt object ----
tdist <- pred_info %>% select("tdist")

# extract indexing in the stars object (note: can calculate col and row indexing
# but do this to check funs)
col_index_stars <- tdist %>%
    slice(1, along="species") %>%
    mutate(id_cell = 1:n(), id_x=NA, id_y = NA) %>%
    select(id_cell, id_x, id_y)

xy_dim <- dim(col_index_stars)
col_index_stars[["id_x"]] <- matrix(rep(1:xy_dim[2], each=xy_dim[1]), 
                                    nrow=xy_dim[2], ncol=xy_dim[1])
col_index_stars[["id_y"]] <- matrix(rep(1:xy_dim[1], xy_dim[2]), 
                                    nrow=xy_dim[2], ncol=xy_dim[1])

col_index_sf <- st_as_sf(col_index_stars, as_points = TRUE)
col_index_dt <- as.data.table(col_index_sf)
rm(col_index_stars, col_index_sf)

# convert tdist to dt object
tdist_sf <- st_as_sf(tdist, as_points = T, na.rm=F)
tdist_dt <- as.data.table(tdist_sf)
tdist_dt[, geometry := NULL]
tdist_dt[,id_cell := col_index_dt$id_cell]
tdist_dt <- tdist_dt[rowSums(!is.na(tdist_dt[,-"id_cell"])) > 0, ]
tdist_dt <- melt(tdist_dt, id.vars="id_cell", variable.name = "species", 
                 variable.factor=T, value.name = "tdist")
tdist_dt <- tdist_dt[!is.na(tdist),]

# manage memory
rm(tdist, tdist_sf)

## repeat for relev ----
relev <- pred_info %>% select("relev")
relev_sf <- st_as_sf(relev, as_points = T, na.rm=F)
relev_dt <- as.data.table(relev_sf)
relev_dt[, geometry := NULL]
relev_dt[,id_cell := col_index_dt$id_cell]
relev_dt <- relev_dt[rowSums(!is.na(relev_dt[,-"id_cell"])) > 0, ]
relev_dt <- melt(relev_dt, id.vars="id_cell", variable.name = "species", 
                 variable.factor=T, value.name = "relev")
relev_dt <- relev_dt[!is.na(relev),]

## merge ----
pred_dt <- merge(tdist_dt, relev_dt, all=T, by=c("id_cell", "species"))
pred_dt[complete.cases(pred_dt),]

## save ----
saveRDS(pred_dt, "outputs/prediction_info_dt.rds")
saveRDS(col_index_dt, "outputs/xy_info_lookup.rds")
