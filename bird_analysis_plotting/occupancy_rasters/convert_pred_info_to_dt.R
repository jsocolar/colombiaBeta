# generate predicted occupancy for a single posterior iteration
# Note: ultimately probably won't bother saving this output, but doing for now
# Note: this is a version that only works with a subset so can run locally 
# Can store about 2 species x point rasters locally (20GB RAM), but not 3
#
# Note: To run within 30GB RAM, there are several rm(list = ls()) calls: do not
# run with unsaved stuff in memory!
rm(list=ls())
library(dplyr); library(stars); library(data.table)
# source("bird_analysis_plotting/get_posterior/get_posterior_z_v6.R")
source("bird_analysis_plotting/occupancy_rasters/helper_functions.R")

# read data
species_names <- unique(readRDS("data/birds.RDS")$species)

pred_info <- readRDS("outputs/prediction_info_trimmed.rds")
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

# repeat for relev
relev <- pred_info %>% select("relev")
relev_sf <- st_as_sf(relev, as_points = T, na.rm=F)
relev_dt <- as.data.table(relev_sf)
relev_dt[, geometry := NULL]
relev_dt[,id_cell := col_index_dt$id_cell]
relev_dt <- relev_dt[rowSums(!is.na(relev_dt[,-"id_cell"])) > 0, ]
relev_dt <- melt(relev_dt, id.vars="id_cell", variable.name = "species", 
                 variable.factor=T, value.name = "relev")
relev_dt <- relev_dt[!is.na(relev),]

pred_dt <- merge(tdist_dt, relev_dt, all=T, by=c("id_cell", "species"))
pred_dt[complete.cases(pred_dt),]

saveRDS(pred_dt, "outputs/prediction_info_dt.rds")
saveRDS(col_index_dt, "outputs/xy_info_lookup.rds")

