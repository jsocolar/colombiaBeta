library(data.table); library(dplyr); library(ggplot2)

#
sample_cols <- function(x) {
    len <- length(unique(x))
    terrain.colors(len)[sample(1:len)]
}

# source("bird_analysis_plotting/occupancy_rasters/helper_functions.R")
pred_dt <- readRDS("outputs/prediction_info_dt.rds")
pred_dt <- pred_dt[!is.na(tdist) & !is.na(relev), ]

N_pasture <- readRDS("outputs/predicted_occupancy_dts/posterior_forest_1.rds")
N_forest <- readRDS("outputs/predicted_occupancy_dts/posterior_pasture_1.rds")

# note: in each cell, a species could be on 48 points (16 * 3)
pred_dt[,`:=`(p_forest = N_forest$N_forest/48,
              p_pasture = N_pasture$N_pasture/48)]
pred_dt[,`:=`(relev = NULL, tdist = NULL)]

# sense check
# source("bird_analysis_plotting/occupancy_rasters/helper_functions.R")
# ggplot(pred_dt[species == "Harpia_harpyja",], 
#        aes(cell_to_x_pos(id_cell), cell_to_y_pos(id_cell), fill=N_pt_forest)) +
#            geom_tile() +
#     coord_equal()

# ggplot(shpfiles, aes(fill=Provincias)) +
#     geom_sf() +
#     geom_sf(data = Colombia, fill=NA, col="black", size=2)
# ggsave("figures/Province_map.png")

# read in spatial lookup datatable
xy_info <- readRDS("outputs/xy_info_lookup.rds")
xy_info <- xy_info[id_cell %in% unique(pred_dt$id_cell),]

# xy_sp <- sf::st_as_sf(xy_info[,c("id_cell", "geometry")]) 
# shpfiles <- shpfiles_crop %>%
#     st_transform(., st_crs(xy_sp))
# xy_sp <- st_join(xy_sp, shpfiles, join = st_intersects)


# Calculate splits ----
## split Colombia into N equal-sized regions
setorder(xy_info, id_x, id_y)
xy_info[,id_split0 := factor(1)]
xy_info[,id_split1 := cut(1:.N, 2, include.lowest = T)]

for(i in 2:14) {
    if(i %% 2 == 0) { # if even
        setorder(xy_info, id_y, id_x)
    } else {
        setorder(xy_info, id_x, id_y) 
    } 
    prev_splits <- names(xy_info)[grepl("id_split", names(xy_info))]
    
    xy_info[,id_split_i := paste0(cut(1:.N, 2, include.lowest = T), .GRP), by=prev_splits]
    setnames(xy_info, "id_split_i", paste0("id_split", i))
}

## plot example regions ----
ggplot(xy_info[seq(1, .N, 4),], aes(id_y, -id_x, col=id_split2)) + 
    geom_tile() +
    coord_equal() +
    guides(colour = "none") +
    scale_colour_manual(values = sample_cols(xy_info$id_split2))
ggsave("figures/split2.png")

ggplot(xy_info[seq(1, .N, 4),], aes(id_y, -id_x, col=id_split8)) + 
    geom_tile() +
    coord_equal() +
    guides(colour = "none") +
    scale_colour_manual(values = sample_cols(xy_info$id_split8))
ggsave("figures/split8.png")

ggplot(xy_info[seq(1, .N, 4),], aes(id_y, -id_x, col=id_split12)) + 
    geom_tile() +
    coord_equal() +
    guides(colour = "none") +
    scale_colour_manual(values = sample_cols(xy_info$id_split12))
ggsave("figures/split12.png")


# calculate regional summaries ----
threshold <- 0.05
### info for loop ----
xy2 <- copy(xy_info)
split_vec <- names(xy_info)[grepl("id_split", names(xy_info))]
summcatch <- vector("list", length(split_vec))
names(summcatch) <- split_vec

### iterate over regions ("splits") ----
for(i in split_vec) {
    print(i)
    # i <- split_vec[1]
    # rename ith split column
    setnames(xy2, i, "id_region")
    # merge ith split column into pred_dt
    pred_dt[xy2, id_region := i.id_region, on="id_cell"]
    
    # calculate number of points in the region (may vary by small amounts due 
    # to total number of cells not dividing exactly) 
    # N_cells <- xy2[, list(N_pts_region = .N * 3 * 16), by="id_region"]

    # summarise the expected number of points occupied in each region, by species
    in_region <- pred_dt[, list(
        in_region = max(p_pasture) > threshold | max(p_forest) > threshold
    ), by = c("species", "id_region")]
    
    # join in_subregion into pred_dt
    pred_dt[in_region, in_region := i.in_region, on = c("species", "id_region")]
    
    # note: filtering for species in region before filtering on whether or not 
    # the within-cell probability exceeds the threshold is redundant. If a species 
    # is present locally (p above threshold) it is by definition in the region 
    # also. Have left it in anyway. 
    point_summary <- pred_dt[
        in_region == TRUE & (p_forest > threshold | p_pasture > threshold), 
        list(
            point_richness_forest = sum(p_forest),
            point_lose_prop = mean(p_forest > p_pasture),
            point_median_logratio = median(log(p_forest/p_pasture))
        ), by = c("id_cell", "id_region")]
    
    # calculate regional summaries
    species_means <- pred_dt[in_region == TRUE, 
                             list(mean_forest = mean(p_forest), 
                                  mean_pasture = mean(p_pasture)), 
                             by = c("id_region", "species")]
    
    # proportion of species losing in a subregion (i.e. regional summary)
    sr_lose_prop <- species_means[, list(sr_lose_prop_regional = 
                                             mean(mean_forest > mean_pasture)), by = "id_region"]
    
    sr_median_logratio <- species_means[, list(sr_median_logratio_regional = 
                                             median(log(mean_forest/mean_pasture))), by = "id_region"]
    
    # calculate the average cell-level summary
    sr_avg_loseprop <- point_summary[, list(
        sr_lose_prop_local = mean(point_lose_prop)), by = "id_region"]
    
    sr_avg_median_logratio <- point_summary[, list(
        sr_lose_prop_local = mean(point_median_logratio)), by = "id_region"]
    
    mean_point_richness_forest <- point_summary[, list(
        point_richness_forest = mean(point_richness_forest)
    ), by = "id_region"]
    
    sr_beta <- species_means[,.N, by="id_region"][mean_point_richness_forest,on="id_region"]
    sr_beta[,beta := N/point_richness_forest]
    
    # join to ensure id_subregions all match (they should, but playing it safe)
    dt_full <- sr_beta[sr_lose_prop, on = "id_region"]
    dt_full <- dt_full[sr_avg_loseprop, on = "id_region"]

    setnames(xy2, "id_region", i)
    summcatch[[i]] <- copy(dt_full)
}

dt_all <- rbindlist(summcatch, idcol = "id")
dt_all[,id := factor(id, levels = split_vec)]

ggplot(dt_all, aes(beta, sr_lose_prop_regional - sr_lose_prop_local)) +
    geom_point(alpha = .2) + 
    stat_smooth()
ggsave("figures/beta_figure_post1.png")

ggplot(dt_all, aes(beta, sr_lose_prop_regional - sr_lose_prop_local)) +
    geom_point(alpha = .2) + 
    facet_wrap(~id) +
    stat_smooth()
ggsave("figures/beta_figure_post1_faceted.png")

saveRDS(dt_all, "outputs/beta_prop_loss_comparisons_post1_thresh05.rds")
