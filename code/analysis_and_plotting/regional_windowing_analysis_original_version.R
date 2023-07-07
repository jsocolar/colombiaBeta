# Script to subset the Colombia study region approx. equal sized subregions before
# calculating regional vs local biodiversity loss statistics.

# housekeeping ----
library(data.table); library(dplyr); library(ggplot2)

sample_cols <- function(x) {
    len <- length(unique(x))
    terrain.colors(len)[sample(1:len)]
}

#source("bird_analysis_plotting/occupancy_rasters/helper_functions.R")

# data ----
pred_dt <- readRDS("outputs/prediction_info_dt.rds")
pred_dt <- pred_dt[!is.na(tdist) & !is.na(relev), ]

# N_forest <- readRDS("outputs/predicted_occupancy_dts/posterior_forest_1.rds")
# N_pasture <- readRDS("outputs/predicted_occupancy_dts/posterior_pasture_1.rds")

post_av <- readRDS("outputs/predicted_occupancy_dts/averaged_posterior_10.rds")

# note: in each cell, a species could be on 48 points (16 * 3)
pred_dt[,`:=`(p_pasture = post_av$N_pasture_update/48,
              p_forest = post_av$N_forest_update/48)]
pred_dt[,`:=`(relev = NULL, tdist = NULL)]

# sense check
# ggplot(pred_dt[species == "Harpia_harpyja",], 
#        aes(cell_to_x_pos(id_cell), cell_to_y_pos(id_cell), fill=N_pt_forest)) +
#            geom_tile() +
#     coord_equal()

# read in spatial lookup datatable
xy_info <- readRDS("outputs/xy_info_lookup.rds")
# remove cell ids that don't have data (i.e. that we don't include in study region)
xy_info <- xy_info[id_cell %in% unique(pred_dt$id_cell),]

# Calculate splits ----
# split Colombia into N equal-sized regions
setorder(xy_info, id_x, id_y)
xy_info[,id_split0 := factor(1)]
xy_info[,id_split1 := cut(1:.N, 2, include.lowest = T)]

# iterate across 14 splits, flipping the ordering from y then x to x then y 
# (and vice versa) each time, before splitting in half. 
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
# ggplot(xy_info[seq(1, .N, 4),], aes(id_y, -id_x, col=id_split2)) + 
#     geom_tile() +
#     coord_equal() +
#     guides(colour = "none") +
#     scale_colour_manual(values = sample_cols(xy_info$id_split2))
# ggsave("figures/split2.png")
# 
# ggplot(xy_info[seq(1, .N, 4),], aes(id_y, -id_x, col=id_split8)) + 
#     geom_tile() +
#     coord_equal() +
#     guides(colour = "none") +
#     scale_colour_manual(values = sample_cols(xy_info$id_split8))
# ggsave("figures/split8.png")
# 
# ggplot(xy_info[seq(1, .N, 4),], aes(id_y, -id_x, col=id_split12)) + 
#     geom_tile() +
#     coord_equal() +
#     guides(colour = "none") +
#     scale_colour_manual(values = sample_cols(xy_info$id_split12))
# ggsave("figures/split12.png")

# calculate regional summaries ----
threshold <- 0.1

# info for loop
xy2 <- copy(xy_info)

## vector of split ids
split_vec <- names(xy_info)[grepl("id_split", names(xy_info))]
## list to catch outputs for each split 
summcatch <- vector("list", length(split_vec))
names(summcatch) <- split_vec

### iterate over regions ("splits")
for(i in split_vec) {
    print(i)
    
    # rename ith split column to "id_region" 
    setnames(xy2, i, "id_region")
    # merge ith split column of xy2 (i.id_region) into pred_dt
    pred_dt[xy2, id_region := i.id_region, on="id_cell"]
    
    # summarise the expected number of points occupied in each region, by species
    in_region <- pred_dt[, list(
        in_region = max(p_pasture) > threshold | max(p_forest) > threshold
    ), by = c("species", "id_region")]
    
    # merge in_region column into pred_dt
    pred_dt[in_region, in_region := i.in_region, on = c("species", "id_region")]
    
    # calculate point level richness and summary statistics (logratio and lose
    # prop) on species exceeding threshold at each point.
    #
    # note: filtering for species in region before filtering on whether or not 
    # the within-cell probability exceeds the threshold is redundant. If a species 
    # is present locally (p above threshold) it is by definition in the region 
    # also. Have left it in anyway. 
    point_summary <- pred_dt[
        in_region == TRUE & (p_forest > threshold | p_pasture > threshold), 
        list(
            point_richness_forest = sum(p_forest),
            point_median_logratio = median(log(p_forest/p_pasture)),
            point_lose_prop = mean(p_forest > p_pasture)
        ), by = c("id_cell", "id_region")]
    
    # calculate regional summaries
    species_means <- pred_dt[in_region == TRUE, 
                             list(mean_forest = mean(p_forest), 
                                  mean_pasture = mean(p_pasture)), 
                             by = c("id_region", "species")]
    
    # proportion of species losing in a subregion (i.e. regional summary)
    sr_lose_prop <- species_means[, list(
        lose_prop_regional = mean(mean_forest > mean_pasture), 
        median_logratio_regional = median(log(mean_forest/mean_pasture))
        ), by = "id_region"]
    
    # calculate the average cell-level summary
    avg_loseprop <- point_summary[, list(
        lose_prop_local = mean(point_lose_prop), 
        median_logratio_local = mean(point_median_logratio)
        ), by = "id_region"]
    
    mean_point_richness_forest <- point_summary[, list(
        point_richness_forest = mean(point_richness_forest)
    ), by = "id_region"]
    
    sr_beta <- species_means[,.N, by="id_region"][mean_point_richness_forest,on="id_region"]
    sr_beta[,beta := N/point_richness_forest]
    
    # join to ensure id_subregions all match (they should, but playing it safe)
    dt_full <- sr_beta[sr_lose_prop, on = "id_region"]
    dt_full <- dt_full[avg_loseprop, on = "id_region"]

    # rename id_region to correct split index
    setnames(xy2, "id_region", i)
    
    # store output
    summcatch[[i]] <- copy(dt_full)
}

# bind all outputs together
dt_all <- rbindlist(summcatch, idcol = "id")
dt_all[,id := factor(id, levels = split_vec)]

# plotting ----
# ggplot(dt_all, aes(beta, lose_prop_regional - lose_prop_local)) +
#     geom_point(alpha = .2) + 
#     stat_smooth()
# ggsave("figures/beta_figure_post1_thresh10.png")
# 
# ggplot(dt_all, aes(beta, median_logratio_regional - median_logratio_local)) +
#     geom_point(alpha = .2) + 
#     stat_smooth()
# 
# # ggsave("figures/beta_figure_post1_thresh10_logratio.png")
# 
# 
# ggplot(dt_all, aes(beta, sr_lose_prop_regional - sr_lose_prop_local)) +
#     geom_point(alpha = .2) + 
#     facet_wrap(~id) +
#     stat_smooth()
# ggsave("figures/beta_figure_post1_faceted_thresh10.png")

# save ----
saveRDS(dt_all, "outputs/beta_prop_loss_comparisons_post1_thresh10.rds")

