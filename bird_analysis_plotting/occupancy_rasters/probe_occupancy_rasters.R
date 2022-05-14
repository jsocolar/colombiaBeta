# probe occupancy rasters

library(ggplot2)
# read files ----
forest_rast_exists <- exists("data/forest_probability_raster_iter1_subset.rds")
pasture_rast_exists <- exists("data/pasture_probability_raster_iter1_subset.rds")

if(!forest_rast_exists | !pasture_rast_exists) {
    source("bird_analysis_plotting/occupancy_rasters/calculate_point_probabilities_local.R")
} else {
    forest_probs <- readRDS("data/forest_probability_raster_iter1_subset.rds")
    pasture_probs <- readRDS("data/pasture_probability_raster_iter1_subset.rds")
}

both <- c(forest_probs, pasture_probs) %>%
    setNames(c("p_forest", "p_pasture")) %>%
    mutate(p_abs_diff = p_forest - p_pasture, 
           p_rel_diff = p_forest/p_pasture)

median_abs_diff <- st_apply(both["p_abs_diff",,,], c("x", "y"), function(x) median(x, na.rm=T))
median_rel_diff <- st_apply(both["p_rel_diff",,,], c("x", "y"), function(x) median(x, na.rm=T))

median_abs_diff_thresh <- st_apply(both["p_abs_diff",,,], c("x", "y"), 
                                   function(x) {x[x<.1] <- NA; median(x, na.rm=T)})
median_rel_diff_thresh <- st_apply(both["p_rel_diff",,,], c("x", "y"), 
                                   function(x) {x[x<.1] <- NA; median(x, na.rm=T)})

all_summ <- c(median_abs_diff, median_rel_diff, 
              median_abs_diff_thresh, median_rel_diff_thresh) %>%
    setNames(c("p_abs_diff", "p_rel_diff", "p_abs_diff_thresh", "p_rel_diff_thresh"))

C_border <- rnaturalearth::ne_countries(country="Colombia", scale=10) %>%
    st_as_sf() %>%
    st_transform(., st_crs(both)) %>%
    st_cast(., "MULTILINESTRING")

# plot
p1 <- ggplot() + 
    geom_stars(data=all_summ[1,,]) +
    theme(legend.position = "bottom", 
          axis.text = element_text(colour="black")) +
    scale_fill_viridis_c(na.value = NA) +
    geom_sf(data=C_border) + xlim(-600000, 650000) + labs(x="", y="")
p2 <- ggplot() + 
    geom_stars(data=all_summ[2,,]) +
    theme(legend.position = "bottom", 
          axis.text = element_text(colour="black")) +
    scale_fill_viridis_c(na.value = NA) +
    geom_sf(data=C_border) + xlim(-600000, 650000) + labs(x="", y="")
p3 <- ggplot() + 
    geom_stars(data=all_summ[3,,]) +
    theme(legend.position = "bottom", 
          axis.text = element_text(colour="black")) +
    scale_fill_viridis_c(na.value = NA) +
    geom_sf(data=C_border) + xlim(-600000, 650000) + labs(x="", y="")
p4 <- ggplot() + 
    geom_stars(data=all_summ[4,,]) +
    theme(legend.position = "bottom", 
          axis.text = element_text(colour="black")) +
    scale_fill_viridis_c(na.value = NA) +
    geom_sf(data=C_border) + xlim(-600000, 650000) + labs(x="", y="")

p_all <- egg::ggarrange(p1, p2, p3, p4)
ggsave("plots_summary_300species.png", p_all)


plot(both["diff",,,1:10])
plot(both["p_forest",,,1:10])
plot(both["p_forest",,,"Vanellus_chilensis"])
plot(both["p_pasture",,,"Vanellus_chilensis"])
plot(both["p_abs_diff",,,"Vanellus_chilensis"])
plot(both["diff",,,"Zonotrichia_capensis"])