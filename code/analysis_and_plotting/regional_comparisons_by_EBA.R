library(sf); library(stars); library(dplyr)

xy_lookup <- read_stars("data/elev_raster/elev_raster/raster_elev_AEA.grd") %>%
    # setNames("elevation") %>%
    mutate(id_cell = 1:n()) %>%
    select(id_cell)

# get bird survey points
birds <- readRDS("outputs/birds.RDS")
point_locs <- birds %>%
    select(site, point, lon, lat) %>%
    unique %>%
    st_as_sf(., coords=c("lon", "lat"), crs="wgs84")

EBA <- read_sf("data/endemic_bird_areas.shp") %>%
    st_make_valid()

CO <- rnaturalearth::ne_countries(country = "Colombia", scale = "small", returnclass = "sf")

CO_high <- rnaturalearth::ne_countries(country = "Colombia", scale = "large", returnclass = "sf") %>%
    st_crop(., CO)

shps_clipped <- st_intersection(EBA, CO_high)
shps_clipped2 <- st_transform(shps_clipped, st_crs(xy_lookup))

plot(shps_clipped %>% select(EBANAME))


library(ggplot2)
to_drop <- c("Chocó", "Darién", "Caribbean", "Costa Central", "Amazon", "white-sand") %>%
    paste0(., collapse="|")

Amazon_rough <- point_locs %>%
    filter(site %in% c("PU", "PS", "SG")) %>%
    st_bbox() %>%
    st_as_sfc() %>%
    st_transform(crs=st_crs(shps_clipped2)) %>%
    st_buffer(100000) %>%
    st_difference(., st_union(shps_clipped2)) %>%
    st_union() %>%
    st_as_sf() %>%
    mutate(label = "Amazon") %>%
    rename(geometry = x) %>%
    st_intersection(CO_high %>% st_transform(st_crs(shps_clipped2)))

shps_clipped3 <- shps_clipped2 %>%
    filter(!grepl(to_drop, EBANAME)) %>%
    mutate(label = case_when(grepl("Central Andean páramo", EBANAME) ~ "Northern Central Andes", 
                             TRUE ~ EBANAME)) %>%
    bind_rows(., Amazon_rough) %>%
    group_by(label) %>%
    summarise(m = unique(label))
plot(shps_clipped3 %>% select(label))

plot_map <- shps_clipped3 %>%
    ggplot() + 
    geom_sf(data=CO_high, fill="grey90") +
    geom_sf(aes(fill=label), col="black", alpha=.8) +
    geom_sf(data=point_locs, aes()) +
    scale_fill_viridis_d() +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"), 
          panel.grid = element_blank())

ggsave("figures/mapped_points_EBA.png", plot = plot_map)

rr <- st_rasterize(shps_clipped3 %>% select(label), xy_lookup %>% mutate(id_cell = NA)) %>%
    mutate(id_cell = as.numeric(xy_lookup$id_cell)) %>%
    rename(label_id = ID)

dd <- as.data.frame(rr) %>%
    as_tibble

EBA_lookup <- shps_clipped3 %>%
    as_tibble %>%
    select(label) %>%
    arrange(label) %>%
    mutate(label_id = 1:n())

cell_lookup <- dd %>%
    filter(!is.na(label_id)) %>%
    left_join(., EBA_lookup)
    
cell_sub <- cell_lookup %>%
    filter(id_cell %in% pred_dt$id_cell)

saveRDS(cell_sub, "outputs/cell_EBA_lookup.rds")

###
library(data.table)
pred_dt <- readRDS("outputs/prediction_info_dt.rds")
pred_dt <- pred_dt[!is.na(tdist) & !is.na(relev), ]
pred_dt[,`:=`(relev = NULL, tdist = NULL)]

# data tables with expected number of occupied points per 2*2km gridcell 
post_av <- readRDS("outputs/predicted_occupancy_dts/averaged_posterior_10.rds")

# note: in each cell, a species could be on 48 points (16 * 3)
pred_dt[,`:=`(p_forest = post_av$N_forest_update/48,
              p_pasture = post_av$N_pasture_update/48)]

cell_sub <- as.data.table(cell_sub)
pred_sub <- pred_dt[cell_sub, on = "id_cell"]

winlose_lookup <- pred_dt[, list(winner = sum(p_pasture) >= sum(p_forest), 
                                 RR = log(sum(p_forest)/sum(p_pasture))), by=species]

summ_EBA <- pred_sub[, list(in_region = any(pmax(p_pasture, p_forest) > .05), 
                            n_pt = .N), by=.(species, label)]
summ_CO <- pred_sub[, list(in_region = any(pmax(p_pasture, p_forest) > .05)), by=.(species)]

win_EBA <- summ_EBA[winlose_lookup, on="species"][in_region == TRUE, mean(winner), by=label]
win_CO <- summ_CO[winlose_lookup, on="species"][in_region == TRUE, mean(winner)]

RR_EBA <- summ_EBA[winlose_lookup, on="species"][in_region == TRUE, mean(RR), by=label]
RR_CO <- summ_CO[winlose_lookup, on="species"][in_region == TRUE, mean(RR)]


p1 <- ggplot(win_EBA, aes(V1, fill=label)) + 
    geom_dotplot(stackgroups = T, binpositions = "all") + 
    geom_vline(aes(xintercept = win_CO), lty="longdash") +
    labs(x = "prop(winners)", title = "(a) Proportion winners") +
    guides(fill="none") +
    scale_fill_viridis_d()

p2 <- ggplot(RR_EBA, aes(V1, fill=label)) + 
    geom_dotplot(stackgroups = T, binpositions = "all") + 
    geom_vline(aes(xintercept = RR_CO), lty="longdash") +
    # scale_y_continuous(limits = c(0, .25)) +
    labs(title = "(b) forest pasture ratio", 
         x = "mean(log(forest/pasture))") +
    # guides(fill="none") +
    scale_fill_viridis_d()

plot_both <- egg::ggarrange(p1, p2)
ggsave("figures/winners&sensitivity_by_region_thresh05.png", plot = plot_both)

drop_list <- c("Cordillera de Mérida", "Darién lowlands", "Central Andean páramo")

mean(RR_EBA[!(EBANAME %in% drop_list),V1/RR_CO])
mean(win_EBA[!(EBANAME %in% drop_list),V1/win_CO])


####
## run multiple iterations ----
# pred_dt <- readRDS("outputs/prediction_info_dt.rds")
# pred_dt <- pred_dt[!is.na(tdist) & !is.na(relev), ]
# pred_dt[,`:=`(relev = NULL, tdist = NULL)]
cell_sub <- as.data.table(cell_sub)
# data tables with expected number of occupied points per 2*2km gridcell 
post_av <- readRDS("outputs/predicted_occupancy_dts/averaged_posterior_10.rds")
catch_list <- c()

for(i in 1:10) {
    print(i)
    forest_i <- paste0("posterior_forest_", i, ".rds")
    pasture_i <- paste0("posterior_pasture_", i, ".rds")
    post_forest <- readRDS(paste0("outputs/predicted_occupancy_dts/", forest_i))
    post_pasture <- readRDS(paste0("outputs/predicted_occupancy_dts/", pasture_i))
    
    # note: in each cell, a species could be on 48 points (16 * 3)
    pred_dt[,`:=`(p_forest = post_forest$N_forest/48,
              p_pasture = post_pasture$N_pasture/48)]
    
    pred_sub <- pred_dt[cell_sub, on = "id_cell"]
    
    winlose_lookup <- pred_dt[, list(winner = sum(p_pasture) >= sum(p_forest), 
                                     RR = log(sum(p_forest)/sum(p_pasture))), by=species]
    
    summ_EBA <- pred_sub[, list(in_region = any(pmax(p_pasture, p_forest) > .05), 
                                n_pt = .N), by=.(species, label)]
    summ_CO <- pred_sub[, list(in_region = any(pmax(p_pasture, p_forest) > .05)), by=.(species)]
    
    win_EBA <- summ_EBA[winlose_lookup, on="species"][in_region == TRUE, mean(winner), by=label]
    win_CO <- summ_CO[winlose_lookup, on="species"][in_region == TRUE, mean(winner)]
    
    RR_EBA <- summ_EBA[winlose_lookup, on="species"][in_region == TRUE, mean(RR), by=label]
    RR_CO <- summ_CO[winlose_lookup, on="species"][in_region == TRUE, mean(RR)]
    
    catch_list[[i]] <- list(RR_region = RR_EBA %>% rename(RR = V1), 
         RR_CO = RR_CO, 
         pwin_region = win_EBA %>% rename(prop_win = V1), 
         pwin_CO = win_CO)
}



RR_CO <- lapply(catch_list, function(x) x$RR_CO) %>%
    unlist


RR_region <- lapply(catch_list, function(x) x$RR_region) %>%
    bind_rows(., .id = "draw") %>%
    mutate(draw = as.integer(draw)) %>%
    bind_rows(., tibble(label = "Colombia", draw = 1:10, RR = RR_CO[draw])) %>%
    mutate(RR_CO = RR_CO[as.integer(draw)]) %>%
    mutate(label = factor(label) %>% relevel(ref = "Colombia"))

pwin_CO <- lapply(catch_list, function(x) x$pwin_CO) %>%
    unlist


pwin_region <- lapply(catch_list, function(x) x$pwin_region) %>%
    bind_rows(., .id = "draw") %>%
    mutate(draw = as.integer(draw)) %>%
    bind_rows(., tibble(label = "Colombia", draw = 1:10, prop_win = pwin_CO[draw])) %>%
    mutate(pwin_CO = pwin_CO[as.integer(draw)]) %>%
    mutate(label = factor(label) %>% relevel(ref = "Colombia"))


ggplot(RR_region, aes(y = RR/RR_CO, x = label)) + 
    ggbeeswarm::geom_beeswarm() +
    coord_flip()

# plotting winner proportion
plot_winners_diff <- ggplot(pwin_region, aes(y = prop_win/pwin_CO, x = label)) + 
    ggbeeswarm::geom_beeswarm() +
    coord_flip() +
    labs(title = "(a) sensitivity", 
         y = "Relative difference", 
         x = "") +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"), 
          panel.grid.minor = element_blank())

plot_winners_abs <- ggplot(pwin_region, aes(y = prop_win, x = label)) + 
    ggbeeswarm::geom_beeswarm() +
    coord_flip() +
    labs(title = "(a) Proportion winners", 
         y = "Proportion winners", 
         x = "") +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"), 
          panel.grid.minor = element_blank())

plot_winners <- egg::ggarrange(plot_winners_abs, 
               plot_winners_diff + theme(axis.text.y = element_blank(), 
                                         plot.title = element_blank()), ncol=2)


## plot sensitivity
plot_RR_diff <- ggplot(RR_region, aes(y = RR/RR_CO, x = label)) + 
    ggbeeswarm::geom_beeswarm() +
    coord_flip() +
    labs(title = "(b) Sensitivity", 
         y = "Relative difference", 
         x = "") +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"), 
          panel.grid.minor = element_blank())

plot_RR_abs <- ggplot(RR_region, aes(y = RR, x = label)) + 
    ggbeeswarm::geom_beeswarm() +
    coord_flip() +
    labs(title = "(b) Sensitivity", 
         y = "Proportion winners", 
         x = "") +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"), 
          panel.grid.minor = element_blank())

plot_RR <- egg::ggarrange(plot_RR_abs, 
                          plot_RR_diff + theme(axis.text.y = element_blank(), 
                                               plot.title = element_blank()), 
                          plot_winners_abs,
                          plot_winners_diff + theme(axis.text.y = element_blank(), 
                                                    plot.title = element_blank()), 
                          ncol=2)

ggsave("figures/sensitivity&prop_winners_by_region.png", plot_RR)
