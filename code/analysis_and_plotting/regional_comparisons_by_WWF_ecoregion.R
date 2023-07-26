# Code to run regional comparisons by WWF ecoregion & plot. 
#
# TODO: split out plotting code from analysis code. 


# packages
library(data.table); library(dplyr); library(ggplot2); library(stars)

# read files
pred_dt <- readRDS("outputs/prediction_info_dt.rds")
pred_dt[,`:=`(elev = NULL, relev = NULL, distance = NULL)]

if(!file.exists("outputs/cell_WWF_lookup.rds")) {
    xy_lookup <- readRDS("outputs/xy_lookup_stars.rds")
    WWF_ecoregions <- readRDS("outputs/WWF_terrestrial_ecoregions.rds") %>%
        mutate(label_id = 1:n())
    WWF_stars <- st_rasterize(WWF_ecoregions %>% select(label = ECO_NAME), 
                              template = xy_lookup %>% mutate(id_cell = NA)) %>%
        mutate(id_cell = as.numeric(xy_lookup$id_cell)) %>%
        rename(label_id = ID)
    
    WWF_df <- as.data.table(WWF_stars) 
    WWF_df <- WWF_df[!is.na(label_id), .(label_id, id_cell)]
    
    WWF_lookup <- WWF_ecoregions %>%
        as_tibble %>%
        select(-geometry) %>%
        rename(label = ECO_NAME)
    
    cell_lookup <- WWF_df %>%
        left_join(., WWF_lookup) %>%
        as.data.table
    
    saveRDS(cell_lookup, "outputs/cell_WWF_lookup.rds")
} else {
    cell_lookup <- readRDS("outputs/cell_WWF_lookup.rds")
}


# calculate sensitivity ----
catch_list <- vector("list", 100)

for(i in 8:100) {
    forest_i <- paste0("posterior_forest_", i, ".rds")
    pasture_i <- paste0("posterior_pasture_", i, ".rds")
    post_forest <- readRDS(paste0("outputs/predicted_occupancy_dts/", forest_i))
    post_pasture <- readRDS(paste0("outputs/predicted_occupancy_dts/", pasture_i))
    
    # note: in each cell, a species could be on 48 points (16 * 3)
    pred_dt[,`:=`(p_forest = post_forest$N_forest/16,
                  p_pasture = post_pasture$N_pasture/16)]
    
    pred_dt[cell_lookup, label := i.label, on = "id_cell"]
    
    winlose_lookup <- pred_dt[, list(winner = sum(p_pasture) >= sum(p_forest), 
                                     RR = log(sum(p_forest)/sum(p_pasture))), 
                              by=species]
    
    summ_region <- pred_dt[, list(in_region = any(pmax(p_pasture, p_forest) > .1), 
                                n_pt = .N), by=.(species, label)]
    summ_CO <- pred_dt[, list(in_region = any(pmax(p_pasture, p_forest) > .1)), 
                       by=.(species)]
    
    # calculate winners
    win_region <- summ_region[winlose_lookup, on="species"][
        in_region == TRUE, .(RR = mean(winner)), by=label
        ]
    
    win_CO <- summ_CO[winlose_lookup, on="species"][
        in_region == TRUE, mean(winner)
        ]
    
    # calculate sensitivity
    RR_region <- summ_region[winlose_lookup, on="species"][
        in_region == TRUE, .(prop_win = mean(RR)), by=label
        ]
    
    RR_CO <- summ_CO[winlose_lookup, on="species"][
        in_region == TRUE, mean(RR)
        ]
    
    
    # store outputs
    catch_list[[i]] <- list(RR_region = RR_region, 
                            RR_CO = RR_CO, 
                            pwin_region = win_region, 
                            pwin_CO = win_CO)
}

catch_old <- catch_list
catch_nulls <- lapply(catch_list, is.null) %>%
    unlist

catch_list <- catch_list[!catch_nulls]

saveRDS(catch_list2, "outputs/temp_catch_list_WWF.rds", compress = FALSE)



# calculate regional gains ----
summ_region <- pred_dt[, list(in_region = any(pmax(p_pasture, p_forest) > .1), 
                              n_pt = .N), by=.(species, label)]

label_vec <- WWF_ecoregions$ECO_NAME

## set up iterating 1000 times
set.seed(101)
region_seqs <- t(replicate(1000, sample(1:13)))


catch_pwin <- region_seqs
catch_pwin[] <- NA
catch_RR <- catch_pwin

for(i in 1:1000) {
    print(i)
    for(j in 1:13) {
        # calculate winners
        in_region <- summ_region[label %in% WWF_lookup$label[region_seqs[i,1:j]]
        ][in_region == TRUE]
        in_region_unique <- in_region[duplicated(species) == FALSE]
        # mean_alpha <- in_region[, .(.N), by=label][,mean(N)]
        # gamma <- in_region_unique[,.N]
        
        temp <- winlose_lookup[in_region_unique, on="species"][,.(prop_win = mean(winner), 
                                                                  RR = mean(RR))]
        
        catch_pwin[i,j] <- temp$prop_win
        catch_RR[i,j] <- temp$RR
        # catch_beta[j] <- mean_alpha/gamma
        
    }
}
    

# plotting ----
RR_df <- reshape2::melt(catch_RR) %>%
    as_tibble %>%
    select(draw = 1, seq_id = 2, RR = 3) 

RR_CO <- RR_df %>%
    filter(seq_id==13) %>%
    slice(1) %>%
    pull(RR)

RR_df2 <- RR_df %>%
    mutate(RR_rel = (RR - min(RR))/(RR_CO-min(RR)))

plot_2_RR <- RR_df2 %>%
    group_by(draw) %>%
    filter(RR_rel > .9) %>%
    slice(1) %>%
    ggplot(aes(seq_id, y = ..count../sum(..count..))) + 
    geom_histogram(boundary=0.5, binwidth = 1, fill="grey90", col="black") +
    scale_x_continuous(breaks = 1:13, limits = c(1, 13)) +
    scale_y_continuous(breaks = seq(0, .5, .1), expand=c(0,0), limits = c(0, .35)) +
    labs(x = "", y = "") +
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.line = element_line(),
          axis.line.y = element_line(),
          panel.grid = element_blank(), 
          panel.border = element_blank())
          # axis.ticks.length.x = unit(3, "pt"), 
          # axis.ticks = element_line(), 
          # axis.line = element_line())

plot_1_RR <- RR_df2 %>%
    group_by(seq_id) %>%
    summarise(mid = mean(RR_rel), 
              lwr = quantile(RR_rel, .05), 
              upr = quantile(RR_rel, .95)) %>%
    ggplot(aes(seq_id, mid, ymin=lwr, ymax=upr)) + 
    geom_line() +
    geom_ribbon(alpha=.1) +
    scale_y_continuous(breaks=seq(0, 1, .2)) +
    scale_x_continuous(breaks=seq(1, 13, 1), limits = c(1, 13)) +
    labs(x = "Number of joins", y = "Sensitivity (scaled)") +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"), 
          panel.grid = element_blank()) +
    geom_hline(yintercept = 1, lty="longdash")

# egg::ggarrange(plot_2, plot_1, heights=c(.3, 1))


## prop_winners
pwin_df <- reshape2::melt(catch_pwin) %>%
    as_tibble %>%
    select(draw = 1, seq_id = 2, pwin = 3) 

pwin_CO <- pwin_df %>%
    filter(seq_id==13) %>%
    slice(1) %>%
    pull(pwin)

pwin_df2 <- pwin_df %>%
    mutate(pwin_rel = (pwin - max(pwin))/(pwin_CO-max(pwin)))


plot_2_pwin <- pwin_df2 %>%
    group_by(draw) %>%
    filter(pwin_rel > .9) %>%
    slice(1) %>%
    ggplot(aes(seq_id, y=..count../sum(..count..))) + 
    geom_histogram(boundary=-0.5, binwidth = 1, fill="grey90", col="black") +
    scale_x_continuous(breaks = 1:13, limits = c(1, 13)) +
    scale_y_continuous(breaks = seq(0, .5, .1), expand=c(0,0), limits = c(0, .35)) +
    labs(x = "", y = "")  +
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.line = element_line(),
          axis.line.y = element_line(),
          panel.grid = element_blank(), 
          panel.border = element_blank())

plot_1_pwin <- pwin_df2 %>%
    group_by(seq_id) %>%
    summarise(mid = mean(pwin_rel), 
              lwr = quantile(pwin_rel, .05), 
              upr = quantile(pwin_rel, .95)) %>%
    ggplot(aes(seq_id, mid, ymin=lwr, ymax=upr)) + 
    geom_line() +
    geom_ribbon(alpha=.1) +
    scale_y_continuous(breaks=seq(0, 1, .2)) +
    scale_x_continuous(breaks=seq(1, 13, 1), limits = c(1, 13)) +
    labs(x = "Number of joins", y = "Proportion winners (scaled)") +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"), 
          panel.grid = element_blank()) +
    geom_hline(yintercept = 1, lty="longdash")
    

plot_all <- egg::ggarrange(plot_2_RR, plot_2_pwin, plot_1_RR, plot_1_pwin, heights=c(.3, 1))
ggsave("figures/accumulation_of_sensitivity_single_iteration_axes.png", plot_all)


## 
# source("code/analysis_and_plotting/helper_functions.R")
# 
# plot(xy_lookup)
# dim(xy_lookup)
# pred_dt[species == "Aburria_aburri"] %>%
#     ggplot(aes(cell_to_x_pos(id_cell, c(712, 517)),
#                cell_to_y_pos(id_cell, c(712, 517)), fill=as.numeric(distance))) + 
#     geom_tile() +
#     coord_equal()
# 
# pred_dt[species == "Aburria_aburri"] %>%
#     ggplot(aes(cell_to_x_pos(id_cell, c(712, 517)),
#                cell_to_y_pos(id_cell, c(712, 517)), fill=p_forest)) + 
#     geom_tile() +
#     coord_equal()
# 
