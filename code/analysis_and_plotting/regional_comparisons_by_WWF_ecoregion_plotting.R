library(sf); set.seed(111)
WWF_ecoregions <- readRDS("outputs/WWF_terrestrial_ecoregions.rds") %>%
    mutate(ECO_NAME = factor(ECO_NAME, levels = sample(ECO_NAME)))

set.seed(101)
col_vec <- sample(viridis::viridis(nrow(WWF_ecoregions)))
names(col_vec) <- levels(WWF_ecoregions$ECO_NAME)

ggplot(WWF_ecoregions) +
    # geom_sf(data=CO, fill="grey90", col="black") +
    geom_sf(aes(fill=ECO_NAME), col="black", alpha=.5) +
    # geom_sf(data=point_locs, pch=21, fill="white", size=2) +
    # scale_fill_manual(values = cols) +
    scale_fill_manual(values = col_vec) +
    scale_x_continuous(breaks=seq(-80, 60, 5)) +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          legend.position = "bottom") +
    labs(fill = "") +
    guides(fill = guide_legend(ncol=3))



RR_CO <- lapply(catch_list, function(x) x$RR_CO) %>%
    unlist

pwin_CO <- lapply(catch_list, function(x) x$pwin_CO) %>%
    unlist

n_draw <- length(catch_list)

RR_region <- lapply(catch_list, function(x) x$RR_region) %>%
    bind_rows(., .id = "draw") %>%
    as_tibble %>%
    rename(RR = prop_win) %>%
    mutate(draw = as.integer(draw)) %>%
    bind_rows(., tibble(label = "Colombia", draw = 1:n_draw, RR = RR_CO[draw])) %>%
    mutate(RR_CO = RR_CO[as.integer(draw)]) %>%
    mutate(rel_diff = exp(RR)/exp(RR_CO)) %>%
    group_by(label) %>%
    mutate(mean_rel_diff = mean(rel_diff)) %>%
    ungroup %>%
    arrange(mean_rel_diff) %>%
    mutate(label_ordered = factor(label, levels = unique(label)),
           label_colour = factor(label, levels = c("Colombia", levels(WWF_ecoregions$ECO_NAME))),
           label_include = !is.na(label_colour) | label == "Colombia") %>%
    arrange(draw)

RR_region %>%
    filter(label != "Colombia") %>%
    summarise(mean(rel_diff), 
              lwr = quantile(rel_diff, .05), 
              upr = quantile(rel_diff, .95))

pwin_region <- lapply(catch_list, function(x) x$pwin_region) %>%
    bind_rows(., .id = "draw") %>%
    as_tibble %>%
    rename(prop_win = RR) %>%
    mutate(draw = as.integer(draw)) %>%
    bind_rows(., tibble(label = "Colombia", draw = 1:n_draw, prop_win = pwin_CO[draw])) %>%
    mutate(prop_win_CO = pwin_CO[as.integer(draw)]) %>%
    mutate(rel_diff = prop_win/prop_win_CO) %>%
    group_by(label) %>%
    mutate(mean_rel_diff = mean(rel_diff)) %>%
    ungroup %>%
    # print(n=20)
    arrange(desc(mean_rel_diff)) %>%
    mutate(label_ordered = factor(label, levels = unique(label)),
           label_colour = factor(label, levels = c("Colombia", levels(WWF_ecoregions$ECO_NAME))),
           label_include = !is.na(label_colour) | label == "Colombia")

levels(RR_region$label_colour)
cols_CO_insert <- c(Colombia = "#FFFFFF", col_vec)
cols2_CO_insert <- c("white", rep("black", 13))

max_rel_RR <- max(RR_region$rel_diff)
plot_RR_rel_diff <- RR_region %>%
    filter(label_include) %>%
    ggplot() +
    ggridges::geom_density_ridges(aes(y = label_ordered, RR, fill=label_colour, 
                                      height = ..density..), stat = "density",
                                  fill=NA, col=NA) +
    ggridges::geom_density_ridges(aes(y = label_ordered, rel_diff,
                                      height = ..density.., 
                                      fill=label_colour, 
                                      col = label_colour),
                                  stat = "density",
                                  alpha = .5) +
    # ggbeeswarm::geom_beeswarm(aes(y = label_ordered, rel_diff))
    scale_fill_manual(values = cols_CO_insert) +
    scale_colour_manual(values = cols2_CO_insert) +
    scale_x_continuous(breaks = seq(0, 3, .2), limits=c(0.2, 8)) +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"),
          axis.text.y = element_blank(),
          panel.grid = element_blank()) +
    guides(fill="none", colour="none") +
    geom_vline(xintercept = 1, lty="longdash") +
    labs(x = "Relative difference",
         y = "") +
    coord_cartesian(xlim=c(0.2, max_rel_RR))

plot_RR_abs <- RR_region %>%
    filter(label_include) %>%
    ggplot(aes(y = label_ordered, RR, fill=label)) +
    ggridges::geom_density_ridges(aes(y = label_ordered, RR, fill=label_colour, 
                                      height = ..density..), stat = "density",
                                  alpha=.5) +
    scale_fill_manual(values = cols_CO_insert) +
    scale_colour_manual(values = cols2_CO_insert) +
    # ggbeeswarm::geom_beeswarm(aes(y = label_ordered, RR)) +
    scale_x_continuous(breaks = log(2^(0:10)), labels=2^(0:10), limits = c(.4, 2.5)) +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"),
          panel.grid = element_blank()) +
    guides(fill="none", colour="none") +
    labs(title = "(a) Sensitivity",
         x = "Sensitivity",
         y = "") 


## proportion winners ----
min_rel_pwin <- min(pwin_region$rel_diff)
plot_pwin_rel_diff <- pwin_region %>%
    filter(label_include) %>%
    ggplot() +
    ggridges::geom_density_ridges(aes(y = label_ordered, prop_win, fill=label_colour, 
                                      height = ..density..), stat = "density",
                                  fill=NA, col=NA) +
    ggridges::geom_density_ridges(aes(y = label_ordered, rel_diff,
                                      height = ..density.., 
                                      fill=label_colour, 
                                      col = label_colour),
                                  stat = "density",
                                  alpha = .5) +
    # ggbeeswarm::geom_beeswarm(aes(y = label_ordered, rel_diff))
    scale_fill_manual(values = cols_CO_insert) +
    scale_colour_manual(values = cols2_CO_insert) +
    scale_x_continuous(breaks = seq(0, 3, .2), limits=c(0.2, 8)) +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"),
          axis.text.y = element_blank(),
          panel.grid = element_blank()) +
    guides(fill="none", colour="none") +
    geom_vline(xintercept = 1, lty="longdash") +
    labs(x = "Relative difference",
         y = "") +
    coord_cartesian(xlim=c(0.95, 1.9))


plot_pwin_abs <- pwin_region %>%
    filter(label_include) %>%
    ggplot(aes(y = label_ordered, prop_win, fill=label)) +
    ggridges::geom_density_ridges(aes(y = label_ordered, prop_win, fill=label_colour, 
                                      height = ..density..), stat = "density",
                                  alpha=.5) +
    scale_fill_manual(values = cols_CO_insert) +
    scale_colour_manual(values = cols2_CO_insert) +
    # ggbeeswarm::geom_beeswarm(aes(y = label_ordered, RR)) +
    scale_x_continuous(limits = c(.18, .46)) +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"),
          panel.grid = element_blank()) +
    guides(fill="none", colour="none") +
    labs(title = "(b) Proportion winners",
         x = "Proportion winners",
         y = "") 

plot_RR <- egg::ggarrange(plot_RR_abs + theme(plot.title = element_text(size=12, face="bold")),
                          plot_RR_rel_diff,
                          plot_pwin_abs + theme(plot.title = element_text(size=12, face="bold")),
                          plot_pwin_rel_diff,
                          ncol=2)

ggsave("figures/sensitivity&prop_winners_by_WWF_region_revised.png", plot_RR, 
       units = "mm", width = 180*.9, height = 210*.9)


plot_RR <- egg::ggarrange(plot_RR_abs + theme(plot.title = element_text(size=12, face="bold")),
                          plot_RR_rel_diff,
                          plot_pwin_abs + theme(plot.title = element_text(size=12, face="bold"), 
                                                axis.text.y = element_blank()),
                          plot_pwin_rel_diff,
                          ncol=4)

ggsave("figures/sensitivity&prop_winners_by_WWF_region_revised2.png", plot_RR, 
       units = "mm",)

