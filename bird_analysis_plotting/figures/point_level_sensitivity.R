# map sensitivity

# housekeeping ----
library(sf); library(ggplot2); library(stars); library(data.table); library(dplyr)

# files 
fnames <- list.files("outputs/posterior_rel_diff/", full.names=T)
xy_lookup <- readRDS("outputs/xy_info_lookup.rds")[,-c("id_x", "id_y")]

# calculate sensitivity summary over posterior
rel_diff_list <- lapply(fnames, readRDS)
rel_diff <- rbindlist(rel_diff_list)
rel_diff <- rel_diff[!is.na(median_log_OR),]
rel_diff <- rel_diff[,list(median_log_RR = mean(median_log_RR), 
                           median_log_OR = mean(median_log_OR), 
                           sd_OR = sd(median_log_OR), 
                           sd_RR = sd(median_log_OR)), by="id_cell"]
rel_diff <- rel_diff[xy_lookup[,.(id_cell, geometry)], on="id_cell"]

rel_diff <- st_as_sf(rel_diff) %>%
    filter(!is.na(median_log_RR))

rd_stars <- st_rasterize(rel_diff, dy=2010, dx=2010)
   
# CO outline for plotting
CO <- rnaturalearth::ne_countries(country = "Colombia", returnclass = "sf", scale = 10) %>%
    st_transform(st_crs(rd_stars))

test <- st_as_sfc(st_bbox(CO)) %>%
    st_buffer(., 50000) %>%
    st_difference(., CO)

# value range for plotting
med_RR_range <- max(abs(range(rel_diff$median_log_RR, na.rm=T)))
med_OR_range <- max(abs(range(rel_diff$median_log_OR, na.rm=T)))

p1 <- ggplot() +
    geom_sf(data = CO, fill="grey70", col=NA) +
    geom_stars(data=rd_stars["median_log_RR",,], downsample = 1, na.rm=T) +
    geom_sf(data = test, fill="grey97", col=NA) +
    geom_sf(data = CO, fill=NA, col="black") +
    theme(axis.text = element_text(colour="black"),
          axis.ticks.length = unit(.2, "cm"),
          panel.border = element_rect(colour="black", fill=NA)) +
    labs(fill = "Median log-ratio", x="", y="") +
    scale_x_continuous(expand=c(0,0), breaks = seq(-100, 100, 5), 
                       limits = c(-660787.5 - 5e4, 675862.5 + 5e4)) +
    scale_y_continuous(expand=c(0,0), limits = c(-924144.7 - 5e4, 931085.3 + 5e4)) +
    scale_fill_gradientn(colours = c(pals::brewer.brbg(201)[1:101],
                                     pals::ocean.curl(n = 201)[100:1]),
                         na.value = NA, 
                         limits = c(-med_RR_range, med_RR_range)) +
    guides(fill = guide_colourbar(barheight = unit( 5 , "cm" ),
                                  ticks.colour = "black",
                                  ticks.linewidth = 2, 
                                  frame.colour = "black", 
                                  frame.linewidth = 1.5, 
                                  label.hjust = 0)) +
    labs(title = "(a) Median log-ratio", fill="")
# p2 <- ggplot() +
#     geom_sf(data = CO, fill="grey70", col=NA) +
#     geom_stars(data=rd_stars["median_log_OR",,], downsample = 1, na.rm=T) +
#     geom_sf(data = test, fill="grey97", col=NA) +
#     geom_sf(data = CO, fill=NA, col="black") +
#     theme(axis.text = element_text(colour="black"),
#           axis.ticks.length = unit(.2, "cm"),
#           panel.border = element_rect(colour="black", fill=NA)) +
#     labs(fill = "Median log-ratio", x="", y="") +
#     scale_x_continuous(expand=c(0,0), breaks = seq(-100, 100, 5), 
#                        limits = c(-660787.5 - 5e4, 675862.5 + 5e4)) +
#     scale_y_continuous(expand=c(0,0), limits = c(-924144.7 - 5e4, 931085.3 + 5e4)) +
#     scale_fill_gradientn(colours = c(pals::brewer.brbg(201)[1:101],
#                                      pals::ocean.curl(n = 201)[100:1]),
#                          na.value = NA, 
#                          limits = c(-med_OR_range, med_OR_range)) +
#     guides(fill = guide_colourbar(barheight = unit( 5 , "cm" ),
#                                   ticks.colour = "black",
#                                   ticks.linewidth = 2, 
#                                   frame.colour = "black", 
#                                   frame.linewidth = 1.5, 
#                                   label.hjust = 0))

p2 <- ggplot() +
    geom_sf(data = CO, fill="grey70", col=NA) +
    geom_stars(data=rd_stars["sd_RR",,], downsample = 1, na.rm=T) +
    geom_sf(data = test, fill="grey97", col=NA) +
    geom_sf(data = CO, fill=NA, col="black") +
    theme(axis.text = element_text(colour="black"),
          axis.ticks.length = unit(.2, "cm"),
          panel.border = element_rect(colour="black", fill=NA)) +
    labs(fill = "Median log-ratio", x="", y="") +
    scale_x_continuous(expand=c(0,0), breaks = seq(-100, 100, 5), 
                       limits = c(-660787.5 - 5e4, 675862.5 + 5e4)) +
    scale_y_continuous(expand=c(0,0), limits = c(-924144.7 - 5e4, 931085.3 + 5e4)) +
    scale_fill_gradientn(colours = c(#pals::brewer.brbg(201)[1:101],
                                     pals::ocean.curl(n = 201)[101:200]),
                         na.value = NA) +
    guides(fill = guide_colourbar(barheight = unit( 5 , "cm" ),
                                  ticks.colour = "black",
                                  ticks.linewidth = 2, 
                                  frame.colour = "black", 
                                  frame.linewidth = 1.5, 
                                  label.hjust = 0)) +
    labs(title = "(b) sd(Median log-ratio)", fill="")


egg::ggarrange(p1, p2, ncol=2)
ggsave("figures/map_temp.png", p1, dpi=400, units="mm", height=116*2, width=107*2)

# plot individual posteriors ----
rd_list <- list()
for(i in 1:12) {
    rd_i <- readRDS(fnames[i])
    rd_i <- rd_i[xy_lookup[,.(id_cell, geometry)], on="id_cell"]
    rd_i[,id_cell := NULL]
    rd_i <- st_as_sf(rd_i) 
    rd_list[[i]] <- st_rasterize(rd_i, dy=2010, dx=2010)
}

rd_all <- rd_list %>%
    do.call(c, .) %>%
    merge %>%
    setNames("median_logratio")

max_abs <- max(abs(c(min(rd_all$median_logratio, na.rm=T), 
                     max(rd_all$median_logratio, na.rm=T))))


plot(rd_all)
