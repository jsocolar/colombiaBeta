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
pred_dt[,`:=`(N_pt_pasture = N_forest$N_pasture,
              N_pt_forest = N_pasture$N_forest)]
pred_dt[,`:=`(relev = NULL, tdist = NULL)]

# source("bird_analysis_plotting/occupancy_rasters/helper_functions.R")
# ggplot(pred_dt[species == "Harpia_harpyja",], 
#        aes(cell_to_x_pos(id_cell), cell_to_y_pos(id_cell), fill=N_pt_forest)) +
#            geom_tile() +
#     coord_equal()

ggplot(shpfiles, aes(fill=Provincias)) +
    geom_sf() +
    geom_sf(data = Colombia, fill=NA, col="black", size=2)
ggsave("figures/Province_map.png")

# read in spatial lookup datatable
xy_info <- readRDS("outputs/xy_info_lookup.rds")
xy_info <- xy_info[id_cell %in% unique(pred_dt$id_cell),]

xy_sp <- sf::st_as_sf(xy_info[,c("id_cell", "geometry")]) 
shpfiles <- shpfiles_crop %>%
    st_transform(., st_crs(xy_sp))
xy_sp <- st_join(xy_sp, shpfiles, join = st_intersects)


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

# splitcols <- c(prev_splits, "id_split14")
# lapply(xy_info, class)
# # xy_info[, lapply(.SD, c(prev_splits, "id_split14")
# xy_info[,(splitcols):= lapply(.SD, as.factor), .SDcols = splitcols]
                 
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
    summ_i <- pred_dt[, list(N_pt_pasture = sum(p_pasture),
                           N_pt_forest = sum(p_forest),
                           logratio = log(sum(p_forest)/sum(p_pasture)),
                           in_region_pasture = max(p_pasture) > .01,
                           in_region_forest = max(p_forest) > .01),
                      by=c("species", "id_region")]

    summ2 <- summ_i[in_region_pasture == TRUE | in_region_forest == TRUE,
                    list(#S_pasture = sum(N_pasture > 0.05)/.N, # "species richness"
                         #S_forest = sum(N_forest > 0.05)/.N,
                         median_logratio = median(logratio),
                         # number of species doing worse in pasture / total number of species in region
                         prop_losers = mean(N_pt_forest > N_pt_pasture), 
                         #prop_extreme_losers = sum(N_pt_forest/N_pt_pasture > 100)/.N,
                         n_species = .N), 
                    by = "id_region"]

    setnames(summ2, "id_region", i)
    setnames(xy2, "id_region", i)
    summcatch[[i]] <- copy(summ2)
}

saveRDS(summcatch, "outputs/summcatch_1pct_threshold.rds")

ggplot(summcatch[[12]], aes(median_logratio, prop_losers)) + geom_point()


s13 <- summcatch$id_split13[xy_info_unique2, on="id_split13"]  
s12 <- summcatch$id_split12[xy_info_unique2, on="id_split12"]  

s14 <- summcatch$id_split14[xy_info, on="id_split14"]  %>%
    left_join(., xy_sp, on = "id_cell")

s4 <- summcatch$id_split4[xy_info_unique2, on="id_split4"] # %>%
    # left_join(., xy_sp, on = "id_cell")


s0 <- summcatch$id_split0[xy_info_unique2, on="id_split0"]  
s4 <- summcatch$id_split4[xy_info, on="id_split4"] 

s6 <- summcatch$id_split6[xy_info, on="id_split6"] 

s14[s4, on = c("id_split14")] %>%
    ggplot(aes(i.n_species/n_species, median_logratio - i.median_logratio, 
               col = Provincias)) +
    geom_point()

s12[s4, on = c("id_split14")] %>%
    ggplot(aes(i.n_species/n_species, prop_losers  - i.prop_losers)) +
    geom_point() +
    facet_wrap(~Provincias, scales="free_x") +
    geom_hline(yintercept = 0) +
    labs(x = "nSp_region/nSp_subregion", 
         y = "prop_losers_subregion - prop_losers_region")

ggsave("figures/s14_s0_comparison_prop_losers.png")

s14[s0, on = c("id_split14")] %>%
    ggplot(aes(i.n_species/n_species, prop_losers  - i.prop_losers)) +
    geom_point() +
    # facet_wrap(~Provincias, scales="free_x") +
    geom_hline(yintercept = 0) +
    labs(x = "nSp_region/nSp_subregion", 
         y = "prop_losers_subregion - prop_losers_region")

ggsave("figures/s14_s0_comparison_prop_losers_nofacet.png")


s14[s0, on = c("id_split14")] %>%
    # filter(i.n_species/n_species > 10) %>%
    ggplot(aes(id_y, -id_x, fill = (i.n_species/n_species > 10))) +
    geom_tile()
    # geom_point() +
    # facet_wrap(~Provincias, scales="free_x") +
    # geom_hline(yintercept = 0) +
    # labs(x = "nSp_region/nSp_subregion", 
    #      y = "prop_losers_subregion - prop_losers_region")



s13[s0, on = c("id_cell")] %>%
    ggplot(aes(id_y, -id_x, fill=i.n_species/n_species)) +
    geom_tile() +
    coord_equal() +
    scale_fill_viridis_c()



    
s13[s0, on = c("id_cell")] %>%
    ggplot(aes(id_y, -id_x, fill=median_logratio)) +
    geom_tile() +
    coord_equal() +
    scale_fill_viridis_c()

s13[, mean(median_logratio)] 


s14 <- summcatch$id_split14[xy_info, on="id_split14"] 
s12 <- summcatch$id_split12[xy_info, on="id_split12"] 
s5 <- summcatch$id_split5[xy_info, on="id_split5"] 

s14[s5, on = c("id_cell")] %>%
    ggplot(aes(id_y, -id_x, fill=i.median_logratio - median_logratio)) +
    geom_tile() +
    coord_equal() +
    scale_fill_gradient2()

s14[s12, mean(i.median_logratio - median_logratio), on = c("id_cell")] 
    
ggplot(aes(i.median_logratio - median_logratio)) +
    geom_histogram()


summcatch$id_split14[xy_info, on="id_split14"] %>%
    ggplot(aes(id_y, -id_x, fill=median_logratio)) +
    geom_tile() +
    coord_equal()


# unpack
region_lookup <- xy_info %>% select(id_split0 : id_split14) %>% unique

catch <- vector("list")
catch2 <- vector("list")
catch3 <- vector("list")

counter <- 0
catch <- list()
for(subregion_i in 14:1) {
    for(region_i in (subregion_i-1):0) {
        
        counter <- counter + 1
        
        regioname_i <- paste0("id_split", region_i)
        subregioname_i <- paste0("id_split", subregion_i)
        
        summ_subregion <- copy(summcatch[[subregioname_i]])
        summ_region <- copy(summcatch[[regioname_i]])
        
        lookup_i <- unique(region_lookup[,.SD, .SDcols = c(regioname_i, subregioname_i)])
        
        summ_region <- summ_region[lookup_i, on=regioname_i]
        summ_subregion <- summ_subregion[summ_region, on = subregioname_i]
        
        # prop losers subregion - prop losers region
        # expect fewer losers in the subregion (i.e. negative)
        # hist(summ_subregion[,prop_losers - i.prop_losers])
        # i. signfies regional estimate, so:
        # +ve diff signifies subregion estimates larger than regional estimates 
        # -ve diff signifies subregion estimates smaller than regional estimates
        catch[[counter]] <- summ_subregion[,list(subregion_scale = subregion_i, 
                                                 id_subregion = get(subregioname_i),
                                                 region_scale = region_i, 
                                                 #S_diff = S_ratio - i.S_ratio, 
                                                 median_lr_diff = median_logratio - i.median_logratio, 
                                                 plosers_diff = prop_losers - i.prop_losers, 
                                                 avg_subregion_species = mean(n_species), 
                                                 region_species = unique(i.n_species)), 
                                           by=eval(regioname_i)]
    }
}


summcatch[[length(summcatch)]]
catch[[1]]
summcatch[[14]][,median_logratio - 2.03] %>% hist

lapply(catch, function(x) mean(x$plosers_diff)) %>% unlist %>% plot

# median_logratio
rbindlist(catch) %>%
    group_by(subregion_scale, region_scale) %>%
    summarise(m = mean(median_lr_diff)) %>%
    ggplot(aes(region_scale, subregion_scale, fill=m)) + 
    geom_tile() +
    geom_text(aes(label = round(m, 2)), col="white") +
    scale_x_continuous(breaks=0:12) +
    scale_y_continuous(breaks=0:12) +
    scale_fill_gradient2()


rbindlist(catch) %>%
    group_by(subregion_scale, region_scale) %>%
    summarise(m = mean(plosers_diff)) %>%
    ggplot(aes(region_scale, subregion_scale, fill=m)) + 
    geom_tile() +
    geom_text(aes(label = round(m, 2)), col="white") +
    scale_x_continuous(breaks=0:12) +
    scale_y_continuous(breaks=0:12) +
    scale_fill_gradient2()




catch_subset <- rbindlist(catch) %>%
    filter(subregion_scale %in% c(4, 8, 14)) 

catch_subset_summ <- catch_subset[, list(median_lr_diff = mean(median_lr_diff), 
                                         lwr = quantile(median_lr_diff, .1), 
                                         upr = quantile(median_lr_diff, .9)), 
                                  by=c("region_scale", "subregion_scale")]    

# mutate(region_scale = paste0("2^-", region_scale)) %>%
ggplot(catch_subset, aes(factor(region_scale), median_lr_diff)) + 
    geom_violin(col="grey50", fill="grey80", scale = "width") +
    geom_point(data=catch_subset_summ, aes(y=median_lr_diff)) +
    geom_linerange(data=catch_subset_summ, aes(ymin = lwr, ymax = upr)) +
    facet_grid(subregion_scale ~ .) + 
    geom_hline(yintercept=0) +
    theme(aspect.ratio = 0.5)  +
    theme_bw() +
    theme(panel.grid=element_blank(), 
          axis.text = element_text(colour = "black"))

ggplot(catch_subset, aes(factor(region_scale), plosers_diff)) + 
    geom_violin(col="grey50", fill="grey80", scale = "width") +
    # geom_point(data=catch_subset_summ, aes(y=median_lr_diff)) +
    # geom_linerange(data=catch_subset_summ, aes(ymin = lwr, ymax = upr)) +
    facet_grid(subregion_scale ~ .) + 
    geom_hline(yintercept=0) +
    theme(aspect.ratio = 0.5)  +
    theme_bw() +
    theme(panel.grid=element_blank(), 
          axis.text = element_text(colour = "black"))
ggsave("figures/plosers_1pct_threshold.png")



catch_subset_summ <- catch_subset[, list(plosers_diff = mean(plosers_diff), 
                                         lwr = quantile(plosers_diff, .1), 
                                         upr = quantile(plosers_diff, .9)), 
                                  by=c("region_scale", "subregion_scale")]  

xy_sp
xy_info[as.data.table(xy_sp), Provincias := i.Provincias, on = "id_cell"]
xy_info_unique <- copy(xy_info[,.SD, .SDcols = eval(names(xy_info)[5:20])])
xy_info_unique2 <- xy_info_unique[, list(Provincias = names(which.max(table(Provincias)))), 
               by = eval(names(xy_info_unique)[-16])]

# mutate(region_scale = paste0("2^-", region_scale)) %>%
ggplot(catch_subset,aes(factor(region_scale), plosers_diff)) + 
    geom_violin(col="grey50", fill="grey80", scale = "width") +
    geom_point(data=catch_subset_summ, aes(y=plosers_diff)) +
    geom_linerange(data=catch_subset_summ, aes(ymin = lwr, ymax = upr)) +
    facet_grid(subregion_scale ~ .) + 
    geom_hline(yintercept=0) +
    theme(aspect.ratio = 0.5)  +
    theme_bw() +
    theme(panel.grid=element_blank(), 
          axis.text = element_text(colour = "black")) +
    labs(y = "prop_losers(subregion) - prop_losers(region)")
ggsave("figures/subregion_14_comparisons.png")



lapply(catch, function(x) setnames(x, 1, "id_region"))

temp <- rbindlist(catch) %>%
    filter(subregion_scale == 14) %>% 
    mutate(id_split14 = id_subregion) %>%
    left_join(., xy_info_unique2 , by="id_split14") 

temp %>%
    ggplot(aes(factor(region_scale), plosers_diff)) + 
    geom_violin(col="grey50", fill="grey80", scale = "width") +
    # geom_point(data=catch_subset_summ, aes(y=plosers_diff)) +
    # geom_linerange(data=catch_subset_summ, aes(ymin = lwr, ymax = upr)) +
    facet_wrap(. ~ Provincias) + 
    geom_hline(yintercept=0) +
    theme(aspect.ratio = 0.5)  +
    theme_bw() +
    theme(panel.grid=element_blank(), 
          axis.text = element_text(colour = "black"))
ggsave("figures/subregion_14_comparisons_faceted.png")

rbindlist(catch) %>%
    filter(subregion_scale == 14) %>%
    left_join(., xy_info_unique2) %>%
    # mutate(region_scale = paste0("2^-", region_scale)) %>%
    ggplot(aes(region_species/avg_subregion_species, median_lr_diff)) +
    geom_point()
    geom_violin() +
    facet_grid(subregion_scale~.) + 
    geom_hline(yintercept=0) +
    theme(aspect.ratio = 0.5) 



# "SR"
rbindlist(catch) %>%
    group_by(subregion_scale, region_scale) %>%
    summarise(m = mean(S_diff)) %>%
    ggplot(aes(region_scale, subregion_scale, fill=m)) + 
    geom_tile() +
    geom_text(aes(label = round(m, 2)), col="white") +
    scale_x_continuous(breaks=0:12) +
    scale_y_continuous(breaks=0:12) +
    scale_fill_gradient2()


rbindlist(catch) %>%
    filter(subregion_scale %in% c( 4, 12, 14)) %>%
    # mutate(region_scale = paste0("2^-", region_scale)) %>%
    ggplot(aes(factor(region_scale), plosers_diff)) + 
    geom_violin() +
    facet_grid(subregion_scale~.) + 
    geom_hline(yintercept=0) +
    theme(aspect.ratio = 0.5) 



plot(unlist(catch2) ~ rep(1:12, each = 4096))
abline(0,0)
lapply(catch, mean)

(means <- sapply(catch3, mean))
(sds <- sapply(catch2, sd))

ggplot() + 
    geom_violin(aes(factor(rep(1:12, each = 4096)), unlist(catch2)))



subregion_samps <- sample(unique(summ_i$id_region), 12)

summ_i[id_region %in% subregion_samps] %>%
    ggplot(aes(N_pasture/16)) + 
    geom_histogram(boundary=0, binwidth=.05, fill="black") +
    geom_histogram(aes(N_forest/16), fill="lightblue", alpha = .8, boundary=0, binwidth=.05)  +
    facet_wrap(~id_region)

summ_i[id_subregion %in% subregion_samps & (N_forest/16 > 0.05 | N_pasture/16 > 0.05)] %>%
    ggplot(aes(N_pasture/16)) + 
    geom_histogram(boundary=0, binwidth=.05, fill="black") +
    geom_histogram(aes(N_forest/16), fill="lightblue", alpha = .8, boundary=0, binwidth=.05)  +
    facet_wrap(~id_subregion)

summ_i[id_subregion %in% subregion_samps & (N_forest/16 > 0.05 | N_pasture/16 > 0.05),
       median(sum(N_forest)/sum(N_pasture)),
       by = id_subregion] 




    
summ2[, exp(logratio_super) - exp(logratio)] %>%
    hist

rbindlist(summcatch, idcol = "id") %>%
    ggplot(aes(logratio_super, logratio)) + geom_point(alpha=.5) + 
    coord_equal() + geom_abline(col="red") + 
    facet_wrap(~id)

summ_all <- rbindlist(summcatch, idcol = "id")

x <- summ_all[, mean(exp(logratio_super) - exp(logratio)), by = "id"]$V1 + 1 %>%
    as.vector %>%
    sum(.)

prod(x)

id         V1
1:  1 0.27664033
2:  2 0.10134633
3:  3 0.11699935
4:  4 0.08646694
5:  5 0.08131091
6:  6 0.08147070
7:  7 0.04993837

rbindlist(summcatch, idcol = "id") %>%
    ggplot(aes(exp(logratio_super) - exp(logratio), y = stat(density))) + 
    geom_histogram() +
    geom_vline(xintercept=0) + 
    # coord_equal() + geom_abline(col="red") + 
    facet_wrap(~id)

    

summ_cell_i <- pred_dt[N_pasture > 1 | N_forest > 1, 
                           list(#N_pasture = (N_pasture/16)/.N, #%/% 0.01, 
                                #N_forest = (N_forest/16)/.N, #%/% 0.01,
                                med_logratio = median(log(N_forest/N_pasture))), 
                           by=c("id_cell", "id_subregion")]
    
    summ_cell_i2 <- summ_cell_i[, list(mean_med_logratio = mean(med_logratio)), 
                                by = c("id_subregion")]
    
    summ_i[summ_cell_i2, mean_med_logratio := i.mean_med_logratio, on = "id_subregion"]
    
    
    summ_i <- pred_dt[,list(N_pasture = sum(N_pasture/16)/.N, #%/% 0.01, 
                            N_forest = sum(N_forest/16)/.N, #%/% 0.01,
                            logratio = log(sum(N_forest)/sum(N_pasture))), 
                      by=c("species", "id_subregion")]
    
    
    
    
    summ_i[]
    summ_i[,sum(N_pasture > .50)/.N]
    summ_i$N_pasture %>% max
    
    sp <- summ_i[, .N, by = c("N_pasture", "id_subregion")]
    sp2 <- sp[N >= N_pasture, max(N_pasture) * 0.01, by = "id_subregion"]
    
    sf <- summ_i[, .N, by = c("N_forest", "id_subregion")]
    sf2 <- sf[N >= N_forest, max(N_forest) * 0.01, by = "id_subregion"]

    mean(sp2$V1)/mean(sf2$V1) 
    0.14/0.18
    
    summ_i[,sum(N_forest > .50)/.N]
    summ_i[,.N, by = cut(N_forest, seq(0, 1e6, ))]
    table(round(summ_i[,"N_forest"], -1)) %>% as.vector()
    
    ggplot(summ_i, aes(N_forest)) + geom_histogram(boundary=0) +
        geom_histogram(aes(N_pasture), fill="red", alpha=.4, boundary=0)
    
    ggplot(summ_i, aes(log(N_forest), log(N_pasture))) + 
        geom_point() + 
        geom_point(data=summ_old, col="blue", alpha=.5) +
        geom_abline(col="red") +
        scale_x_continuous(limits = c(log(0.01), log(1))) + #, breaks = 2^(0:-4)) +
        scale_y_continuous(limits = c(log(0.01), log(1))) + #, breaks = 2^(0:-4)) +
        coord_equal() 
    # geom_hline(yintercept = 0.01, col="red", lty = "longdash", limits = c(0.01)) +
    # geom_vline(xintercept = 0.01, col="red", lty = "longdash", limits = c(0.01))
    cut(c(0.1, 0.1, 0.2, .5), seq(0, 1, .1))
    c(0.1, 0.1, 0.2, .5) %/% 0.01 * 0.01
    
    table(c(0.1, 0.1, 0.2, .5) %/% 0.1)
    
    # summ_old <- copy(summ_i)
    
    ggsave("temp_plot.png")
    
    geom_histogram(boundary=0) +
        geom_histogram(aes(N_pasture), fill="red", alpha=.4, boundary=0)
    
    
    summ_catch[[i]] <- summ_i[,list(med_lr = median(logratio)), by="id_subregion"]
}


lapply(summ_catch, function(x) x[,list(med_lr = mean(logratio)), by="id_subregion"]$med_lr %>% median) %>% unlist

plot(lapply(summ_catch, function(x) median(x$med_lr)) %>% unlist ~ res_seq)
library(ggplot2)

allcatch <- rbindlist(summ_catch, idcol = "id")

ggplot(allcatch) + 
    geom_histogram(aes(x=med_lr)) +
    facet_grid(id~.)

xy_info <- xy_info %>%
    # calculate 20 km subregions
    dplyr::mutate(id_subregion_x = ceiling(id_x/200), 
                  id_subregion_y = ceiling(id_y/200), 
                  id_subregion = xy_to_cell(id_subregion_x, id_subregion_y, c(max(id_subregion_x), max(id_subregion_y))))

pred_dt[xy_info, id_subregion := i.id_subregion, on="id_cell"]

summ2 <- pred_dt[N_forest > 1/3 | N_pasture > 1/3, list(N_pasture = sum(N_pasture), 
                       N_forest = sum(N_forest),
                       ratio = sum(N_forest)/sum(N_pasture)), by=c("species", "id_subregion")]

summ3 <- summ2[,list(med_lr = mean(log(ratio))), by="id_subregion"]

hist(summ3$med_lr)
mean(log(summ$N_forest/summ$N_pasture))
median(summ3$med_lr)



mean(log(summ$N_forest/summ$N_pasture))
