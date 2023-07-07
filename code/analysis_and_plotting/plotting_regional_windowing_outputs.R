# library(data.table); library(dplyr); library(ggplot2)
# 
# sample_cols <- function(x) {
#     len <- length(unique(x))
#     terrain.colors(len)[sample(1:len)]
# }
# 
# #source("bird_analysis_plotting/occupancy_rasters/helper_functions.R")
# 
# # data ----
# pred_dt <- readRDS("outputs/prediction_info_dt.rds")
# pred_dt <- pred_dt[!is.na(tdist) & !is.na(relev), ]
# 
# # N_forest <- readRDS("outputs/predicted_occupancy_dts/posterior_forest_1.rds")
# # N_pasture <- readRDS("outputs/predicted_occupancy_dts/posterior_pasture_1.rds")
# 
# post_av <- readRDS("outputs/predicted_occupancy_dts/averaged_posterior_10.rds")
# 
# # note: in each cell, a species could be on 48 points (16 * 3)
# pred_dt[,`:=`(p_pasture = post_av$N_pasture_update/48,
#               p_forest = post_av$N_forest_update/48)]
# pred_dt[,`:=`(relev = NULL, tdist = NULL)]



pred_summ <- pred_dt[,list(N_pasture = sum(p_pasture*48), N_forest = sum(p_forest*48)), by="species"]

ggplot(pred_summ, aes((N_pasture), (N_forest))) + 
    geom_point(alpha=.6) +
    geom_abline(col="red") +
    scale_x_continuous(trans="log", breaks = 10^(-10:10)) + 
    scale_y_continuous(trans="log", breaks = 10^(-100:10)) + 
    coord_equal() + 
    stat_smooth() +
    labs(x = "Number of pasture points", 
         y = "Number of forest points") + 
    # theme_bw() +
    theme(axis.text = element_text(colour = "black"), 
          panel.border = element_rect(colour="black", fill=NA))
ggsave("figures/predicted_number_of_points.png", units="mm", height=120, width=160)


ggplot(pred_summ, aes((N_forest), N_pasture)) + 
    geom_point(alpha=.6) +
    geom_abline(col="red") +
    scale_x_continuous(trans="log", breaks = 10^(-10:10)) + 
    scale_y_continuous(trans="log", breaks = 10^(-100:10)) + 
    coord_equal() + 
    stat_smooth() +
    labs(x = "Number of pasture points", 
         y = "Number of forest points") + 
    # theme_bw() +
    theme(axis.text = element_text(colour = "black"), 
          panel.border = element_rect(colour="black", fill=NA))

pred_summ[, `:=`(f_p = cut(N_pasture, seq(0, 1e7, 1000)), 
                 f_f = cut(N_forest, seq(0, 1e7, 1000)))]
ps2 <- pred_summ[,list(N_pasture = mean(N_pasture), N_forest = mean(N_forest)), by=.(f_p, f_f)]


ggplot(pred_summ, aes((N_forest), N_pasture)) + 
    geom_point(alpha=.6) +
    geom_abline(col="red") +
    scale_x_continuous(trans="log", breaks = 10^(-10:10)) + 
    scale_y_continuous(trans="log", breaks = 10^(-100:10)) + 
    coord_equal() + 
    stat_smooth() +
    # geom_point(data=ps2, col="red") +
    labs(x = "Number of pasture points", 
         y = "Number of forest points") + 
    # theme_bw() +
    theme(axis.text = element_text(colour = "black"), 
          panel.border = element_rect(colour="black", fill=NA))




df <- readRDS("outputs/beta_prop_loss_comparisons_post1_thresh10.rds")
df
# df[,list(beta = mean(beta), 
#          eff = mean(median_logratio_regional - median_logratio_local)), by=id] %>%

ggplot(df, aes(beta, exp(median_logratio_regional - median_logratio_local), col=N)) +
    geom_point(alpha = .5) +
    geom_hline(yintercept = 1, lty="longdash") + 
    theme(axis.text = element_text(colour = "black"), 
          strip.text = element_text(hjust =0, face="bold"), 
          strip.background = element_rect(colour=NA),
          panel.border = element_rect(colour="black", fill=NA)) +
    labs(y = bquote(medianLR[regional]-medianLR[point]), 
         x = "Beta") +
    #scale_y_continuous(trans="exp") +
    scale_colour_viridis_c() + 
    facet_wrap(~cut(N, seq(0, 1600, 100)))

df_sub <- df %>%
    group_by(id) %>%
    mutate(row_id = 1:n()) %>%
    filter(row_id %in% if(max(row_id) > 200) sample(row_id, 200) else row_id)

ggplot(df, aes(id, exp(median_logratio_regional - median_logratio_local), col=N)) +
    geom_violin() + 
    geom_jitter(data=df_sub, height=0, width=.1, alpha=.2) +
    # geom_hline(yintercept = 0, lty="longdash") + 
    theme(axis.text = element_text(colour = "black"), 
          strip.text = element_text(hjust =0, face="bold"), 
          strip.background = element_rect(colour=NA),
          panel.border = element_rect(colour="black", fill=NA)) +
    labs(y = bquote(medianLR[regional]-medianLR[point]), 
         x = "Beta") +
    #scale_y_continuous(trans="exp") +
    scale_colour_viridis_c() + facet_wrap(~cut(point_richness_forest, 6)) 



# plot local vs regional losses ----
ggplot(df, mapping = aes(id, exp(median_logratio_local - median_logratio_regional))) +
    geom_blank() +
    stat_summary(data = df[id %in% paste0("id_split", 0:2)],
                 fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid") +
    geom_boxplot(data = df[!(id %in% paste0("id_split", 0:2))], fill="grey95") +
    geom_point(data = df[id %in% paste0("id_split", 1:2)]) + 
    geom_hline(yintercept = 1, lty="longdash") +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"),
          panel.grid = element_blank()) +
    scale_y_continuous(trans="log", breaks=2^(-5:5)) +
    scale_x_discrete(labels=c(bquote(s[0]), 
                              bquote(s[1]),
                              bquote(s[2]),
                              bquote(s[3]),
                              bquote(s[4]),
                              bquote(s[5]),
                              bquote(s[6]),
                              bquote(s[7]),
                              bquote(s[8]),
                              bquote(s[9]),
                              bquote(s[10]),
                              bquote(s[11]),
                              bquote(s[12]),
                              bquote(s[13]),
                              bquote(s[14])
    )) +
    labs(x = "Region scale", 
         y = bquote(medianLR[regional]-medianLR[point]))

ggplot(df, mapping = aes(id, exp(median_logratio_local - median_logratio_regional))) +
    geom_blank() +
    stat_summary(data = df[id %in% paste0("id_split", 0:2)],
                 fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid") +
    geom_boxplot(data = df[!(id %in% paste0("id_split", 0:2))], fill="grey95") +
    geom_point(data = df[id %in% paste0("id_split", 1:2)]) + 
    geom_hline(yintercept = 1, lty="longdash") +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"),
          panel.grid = element_blank()) +
    scale_y_continuous(trans="log", breaks=2^(-5:5)) +
    scale_x_discrete(labels=c(bquote(s[0]), 
                              bquote(s[1]),
                              bquote(s[2]),
                              bquote(s[3]),
                              bquote(s[4]),
                              bquote(s[5]),
                              bquote(s[6]),
                              bquote(s[7]),
                              bquote(s[8]),
                              bquote(s[9]),
                              bquote(s[10]),
                              bquote(s[11]),
                              bquote(s[12]),
                              bquote(s[13]),
                              bquote(s[14])
    )) +
    labs(x = "Region scale", 
         y = bquote(medianLR[point]-medianLR[regional]))
ggsave("figures/scaling_region_size_flipped.png")


ggplot(df, mapping = aes(id, (1-lose_prop_regional) - (1-lose_prop_local))) +
    geom_blank() +
    stat_summary(data = df[id %in% paste0("id_split", 0:2)],
                 fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid") +
    geom_boxplot(data = df[!(id %in% paste0("id_split", 0:2))], fill="grey95") +
    geom_point(data = df[id %in% paste0("id_split", 1:2)]) + 
    geom_hline(yintercept = 0, lty="longdash") +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"),
          panel.grid = element_blank()) +
    scale_y_continuous() +
    scale_x_discrete(labels=c(bquote(s[0]), 
                              bquote(s[1]),
                              bquote(s[2]),
                              bquote(s[3]),
                              bquote(s[4]),
                              bquote(s[5]),
                              bquote(s[6]),
                              bquote(s[7]),
                              bquote(s[8]),
                              bquote(s[9]),
                              bquote(s[10]),
                              bquote(s[11]),
                              bquote(s[12]),
                              bquote(s[13]),
                              bquote(s[14])
    )) +
    labs(x = "Region scale", 
         y = bquote(medianLR[regional]-medianLR[point]))

ggsave("figures/scaling_region_size.png")



# p1 <- 
ggplot(df, aes(beta, exp(median_logratio_regional - median_logratio_local))) +
    geom_point(alpha = .1) +
    geom_hline(yintercept = 0, lty="longdash") + 
    scale_y_continuous(trans="log", breaks = 2^(-5:5)) +
    theme_bw() +
    geom_hline(yintercept = 1, lty = "longdash") +
    theme(axis.text = element_text(colour = "black"), 
          panel.grid = element_blank(),
          strip.text = element_text(hjust =0, face="bold"), 
          strip.background = element_rect(colour=NA),
          panel.border = element_rect(colour="black", fill=NA)) +
    labs(y = bquote(medianLR[regional]-medianLR[local]), 
         x = "Beta")
ggsave("figures/LR_beta.png")

ggplot(df, aes(beta, exp(median_logratio_local - median_logratio_regional))) +
    geom_point(alpha = .1) +
    geom_hline(yintercept = 0, lty="longdash") + 
    scale_y_continuous(trans="log", breaks = 2^(-5:5)) +
    theme_bw() +
    geom_hline(yintercept = 1, lty = "longdash") +
    theme(axis.text = element_text(colour = "black"), 
          panel.grid = element_blank(),
          strip.text = element_text(hjust =0, face="bold"), 
          strip.background = element_rect(colour=NA),
          panel.border = element_rect(colour="black", fill=NA)) +
    labs(y = bquote(medianLR[local]-medianLR[regional]), 
         x = "Beta")
ggsave("figures/LR_beta_flipped.png")



####
source("code/analysis_and_plotting/helper_functions.R")
ggplot(xy_info[seq(1, .N, 4),], aes(id_y, -id_x, fill=exp(median_logratio_local - median_logratio_regional))) +
    geom_tile() +
    coord_equal() +
    guides(colour = "none") +
    scale_fill_viridis_c()
# ggsave("figures/split2.png")

df_sub <- copy(df[id == "id_split14",])
xy_sub <- copy(xy_info)
xy_sub[,id_region := id_split14]

df_sub[xy_sub, on="id_region"] %>%
    ggplot(aes(id_y, -id_x, fill = exp(median_logratio_local - median_logratio_regional))) +
    geom_tile() +
    coord_equal() +
    guides(colour = "none") +
    scale_fill_gradient2(midpoint=1)


df_sub <- copy(df[id == "id_split6",])
xy_sub <- copy(xy_info)
xy_sub[,id_region := id_split6]

df_sub[xy_sub, on="id_region"] %>%
    ggplot(aes(id_y, -id_x, fill = exp(median_logratio_local - median_logratio_regional))) +
    geom_tile() +
    coord_equal() +
    guides(colour = "none") +
    scale_fill_gradient2(midpoint=1) + 
    theme_






df[is.infinite(log(exp(exp(median_logratio_regional - median_logratio_local)))),]

ggplot(df, aes(point_richness_forest, median_logratio_regional - median_logratio_local)) +
    geom_point(alpha = .1) +
    # geom_hline(yintercept = 0, lty="longdash") + 
    theme(axis.text = element_text(colour = "black"), 
          strip.text = element_text(hjust =0, face="bold"), 
          strip.background = element_rect(colour=NA),
          panel.border = element_rect(colour="black", fill=NA)) +
    # labs(y = bquote(medianLR[regional]-medianLR[point]), 
         # x = "Beta") +
    scale_colour_viridis_c() +
    facet_wrap(~id)
    facet_wrap(~cut(N, seq(0, 1800, 200)))



ggplot(df, aes(beta, median_logratio_regional - median_logratio_local)) +
    geom_point(alpha = .1) +
    geom_hline(yintercept = 0, lty="longdash") + 
    theme(axis.text = element_text(colour = "black"), 
          strip.text = element_text(hjust =0, face="bold"), 
          strip.background = element_rect(colour=NA),
          panel.border = element_rect(colour="black", fill=NA)) +
    labs(y = bquote(medianLR[regional]-medianLR[point]), 
         x = "Beta") + 
    scale_colour_viridis_c() +
    facet_wrap(~cut(1/point_richness_forest, 6))

p1
id_subset <- levels(df$id)[c(seq(1, 14, 2), 14, 15)]

p2 <- ggplot(df[id %in% id_subset,], aes(beta, median_logratio_regional - median_logratio_local)) +
    geom_point(alpha = .1) + 
    facet_wrap(~id) + 
    geom_hline(yintercept = 0, lty="longdash") + 
    theme(axis.text = element_text(colour = "black"), 
          strip.text = element_text(hjust =0, face="bold"), 
          strip.background = element_rect(colour=NA, fill=NA),
          panel.border = element_rect(colour="black", fill=NA)) +
    labs(y = bquote(medianLR[regional]-medianLR[point]), 
         x = "Beta")

p_both <- egg::ggarrange(p1, p2, ncol=1)
ggsave("figures/regional_loss.png", plot = p_both, units="mm", width=100*1.2, height=162*1.2)


p1 <- ggplot(df, aes(beta, lose_prop_regional  - lose_prop_local)) +
    geom_point(alpha = .1) +
    geom_hline(yintercept = 0, lty="longdash") + 
    theme(axis.text = element_text(colour = "black"), 
          strip.text = element_text(hjust =0, face="bold"), 
          strip.background = element_rect(colour=NA),
          panel.border = element_rect(colour="black", fill=NA))

p2 <- ggplot(df[id %in% id_subset,], aes(beta, lose_prop_regional  - lose_prop_local)) +
    geom_point(alpha = .1) + 
    facet_wrap(~id) + 
    geom_hline(yintercept = 0, lty="longdash") + 
    theme(axis.text = element_text(colour = "black"), 
          strip.text = element_text(hjust =0, face="bold"), 
          strip.background = element_rect(colour=NA, fill=NA),
          panel.border = element_rect(colour="black", fill=NA))

p_both <- egg::ggarrange(p1, p2, ncol=1)
ggsave("figures/regional_loss_losers.png", plot = p_both, units="mm", width=100*1.2, height=162*1.2)



hist(df[id == "id_split14", N])
