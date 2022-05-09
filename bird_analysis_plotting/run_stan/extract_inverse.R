library("cmdstanr")
v6_2 <- read_cmdstan_csv("/Users/jacobsocolar/Downloads/occupancy_v6_threads-202102191542-1-557a60.csv")

inv_metric <- v6_2$inv_metric[[1]]

saveRDS(inv_metric, "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/inv_metric_6_2.RDS")
