# Script to keep track of how many mystery recordings in different categories I have

setwd('/Users/jacobsocolar/Dropbox/Work/Colombia/Colombia_recordings')
sites <- list.files()

low <- 0
med <- 0
high <- 0
veryhigh <- 0
for(i in 1:length(sites)){
  files2 <- list.files(sites[i])
  low <- low + sum(grepl('_low', files2))
  med <- med + sum(grepl('_med', files2))
  high <- high + sum(grepl('_high', files2))
  veryhigh <- veryhigh + sum(grepl('_veryhigh', files2))
  
}
