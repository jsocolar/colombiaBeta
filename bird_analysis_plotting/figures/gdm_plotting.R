library(gdm)
dev.off()
birds_raw <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/gdm_outputs/birds/gdms_raw_raup.RDS")
birds_modeled <- readRDS("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/Analysis/gdm_outputs/birds/gdms_modeled.RDS")

ispline_q_extract <- function(model_list, predictor, p_quantile){
  spline_list <- lapply(model_list, isplineExtract)
  values_matrix <- do.call(cbind, lapply(spline_list, function(x){x$y[,predictor]}))
  quantile <- apply(values_matrix, 1, function(x){quantile(x, p_quantile)})
  return(quantile)
}

##### Raw data ####

pdf(file = "/Users/jacobsocolar/Dropbox/Work/Colombia/Writing/Figures/gdm_raw.PDF",
    width = 6, height = 4)
par(mfrow=c(2,3))
par(mar=c(3,2,0.5,0.5))

# Geographic distance
pred <- "Geographic"
ymax <- 1.05*max(max(ispline_q_extract(birds_raw$forest_bb, pred, .95)), max(ispline_q_extract(birds_raw$pasture_bb, pred, .95)))
plot(isplineExtract(birds_raw$forest_obs)$x[,pred]*111, isplineExtract(birds_raw$forest_obs)$y[,pred],
     type="n", xlab="km", ylim = c(0,9), ylab="", mgp=c(1.9,.7,0))
title("Geographic distance", line = -1.5, cex.main=1)
x <- isplineExtract(birds_raw$forest_obs)$x[,pred]*111
xp <- isplineExtract(birds_raw$pasture_obs)$x[,pred]*111
for(i in 1:9){
  p_lower <- .5 - .05*i
  p_upper <- .5 + .05*i
  q_lower_f <- ispline_q_extract(birds_raw$forest_bb, pred, p_lower)
  q_upper_f <- ispline_q_extract(birds_raw$forest_bb, pred, p_upper)
  polygon(c(x, rev(x)), c(q_lower_f, rev(q_upper_f)), col = rgb(red=0.2, green=0.2, blue=1.0, alpha=0.2), border = F)
  q_lower_p <- ispline_q_extract(birds_raw$pasture_bb, pred, p_lower)
  q_upper_p <- ispline_q_extract(birds_raw$pasture_bb, pred, p_upper)
  polygon(c(xp, rev(xp)), c(q_lower_p, rev(q_upper_p)), col = rgb(red=0.8, green=0.14, blue=0.15, alpha=0.2), border = F)
}
lines(isplineExtract(birds_raw$forest_obs)$x[,pred]*111, isplineExtract(birds_raw$forest_obs)$y[,pred], 
      type="l", col = "black", lwd=1)
lines(isplineExtract(birds_raw$pasture_obs)$x[,pred]*111, isplineExtract(birds_raw$pasture_obs)$y[,pred],
      type="l", col = "black", lwd=1, lty=2)

# Elevation
pred <- "elev_ALOS"
ymax <- 1.05*max(max(ispline_q_extract(birds_raw$forest_bb, pred, .95)), max(ispline_q_extract(birds_raw$pasture_bb, pred, .95)))
plot(isplineExtract(birds_raw$forest_obs)$x[,pred], isplineExtract(birds_raw$forest_obs)$y[,pred],
     type="n", xlab="masl", ylim = c(0,9), ylab = "", yaxt="n", mgp=c(1.9,.7,0))
title("Elevation", line = -1.5, cex.main=1)
x <- isplineExtract(birds_raw$forest_obs)$x[,pred]
xp <- isplineExtract(birds_raw$pasture_obs)$x[,pred]
for(i in 1:9){
  p_lower <- .5 - .05*i
  p_upper <- .5 + .05*i
  q_lower_f <- ispline_q_extract(birds_raw$forest_bb, pred, p_lower)
  q_upper_f <- ispline_q_extract(birds_raw$forest_bb, pred, p_upper)
  polygon(c(x, rev(x)), c(q_lower_f, rev(q_upper_f)), col = rgb(red=0.2, green=0.2, blue=1.0, alpha=0.2), border = F)
  q_lower_p <- ispline_q_extract(birds_raw$pasture_bb, pred, p_lower)
  q_upper_p <- ispline_q_extract(birds_raw$pasture_bb, pred, p_upper)
  polygon(c(xp, rev(xp)), c(q_lower_p, rev(q_upper_p)), col = rgb(red=0.8, green=0.14, blue=0.15, alpha=0.2), border = F)
}
lines(isplineExtract(birds_raw$forest_obs)$x[,pred], isplineExtract(birds_raw$forest_obs)$y[,pred], 
      type="l", col = "black", lwd=1)
lines(isplineExtract(birds_raw$pasture_obs)$x[,pred], isplineExtract(birds_raw$pasture_obs)$y[,pred],
      type="l", col = "black", lwd=1, lty=2)

# Dummy plot to fill layout screen
plot(1,1, type='n', axes=F, xlab="", ylab="")

# Mountain barrier
pred <- "matrix_1"
ymax <- 1.05*max(max(ispline_q_extract(birds_raw$forest_bb, pred, .95)), max(ispline_q_extract(birds_raw$pasture_bb, pred, .95)))
plot(isplineExtract(birds_raw$forest_obs)$x[,pred], isplineExtract(birds_raw$forest_obs)$y[,pred],
     type="n", xlab="m below 4100", ylab="", ylim = c(0,1.1), mgp=c(1.9,.7,0))
title("Mountain barrier", line = -1.5, cex.main=1)
x <- isplineExtract(birds_raw$forest_obs)$x[,pred]
xp <- isplineExtract(birds_raw$pasture_obs)$x[,pred]
for(i in 1:9){
  p_lower <- .5 - .05*i
  p_upper <- .5 + .05*i
  q_lower_f <- ispline_q_extract(birds_raw$forest_bb, pred, p_lower)
  q_upper_f <- ispline_q_extract(birds_raw$forest_bb, pred, p_upper)
  polygon(c(x, rev(x)), c(q_lower_f, rev(q_upper_f)), col = rgb(red=0.2, green=0.2, blue=1.0, alpha=0.2), border = F)
  q_lower_p <- ispline_q_extract(birds_raw$pasture_bb, pred, p_lower)
  q_upper_p <- ispline_q_extract(birds_raw$pasture_bb, pred, p_upper)
  polygon(c(xp, rev(xp)), c(q_lower_p, rev(q_upper_p)), col = rgb(red=0.8, green=0.14, blue=0.15, alpha=0.2), border = F)
}
lines(isplineExtract(birds_raw$forest_obs)$x[,pred], isplineExtract(birds_raw$forest_obs)$y[,pred], 
      type="l", col = "black", lwd=1)
lines(isplineExtract(birds_raw$pasture_obs)$x[,pred], isplineExtract(birds_raw$pasture_obs)$y[,pred],
      type="l", col = "black", lwd=1, lty=2)

# Valley barrier
pred <- "matrix_2"
ymax <- 1.05*max(max(ispline_q_extract(birds_raw$forest_bb, pred, .95)), max(ispline_q_extract(birds_raw$pasture_bb, pred, .95)))
plot(isplineExtract(birds_raw$forest_obs)$x[,pred], isplineExtract(birds_raw$forest_obs)$y[,pred],
     type="n", xlab="masl", ylab = "", yaxt ='n', ylim = c(0,1.1), mgp=c(1.9,.7,0))
title("Valley barrier", line = -1.5, cex.main=1)
x <- isplineExtract(birds_raw$forest_obs)$x[,pred]
xp <- isplineExtract(birds_raw$pasture_obs)$x[,pred]
for(i in 1:9){
  p_lower <- .5 - .05*i
  p_upper <- .5 + .05*i
  q_lower_f <- ispline_q_extract(birds_raw$forest_bb, pred, p_lower)
  q_upper_f <- ispline_q_extract(birds_raw$forest_bb, pred, p_upper)
  polygon(c(x, rev(x)), c(q_lower_f, rev(q_upper_f)), col = rgb(red=0.2, green=0.2, blue=1.0, alpha=0.2), border = F)
  q_lower_p <- ispline_q_extract(birds_raw$pasture_bb, pred, p_lower)
  q_upper_p <- ispline_q_extract(birds_raw$pasture_bb, pred, p_upper)
  polygon(c(xp, rev(xp)), c(q_lower_p, rev(q_upper_p)), col = rgb(red=0.8, green=0.14, blue=0.15, alpha=0.2), border = F)
}
lines(isplineExtract(birds_raw$forest_obs)$x[,pred], isplineExtract(birds_raw$forest_obs)$y[,pred], 
      type="l", col = "black", lwd=1)
lines(isplineExtract(birds_raw$pasture_obs)$x[,pred], isplineExtract(birds_raw$pasture_obs)$y[,pred],
      type="l", col = "black", lwd=1, lty=2)

# Precipitation
pred <- "precip_ceccherini"
ymax <- 1.05*max(max(ispline_q_extract(birds_raw$forest_bb, pred, .95)), max(ispline_q_extract(birds_raw$pasture_bb, pred, .95)))
plot(isplineExtract(birds_raw$forest_obs)$x[,pred], isplineExtract(birds_raw$forest_obs)$y[,pred],
     type="n", xlab="mm", ylab = "", yaxt='n', ylim = c(0,1.1), mgp=c(1.9,.7,0))
title("Annual precipitation", line = -1.5, cex.main=1)
x <- isplineExtract(birds_raw$forest_obs)$x[,pred]
xp <- isplineExtract(birds_raw$pasture_obs)$x[,pred]
for(i in 1:9){
  p_lower <- .5 - .05*i
  p_upper <- .5 + .05*i
  q_lower_f <- ispline_q_extract(birds_raw$forest_bb, pred, p_lower)
  q_upper_f <- ispline_q_extract(birds_raw$forest_bb, pred, p_upper)
  polygon(c(x, rev(x)), c(q_lower_f, rev(q_upper_f)), col = rgb(red=0.2, green=0.2, blue=1.0, alpha=0.2), border = F)
  q_lower_p <- ispline_q_extract(birds_raw$pasture_bb, pred, p_lower)
  q_upper_p <- ispline_q_extract(birds_raw$pasture_bb, pred, p_upper)
  polygon(c(xp, rev(xp)), c(q_lower_p, rev(q_upper_p)), col = rgb(red=0.8, green=0.14, blue=0.15, alpha=0.2), border = F)
}
lines(isplineExtract(birds_raw$forest_obs)$x[,pred], isplineExtract(birds_raw$forest_obs)$y[,pred], 
      type="l", col = "black", lwd=1)
lines(isplineExtract(birds_raw$pasture_obs)$x[,pred], isplineExtract(birds_raw$pasture_obs)$y[,pred],
      type="l", col = "black", lwd=1, lty=2)

dev.off()

##### Modeled data ####
# Geographic distance
pred <- "Geographic"
ymax <- 1.05*max(max(ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, .95)), max(ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, .95)))
plot(isplineExtract(birds_raw$forest_obs)$x[,pred]*111, isplineExtract(birds_raw$forest_obs)$y[,pred],
     type="n", xlab="Geographic distance (km)", ylab="Partial ecological distance",
     ylim = c(0,ymax))
x <- isplineExtract(birds_raw$forest_obs)$x[,pred]*111
xp <- isplineExtract(birds_raw$pasture_obs)$x[,pred]*111
for(i in 1:9){
  p_lower <- .5 - .05*i
  p_upper <- .5 + .05*i
  q_lower_f <- ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, p_lower)
  q_upper_f <- ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, p_upper)
  polygon(c(x, rev(x)), c(q_lower_f, rev(q_upper_f)), col = rgb(red=0.2, green=0.2, blue=1.0, alpha=0.2), border = F)
  q_lower_p <- ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, p_lower)
  q_upper_p <- ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, p_upper)
  polygon(c(xp, rev(xp)), c(q_lower_p, rev(q_upper_p)), col = rgb(red=0.8, green=0.14, blue=0.15, alpha=0.2), border = F)
}
lines(isplineExtract(birds_raw$forest_obs)$x[,pred]*111, ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, .5), 
      type="l", col = "black", lwd=1)
lines(isplineExtract(birds_raw$pasture_obs)$x[,pred]*111, ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, .5),
      type="l", col = "black", lwd=1, lty=2)

# Elevation
pred <- "elev_ALOS"
ymax <- 1.05*max(max(ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, .95)), max(ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, .95)))
plot(isplineExtract(birds_raw$forest_obs)$x[,pred], isplineExtract(birds_raw$forest_obs)$y[,pred],
     type="n", xlab="Elevation (m)", ylab="Partial ecological distance",
     ylim = c(0,ymax))
x <- isplineExtract(birds_raw$forest_obs)$x[,pred]
xp <- isplineExtract(birds_raw$pasture_obs)$x[,pred]
for(i in 1:9){
  p_lower <- .5 - .05*i
  p_upper <- .5 + .05*i
  q_lower_f <- ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, p_lower)
  q_upper_f <- ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, p_upper)
  polygon(c(x, rev(x)), c(q_lower_f, rev(q_upper_f)), col = rgb(red=0.2, green=0.2, blue=1.0, alpha=0.2), border = F)
  q_lower_p <- ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, p_lower)
  q_upper_p <- ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, p_upper)
  polygon(c(xp, rev(xp)), c(q_lower_p, rev(q_upper_p)), col = rgb(red=0.8, green=0.14, blue=0.15, alpha=0.2), border = F)
}
lines(isplineExtract(birds_raw$forest_obs)$x[,pred], ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, .5), 
      type="l", col = "black", lwd=1)
lines(isplineExtract(birds_raw$pasture_obs)$x[,pred], ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, .5),
      type="l", col = "black", lwd=1, lty=2)

# Precipitation
pred <- "precip_ceccherini"
ymax <- 1.05*max(max(ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, .95)), max(ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, .95)))
plot(isplineExtract(birds_raw$forest_obs)$x[,pred], isplineExtract(birds_raw$forest_obs)$y[,pred],
     type="n", xlab="Annual precipitation (mm)", ylab="Partial ecological distance",
     ylim = c(0,ymax))
x <- isplineExtract(birds_raw$forest_obs)$x[,pred]
xp <- isplineExtract(birds_raw$pasture_obs)$x[,pred]
for(i in 1:9){
  p_lower <- .5 - .05*i
  p_upper <- .5 + .05*i
  q_lower_f <- ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, p_lower)
  q_upper_f <- ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, p_upper)
  polygon(c(x, rev(x)), c(q_lower_f, rev(q_upper_f)), col = rgb(red=0.2, green=0.2, blue=1.0, alpha=0.2), border = F)
  q_lower_p <- ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, p_lower)
  q_upper_p <- ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, p_upper)
  polygon(c(xp, rev(xp)), c(q_lower_p, rev(q_upper_p)), col = rgb(red=0.8, green=0.14, blue=0.15, alpha=0.2), border = F)
}
lines(isplineExtract(birds_raw$forest_obs)$x[,pred], ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, .5), 
      type="l", col = "black", lwd=1)
lines(isplineExtract(birds_raw$pasture_obs)$x[,pred], ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, .5),
      type="l", col = "black", lwd=1, lty=2)

# Mountain barrier
pred <- "matrix_1"
ymax <- 1.05*max(max(ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, .95)), max(ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, .95)))
plot(isplineExtract(birds_raw$forest_obs)$x[,pred], isplineExtract(birds_raw$forest_obs)$y[,pred],
     type="n", xlab="Mountain barrier (m below 4100)", ylab="Partial ecological distance",
     ylim = c(0,ymax))
x <- isplineExtract(birds_raw$forest_obs)$x[,pred]
xp <- isplineExtract(birds_raw$pasture_obs)$x[,pred]
for(i in 1:9){
  p_lower <- .5 - .05*i
  p_upper <- .5 + .05*i
  q_lower_f <- ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, p_lower)
  q_upper_f <- ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, p_upper)
  polygon(c(x, rev(x)), c(q_lower_f, rev(q_upper_f)), col = rgb(red=0.2, green=0.2, blue=1.0, alpha=0.2), border = F)
  q_lower_p <- ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, p_lower)
  q_upper_p <- ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, p_upper)
  polygon(c(xp, rev(xp)), c(q_lower_p, rev(q_upper_p)), col = rgb(red=0.8, green=0.14, blue=0.15, alpha=0.2), border = F)
}
lines(isplineExtract(birds_raw$forest_obs)$x[,pred], ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, .5), 
      type="l", col = "black", lwd=1)
lines(isplineExtract(birds_raw$pasture_obs)$x[,pred], ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, .5),
      type="l", col = "black", lwd=1, lty=2)

# Valley barrier
pred <- "matrix_2"
ymax <- 1.05*max(max(ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, .95)), max(ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, .95)))
plot(isplineExtract(birds_raw$forest_obs)$x[,pred], isplineExtract(birds_raw$forest_obs)$y[,pred],
     type="n", xlab="Valley barrier (m above 0)", ylab="Partial ecological distance",
     ylim = c(0,ymax))
x <- isplineExtract(birds_raw$forest_obs)$x[,pred]
xp <- isplineExtract(birds_raw$pasture_obs)$x[,pred]
for(i in 1:9){
  p_lower <- .5 - .05*i
  p_upper <- .5 + .05*i
  q_lower_f <- ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, p_lower)
  q_upper_f <- ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, p_upper)
  polygon(c(x, rev(x)), c(q_lower_f, rev(q_upper_f)), col = rgb(red=0.2, green=0.2, blue=1.0, alpha=0.2), border = F)
  q_lower_p <- ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, p_lower)
  q_upper_p <- ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, p_upper)
  polygon(c(xp, rev(xp)), c(q_lower_p, rev(q_upper_p)), col = rgb(red=0.8, green=0.14, blue=0.15, alpha=0.2), border = F)
}
lines(isplineExtract(birds_raw$forest_obs)$x[,pred], ispline_q_extract(birds_modeled$forest_gdm_rep_bb, pred, .5), 
      type="l", col = "black", lwd=1)
lines(isplineExtract(birds_raw$pasture_obs)$x[,pred], ispline_q_extract(birds_modeled$pasture_gdm_rep_bb, pred, .5),
      type="l", col = "black", lwd=1, lty=2)

