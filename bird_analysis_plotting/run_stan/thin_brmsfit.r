fm3 <- readRDS("/home/ec2-user/outputs/fm3.RDS")

thin_brms <- function(brmsfit, thin) {
    brmsfit <- fm3
  sim <- brmsfit$fit@sim 
  for(chain in 1:sim$chains) {
    sim$samples[[chain]] <- sim$samples[[chain]][seq(sim$warmup+1, sim$iter, thin), ]
    attributes(sim$samples[[chain]])$sampler_params <-
      attributes(sim$samples[[chain]])$sampler_params[seq(sim$warmup+1, sim$iter, thin), ]
  }
  
  # Update the meta-info
   sim$warmup2 <- rep(0, sim$chains) 
   sim$n_save <- length(seq(sim$warmup+1, sim$iter, thin)) |>
     rep(4)
   sim$thin <- thin
  
#   stan_args <- brmsfit$fit@stan_args[[1]]
#   stan_args$warmup <- 0
#   stan_args$iter <- sim$iter
#   stan_args$thin <- thin
  
  brmsfit$fit@sim <- sim
#  brmsfit$fit@stan_args <- list(stan_args)
  
  brmsfit
}

fm3_thin_300 <- thin_brms(fm3, 300)
summary(fm3_thin_300)

fm3_thin_100 <- thin_brms(fm3, 100)
summary(fm3_thin_100)

fm3_thin_10 <- thin_brms(fm3, 10)
summary(fm3_thin_10)

saveRDS(fm3_thin_300, "/home/ec2-user/outputs/fm3_thin_300.RDS")
saveRDS(fm3_thin_100, "/home/ec2-user/outputs/fm3_thin_100.RDS")
saveRDS(fm3_thin_10, "/home/ec2-user/outputs/fm3_thin_10.RDS")
