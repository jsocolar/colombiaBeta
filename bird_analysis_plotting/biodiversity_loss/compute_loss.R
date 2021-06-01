# Functions to compute the loss of conservation value over all cells (get_avg_cell_logratios)
# or over an arbitrary collection of cells (get_regional_logratios)

# A helper function
avg_logratio <- function(fp_pp_i, cutoff_use){
  n_col <- length(fp_pp_i)/2
  if(n_col != round(n_col)){stop("dimension error in fp_pp_i: expected even number of columns")}
  fp_i <- fp_pp_i[1:n_col]
  pp_i <- fp_pp_i[(n_col+1):(2*n_col)]
  incl <- ((fp_i) >= cutoff_use) | ((pp_i) >= cutoff_use)
  fp_inc <- fp_i[incl]
  pp_inc <- pp_i[incl]
  a_l <- mean(log(fp_inc/pp_inc))
  m_l <- median(log(fp_inc/pp_inc))
  n <- sum(incl)
  return(list(avg_logratio = a_l, med_logratio = m_l, n = n))
}

get_avg_cell_logratios <- function(forest_probs, pasture_probs, cutoff_type, cutoff){
  # forest_probs, pasture_probs: a dataframe with columns for cell_id, x, y, and then one column of psi for each species
  # cutoff_type: either "absolute" or "relative"
  # cutoff for including a species in a cell.  
  #   If cutoff_type = "absolute" then a species must be present in EITHER forest OR pasture with probability >= cutoff. 
  #   If cutoff_type = "relative" then a species must be present in EITHER forest OR pasture with probability >=
  #     cutoff * maxprob, where maxprob is the maximum probability in EITHER forest OR pasture attained anywhere in Colombia.
  if(!all.equal(dim(forest_probs), dim(pasture_probs))){stop("forest_probs and pasture_probs have different dimensions")}
  if(cutoff <= 0){stop("cutoff must be greater than zero and less than one")}
  if(cutoff >= 1){stop("cutoff must be greater than zero and less than one")}
  fp <- forest_probs[, 5:ncol(forest_probs)]
  pp <- pasture_probs[, 5:ncol(pasture_probs)]
  fp_pp <- cbind(fp, pp)
  
  if(cutoff_type == "relative"){
    ap <- rbind(fp, pp)
    maxprobs <- apply(ap, 2, max)
    cutoff_use <- cutoff*maxprobs
  }else if(cutoff_type == "absolute"){
    cutoff_use <- rep(cutoff, ncol(fp))
  }else{stop("cutoff type must be one of 'relative' or 'absolute'")}
  
  aln <- apply(fp_pp, 1, avg_logratio, cutoff_use = cutoff_use)
  
  output <- as.data.frame(do.call(rbind, aln))
  output$avg_logratio <- as.numeric(output$avg_logratio)
  output$med_logratio <- as.numeric(output$med_logratio)
  output$n <- as.integer(output$n)
  return(output)
}


get_regional_logratios <- function(forest_probs, pasture_probs, cutoff_type, cutoff, cell_positions = NULL){
  if(!all.equal(dim(forest_probs), dim(pasture_probs))){stop("forest_probs and pasture_probs have different dimensions")}
  if(cutoff <= 0){stop("cutoff must be greater than zero and less than one")}
  if(cutoff >= 1){stop("cutoff must be greater than zero and less than one")}
  
  fp <- forest_probs[, 5:ncol(forest_probs)]
  pp <- pasture_probs[, 5:ncol(pasture_probs)]
  if(!is.null(cell_positions)){
    fp <- fp[cell_positions, ]
    pp <- pp[cell_positions, ]
  }
  
  if(cutoff_type == "relative"){
    ap <- rbind(fp, pp)
    maxprobs <- apply(ap, 2, max)
    cutoff_use <- cutoff*maxprobs
  }else if(cutoff_type == "absolute"){
    cutoff_use <- rep(cutoff, ncol(fp))
  }else{stop("cutoff type must be one of 'relative' or 'absolute'")}
  
  fp_max <- apply(fp, 2, max)
  pp_max <- apply(pp, 2, max)
  incl <- ((fp_max) >= cutoff_use) | ((pp_max) >= cutoff_use)
  
  fp_occ <- colSums(fp)[incl]
  pp_occ <- colSums(pp)[incl]
  
  avg_logratio <- mean(log(fp_occ/pp_occ))
  med_logratio <- median(log(fp_occ/pp_occ))
  return(list(avg_logratio=avg_logratio, med_logratio=med_logratio, n = sum(incl)))
}


get_sample_percent_decline <- function(forest_samples, pasture_samples, cutoff, cell_positions = NULL) {
  if(!all.equal(dim(forest_samples), dim(pasture_samples))){stop("forest_samples and pasture_samples have different dimensions")}
  if(cutoff < 1){stop("cutoff must be at least one")}
  fs <- forest_samples[, 5:ncol(forest_samples)]
  ps <- pasture_samples[, 5:ncol(pasture_samples)]
  if(!is.null(cell_positions)){
    fs <- fs[cell_positions, ]
    ps <- ps[cell_positions, ]
  }
  f_plus_p <- fs + ps
  fs <- fs[, colSums(f_plus_p) >= cutoff]
  ps <- ps[, colSums(f_plus_p) >= cutoff]
  ratios_pointwise <- fs/ps
  decline_frac_pointwise <- apply(ratios_pointwise, 1, function(x){sum(x > 1, na.rm = T)/sum(!is.na(x))})
  ratios_total <- colSums(fs)/colSums(ps)
  nsp <- ncol(fs)
  decline_frac_total <- sum(ratios_total > 1)/nsp
  return(list(pointwise = decline_frac_pointwise, total = decline_frac_total, nsp = nsp))
}
