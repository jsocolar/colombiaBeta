// 
//
// To do: write a function that accepts varying number of site visits

functions{
    real partial_sum(int[,] det_slice, 
                     int start, int end, 
                     vector b0, 
                     vector b1, 
                     vector b4,
                     vector d0, 
                     real d2, 
                     vector d5,
                     vector pt_cov1, 
                     vector sp_cov1,
                     row_vector[] vis_cov1,
                     int[] id_sp, 
                     int[] id_sp_cl, 
                     int[] Q) {
        // indexing vars
        int len = end - start;
        int r0 = start - 1;
        
        vector[len] lp;
        vector[len] logit_psi;
        row_vector[4] logit_theta[len];
        for (r in 1:len) {
            // calculate psi & theta
            logit_psi[r] = b0[id_sp[r0+r]] + b1[id_sp_cl[r0+r]] + b4[id_sp[r0+r]]*pt_cov1[r0+r];
            logit_theta[r] = d0[id_sp[r0+r]] + d2*sp_cov1[r0+r] + d5[id_sp[r0+r]]*vis_cov1[r0+r];
            // likelihood
            if (Q[r0 + r] == 1) 
                lp[r] = log_inv_logit(logit_psi[r]) +
                    bernoulli_logit_lpmf(det_slice[r0 + r] | logit_theta[r]);
            else lp[r] = log_sum_exp(
                log_inv_logit(logit_psi[r]) +
                    log1m_inv_logit(logit_theta[r, 1]) +
                    log1m_inv_logit(logit_theta[r, 2]) +
                    log1m_inv_logit(logit_theta[r, 3]) +
                    log1m_inv_logit(logit_theta[r, 4]),
                log1m_inv_logit(logit_psi[r]));
        } 
        return sum(lp);
    }
}
data {
    // dimensions
    int<lower=1> n_visit; //fixed number of visits
    int<lower=1> n_species; //number of species
    int<lower=1> n_sp_cl; //number of unique species:cluster 
    int<lower=1> n_sp_pt; //number of unique species:point
    int<lower=1> n_pt; //number of unique points
    int<lower=1> n_visit_max; // currently redundant, but will set function later
    int<lower=1> n_tot; // nrows in df
    int<lower=1> n_1; // number of Q=1 points, 
    
    // indexing variables
    int<lower=1> id_pt[n_tot];
    int<lower=1> id_cl[n_tot];
    int<lower=1> id_sp[n_tot]; 
    int<lower=1> id_sp_cl[n_tot]; // species:cluster
    int<lower=0, upper=1> Q[n_tot]; // species:cluster
    
    // data & covariates
    vector[n_tot] sp_cov1; //species covariate 1
    vector[n_tot] sp_cov2; //species covariate 2
    vector[n_tot] pt_cov1; //point covariate 1
    vector[n_tot] pt_cov2; //point covariate 2
    row_vector[n_visit_max] vis_cov1[n_tot]; //visit covariate 1
    int<lower=0, upper=1> det_data[n_tot, n_visit_max]; // detection history
    
    // grainsize for reduce_sum 
    int<lower=1> grainsize;
}
parameters {
    real mu_b0;
    real<lower=0> sigma_b0;
    vector[n_species] b0_raw;
    
    real<lower=0> sigma_b1;
    vector[n_sp_cl] b1_raw;
    
    real b2;
    real b3;
    
    real mu_b4;
    real<lower=0> sigma_b4;
    vector[n_species] b4_raw;
    
    real mu_b5;
    real<lower=0> sigma_b5;
    vector[n_species] b5_raw;
    
    real mu_b6;
    real<lower=0> sigma_b6;
    vector[n_species] b6_raw;
    
    real mu_d0;
    real<lower=0> sigma_d0;
    vector[n_species] d0_raw;
    
    real<lower=0> sigma_d1;
    vector[n_sp_pt] d1_raw; 
    
    real d2;
    real d3;
    real mu_d4;
    real<lower=0> sigma_d4;
    vector[n_species] d4_raw;
    
    real mu_d5;
    real<lower=0> sigma_d5;
    vector[n_species] d5_raw;
}
transformed parameters{
    // occupancy
    vector[n_species] b0 = mu_b0 + b0_raw * sigma_b0;
    vector[n_sp_cl] b1 = b1_raw * sigma_b1;
    vector[n_species] b4 = mu_b4 + b4_raw * sigma_b4;
    vector[n_species] b5 = mu_b5 + b5_raw * sigma_b5;
    vector[n_species] b6 = mu_b6 + b6_raw * sigma_b6;
    
    // detection
    vector[n_species] d0 = mu_d0 + d0_raw * sigma_d0;
    vector[n_sp_pt] d1 = d1_raw * sigma_d1;
    vector[n_species] d4 = mu_d4 + d4_raw * sigma_d4;
    vector[n_species] d5 = mu_d5 + d5_raw * sigma_d5;
}
model {
    //Likelihood
    target += partial_sum(det_data, 1, grainsize, b0, b1, b4, d0, d2, d5, pt_cov1, sp_cov1, 
        vis_cov1, id_sp, id_sp_cl, Q);
    target += partial_sum(det_data, grainsize+1, n_tot, b0, b1, b4, d0, d2, d5, pt_cov1, 
        sp_cov1, vis_cov1, id_sp, id_sp_cl, Q);
    // Hyper-priors:
    mu_b0 ~ normal(0,10);
    b2 ~ normal(0,10);
    b3 ~ normal(0,10);
    mu_b4 ~ normal(0,10);
    mu_b5 ~ normal(0,10);
    mu_b6 ~ normal(0,10);
    
    mu_d0 ~ normal(0,10);
    d2 ~ normal(0,10);
    d3 ~ normal(0,10);
    mu_d4 ~ normal(0,10);
    mu_d5 ~ normal(0,10);
    
    sigma_b0 ~ normal(0,10);
    sigma_b1 ~ normal(0,10);
    sigma_b4 ~ normal(0,10);
    sigma_b5 ~ normal(0,10);
    sigma_b6 ~ normal(0,10);
    
    sigma_d0 ~ normal(0,10);
    sigma_d1 ~ normal(0,10);
    sigma_d4 ~ normal(0,10);
    sigma_d5 ~ normal(0,10);
    
    //Random Effects
    b0_raw ~ normal(0, 1);
    to_vector(b1_raw) ~ normal(0, 1);
    b4_raw ~ normal(0, 1);
    b5_raw ~ normal(0, 1);
    b6_raw ~ normal(0, 1);
    
    d0_raw ~ normal(0, 1);
    to_vector(d1_raw) ~ normal(0, 1);
    d4_raw ~ normal(0, 1);
    d5_raw ~ normal(0, 1);
}
