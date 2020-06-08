// Reformulated model to accept ragged datastructure (varying number of points 
// within species). As written, data has to be sorted by Q, 1 first, then 0s, 
// n_1 indexes the length of 1s. This might be redundant, or might marginally
// increase speed. 
//
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
    int<lower =1> id_sp[n_tot]; 
    int<lower =1> id_sp_cl[n_tot]; // species:cluster
    
    // data & covariates
    // note: det_data has to be passed sorted by Q
    vector[n_tot] sp_cov1; //species covariate 1
    vector[n_tot] sp_cov2; //species covariate 2
    vector[n_tot] pt_cov1; //point covariate 1
    vector[n_tot] pt_cov2; //point covariate 2
    //matrix[n_tot, n_visit_max] vis_cov1; //visit covariate 1
    row_vector[n_visit_max] vis_cov1[n_tot]; //visit covariate 1
    int<lower=0, upper=1> det_data[n_tot, n_visit_max]; // detection history
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
    // psi & theta
    real logit_psi[n_tot];
    //matrix[n_tot, n_visit_max] logit_theta;
    row_vector[n_visit_max] logit_theta[n_tot];
    vector[n_tot] log_prob_increment;
    // fill psi and theta
    for (r in 1:n_tot) {
        logit_psi[r] = b0[id_sp[r]] + b1[id_sp_cl[r]] + b4[id_sp[r]]*pt_cov1[r];
        logit_theta[r] = d0[id_sp[r]] + d2*sp_cov1[r] + d5[id_sp[r]]*vis_cov1[r];
    }
    
    // need to write a function that will do psi * det1 * det2 ... on a vector. 
    // Leaving in loop for now
    for(r in 1:n_1) {
       log_prob_increment[r] = log_inv_logit(logit_psi[r]) +
            bernoulli_logit_lpmf(det_data[r] | logit_theta[r]);
    }
     for(r2 in (n_1+1):n_tot) {
        log_prob_increment[r2] = log_sum_exp(
            log_inv_logit(logit_psi[r2]) +
                log1m_inv_logit(logit_theta[r2, 1]) +
                log1m_inv_logit(logit_theta[r2, 2]) +
                log1m_inv_logit(logit_theta[r2, 3]) +
                log1m_inv_logit(logit_theta[r2, 4]),
            log1m_inv_logit(logit_psi[r2]));
    }
}
model {
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
    
    //Likelihood (data level)
    target += sum(log_prob_increment);
}
