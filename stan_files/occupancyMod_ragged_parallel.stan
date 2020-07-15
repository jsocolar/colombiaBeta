// 
//
// To do: write a function that accepts varying number of site visits

functions{
    real partial_sum(int ns,               // number of species
                    int ng,                // number of genera
                    int np,                // number of points
                    int nc,                // number of clusters
                    int nsc,               // number of species-clusters
                    
                    int[] id_sp,           // integer IDs for species
                    int[] id_gen,          // integer IDs for genera
                    int[] id_sp_cl,        // integer IDs for species-by-cluster
                    int[] Q,               // Q-vector (0/1 to indicate whether species-point has a detection)
                    int[] nv,              // number of visits to the given species-point (= number of visits to the point)
                    
                    
                    
                    
                    
                    int[,] det_slice,      // detection slice: a slice of the detection array, 
                                                // where rows are species-points and columns are 
                                                // visits
                     int start,             // the starting row of the detection slice
                     int end,               // the ending row of the detection slice
                     
                     // occupancy terms
                     vector b0_cs,          // vector of intercepts by cluster-species
                     vector b0_s,           // vector of intercepts by species
                     vector b0_g,           // vector of intercepts by genus
                     vector b1_l_s,         // vector of slopes for elevation (linear) by species
                     vector b1_q_s,         // vector of slopes for elevation (squared) by species
                     vector b2_s,           // vector of slopes for pasture by species
                     vector b2_g,           // vector of slopes for pasture by genus
                     
                     // detection terms
                     vector d0_s,           // vector of intercepts by species
                     vector d0_g,           // vector of intercepts by genus
                     vector d1_s,           // vector of slopes for pasture by species
                     vector d2_s,           // vector of slopes for time-of-day by species
                     
                     row_vector[nsp] scaled_elev[npt],   // array of scaled elevations: rows are species, columns are points.
                     row_vector[]
                     vector pt_cov1,        // values of point-associated covariate 1
                     vector sp_cov1,        // values of species-associated covariate 1
                     row_vector[] vis_cov1 // values of visit-associated covariate 1
                     ) {            // number of visits in a given row
                     
        // indexing variables
        int len = end - start;
        int r0 = start - 1;
        
        // variables for computation
        vector[len] lp;                     // a vector to store the log-probability contribution from 
                                                // each row of det_slice
        real logit_psi;                     // a place to hold the value of logit_psi for a given row.
                                                // not a vector; overwritten at each iteration of the loop.
        vector[4] logit_theta;              // a place to hold values of logit_theta for a given row.
                                                // overwritten at each iteration of the loop.
        real theta_term;                    // a place to compute the sum of log1m_inv_logit(logit_theta) for 
                                                // a given row. Overwritten at each iteration of the loop.
                                                
        // The overall idea here is that we will let theta always be of length 4, with trailing zeros on det_data
        // and the visit covariates for any row with < 4 visits.  Then, below, we use only the relevant terms of 
        // theta whenever nv < 4.
        for (r in 1:len) {
            // calculate psi & theta
            logit_psi = b0[id_sp[r0+r]] + b1[id_sp_cl[r0+r]] + b4[id_sp[r0+r]]*pt_cov1[r0+r];
            logit_theta = d0[id_sp[r0+r]] + d2*sp_cov1[r0+r] + d5[id_sp[r0+r]]*vis_cov1[r0+r];  
            // Where nv < 4, logit_theta will have some trailing dummy terms.
            
            // likelihood
            if (Q[r0 + r] == 1) {
                lp[r] = log_inv_logit(logit_psi) +
                    bernoulli_logit_lpmf(det_slice[r0 + r, 1:nv[r0 + r]] | logit_theta[1:nv[r0 + r]]);
                // First term is probability of occupancy; second term is probabily of observed detection
                // data given occupancy.
            } else {
                theta_term = 0;
                for(v in 1:nv[r0 + r]) {
                    theta_term += log1m_inv_logit(logit_theta[v]);
                }
                // theta_term is the sum of log1m_inv_logit(logit_theta) for all the relevant terms of theta
                lp[r] = log_sum_exp(log_inv_logit(logit_psi) + theta_term, log1m_inv_logit(logit_psi[r]));
            }
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
    target += reduce_sum(partial_sum, det_data, grainsize, b0, b1, b4, d0, d2, d5, pt_cov1, sp_cov1, 
        vis_cov1, id_sp, id_sp_cl, Q);

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
