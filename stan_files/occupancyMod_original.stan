// Verbatim copy of model written in `occupancy simulation 2.R` file, save for 
// minor indentation tweaks. Note: simulation model has been tweaked so simulation
// exactly matches model. 

data {
    int<lower=1> n_point; //number of sites
    int<lower=1> n_visit; //fixed number of visits
    int<lower=1> n_species; //number of species
    int<lower=1> n_cluster; //number of clusters
    int<lower=0, upper=1> det_data[n_point, n_visit, n_species]; //detection history
    int<lower=0, upper=1> Q[n_point, n_species]; //at least one detection
    int<lower=0, upper=n_cluster> clusterID[n_point]; //cluster identifier
    vector[n_species] sp_cov1; //species covariate 1
    vector[n_species] sp_cov2; //species covariate 2
    vector[n_point] pt_cov1; //point covariate 1
    vector[n_point] pt_cov2; //point covariate 2
    matrix[n_point, n_visit] vis_cov1; //visit covariate 1
}
parameters {
    real mu_b0;
    real<lower=0> sigma_b0;
    vector[n_species] b0_raw;
    
    real<lower=0> sigma_b1;
    matrix[n_species, n_cluster] b1_raw;
    
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
    matrix[n_species, n_point] d1_raw;
    
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
    vector[n_species] b0 = mu_b0 + b0_raw * sigma_b0;
    matrix[n_species, n_cluster] b1 = b1_raw * sigma_b1;
    vector[n_species] b4 = mu_b4 + b4_raw * sigma_b4;
    vector[n_species] b5 = mu_b5 + b5_raw * sigma_b5;
    vector[n_species] b6 = mu_b6 + b6_raw * sigma_b6;
    
    vector[n_species] d0 = mu_d0 + d0_raw * sigma_d0;
    matrix[n_species, n_point] d1 = d1_raw * sigma_d1;
    vector[n_species] d4 = mu_d4 + d4_raw * sigma_d4;
    vector[n_species] d5 = mu_d5 + d5_raw * sigma_d5;
    real logit_psi[n_point, n_species];
    real logit_theta[n_point, n_visit, n_species];
    matrix[n_point, n_species] log_prob_increment;
    for(i in 1:n_point){
        for(k in 1:n_species){
            logit_psi[i,k] = b0[k] + b1[k, clusterID[i]] + b4[k]*pt_cov1[i]; 
        }
    }
    for(i in 1:n_point){
        for(j in 1:n_visit){
            for(k in 1:n_species){
                logit_theta[i,j,k] = d0[k] + d2*sp_cov1[k] + d5[k]*vis_cov1[i,j];
            }
        }
    }
    
    for(i in 1:n_point){
        for(k in 1:n_species){
            if(Q[i,k] == 1)
                log_prob_increment[i,k] = log_inv_logit(logit_psi[i,k]) + 
                    bernoulli_logit_lpmf(det_data[i,1,k] | logit_theta[i,1,k]) + 
                    bernoulli_logit_lpmf(det_data[i,2,k] | logit_theta[i,2,k]) +
                    bernoulli_logit_lpmf(det_data[i,3,k] | logit_theta[i,3,k]) +
                    bernoulli_logit_lpmf(det_data[i,4,k] | logit_theta[i,4,k]);
            else
                log_prob_increment[i,k] = log_sum_exp(log_inv_logit(logit_psi[i,k]) + 
                        log1m_inv_logit(logit_theta[i,1,k]) + 
                        log1m_inv_logit(logit_theta[i,2,k]) + 
                        log1m_inv_logit(logit_theta[i,3,k]) + 
                        log1m_inv_logit(logit_theta[i,4,k]), 
                    log1m_inv_logit(logit_psi[i,k]));
        }
    }
}
model {
    //Hyper-priors:
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
