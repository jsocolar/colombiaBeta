// This is a Stan model for the full Colombia bird dataset.
// 
//

functions{
    real partial_sum(

// numbers of random effect levels
        int n_sp,               // number of species
        int n_fam,              // number of families
        int n_sp_cl,            // number of species-clusters
        
// Model terms             
    // Occupancy
        // Spatial effects
        vector b0_sp_cl,       // intercepts by species-cluster
        
        // Taxonomic effects
        vector b0_sp,          // intercepts by species
        vector b0_fam,         // intercepts by family
        
        // Elevation effects             
        vector b1_relev_sp,    // slopes for elevation by species
        vector b1_relev2_sp,   // slopes for elevation^2 by species
        // add new terms for lowland species?
        
        // Pasture effects             
        vector b2_pasture_sp,  // slopes for pasture by species
        vector b2_pasture_fam, // slopes for pasture by family            

        // Trait effects
        real b3_eastOnly, 
        real b3_westOnly,
        real b3_snsmOnly,
        real b3_notWandes,
        real b3_notEandes,
        
        real b3_elevMedian,
        real b3_elevBreadth,
        
        real b3_forestPresent,         // universal slopes for habitat
        real b3_forestSpecialist,
        real b3_tfSpecialist,
        real b3_dryForestPresent,
        real b3_floodDrySpecialist,
        real b3_floodSpecialist,
        real b3_aridPresent,
        
        real b3_migratory,
        real b3_mass,
        
        real b3_dietInvert,              // universal slopes for diet (reference category == omnivore)
        real b3_dietCarn,
        real b3_dietFruitNect,
        real b3_dietGran,
        
        real b3_x_elevMedian_forestPresent,
        real b3_x_elevMedian_forestSpecialist,
        
        // pasture-x-trait interactions
        real b4_eastOnly, 
        real b4_westOnly,
        real b4_snsmOnly,
        real b4_notWandes,
        real b4_notEandes,
        
        real b4_elevMedian,
        real b4_elevBreadth,
        
        real b4_forestPresent,
        real b4_forestSpecialist,
        real b4_tfSpecialist,
        real b4_dryForestPresent,
        real b4_floodDrySpecialist,
        real b4_floodSpecialist,
        real b4_aridPresent,
        
        real b4_migratory,
        real b4_mass,
        
        real b4_dietInvert,
        real b4_dietCarn,
        real b4_dietFruitNect,
        real b4_dietGran,
        
        real b4_x_elevMedian_forestPresent,
        real b4_x_elevMedian_forestSpecialist,
                     
    // Detection
        // Taxonomic effects
        vector d0_sp,          // intercepts by species
        vector d0_fam,         // intercepts by family
        
        // Pasture effects
        vector d1_pasture_sp,  // slopes for pasture by species
        vector d1_pasture_fam, // slopes by family
        
        // Nuisance effects
        vector d2_time_sp,     // slopes for time-of-day by species
        real d2_obsSM,         // fixed slope for observer Simon
        real d2_obsDE,         // fixed slope for observer David
        real d2_obsJG,         // fixed slope for observer James
        
        // Trait effects
        real d3_mass,
        real d3_elevMedian,
        real d3_migratory,
                     
// Data
    // Integer IDs for random effect grouping terms
        int[] id_sp,           // species
        int[] id_fam,          // families
        int[] id_sp_cl,        // species-cluster
        
    // Q-vector
        int[] Q,               // 0/1 to indicate whether species-point has a detection
        
    // Number of visits
        int[] nv,              // number of visits to the given species-point (= number of visits to the point)
    
    // Elevations
        vector relev,    
        vector relev2, 

        
    // Pasture             
        vector pasture, 

    // Traits
        real[] eastOnly, 
        real[] westOnly,
        real[] snsmOnly,
        real[] notWandes,
        real[] notEandes,
        
        real[] elevMedian,
        real[] elevBreadth,
        
        real[] forestPresent, 
        real[] forestSpecialist,
        real[] tfSpecialist,
        real[] dryForestPresent,
        real[] floodDrySpecialist,
        real[] floodSpecialist,
        real[] aridPresent,
        
        real[] migratory,
        real[] mass,
        
        int[] dietInvert,
        int[] dietCarn,
        int[] dietFruitNect,
        int[] dietGran,
        
        row_vector[] time,     
        row_vector[] obsSM,
        row_vector[] obsDE,
        row_vector[] obsJG,


// Data slicing and indexing
    // Detection slice            
        int[,] det_slice,      // a slice of the detection array where rows are species-points and columns are visits
        
    // cutpoints for slicing                       
        int start,             // the starting row of the detection slice
        int end              // the ending row of the detection slice   
                     
    ){   // End of function inputs, beginning of computational part of function.
// indexing variables
        int len = end - start;
        int r0 = start - 1;
        
// variables for computation
        vector[len] lp;             // a vector to store the log-probability contribution from each row of det_slice
        real logit_psi;             // container for the value of logit_psi for a given row. Overwritten at each iteration of the loop.
        row_vector[4] logit_theta;      // container values of logit_theta for a given row. Overwritten at each iteration of the loop.
                                        // theta will always have length 4. For points with < 4 visits, the trailing elements of theta
                                        // will be nonsense, computed from trailing 0s on the visit covariates.
                                        // det_data will have trailing -99s, thus ensuring that none of these nonexistent visits
                                        // accidentally slips into analysis.
        
        for (r in 1:len) {
            // calculate psi & theta
            logit_psi = 
                b0_sp_cl[id_sp_cl[r0+r]] + 
                b0_sp[id_sp[r0+r]] + b0_fam[id_fam[r0+r]] +
                
                b1_relev_sp[id_sp[r0+r]]*relev[r0+r] + b1_relev2_sp[id_sp[r0+r]]*relev2[r0+r] +
                
                b2_pasture_sp[id_sp[r0+r]]*pasture[r0+r] + b2_pasture_fam[id_fam[r0+r]]*pasture[r0+r] +
                
                b3_eastOnly*eastOnly[r0+r] + b3_westOnly*westOnly[r0+r] + b3_snsmOnly*snsmOnly[r0+r] + b3_notWandes*notWandes[r0+r] + b3_notEandes*notEandes[r0+r] +
                b3_elevMedian*elevMedian[r0+r] + b3_elevBreadth*elevBreadth[r0+r] + 
                b3_forestPresent*forestPresent[r0+r] + b3_forestSpecialist*forestSpecialist[r0+r] + b3_tfSpecialist*tfSpecialist[r0+r] + 
                    b3_dryForestPresent*dryForestPresent[r0+r] + b3_floodDrySpecialist*floodDrySpecialist[r0+r] + b3_floodSpecialist*floodSpecialist[r0+r] +
                    b3_aridPresent*aridPresent[r0+r] +
                b3_migratory*migratory[r0+r] + b3_mass*mass[r0+r] + 
                b3_dietInvert*dietInvert[r0+r] + b3_dietCarn*dietCarn[r0+r] + b3_dietFruitNect*dietFruitNect[r0+r] + b3_dietGran*dietGran[r0+r] +
                b3_x_elevMedian_forestPresent*elevMedian[r0+r]*forestPresent[r0+r] + b3_x_elevMedian_forestSpecialist*elevMedian[r0+r]*forestSpecialist[r0+r] +
                
                
                b4_eastOnly*eastOnly[r0+r]*pasture[r0+r] + b4_westOnly*westOnly[r0+r]*pasture[r0+r] + b4_snsmOnly*snsmOnly[r0+r]*pasture[r0+r] + 
                    b4_notWandes*notWandes[r0+r]*pasture[r0+r] + b4_notEandes*notEandes[r0+r]*pasture[r0+r] +
                b4_elevMedian*elevMedian[r0+r]*pasture[r0+r] + b4_elevBreadth*elevBreadth[r0+r]*pasture[r0+r] + 
                b4_forestPresent*forestPresent[r0+r]*pasture[r0+r] + b4_forestSpecialist*forestSpecialist[r0+r]*pasture[r0+r] + 
                    b4_tfSpecialist*tfSpecialist[r0+r]*pasture[r0+r] + b4_dryForestPresent*dryForestPresent[r0+r]*pasture[r0+r] + 
                    b4_floodDrySpecialist*floodDrySpecialist[r0+r]*pasture[r0+r] + b4_floodSpecialist*floodSpecialist[r0+r]*pasture[r0+r] +
                    b4_aridPresent*aridPresent[r0+r]*pasture[r0+r] +
                b4_migratory*migratory[r0+r]*pasture[r0+r] + b4_mass*mass[r0+r]*pasture[r0+r] + 
                b4_dietInvert*dietInvert[r0+r]*pasture[r0+r] + b4_dietCarn*dietCarn[r0+r]*pasture[r0+r] + b4_dietFruitNect*dietFruitNect[r0+r]*pasture[r0+r] + 
                    b4_dietGran*dietGran[r0+r]*pasture[r0+r] +
                b4_x_elevMedian_forestPresent*elevMedian[r0+r]*forestPresent[r0+r]*pasture[r0+r] + 
                    b4_x_elevMedian_forestSpecialist*elevMedian[r0+r]*forestSpecialist[r0+r]*pasture[r0+r];
            
            // calculate logit_theta.  Note that this is vectorized because the time and observer covariates are row vectors
            // where nv < 4, logit_theta will have some trailing dummy terms.
            logit_theta = 
                d0_sp[id_sp[r0+r]] + d0_fam[id_fam[r0+r]] +
                d1_pasture_sp[id_sp[r0+r]]*pasture[r0+r] + d1_pasture_fam[id_fam[r0+r]]*pasture[r0+r] +
                d2_time_sp[id_sp[r0+r]]*time[r0+r] + d2_obsSM*obsSM[r0+r] + d2_obsDE*obsDE[r0+r] + d2_obsJG*obsJG[r0+r] +
                d3_mass*mass[r0+r] + d3_elevMedian*elevMedian[r0+r] + d3_migratory*migratory[r0+r];
            
            // likelihood
            if (Q[r0 + r] == 1) {
                lp[r] = log_inv_logit(logit_psi) +  // likelihood of occupancy
                    bernoulli_logit_lpmf(det_slice[r0 + r, 1:nv[r0 + r]] | logit_theta[1:nv[r0 + r]]);   // likelihood of observed detection history given occupancy
            } else {
                lp[r] = log_sum_exp(
                            log_inv_logit(logit_psi) + // likelihood of occupancy
                                sum(log1m_inv_logit(logit_theta[1:nv[r0 + r]])), // likelihood of all-zero detection history given occupancy
                            log1m_inv_logit(logit_psi)); // likelihood of non-occupancy
            }
        } 
    return sum(lp); // perhaps instead of storing lp as vector and summing it would be better to just do lp += ?
    } // end of function
} // Close the functions block

data {
    // grainsize for reduce_sum 
        int<lower=1> grainsize;  
    
    // dimensions
        int<lower=1> n_sp; // number of species
        int<lower=1> n_fam; // number of families
        int<lower=1> n_sp_cl; // number of unique species:cluster 
        int<lower=1> n_tot; // number of total rows (number of unique species:point)
        int<lower=1> n_visit_max; // maximum number of visits to a point
    
    // Detection array
        int<lower=0, upper=1> det_data[n_tot, n_visit_max]; // detection history
    
    // Q
        int<lower=0, upper=1> Q[n_tot];
    
    // number of visits
        int<lower=1, upper=n_visit_max> nv[n_tot];
    
    // Covariates
        // Random effect levels
        int<lower=1> id_sp[n_tot];
        int<lower=1> id_fam[n_tot];
        int<lower=1> id_sp_cl[n_tot];
        
        // Elevation
        vector[n_tot] relev;
        vector[n_tot] relev2;
        
        // Pasture
        vector[n_tot] pasture;
        
        // Traits
        int<lower=0, upper=1> eastOnly[n_tot];
        int<lower=0, upper=1> westOnly[n_tot];
        int<lower=0, upper=1> snsmOnly[n_tot];
        int<lower=0, upper=1> notWandes[n_tot];
        int<lower=0, upper=1> notEandes[n_tot];

        vector[n_tot] elevMedian;
        vector[n_tot] elevBreadth;
        
        int<lower=0, upper=1> forestPresent[n_tot];
        int<lower=0, upper=1> forestSpecialist[n_tot];
        int<lower=0, upper=1> tfSpecialist[n_tot];
        int<lower=0, upper=1> dryForestPresent[n_tot];
        int<lower=0, upper=1> floodDrySpecialist[n_tot];        
        int<lower=0, upper=1> floodSpecialist[n_tot];        
        int<lower=0, upper=1> aridPresent[n_tot];        

        int<lower=0, upper=1> migratory[n_tot];        
        vector[n_tot] mass;
        
        int<lower=0, upper=1> dietInvert[n_tot];
        int<lower=0, upper=1> dietCarn[n_tot];
        int<lower=0, upper=1> dietFruitNect[n_tot];
        int<lower=0, upper=1> dietGran[n_tot];
        
        row_vector[n_visit_max] time[n_tot];
        row_vector[n_visit_max] obsSM[n_tot];
        row_vector[n_visit_max] obsDE[n_tot];
        row_vector[n_visit_max] obsJG[n_tot];
}

parameters {
// Occupancy

    // Intercepts
        real mu_b0;
    
        real<lower=0> sigma_b0_sp_cl;
        vector[n_sp_cl] b0_sp_cl_raw;
    
        real<lower=0> sigma_b0_sp;
        vector[n_sp] b0_sp_raw;
    
        real<lower=0> sigma_b0_fam;
        vector[n_fam] b0_fam_raw;
    
    // Slopes
        // Elevation effects
            real mu_b1_relev;
            real<lower=0> sigma_b1_relev_sp;
            vector[n_sp] b1_relev_sp_raw;
            
            real<upper=0> mu_b1_relev2;   // constrain quadratic effect to be negative
            real<lower=0> sigma_b1_relev2_sp;
            vector[n_sp] b1_relev2_sp_raw;
    
        // Pasture effects
            real mu_b2_pasture;
            
            real<lower=0> sigma_b2_pasture_sp;
            vector[n_sp] b2_pasture_sp_raw;    
            
            real<lower=0> sigma_b2_pasture_fam;
            vector[n_fam] b2_pasture_fam_raw;    
        
        // Trait effects
            real b3_eastOnly;
            real b3_westOnly;
            real b3_snsmOnly;
            real b3_notWandes;
            real b3_notEandes;
        
            real b3_elevMedian;
            real b3_elevBreadth;
        
            real b3_forestPresent;
            real b3_forestSpecialist;
            real b3_tfSpecialist;
            real b3_dryForestPresent;
            real b3_floodDrySpecialist;
            real b3_floodSpecialist;
            real b3_aridPresent;
        
            real b3_migratory;
            real b3_mass;
        
            real b3_dietInvert;
            real b3_dietCarn;
            real b3_dietFruitNect;
            real b3_dietGran;
        
            real b3_x_elevMedian_forestPresent;
            real b3_x_elevMedian_forestSpecialist;
        
        // pasture-x-trait interactions
            real b4_eastOnly; 
            real b4_westOnly;
            real b4_snsmOnly;
            real b4_notWandes;
            real b4_notEandes;
        
            real b4_elevMedian;
            real b4_elevBreadth;
        
            real b4_forestPresent;
            real b4_forestSpecialist;
            real b4_tfSpecialist;
            real b4_dryForestPresent;
            real b4_floodDrySpecialist;
            real b4_floodSpecialist;
            real b4_aridPresent;
        
            real b4_migratory;
            real b4_mass;
        
            real b4_dietInvert;
            real b4_dietCarn;
            real b4_dietFruitNect;
            real b4_dietGran;
        
            real b4_x_elevMedian_forestPresent;
            real b4_x_elevMedian_forestSpecialist;

// Detection

    // Intercepts
        real mu_d0;
    
        real<lower=0> sigma_d0_sp;
        vector[n_sp] d0_sp_raw;
    
        real<lower=0> sigma_d0_fam;
        vector[n_fam] d0_fam_raw;
    
    // Slopes
        
        // Pasture effects
            real mu_d1_pasture;
            
            real<lower=0> sigma_d1_pasture_sp;
            vector[n_sp] d1_pasture_sp_raw;
    
            real<lower=0> sigma_d1_pasture_fam;
            vector[n_fam] d1_pasture_fam_raw;
        
        // Nuisance effects
            real mu_d2_time;
            real<lower=0> sigma_d2_time_sp;
            vector[n_sp] d2_time_sp_raw;
            

            real d2_obsSM;       
            real d2_obsDE;       
            real d2_obsJG;      
        
        // Trait effects
            real d3_mass;
            real d3_elevMedian;
            real d3_migratory;
}
transformed parameters{
    // occupancy
    vector[n_sp] b0_sp = mu_b0 + b0_sp_raw * sigma_b0_sp;
    vector[n_fam] b0_fam = b0_fam_raw * sigma_b0_fam;
    vector[n_sp_cl] b0_sp_cl = b0_sp_cl_raw * sigma_b0_sp_cl;
    
    vector[n_sp] b1_relev_sp = mu_b1_relev + b1_relev_sp_raw * sigma_b1_relev_sp;
    vector[n_sp] b1_relev2_sp = mu_b1_relev2 + b1_relev2_sp_raw * sigma_b1_relev2_sp;

    vector[n_sp] b2_pasture_sp = mu_b2_pasture + b2_pasture_sp_raw * sigma_b2_pasture_sp;
    vector[n_fam] b2_pasture_fam = b2_pasture_fam_raw * sigma_b2_pasture_fam;

    // detection
    vector[n_sp] d0_sp = mu_d0 + d0_sp_raw * sigma_d0_sp;
    vector[n_fam] d0_fam = d0_fam_raw * sigma_d0_fam;

    vector[n_sp] d1_pasture_sp = mu_d1_pasture + d1_pasture_sp_raw * sigma_d1_pasture_sp;
    vector[n_fam] d1_pasture_fam = d1_pasture_fam_raw * sigma_d1_pasture_fam;

    vector[n_sp] d2_time_sp = mu_d2_time + d2_time_sp_raw * sigma_d2_time_sp;
}

model {
    //Likelihood
    target += reduce_sum(
        // partial_sum function
        partial_sum, 
        
        // variable sizes
        n_sp, n_fam, n_sp_cl, b0_sp_cl, 
        
        // parameters
        b0_sp, b0_fam, b1_relev_sp, b1_relev2_sp, b2_pasture_sp, b2_pasture_fam,
        b3_eastOnly, b3_westOnly, b3_snsmOnly, b3_notWandes, b3_notEandes, b3_elevMedian, b3_elevBreadth,
        b3_forestPresent, b3_forestSpecialist, b3_tfSpecialist, b3_dryForestPresent, b3_floodDrySpecialist,
        b3_floodSpecialist, b3_aridPresent, b3_migratory, b3_mass, b3_dietInvert, b3_dietCarn, b3_dietFruitNect,
        b3_dietGran, b3_x_elevMedian_forestPresent, b3_x_elevMedian_forestSpecialist, b4_eastOnly, b4_westOnly,
        b4_snsmOnly, b4_notWandes, b4_notEandes,b4_elevMedian, b4_elevBreadth, b4_forestPresent, b4_forestSpecialist,
        b4_tfSpecialist, b4_dryForestPresent, b4_floodDrySpecialist, b4_floodSpecialist, b4_aridPresent, b4_migratory,
        b4_mass, b4_dietInvert, b4_dietCarn, b4_dietFruitNect, b4_dietGran, b4_x_elevMedian_forestPresent,
        b4_x_elevMedian_forestSpecialist,
        d0_sp, d0_fam, d1_pasture_sp, d1_pasture_fam, d2_time_sp, d2_obsSM, d2_obsDE, d2_obsJG, d3_mass, d3_elevMedian,
        d3_migratory,
        
        // Data
        id_sp, id_fam, id_sp_cl, Q, nv, relev, relev2, pasture, eastOnly, westOnly, snsmOnly, notWandes, notEandes,
        elevMedian, elevBreadth, forestPresent, forestSpecialist, tfSpecialist, dryForestPresent, floodDrySpecialist,
        floodSpecialist, aridPresent, migratory, mass, dietInvert, dietCarn, dietFruitNect, dietGran, time, obsSM,
        obsDE, obsJG,
        
        det_data);


    // Priors
    mu_b0 ~ student_t(7.763, 0, 1.566);  // This is Dorazio's suggested prior, which is approximately uniform on the probability scale between 0.1 and 0.99. See also Northrup & Gerber 2018
    
    sigma_b0_sp_cl ~ normal(0, 2);
    b0_sp_cl_raw ~ normal(0, 1);
    
    sigma_b0_sp ~ normal(0, 2);
    b0_sp_raw ~ normal(0, 1);
    
    sigma_b0_fam ~ normal(0, 2);
    b0_fam_raw ~ normal(0, 1);
    
    mu_b1_relev ~ normal(0, 1);   // meant to be somewhat informative, as overall elevation relationships should be pretty flat
    sigma_b1_relev_sp ~ normal(0, 2);
    b1_relev_sp_raw ~ normal(0, 1);
    
    mu_b1_relev2 ~ normal(0, 3);   // this is a half-normal prior as mu_b1_relev2 is constrained to be negative in the parameter declaration
    sigma_b1_relev2_sp ~ normal(0, 2);
    b1_relev_sp_raw ~ normal(0, 1);

    mu_b2_pasture ~ normal(0, 2);
    
    sigma_b2_pasture_sp ~ normal(0, 2);
    b2_pasture_sp_raw ~ normal(0, 1);
    
    sigma_b2_pasture_fam ~ normal(0, 2);
    b2_pasture_fam_raw ~ normal(0, 1);
    
    b3_eastOnly ~ normal(0, 2);
    b3_westOnly ~ normal(0, 2);
    b3_snsmOnly ~ normal(0, 2);
    b3_notWandes ~ normal(0, 2);
    b3_notEandes ~ normal(0, 2);
    
    b3_elevMedian ~ normal(0, 2);
    b3_elevBreadth ~ normal(0, 2);
    b3_forestPresent ~ normal(0, 2);
    b3_forestSpecialist ~ normal(0, 2);
    b3_tfSpecialist ~ normal(0, 2);
    b3_dryForestPresent ~ normal(0, 2);
    b3_floodDrySpecialist ~ normal(0, 2);
    b3_floodSpecialist ~ normal(0, 2);
    b3_aridPresent ~ normal(0, 2);
        
    b3_migratory ~ normal(0, 2);
    b3_mass ~ normal(0, 2);
    
    b3_dietInvert ~ normal(0, 2);
    b3_dietCarn ~ normal(0, 2);
    b3_dietFruitNect ~ normal(0, 2);
    b3_dietGran ~ normal(0, 2);
    
    b3_x_elevMedian_forestPresent ~ normal(0, 2);
    b3_x_elevMedian_forestSpecialist ~ normal(0, 2);
    
    
    b4_eastOnly ~ normal(0, 2);
    b4_westOnly ~ normal(0, 2);
    b4_snsmOnly ~ normal(0, 2);
    b4_notWandes ~ normal(0, 2);
    b4_notEandes ~ normal(0, 2);
    
    b4_elevMedian ~ normal(0, 2);
    b4_elevBreadth ~ normal(0, 2);
    b4_forestPresent ~ normal(0, 2);
    b4_forestSpecialist ~ normal(0, 2);
    b4_tfSpecialist ~ normal(0, 2);
    b4_dryForestPresent ~ normal(0, 2);
    b4_floodDrySpecialist ~ normal(0, 2);
    b4_floodSpecialist ~ normal(0, 2);
    b4_aridPresent ~ normal(0, 2);
        
    b4_migratory ~ normal(0, 2);
    b4_mass ~ normal(0, 2);
    
    b4_dietInvert ~ normal(0, 2);
    b4_dietCarn ~ normal(0, 2);
    b4_dietFruitNect ~ normal(0, 2);
    b4_dietGran ~ normal(0, 2);
    
    b4_x_elevMedian_forestPresent ~ normal(0, 2);
    b4_x_elevMedian_forestSpecialist ~ normal(0, 2);

// Detection
    mu_d0 ~ student_t(7.763, 0, 1.566);  // This is Dorazio's suggested prior, which is approximately uniform on the probability scale between 0.1 and 0.99. See also Northrup & Gerber 2018
    
    sigma_d0_sp ~ normal(0, 2);
    d0_sp_raw ~ normal(0, 1);
    
    sigma_d0_fam ~ normal(0, 2);
    d0_fam_raw ~ normal(0, 1);

    mu_d1_pasture ~ normal(0, 2);
    
    sigma_d1_pasture_sp ~ normal(0, 2);
    d1_pasture_sp_raw ~ normal(0, 1);
    
    sigma_d1_pasture_fam ~ normal(0, 2);
    d1_pasture_fam_raw ~ normal(0, 2);
    
    mu_d2_time ~ normal(0, 2);
    sigma_d2_time_sp ~ normal(0, 2);
    d2_time_sp_raw ~ normal(0, 1);

    d2_obs_SM ~ normal(0, 1); // intentionally somewhat informative
    d2_obs_DE ~ normal(0, 1);
    d2_obs_JG ~ normal(0, 1);        

    d3_mass ~ normal(0, 2);
    d3_elevMedian ~ normal(0, 2);
    d3_migratory ~ normal(0, 2);
}
