// This is a Stan model for the full Colombia bird dataset, version 7.0
// Changes:   Implementing the redundancy trick

functions{
  real partial_sum(
    // Function arguments:
      // Data slicing and indexing
        // Detection slice            
          int[,] det_slice,      // slice of detection array. 
        // cutpoints for slicing                       
          int start,             // the starting row of the detection slice
          int end,              // the ending row of the detection slice   
      // numbers of random effect levels
        int n_spCl,             // number of species-clusters
        int n_spSr,             // number of species-subregions
        int n_sp,               // number of species
        int n_fam,              // number of families
        int n_spObs,            // number of species-observers
      // Parameters             
        // Occupancy
          // Intercept
            real mu_b0,
            // Spatial effects (zero-centered)
              vector b0_spCl,       // intercepts by species-cluster
              vector b0_spSr,       // intercepts by species-subregion
            // Taxonomic effects (zero-centered)
              vector b0_sp,          // intercepts by species
              vector b0_fam,         // intercepts by family
          // Slopes
            // Elevation effects
              real mu_b1_relev, 
              vector b1_relev_sp,    // slopes by species (zero-centered)
              real mu_b1_relev2,
              vector b1_relev2_sp,   // slopes by species (zero-centered)
              vector b1_lowland_unique,          // modified logit-quadratic lowland spp: intercept
              vector b1_x_lowland_relev_unique,  //    linear term
              vector b1_x_lowland_relev2_unique, //    quadratic term
          // Pasture effects
            real mu_b2_pasture,
            // Taxonomic effects (zero-centered)
            vector b2_pasture_sp,  // slopes by species
            vector b2_pasture_fam, // slopes by family 
          // Trait effects; covariates vary by species and all slopes are universal
            // Biogeography
              // Barriers
                vector b3_mountain_barrier_unique,
                vector b3_valley_barrier_unique,
              // Elevations
                vector b3_elevMedian_unique,
                vector b3_elevBreadth_unique,
              // Habitats
                vector b3_forestPresent_unique,
                vector b3_forestSpecialist_unique,
                vector b3_tfSpecialist_unique,
                vector b3_dryForestPresent_unique,
                vector b3_floodDrySpecialist_unique,
                vector b3_aridPresent_unique,
              // Migratory
                vector b3_migratory_unique,
              // Interactions
                real b3_x_elevMedian_forestPresent,
                real b3_x_elevMedian_forestSpecialist,
            // Functional traits
              vector b3_mass_unique,
              // diet (effects-coded; reference category = omnivore)
                vector b3_dietInvert_unique,
                vector b3_dietCarn_unique,
                vector b3_dietFruitNect_unique,
                vector b3_dietGran_unique,
          // pasture-x-trait interactions (same traits as b3).
            // Biogeography
              // Barriers
                vector b4_mountain_barrier_unique, 
                vector b4_valley_barrier_unique,
              // Elevations
                real b4_elevMedian,
                real b4_elevBreadth,
              // Habitats
                vector b4_forestPresent_unique,
                vector b4_forestSpecialist_unique,
                vector b4_tfSpecialist_unique,
                vector b4_dryForestPresent_unique,
                vector b4_floodDrySpecialist_unique,
                vector b4_aridPresent_unique,
              // Migratory
                vector b4_migratory_unique,
              // Interactions
                real b4_x_elevMedian_forestPresent,
                real b4_x_elevMedian_forestSpecialist,
            // Functional traits
              real b4_mass,
              // diet (effects coded; reference category = omnivore)
                vector b4_dietInvert_unique,
                vector b4_dietCarn_unique,
                vector b4_dietFruitNect_unique,
                vector b4_dietGran_unique,
          // Range maps
            real mu_b5_distance_to_range,
            vector b5_distance_to_range_sp, // species effects (zero-centered)
        // Detection
          // Intercept
            real mu_d0,
            // Taxonomic effects (zero-centered)
              vector d0_sp,          // species
              vector d0_fam,         // family
            // Species-observer effect (zero-centered)
              row_vector d0_spObs,
          // Slopes
            // Pasture effects
              real mu_d1_pasture,
              // Taxonomic effects (zero-centered)
                vector d1_pasture_sp,  // species
                vector d1_pasture_fam, // family
            // Nuisance effects
              // Time
                real mu_d2_time,
                vector d2_time_sp,     // species effects (zero-centered)
              // Observer
                real d2_obsSM,         // fixed slope for observer Mills
                real d2_obsDE,         // fixed slope for observer Edwards
                real d2_obsJG,         // fixed slope for observer Gilroy
              // Trait effects (universal slopes)
                vector d3_mass_unique,
                vector d3_elevMedian_unique,
                vector d3_migratory_unique,
                vector d3_dietCarn_unique,
                real d3_x_time_elevMedian,
      // Data
        // Integer IDs for random effect grouping terms
          int[] id_spCl,         // species-cluster
          int[] id_spSr,         // subregion
          int[] id_sp,           // species
          int[] id_fam,          // families
          int[,] id_spObs,       // species-observer
        // Q-vector (binary indicator of whether a species-point has any detections)
          int[] Q,               
        // Number of visits
          int[] nv,
        // Elevations
          vector relev,    
          vector relev2, 
        // Lowland occurrence
          int[] lowland,
        // Pasture             
          vector pasture, 
        // Biogeography
          // Barriers
            int[] mountain_barrier, 
            int[] valley_barrier,
          // // Elevations
            int[] elevMedian_id,
            int[] elevBreadth_id,
          // Habitats
            int[] forestPresent,
            int[] forestSpecialist,
            int[] tfSpecialist,
            int[] dryForestPresent,
            int[] floodDrySpecialist,
            int[] aridPresent,
          // Migratory
            int[] migratory,
        // Traits
            int[] mass_id,
            // diet (effects-coded; reference category = omnivore)
              int[] dietInvert,
              int[] dietCarn,
              int[] dietFruitNect,
              int[] dietGran,
        // Distance to range
          vector distance_to_range,
      // Detection nuisance covariates
        matrix time,     
        matrix obsSM,
        matrix obsDE,
        matrix obsJG,
      // pre-computed interactions
        int[] lowland_x_relev,
        int[] lowland_x_relev2,
        vector elevMedian_x_forestPresent,
        vector elevMedian_x_forestSpecialist,
        int[] mountainBarrier_x_pasture,
        int[] valleyBarrier_x_pasture,
        vector elevMedian_x_pasture,
        vector elevBreadth_x_pasture,
        int[] forestPresent_x_pasture,
        int[] forestSpecialist_x_pasture,
        int[] tfSpecialist_x_pasture,
        int[] dryForestPresent_x_pasture,
        int[] floodDrySpecialist_x_pasture,
        int[] aridPresent_x_pasture,
        int[] migratory_x_pasture,
        vector elevMedian_x_forestPresent_x_pasture,
        vector elevMedian_x_forestSpecialist_x_pasture,
        vector mass_x_pasture,
        int[] dietInvert_x_pasture,
        int[] dietCarn_x_pasture,
        int[] dietFruitNect_x_pasture,
        int[] dietGran_x_pasture,
        matrix time_x_elev
  ){   // End function arguments, begin computation
    // indexing variables
      int len = 1 + end - start;
      int r0 = start - 1;
    // slices of Q and nv
      int Q_slice[len] = Q[start:end];
      int nv_slice[len] = nv[start:end];
      int id_spObs_slice[len,4] = id_spObs[start:end,];
    // variables for computation
      vector[len] lp;             // container for row-wise log-probability increments
      vector[len] logit_psi;      // rowise logit_psi
      row_vector[4] logit_theta;  // single-row logit_theta
            // For points with < n_visit_max (4) visits, trailing elements of 
              // logit_theta are nonsense, computed from trailing 0s on the visit
              //  covariates. det_data has trailing -1s, ensuring that nonexistent 
              //  visits cannot accidentally slip into analysis.
      vector[len] logit_theta_vector; // vectorizable part of logit_theta
      matrix[len, 4] logit_theta_matrix;  // visit-specific logit_theta terms
    // Computation:
      logit_psi =
        // Intercepts
          mu_b0 + 
          // Spatial
          b0_spCl[id_spCl[start:end]] + b0_spSr[id_spSr[start:end]] +
          // Taxonomic
          b0_sp[id_sp[start:end]] + b0_fam[id_fam[start:end]] +
        // Slopes
          // Elevation
            // Linear
              (mu_b1_relev + b1_relev_sp[id_sp[start:end]]) .* relev[start:end] + 
            // Quadratic 
              (mu_b1_relev2 + b1_relev2_sp[id_sp[start:end]]) .* relev2[start:end] +
            // Lowland
              b1_lowland_unique[lowland[start:end]] + b1_x_lowland_relev_unique[lowland_x_relev[start:end]] + 
              b1_x_lowland_relev2_unique[lowland_x_relev2[start:end]] +
          // Pasture
            (mu_b2_pasture + b2_pasture_sp[id_sp[start:end]] + b2_pasture_fam[id_fam[start:end]]) .* pasture[start:end] +
          // Biogeography
            // Barriers
              b3_mountain_barrier_unique[mountain_barrier[start:end]] + b3_valley_barrier_unique[valley_barrier[start:end]] +
            // Elevations
              b3_elevMedian_unique[elevMedian_id[start:end]] + b3_elevBreadth_unique[elevBreadth_id[start:end]] + 
            // Habitats
              b3_forestPresent_unique[forestPresent[start:end]] + b3_forestSpecialist_unique[forestSpecialist[start:end]] + 
              b3_tfSpecialist_unique[tfSpecialist[start:end]] + b3_dryForestPresent_unique[dryForestPresent[start:end]] + 
              b3_floodDrySpecialist_unique[floodDrySpecialist[start:end]] + b3_aridPresent_unique[aridPresent[start:end]] +
            // Migratory
              b3_migratory_unique[migratory[start:end]] + 
            // Interactions
              b3_x_elevMedian_forestPresent*elevMedian_x_forestPresent[start:end] + 
              b3_x_elevMedian_forestSpecialist*elevMedian_x_forestSpecialist[start:end] +
          // Traits
            b3_mass_unique[mass_id[start:end]] + 
            // Diet
              b3_dietInvert_unique[dietInvert[start:end]] + b3_dietCarn_unique[dietCarn[start:end]] + 
              b3_dietFruitNect_unique[dietFruitNect[start:end]] + b3_dietGran_unique[dietGran[start:end]] +
          // Pasture interactions
            // Biogeography
              // Barriers
                b4_mountain_barrier_unique[mountainBarrier_x_pasture[start:end]] + b4_valley_barrier_unique[valleyBarrier_x_pasture[start:end]] +
              // Elevations
                b4_elevMedian*elevMedian_x_pasture[start:end] + b4_elevBreadth*elevBreadth_x_pasture[start:end] + 
              // Habitats
                b4_forestPresent_unique[forestPresent_x_pasture[start:end]] + b4_forestSpecialist_unique[forestSpecialist_x_pasture[start:end]] + 
                b4_tfSpecialist_unique[tfSpecialist_x_pasture[start:end]] + b4_dryForestPresent_unique[dryForestPresent_x_pasture[start:end]] + 
                b4_floodDrySpecialist_unique[floodDrySpecialist_x_pasture[start:end]] + b4_aridPresent_unique[aridPresent_x_pasture[start:end]] +
              // Migratory
                b4_migratory_unique[migratory_x_pasture[start:end]] + 
              // Interactions
                b4_x_elevMedian_forestPresent*elevMedian_x_forestPresent_x_pasture[start:end] + 
                b4_x_elevMedian_forestSpecialist*elevMedian_x_forestSpecialist_x_pasture[start:end] +
              // Traits
                b4_mass*mass_x_pasture[start:end] + 
                // Diet
                  b4_dietInvert_unique[dietInvert_x_pasture[start:end]] + b4_dietCarn_unique[dietCarn_x_pasture[start:end]] + 
                  b4_dietFruitNect_unique[dietFruitNect_x_pasture[start:end]] + b4_dietGran_unique[dietGran_x_pasture[start:end]] +
          // Range maps
            (mu_b5_distance_to_range + b5_distance_to_range_sp[id_sp[start:end]]) .* distance_to_range[start:end];
      
      logit_theta_vector = 
        // Intercepts
          mu_d0 +
          // Taxonomic
            d0_sp[id_sp[start:end]] + d0_fam[id_fam[start:end]] + 
          // Species-observer: see loop
        // Slopes
          // Pasture
            (mu_d1_pasture + d1_pasture_sp[id_sp[start:end]] + d1_pasture_fam[id_fam[start:end]]) .* pasture[start:end] +
          // Nuisance
            // Time: see matrix part
            // Observer: see matrix part
            // Trait
              d3_mass_unique[mass_id[start:end]] + d3_elevMedian_unique[elevMedian_id[start:end]] + 
              d3_migratory_unique[migratory[start:end]] + 
              d3_dietCarn_unique[dietCarn[start:end]];
        
      logit_theta_matrix =
        // Time
          diag_pre_multiply(mu_d2_time + d2_time_sp[id_sp[start:end]], time[start:end,]) + 
        // Observer
          d2_obsSM*obsSM[start:end,] + d2_obsDE*obsDE[start:end,] + d2_obsJG*obsJG[start:end,] +
        // Interaction
          d3_x_time_elevMedian*time_x_elev[start:end,];
        
      for (r in 1:len) {  // loop over species-points
        logit_theta = logit_theta_vector[r] + logit_theta_matrix[r,1:nv_slice[r]] + d0_spObs[id_spObs_slice[r,1:nv_slice[r]]];
        // likelihood
        if (Q_slice[r] == 1) {
          lp[r] = log_inv_logit(logit_psi[r]) +  // likelihood of occupancy
          bernoulli_logit_lpmf(det_slice[r, 1:nv_slice[r]] | logit_theta);   // conditional likelihood of observed history
        } else {
          lp[r] = log_sum_exp(
            log_inv_logit(logit_psi[r]) + // likelihood of occupancy
            sum(log1m_inv_logit(logit_theta)), // conditional likelihood of all-zero detection history
            log1m_inv_logit(logit_psi[r])); // likelihood of non-occupancy
        }
      } 
    return sum(lp);
  } // end of function
} // Close the functions block

data {
  // grainsize for reduce_sum 
    int<lower=1> grainsize;  
  // dimensions
    int<lower=1> n_spCl; // number of unique species:cluster 
    int<lower=1> n_spSr; // number of unique species:subregion
    int<lower=1> n_sp; // number of species
    int<lower=1> n_fam; // number of families
    int<lower=1> n_spObs; // number of unique species:observer
    int<lower=1> n_tot; // number of total rows (number of unique species:point)
    int<lower=1> n_visit_max; // maximum number of visits to a point
    
    int<lower=1> n_elevMedian;
    int<lower=1> n_elevBreadth;
    int<lower=1> n_mass;
    int<lower=1> n_lowland_x_relev;
    int<lower=1> n_lowland_x_relev2;
  // Detection array
    int<lower=-1, upper=1> det_data[n_tot, n_visit_max]; // detection history. -1s encode missing values
  // Q
    int<lower=0, upper=1> Q[n_tot];
  // number of visits
    int<lower=1, upper=n_visit_max> nv[n_tot];
  // Covariates
    // Random effect levels
      int<lower=1> id_spCl[n_tot];
      int<lower=1> id_spSr[n_tot];
      int<lower=1> id_sp[n_tot];
      int<lower=1> id_fam[n_tot];
      int<lower=0> id_spObs[n_tot, n_visit_max]; // Nonexistent visits are coded as zero.
    // Elevation and lowland
      vector[n_tot] relev;
      vector[n_tot] relev2;
      int lowland[n_tot];
    // Pasture
      vector[n_tot] pasture;
    // Biogeogaphy
      // Barriers
        int<lower=1,upper=2> mountain_barrier[n_tot]; // a vector of 1s and 2s; 1 if no mountain barrier, 2 if mountain barrier.
        int<lower=1,upper=2> valley_barrier[n_tot];
      // Elevations
        vector[n_elevMedian] elevMedian_u;
        int<lower=1> elevMedian_id[n_tot];
        
        vector[n_elevBreadth] elevBreadth_u;
        int<lower=1> elevBreadth_id[n_tot];
      // Habitats
        int<lower=1,upper=2> forestPresent[n_tot];
        int<lower=1,upper=2> forestSpecialist[n_tot];
        int<lower=1,upper=2> tfSpecialist[n_tot];
        int<lower=1,upper=2> dryForestPresent[n_tot];
        int<lower=1,upper=2> floodDrySpecialist[n_tot];        
        int<lower=1,upper=2> aridPresent[n_tot];
      // Migratory
        int<lower=1,upper=2> migratory[n_tot];
    // Traits
      vector[n_mass] mass_u;
      int<lower=1> mass_id[n_tot];
      // Diet (effects-coded; reference category = omnivore)
        int<lower=1,upper=2> dietInvert[n_tot];
        int<lower=1,upper=2> dietCarn[n_tot];
        int<lower=1,upper=2> dietFruitNect[n_tot];
        int<lower=1,upper=2> dietGran[n_tot];
    // Distance to range
      vector[n_tot] distance_to_range;
    // Nuisance detection covariates
      matrix[n_tot, n_visit_max] time;
      matrix[n_tot, n_visit_max] obsSM;
      matrix[n_tot, n_visit_max] obsDE;
      matrix[n_tot, n_visit_max] obsJG;
      
    // pre-computed interactions
      int<lower=1,upper=2> mountainBarrier_x_pasture[n_tot];
      int<lower=1,upper=2> valleyBarrier_x_pasture[n_tot];
      vector[n_tot] elevMedian_x_forestPresent;
      vector[n_tot] elevMedian_x_forestSpecialist;
      int<lower=1,upper=2> forestPresent_x_pasture[n_tot];
      int<lower=1,upper=2> forestSpecialist_x_pasture[n_tot];
      int<lower=1,upper=2> tfSpecialist_x_pasture[n_tot];
      int<lower=1,upper=2> dryForestPresent_x_pasture[n_tot];
      int<lower=1,upper=2> floodDrySpecialist_x_pasture[n_tot];
      int<lower=1,upper=2> aridPresent_x_pasture[n_tot];
      int<lower=1,upper=2> migratory_x_pasture[n_tot];
      int<lower=1,upper=2> dietInvert_x_pasture[n_tot];
      int<lower=1,upper=2> dietCarn_x_pasture[n_tot];
      int<lower=1,upper=2> dietFruitNect_x_pasture[n_tot];
      int<lower=1,upper=2> dietGran_x_pasture[n_tot];
      vector[n_tot] elevMedian_x_forestPresent_x_pasture;
      vector[n_tot] elevMedian_x_forestSpecialist_x_pasture;
      vector[n_tot] elevMedian_x_pasture;
      vector[n_tot] elevBreadth_x_pasture;
      vector[n_tot] mass_x_pasture;
      matrix[n_tot, n_visit_max] time_x_elev;
      vector[n_lowland_x_relev] lowland_x_relev_u;
      int<lower=1> lowland_x_relev_id[n_tot];
      vector[n_lowland_x_relev2] lowland_x_relev2_u;
      int<lower=1> lowland_x_relev2_id[n_tot];
} // Close the data block

parameters {
  // Occupancy
    // Intercepts
      // Overall
        real mu_b0;
      // zero-centered taxonomic
        real<lower=0> sigma_b0_taxonomic;
        real<lower=0,upper=1> p_b0_taxonomic_sp;
        vector[n_sp] b0_sp_raw;
        vector[n_fam] b0_fam_raw;
      // zero-centered spatial
        real<lower=0> sigma_b0_spatial;
        real<lower=0,upper=1> p_b0_spatial_cl;
        vector[n_spCl] b0_spCl_raw;
        vector[n_spSr] b0_spSr_raw;
    // Slopes
      // zero-centered relev by species
        real mu_b1_relev;
        real<lower=0> sigma_b1_relev_sp;
        vector[n_sp] b1_relev_sp_raw;
      // zero-centered relev2 by species
        real mu_b1_relev2;
        real<lower=0> sigma_b1_relev2_sp;
        vector[n_sp] b1_relev2_sp_raw;
      // lowland interactions (universal)
        real b1_lowland;
        real b1_x_lowland_relev;
        real b1_x_lowland_relev2;
      // Pasture
        real mu_b2_pasture;
        // zero-centered taxonomic
          real<lower=0> sigma_b2_taxonomic;
          real<lower=0,upper=1> p_b2_taxonomic_sp;
          vector[n_sp] b2_pasture_sp_raw;
          vector[n_fam] b2_pasture_fam_raw;
      // Trait effects
        // Biogeography
          // Barriers
            real b3_mountain_barrier;
            real b3_valley_barrier;
          // Elevations
            real b3_elevMedian;
            real b3_elevBreadth;
          // Habitat
            real b3_forestPresent;
            real b3_forestSpecialist;
            real b3_tfSpecialist;
            real b3_dryForestPresent;
            real b3_floodDrySpecialist;
            real b3_aridPresent;
          // Migration
            real b3_migratory;
          // Interactions
            real b3_x_elevMedian_forestPresent;
            real b3_x_elevMedian_forestSpecialist;
        // Functional traits
          real b3_mass;
          real b3_dietInvert;
          real b3_dietCarn;
          real b3_dietFruitNect;
          real b3_dietGran;
      // Pasture-x-trait interactions
        // Biogeography
          // Barriers
            real b4_mountain_barrier; 
            real b4_valley_barrier;
          // Elevations
            real b4_elevMedian;
            real b4_elevBreadth;
          // Habitat
            real b4_forestPresent;
            real b4_forestSpecialist;
            real b4_tfSpecialist;
            real b4_dryForestPresent;
            real b4_floodDrySpecialist;
            real b4_aridPresent;
          // Migration
            real b4_migratory;
          // Interactions
            real b4_x_elevMedian_forestPresent;
            real b4_x_elevMedian_forestSpecialist;
        // Functional traits
          real b4_mass;
          real b4_dietInvert;
          real b4_dietCarn;
          real b4_dietFruitNect;
          real b4_dietGran;
      // Range maps
        real mu_b5_distance_to_range;
        real<lower=0> sigma_b5_distance_to_range_sp;
        vector[n_sp] b5_distance_to_range_sp_raw;
  // Detection
    // Intercepts
      real mu_d0;
      // zero-centered taxonomic
        real<lower=0> sigma_d0_taxonomic;
        real<lower=0,upper=1> p_d0_taxonomic_sp;
        vector[n_sp] d0_sp_raw;
        vector[n_fam] d0_fam_raw;
      // zero-centered species/observer
        real<lower=0> sigma_d0_spObs;
        row_vector[n_spObs] d0_spObs_raw;
    // Slopes
      // Pasture
        real mu_d1_pasture;
        // zero-centered taxonomic
          real<lower=0> sigma_d1_taxonomic;
          real<lower=0,upper=1> p_d1_taxonomic_sp;
          vector[n_sp] d1_pasture_sp_raw;
          vector[n_sp] d1_pasture_fam_raw;
      // zero-centered time by species
        real mu_d2_time;
        real<lower=0> sigma_d2_time_sp;
        vector[n_sp] d2_time_sp_raw;
      // observer
        real d2_obsSM;       
        real d2_obsDE;       
        real d2_obsJG;      
      // Trait effects
        real d3_mass;
        real d3_elevMedian;
        real d3_migratory;
        real d3_dietCarn;
        real d3_x_time_elevMedian;
} // End parameters block

transformed parameters {
  // Standard deviations for effects with combined variance terms
    real<lower=0> sigma_b0_spCl = sqrt(sigma_b0_spatial^2 * p_b0_spatial_cl);
    real<lower=0> sigma_b0_spSr = sqrt(sigma_b0_spatial^2 * (1-p_b0_spatial_cl));
    real<lower=0> sigma_b0_sp = sqrt(sigma_b0_taxonomic^2 * p_b0_taxonomic_sp);
    real<lower=0> sigma_b0_fam = sqrt(sigma_b0_taxonomic^2 * (1-p_b0_taxonomic_sp));
    real<lower=0> sigma_b2_pasture_sp = sqrt(sigma_b2_taxonomic^2 * p_b2_taxonomic_sp);
    real<lower=0> sigma_b2_pasture_fam = sqrt(sigma_b2_taxonomic^2 * (1-p_b2_taxonomic_sp));
    real<lower=0> sigma_d0_sp = sqrt(sigma_d0_taxonomic^2 * p_d0_taxonomic_sp);
    real<lower=0> sigma_d0_fam = sqrt(sigma_d0_taxonomic^2 * (1-p_d0_taxonomic_sp));
    real<lower=0> sigma_d1_pasture_sp = sqrt(sigma_d1_taxonomic^2 * p_d1_taxonomic_sp);
    real<lower=0> sigma_d1_pasture_fam = sqrt(sigma_d1_taxonomic^2 * (1-p_d1_taxonomic_sp));
}

model {
  // Redundancy trick
    vector[2] binary_contrasts = [-1, 1]';
    
    vector[2] b3_mountain_barrier_unique = b3_mountain_barrier * binary_contrasts;
    vector[2] b3_valley_barrier_unique = b3_valley_barrier * binary_contrasts;
    vector[2] b3_forestPresent_unique = b3_forestPresent * binary_contrasts;
    vector[2] b3_forestSpecialist_unique = b3_forestSpecialist * binary_contrasts;
    vector[2] b3_tfSpecialist_unique = b3_tfSpecialist * binary_contrasts;
    vector[2] b3_dryForestPresent_unique = b3_dryForestPresent * binary_contrasts;
    vector[2] b3_floodDrySpecialist_unique = b3_floodDrySpecialist * binary_contrasts;
    vector[2] b3_aridPresent_unique = b3_aridPresent * binary_contrasts;
    vector[2] b3_migratory_unique = b3_migratory * binary_contrasts;
    vector[2] b3_dietInvert_unique = b3_dietInvert * binary_contrasts;
    vector[2] b3_dietCarn_unique = b3_dietCarn * binary_contrasts;
    vector[2] b3_dietFruitNect_unique = b3_dietFruitNect * binary_contrasts;
    vector[2] b3_dietGran_unique = b3_dietGran * binary_contrasts;
    
    vector[2] b4_mountain_barrier_unique = b4_mountain_barrier * binary_contrasts;
    vector[2] b4_valley_barrier_unique = b4_valley_barrier * binary_contrasts;
    vector[2] b4_forestPresent_unique = b4_forestPresent * binary_contrasts;
    vector[2] b4_forestSpecialist_unique = b4_forestSpecialist * binary_contrasts;
    vector[2] b4_tfSpecialist_unique = b4_tfSpecialist * binary_contrasts;
    vector[2] b4_dryForestPresent_unique = b4_dryForestPresent * binary_contrasts;
    vector[2] b4_floodDrySpecialist_unique = b4_floodDrySpecialist * binary_contrasts;
    vector[2] b4_aridPresent_unique = b4_aridPresent * binary_contrasts;
    vector[2] b4_migratory_unique = b4_migratory * binary_contrasts;
    vector[2] b4_dietInvert_unique = b4_dietInvert * binary_contrasts;
    vector[2] b4_dietCarn_unique = b4_dietCarn * binary_contrasts;
    vector[2] b4_dietFruitNect_unique = b4_dietFruitNect * binary_contrasts;
    vector[2] b4_dietGran_unique = b4_dietGran * binary_contrasts;
    
    vector[n_elevMedian] b3_elevMedian_unique = b3_elevMedian * elevMedian_u;
    vector[n_elevBreadth] b3_elevBreadth_unique = b3_elevBreadth * elevBreadth_u;
    vector[n_mass] b3_mass_unique = b3_mass * mass_u;
   
    vector[2] d3_migratory_unique = d3_migratory * binary_contrasts;
    vector[2] d3_dietCarn_unique = d3_dietCarn * binary_contrasts;
    
    vector[n_mass] d3_mass_unique = d3_mass * mass_u;
    vector[n_elevMedian] d3_elevMedian_unique = d3_elevMedian * elevMedian_u;
    
    vector[2] b1_lowland_unique = b1_lowland * binary_contrasts;
    vector[n_lowland_x_relev] b1_x_lowland_relev_unique = b1_x_lowland_relev * lowland_x_relev_u;
    vector[n_lowland_x_relev2] b1_x_lowland_relev2_unique = b1_x_lowland_relev2 * lowland_x_relev2_u;

  // Manual non-centering to avoid sticky boundaries
    vector[n_sp] b0_sp = b0_sp_raw * sigma_b0_sp;
    vector[n_fam] b0_fam = b0_fam_raw * sigma_b0_fam;
    vector[n_spCl] b0_spCl = b0_spCl_raw * sigma_b0_spCl;
    vector[n_spSr] b0_spSr = b0_spSr_raw * sigma_b0_spSr;
    vector[n_sp] b1_relev_sp = b1_relev_sp_raw * sigma_b1_relev_sp;
    vector[n_sp] b1_relev2_sp = b1_relev2_sp_raw * sigma_b1_relev2_sp;
    vector[n_sp] b2_pasture_sp = b2_pasture_sp_raw * sigma_b2_pasture_sp;
    vector[n_fam] b2_pasture_fam = b2_pasture_fam_raw * sigma_b2_pasture_fam;
    vector[n_sp] b5_distance_to_range_sp = b5_distance_to_range_sp_raw * sigma_b5_distance_to_range_sp;
    vector[n_sp] d0_sp = d0_sp_raw * sigma_d0_sp;
    vector[n_fam] d0_fam = d0_fam_raw * sigma_d0_fam;
    row_vector[n_spObs] d0_spObs = d0_spObs_raw * sigma_d0_spObs;
    vector[n_sp] d1_pasture_sp = d1_pasture_sp_raw * sigma_d1_pasture_sp;
    vector[n_fam] d1_pasture_fam = d1_pasture_fam_raw * sigma_d1_pasture_fam;
    vector[n_sp] d2_time_sp = d2_time_sp_raw * sigma_d2_time_sp;
    profile("priors"){
  // Priors and Jacobian adjustments
    // Occupancy
      // Intercepts
        mu_b0 ~ normal(-7, 2.5);
        // Spatial
          sigma_b0_spatial ~ normal(0, 3);
          // implicit uniform prior on p_b0_spatial_cl
          b0_spCl_raw ~ std_normal();
          b0_spSr_raw ~ std_normal();
        // Taxonomic
          sigma_b0_taxonomic ~ normal(0, 2);
          // implicit uniform prior on p_b0_taxonomic_sp
          b0_sp_raw ~ std_normal();
          b0_fam_raw ~ std_normal();
      // Slopes
        // Elevations
          // Linear
            mu_b1_relev ~ normal(0, 5);
            sigma_b1_relev_sp ~ normal(0, 2);
            b1_relev_sp_raw ~ std_normal();
          // Quadratic
            mu_b1_relev2 ~ normal(0, 5);
            sigma_b1_relev2_sp ~ normal(0, 2);
            b1_relev2_sp_raw ~ std_normal();
          // Lowland
            b1_lowland ~ normal(0, 1);
            b1_x_lowland_relev ~ normal(0, 5);
            b1_x_lowland_relev2 ~ normal(0, 5);
        // Pasture
          mu_b2_pasture ~ normal(0, 1);
          // Taxonomic
            sigma_b2_taxonomic ~ normal(0, 1);
            // implicit uniform prior on p_b2_taxonomic_sp
            b2_pasture_sp_raw ~ std_normal();
            b2_pasture_fam_raw ~ std_normal();
        // Biogeographic
          b3_mountain_barrier ~ normal(0, 1);
          b3_valley_barrier ~ normal(0, 1);
          b3_elevMedian ~ normal(0, 1);
          b3_elevBreadth ~ normal(0, 1);
          b3_forestPresent ~ normal(0, 1);
          b3_forestSpecialist ~ normal(0, 1);
          b3_tfSpecialist ~ normal(0, 1);
          b3_dryForestPresent ~ normal(0, 1);
          b3_floodDrySpecialist ~ normal(0, 1);
          b3_aridPresent ~ normal(0, 1);
          b3_migratory ~ normal(0, 1);
          b3_x_elevMedian_forestPresent ~ normal(0, 1);
          b3_x_elevMedian_forestSpecialist ~ normal(0, 1);
        // Functional
          b3_mass ~ normal(0, .5);
          b3_dietInvert ~ normal(0, .5);
          b3_dietCarn ~ normal(0, .5);
          b3_dietFruitNect ~ normal(0, .5);
          b3_dietGran ~ normal(0, .5);
        // Pasture interactions
          // Biogeographic
            b4_mountain_barrier ~ normal(0, 1);
            b4_valley_barrier ~ normal(0, 1);
            b4_elevMedian ~ normal(0, 1);
            b4_elevBreadth ~ normal(0, 1);
            b4_forestPresent ~ normal(0, 1);
            b4_forestSpecialist ~ normal(0, 1);
            b4_tfSpecialist ~ normal(0, 1);
            b4_dryForestPresent ~ normal(0, 1);
            b4_floodDrySpecialist ~ normal(0, 1);
            b4_aridPresent ~ normal(0, 1);
            b4_migratory ~ normal(0, 1);
            b4_x_elevMedian_forestPresent ~ normal(0, 1);
            b4_x_elevMedian_forestSpecialist ~ normal(0, 1);
          // Functional
            b4_mass ~ normal(0, .5);
            b4_dietInvert ~ normal(0, .5);
            b4_dietCarn ~ normal(0, .5);
            b4_dietFruitNect ~ normal(0, .5);
            b4_dietGran ~ normal(0, .5);
        // Distance-to-range
          mu_b5_distance_to_range ~ normal(0, 5);
          sigma_b5_distance_to_range_sp ~ normal(0, 2);
          b5_distance_to_range_sp_raw ~ std_normal();
    // Detection
      // Intercepts
        mu_d0 ~ normal(-3,1);
        // Taxonomic
          sigma_d0_taxonomic ~ normal(0, 2);
          // Implicit uniform prior on p_d0_taxonomic
          d0_sp_raw ~ std_normal();
          d0_fam_raw ~ std_normal();
        // Observer
          sigma_d0_spObs ~ normal(0, 1);
          d0_spObs_raw ~ std_normal();
      // Slopes
        // Pasture
          mu_d1_pasture ~ normal(0,.75);
          // Taxonomic
            sigma_d1_taxonomic ~ normal(0,1);
            // Implicit uniform prior on p_d1_taxonomic
            d1_pasture_sp_raw ~ std_normal();
            d1_pasture_fam_raw ~ std_normal();
        // Time
          mu_d2_time ~ normal(0, .5);
          sigma_d2_time_sp ~ normal(0, 1);
          d2_time_sp_raw ~ std_normal();
        // Observer
          d2_obsSM ~ normal(0, .25);
          d2_obsDE ~ normal(0, .25);
          d2_obsJG ~ normal(0, .25);        
        // Traits
          d3_mass ~ normal(0, .5);
          d3_elevMedian ~ normal(0, .5);
          d3_migratory ~ normal(0, 1);
          d3_dietCarn ~ normal(0, .5);
          d3_x_time_elevMedian ~ normal(0,.5);
  }
  profile("reduce_sum"){
  // Likelihood computed via reduce_sum
    target += reduce_sum_static(
      // partial_sum function
        partial_sum, 
      // surveys
        det_data,
      // grainsize
        grainsize,
      // variable sizes
        n_spCl, n_spSr, n_sp, n_fam, n_spObs,
      // parameters: 
        // Occupancy intercepts
          mu_b0, b0_spCl, b0_spSr, b0_sp, b0_fam, 
        // Occupancy elevation
          mu_b1_relev, b1_relev_sp, mu_b1_relev2, b1_relev2_sp, 
          b1_lowland_unique, b1_x_lowland_relev_unique, b1_x_lowland_relev2_unique,
        // Occupancy pasture
          mu_b2_pasture, b2_pasture_sp, b2_pasture_fam,
        // Occupancy traits
          b3_mountain_barrier_unique, b3_valley_barrier_unique,
          b3_elevMedian_unique, b3_elevBreadth_unique,
          b3_forestPresent_unique, b3_forestSpecialist_unique, b3_tfSpecialist_unique, 
          b3_dryForestPresent_unique, b3_floodDrySpecialist_unique, 
          b3_aridPresent_unique, 
          b3_migratory_unique, 
          b3_x_elevMedian_forestPresent, b3_x_elevMedian_forestSpecialist,
          b3_mass_unique, 
          b3_dietInvert_unique, b3_dietCarn_unique, b3_dietFruitNect_unique,
          b3_dietGran_unique,  
          b4_mountain_barrier_unique, b4_valley_barrier_unique,
          b4_elevMedian, b4_elevBreadth, 
          b4_forestPresent_unique, b4_forestSpecialist_unique, b4_tfSpecialist_unique, b4_dryForestPresent_unique, b4_floodDrySpecialist_unique, b4_aridPresent_unique, 
          b4_migratory_unique,
          b4_x_elevMedian_forestPresent, b4_x_elevMedian_forestSpecialist,
          b4_mass, b4_dietInvert_unique, b4_dietCarn_unique, b4_dietFruitNect_unique, b4_dietGran_unique, 
        // Occupancy distance-to-range
          mu_b5_distance_to_range, b5_distance_to_range_sp,
        // Detection intercept  
          mu_d0, d0_sp, d0_fam, d0_spObs,
        // Detection pasture
          mu_d1_pasture, d1_pasture_sp, d1_pasture_fam, 
        // Detection nuisance
          mu_d2_time, d2_time_sp, d2_obsSM, d2_obsDE, d2_obsJG, 
        // Detection traits  
          d3_mass_unique, d3_elevMedian_unique, d3_migratory_unique, d3_dietCarn_unique, d3_x_time_elevMedian,
      // Data
        // random effect levels
          id_spCl, id_spSr, id_sp, id_fam, id_spObs,
        // Q and nv
          Q, nv, 
        // covariates
          relev, relev2, lowland_id, pasture, mountain_barrier, valley_barrier,
          elevMedian_id, elevBreadth_id, 
          forestPresent, forestSpecialist, tfSpecialist, dryForestPresent, floodDrySpecialist, aridPresent, migratory, 
          mass_id, dietInvert, dietCarn, dietFruitNect, dietGran, 
          distance_to_range,
          time, obsSM, obsDE, obsJG,
        // pre-computed interactions
          lowland_x_relev_id, lowland_x_relev2_id, elevMedian_x_forestPresent, elevMedian_x_forestSpecialist, mountainBarrier_x_pasture,
          valleyBarrier_x_pasture, elevMedian_x_pasture, elevBreadth_x_pasture, forestPresent_x_pasture, 
          forestSpecialist_x_pasture, tfSpecialist_x_pasture, dryForestPresent_x_pasture, floodDrySpecialist_x_pasture,
          aridPresent_x_pasture, 
          migratory_x_pasture, 
          elevMedian_x_forestPresent_x_pasture, 
          elevMedian_x_forestSpecialist_x_pasture, 
          mass_x_pasture, 
          dietInvert_x_pasture, dietCarn_x_pasture, dietFruitNect_x_pasture, dietGran_x_pasture, time_x_elev
    ); // end reduce_sum call
  }
} // end model block

