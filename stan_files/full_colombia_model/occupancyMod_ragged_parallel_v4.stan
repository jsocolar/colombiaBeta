// This is a Stan model for the full Colombia bird dataset.

functions{
  real partial_sum(
    
    // Function arguments:
      
      // Data slicing and indexing
    // Detection slice            
    int[,] det_slice,      // a slice of the detection array where rows are species-points and columns are visits
    
    // cutpoints for slicing                       
    int start,             // the starting row of the detection slice
    int end,              // the ending row of the detection slice   
    
    // numbers of random effect levels
    int n_spCl,             // number of species-clusters
    int n_sp,               // number of species
    int n_fam,              // number of families
    
    // Parameters             
    // Occupancy
    // Spatial effects
    vector b0_spCl,       // intercepts by species-cluster (zero-centered)
    
    // Taxonomic effects
    vector b0_sp,          // intercepts by species (includes the overall intercept)
    vector b0_fam,         // intercepts by family (zero-centered)
    
    // Elevation effects             
    vector b1_relev_sp,    // slopes for elevation by species
    vector b1_relev2_sp,   // slopes for elevation^2 by species
    
    real b1_lowland,             // modified logit-quadratic for species with minima of 0: intercept
    real b1_x_lowland_relev,     // modified logit-quadratic for species with minima of 0: linear term
    real b1_x_lowland_relev2,    // modified logit-quadratic for species with minima of 0: quadratic term
    
    // Pasture effects             
    vector b2_pasture_sp,  // slopes for pasture by species (includes the overall pasture effect)
    vector b2_pasture_fam, // slopes for pasture by family (zero-centered)         
    
    // Trait effects; covariates vary by species and all slopes are universal
    // Biogeography
    real b3_eastOnly, 
    real b3_westOnly,
    real b3_snsmOnly,
    real b3_notWandes,
    real b3_notEandes,
    
    // Elevations
    real b3_elevMedian,
    real b3_elevBreadth,
    
    // Habitats
    real b3_forestPresent,
    real b3_forestSpecialist,
    real b3_tfSpecialist,
    real b3_dryForestPresent,
    real b3_floodDrySpecialist,
    real b3_floodSpecialist,
    real b3_aridPresent,
    
    // Functional traits
    real b3_migratory,
    real b3_mass,
    
    real b3_dietInvert,              // diet (reference category == omnivore)
    real b3_dietCarn,
    real b3_dietFruitNect,
    real b3_dietGran,
    
    // Interactions
    real b3_x_elevMedian_forestPresent,
    real b3_x_elevMedian_forestSpecialist,
    
    // pasture-x-trait interactions: covariates vary by species and all slopes are universal. Identical to the b3 terms, but multipled by the pasture covariate.
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
    vector d0_sp,          // intercepts by species (includes the overall intercept)
    vector d0_fam,         // intercepts by family (zero-centered)
    
    // Pasture effects
    vector d1_pasture_sp,  // slopes for pasture by species
    vector d1_pasture_fam, // slopes by family
    
    // Nuisance effects
    vector d2_time_sp,     // slopes for time-of-day by species
    real d2_obsSM,         // fixed slope for observer Simon
    real d2_obsDE,         // fixed slope for observer David
    real d2_obsJG,         // fixed slope for observer James
    
    // Trait effects (universal slopes)
    real d3_mass,
    real d3_elevMedian,
    real d3_migratory,
    real d3_dietCarn,
    real d3_x_time_elevMedian,
    
    // Data
    // Integer IDs for random effect grouping terms
    int[] id_spCl,         // species-cluster
    int[] id_sp,           // species
    int[] id_fam,          // families
    
    // Q-vector
    int[] Q,               // 0/1 to indicate whether species-point has a detection
    
    // Number of visits
    int[] nv,              // number of visits to the given species-point (= number of visits to the point)
    
    // Elevations
    vector relev,    
    vector relev2, 
    // Is species lowland (a species trait, but interacts with elevations)
    vector lowland,
    
    // Pasture             
    vector pasture, 
    
    // Traits
    vector eastOnly, 
    vector westOnly,
    vector snsmOnly,
    vector notWandes,
    vector notEandes,
    
    vector elevMedian,
    vector elevBreadth,
    
    vector forestPresent, 
    vector forestSpecialist,
    vector tfSpecialist,
    vector dryForestPresent,
    vector floodDrySpecialist,
    vector floodSpecialist,
    vector aridPresent,
    
    vector migratory,
    vector mass,
    
    vector dietInvert,
    vector dietCarn,
    vector dietFruitNect,
    vector dietGran,
    
    // Detection nuisance covariates
    row_vector[] time,     
    row_vector[] obsSM,
    row_vector[] obsDE,
    row_vector[] obsJG
  ){   // End function arguments, begin computation
    
    // indexing variables
    int len = 1 + end - start;
    int r0 = start - 1;
    
    // variables for computation
    vector[len] lp;             // a vector to store the log-probability contribution from each row of det_slice
    real logit_psi;             // container for the value of logit_psi for a given row. Overwritten at each iteration of the loop.
    row_vector[4] logit_theta;      // container values of logit_theta for a given row. Overwritten at each iteration of the loop.
    // logit_theta will always have length 4. For points with < 4 visits, the trailing elements of theta
    // will be nonsense, computed from trailing 0s on the visit covariates.
    // det_data will have trailing -99s, thus ensuring that none of these nonexistent visits
    // accidentally slips into analysis.
    
    // Computation:
      for (r in 1:len) {  // loop over species-points
        // calculate logit_psi for that species-point
        logit_psi = 
          b0_spCl[id_spCl[r0+r]] + 
          b0_sp[id_sp[r0+r]] + b0_fam[id_fam[r0+r]] +
          
          b1_relev_sp[id_sp[r0+r]]*relev[r0+r] + b1_relev2_sp[id_sp[r0+r]]*relev2[r0+r] +
          b1_lowland*lowland[r0+r] + b1_x_lowland_relev*lowland[r0+r]*relev[r0+r] + b1_x_lowland_relev2*lowland[r0+r]*relev2[r0+r] +
          
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
        
        // calculate logit_theta for that species_point.  Note that this is vectorized because the time and observer covariates are row vectors
        // where nv < 4, logit_theta will have some trailing dummy terms.
        logit_theta = 
          d0_sp[id_sp[r0+r]] + d0_fam[id_fam[r0+r]] +
          d1_pasture_sp[id_sp[r0+r]]*pasture[r0+r] + d1_pasture_fam[id_fam[r0+r]]*pasture[r0+r] +
          d2_time_sp[id_sp[r0+r]]*time[r0+r] + d2_obsSM*obsSM[r0+r] + d2_obsDE*obsDE[r0+r] + d2_obsJG*obsJG[r0+r] +
          d3_mass*mass[r0+r] + d3_elevMedian*elevMedian[r0+r] + d3_migratory*migratory[r0+r] + d3_dietCarn*dietCarn[r0+r] + d3_x_time_elevMedian*time[r0+r]*elevMedian[r0+r];
        
        // likelihood
        if (Q[r0 + r] == 1) {
          lp[r] = log_inv_logit(logit_psi) +  // likelihood of occupancy
          bernoulli_logit_lpmf(det_slice[r, 1:nv[r0 + r]] | logit_theta[1:nv[r0 + r]]);   // likelihood of observed detection history given occupancy
        } else {
          lp[r] = log_sum_exp(
            log_inv_logit(logit_psi) + // likelihood of occupancy
            sum(log1m_inv_logit(logit_theta[1:nv[r0 + r]])), // likelihood of all-zero detection history given occupancy
            log1m_inv_logit(logit_psi)); // likelihood of non-occupancy
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
  int<lower=1> n_sp; // number of species
  int<lower=1> n_fam; // number of families
  int<lower=1> n_tot; // number of total rows (number of unique species:point)
  int<lower=1> n_visit_max; // maximum number of visits to a point
  
  // Detection array
  int<lower=-1, upper=1> det_data[n_tot, n_visit_max]; // detection history. -1s encode missing values
  
  // Q
  int<lower=0, upper=1> Q[n_tot];
  
  // number of visits
  int<lower=1, upper=n_visit_max> nv[n_tot];
  
  // Covariates
  // Random effect levels
  int<lower=1> id_spCl[n_tot];
  int<lower=1> id_sp[n_tot];
  int<lower=1> id_fam[n_tot];
  
  // Elevation and lowland trait
  vector[n_tot] relev;
  vector[n_tot] relev2;
  vector[n_tot] lowland;
  
  // Pasture
  vector[n_tot] pasture;
  
  // Traits
  vector[n_tot] eastOnly;
  vector[n_tot] westOnly;
  vector[n_tot] snsmOnly;
  vector[n_tot] notWandes;
  vector[n_tot] notEandes;
  
  vector[n_tot] elevMedian;
  vector[n_tot] elevBreadth;
  
  vector[n_tot] forestPresent;
  vector[n_tot] forestSpecialist;
  vector[n_tot] tfSpecialist;
  vector[n_tot] dryForestPresent;
  vector[n_tot] floodDrySpecialist;        
  vector[n_tot] floodSpecialist;        
  vector[n_tot] aridPresent;        
  
  vector[n_tot] migratory;        
  vector[n_tot] mass;
  
  vector[n_tot] dietInvert;
  vector[n_tot] dietCarn;
  vector[n_tot] dietFruitNect;
  vector[n_tot] dietGran;
  
  // Nuisance detection covariates
  row_vector[n_visit_max] time[n_tot];
  row_vector[n_visit_max] obsSM[n_tot];
  row_vector[n_visit_max] obsDE[n_tot];
  row_vector[n_visit_max] obsJG[n_tot];
  
  
  // offsets and multipliers
  real mu_b0_off; 
  real mu_b0_mult;
  
  real log_sigma_b0_sp_off;
  real log_sigma_b0_sp_mult;
  
  vector[n_sp] b0_sp_off;
  vector[n_sp] b0_sp_mult;
  
  real log_sigma_b0_fam_off;
  real log_sigma_b0_fam_mult;
  
  vector[n_fam] b0_fam_off;
  vector[n_fam] b0_fam_mult;
  
  real mu_b1_relev_off; 
  real mu_b1_relev_mult;
  
  real log_sigma_b1_relev_sp_off;
  real log_sigma_b1_relev_sp_mult;
  
  vector[n_sp] b1_relev_sp_off;
  vector[n_sp] b1_relev_sp_mult;
  
  real mu_b1_relev2_off;
  real mu_b1_relev2_mult;
  
  real log_sigma_b1_relev2_sp_off;
  real log_sigma_b1_relev2_sp_mult;
  
  vector[n_sp] b1_relev2_sp_off;
  vector[n_sp] b1_relev2_sp_mult;
  
  real b1_lowland_off;
  real b1_lowland_mult;
  real b1_x_lowland_relev_off;
  real b1_x_lowland_relev_mult;
  real b1_x_lowland_relev2_off;
  real b1_x_lowland_relev2_mult;
  
  real mu_b2_pasture_off;
  real mu_b2_pasture_mult;
  
  real log_sigma_b2_pasture_sp_off;
  real log_sigma_b2_pasture_sp_mult;
  
  vector[n_sp] b2_pasture_sp_off;
  vector[n_sp] b2_pasture_sp_mult;
  
  real log_sigma_b2_pasture_fam_off;
  real log_sigma_b2_pasture_fam_mult;
  
  vector[n_fam] b2_pasture_fam_off;
  vector[n_fam] b2_pasture_fam_mult;
  
  real b3_eastOnly_off;
  real b3_eastOnly_mult;
  
  real b3_westOnly_off;
  real b3_westOnly_mult;
  
  real b3_snsmOnly_off;
  real b3_snsmOnly_mult;
  
  real b3_notWandes_off;
  real b3_notWandes_mult;
  
  real b3_notEandes_off;
  real b3_notEandes_mult;
  
  real b3_elevMedian_off;
  real b3_elevMedian_mult;
  
  real b3_elevBreadth_off;
  real b3_elevBreadth_mult;
  
  real b3_forestPresent_off;
  real b3_forestPresent_mult;
  
  real b3_forestSpecialist_off;
  real b3_forestSpecialist_mult;
  
  real b3_tfSpecialist_off;
  real b3_tfSpecialist_mult;
  
  real b3_dryForestPresent_off;
  real b3_dryForestPresent_mult;
  
  real b3_floodDrySpecialist_off;
  real b3_floodDrySpecialist_mult;
  
  real b3_floodSpecialist_off;
  real b3_floodSpecialist_mult;
  
  real b3_aridPresent_off; 
  real b3_aridPresent_mult;
  
  real b3_migratory_off;
  real b3_migratory_mult;
  
  real b3_mass_off;
  real b3_mass_mult;
  
  real b3_dietInvert_off;
  real b3_dietInvert_mult;
  
  real b3_dietCarn_off;
  real b3_dietCarn_mult;
  
  real b3_dietFruitNect_off;
  real b3_dietFruitNect_mult;
  
  real b3_dietGran_off;
  real b3_dietGran_mult;
  
  real b3_x_elevMedian_forestPresent_off;
  real b3_x_elevMedian_forestPresent_mult;
  
  real b3_x_elevMedian_forestSpecialist_off; 
  real b3_x_elevMedian_forestSpecialist_mult;
  
  real b4_eastOnly_off; 
  real b4_eastOnly_mult;
  
  real b4_westOnly_off;
  real b4_westOnly_mult;
  
  real b4_snsmOnly_off;
  real b4_snsmOnly_mult;
  
  real b4_notWandes_off; 
  real b4_notWandes_mult;
  
  real b4_notEandes_off;
  real b4_notEandes_mult;
  
  real b4_elevMedian_off;
  real b4_elevMedian_mult;
  
  real b4_elevBreadth_off;
  real b4_elevBreadth_mult;
  
  real b4_forestPresent_off;
  real b4_forestPresent_mult;
  
  real b4_forestSpecialist_off;
  real b4_forestSpecialist_mult;
  
  real b4_tfSpecialist_off;
  real b4_tfSpecialist_mult;
  
  real b4_dryForestPresent_off;
  real b4_dryForestPresent_mult;
  
  real b4_floodDrySpecialist_off;
  real b4_floodDrySpecialist_mult;
  
  real b4_floodSpecialist_off;
  real b4_floodSpecialist_mult;
  
  real b4_aridPresent_off;
  real b4_aridPresent_mult;
  
  real b4_migratory_off;
  real b4_migratory_mult;
  
  real b4_mass_off;
  real b4_mass_mult;
  
  real b4_dietInvert_off;
  real b4_dietInvert_mult;
  
  real b4_dietCarn_off;
  real b4_dietCarn_mult;
  
  real b4_dietFruitNect_off;
  real b4_dietFruitNect_mult;
  
  real b4_dietGran_off;
  real b4_dietGran_mult;
  
  real b4_x_elevMedian_forestPresent_off;
  real b4_x_elevMedian_forestPresent_mult;
  
  real b4_x_elevMedian_forestSpecialist_off;
  real b4_x_elevMedian_forestSpecialist_mult;
  
  real mu_d0_off;
  real mu_d0_mult;
  
  real log_sigma_d0_sp_off;
  real log_sigma_d0_sp_mult;
  
  vector[n_sp] d0_sp_off;
  vector[n_sp] d0_sp_mult;
  
  real log_sigma_d0_fam_off;
  real log_sigma_d0_fam_mult;
  
  vector[n_fam] d0_fam_off;
  vector[n_fam] d0_fam_mult;
  
  real mu_d1_pasture_off;
  real mu_d1_pasture_mult;
  
  real log_sigma_d1_pasture_sp_off;
  real log_sigma_d1_pasture_sp_mult;
  
  vector[n_sp] d1_pasture_sp_off;
  vector[n_sp] d1_pasture_sp_mult;
  
  real log_sigma_d1_pasture_fam_off;
  real log_sigma_d1_pasture_fam_mult;
  
  vector[n_fam] d1_pasture_fam_off;
  vector[n_fam] d1_pasture_fam_mult;
  
  real mu_d2_time_off;
  real mu_d2_time_mult;
  
  real log_sigma_d2_time_sp_off;
  real log_sigma_d2_time_sp_mult;
  
  vector[n_sp] d2_time_sp_off;
  vector[n_sp] d2_time_sp_mult;
  
  real d2_obsSM_off;
  real d2_obsSM_mult;
  
  real d2_obsDE_off;
  real d2_obsDE_mult;
  
  real d2_obsJG_off;
  real d2_obsJG_mult;      
  
  real d3_mass_off;
  real d3_mass_mult;
  
  real d3_elevMedian_off;
  real d3_elevMedian_mult;
  
  real d3_migratory_off; 
  real d3_migratory_mult;
  
  real d3_dietCarn_off;
  real d3_dietCarn_mult;
  
  real d3_x_time_elevMedian_off;
  real d3_x_time_elevMedian_mult;
} // Close the data block

parameters {
  // Occupancy
  // Intercepts
  real<offset=mu_b0_off, multiplier=mu_b0_mult> mu_b0;
  
  real<lower=0> sigma_b0_spCl;
  vector<multiplier=sigma_b0_spCl>[n_spCl] b0_spCl;
  
  real<offset=log_sigma_b0_sp_off, multiplier=log_sigma_b0_sp_mult> log_sigma_b0_sp;
  vector<offset=mu_b0, multiplier=exp(log_sigma_b0_sp)>[n_sp] b0_sp;
  
  real<offset=log_sigma_b0_fam_off, multiplier=log_sigma_b0_fam_mult> log_sigma_b0_fam;
  vector<offset=0, multiplier=exp(log_sigma_b0_fam)>[n_fam] b0_fam;
  
  // Slopes
  // Elevation effects
  real<offset=mu_b1_relev_off, multiplier=mu_b1_relev_mult> mu_b1_relev;
  real<offset=log_sigma_b1_relev_sp_off, multiplier=log_sigma_b1_relev_sp_mult> log_sigma_b1_relev_sp;
  vector<offset=b1_relev_sp_off, multiplier=b1_relev_sp_mult>[n_sp] b1_relev_sp;
  
  real<offset=mu_b1_relev2_off, multiplier=mu_b1_relev2_mult> mu_b1_relev2;
  real<offset=log_sigma_b1_relev2_sp_off, multiplier=log_sigma_b1_relev2_sp_mult> log_sigma_b1_relev2_sp;
  vector<offset=b1_relev2_sp_off, multiplier=b1_relev2_sp_mult>[n_sp] b1_relev2_sp;
  
  real<offset=b1_lowland_off, multiplier=b1_lowland_mult> b1_lowland;
  real<offset=b1_x_lowland_relev_off, multiplier=b1_x_lowland_relev_mult> b1_x_lowland_relev;
  real<offset=b1_x_lowland_relev2_off, multiplier=b1_x_lowland_relev2_mult> b1_x_lowland_relev2;
  
  // Pasture effects
  real<offset=mu_b2_pasture_off, multiplier=mu_b2_pasture_mult> mu_b2_pasture;
  real<offset=log_sigma_b2_pasture_sp_off, multiplier=log_sigma_b2_pasture_sp_mult> log_sigma_b2_pasture_sp;
  vector<offset=mu_b2_pasture, multiplier=exp(log_sigma_b2_pasture_sp)>[n_sp] b2_pasture_sp;    
  
  real<offset=log_sigma_b2_pasture_fam_off, multiplier=log_sigma_b2_pasture_fam_mult> log_sigma_b2_pasture_fam;
  vector<offset=0, multiplier=exp(log_sigma_b2_pasture_fam)>[n_fam] b2_pasture_fam;    
  
  // Trait effects
  real<offset=b3_eastOnly_off, multiplier=b3_eastOnly_mult> b3_eastOnly;
  real<offset=b3_westOnly_off, multiplier=b3_westOnly_mult> b3_westOnly;
  real<offset=b3_snsmOnly_off, multiplier=b3_snsmOnly_mult> b3_snsmOnly;
  real<offset=b3_notWandes_off, multiplier=b3_notWandes_mult> b3_notWandes;
  real<offset=b3_notEandes_off, multiplier=b3_notEandes_mult> b3_notEandes;
  
  real<offset=b3_elevMedian_off, multiplier=b3_elevMedian_mult> b3_elevMedian;
  real<offset=b3_elevBreadth_off, multiplier=b3_elevBreadth_mult> b3_elevBreadth;
  
  real<offset=b3_forestPresent_off, multiplier=b3_forestPresent_mult> b3_forestPresent;
  real<offset=b3_forestSpecialist_off, multiplier=b3_forestSpecialist_mult> b3_forestSpecialist;
  real<offset=b3_tfSpecialist_off, multiplier=b3_tfSpecialist_mult> b3_tfSpecialist;
  real<offset=b3_dryForestPresent_off, multiplier=b3_dryForestPresent_mult> b3_dryForestPresent;
  real<offset=b3_floodDrySpecialist_off, multiplier=b3_floodDrySpecialist_mult> b3_floodDrySpecialist;
  real<offset=b3_floodSpecialist_off, multiplier=b3_floodSpecialist_mult> b3_floodSpecialist;
  real<offset=b3_aridPresent_off, multiplier=b3_aridPresent_mult> b3_aridPresent;
  
  real<offset=b3_migratory_off, multiplier=b3_migratory_mult> b3_migratory;
  real<offset=b3_mass_off, multiplier=b3_mass_mult> b3_mass;
  
  real<offset=b3_dietInvert_off, multiplier=b3_dietInvert_mult> b3_dietInvert;
  real<offset=b3_dietCarn_off, multiplier=b3_dietCarn_mult> b3_dietCarn;
  real<offset=b3_dietFruitNect_off, multiplier=b3_dietFruitNect_mult> b3_dietFruitNect;
  real<offset=b3_dietGran_off, multiplier=b3_dietGran_mult> b3_dietGran;
  
  real<offset=b3_x_elevMedian_forestPresent_off, multiplier=b3_x_elevMedian_forestPresent_mult> b3_x_elevMedian_forestPresent;
  real<offset=b3_x_elevMedian_forestSpecialist_off, multiplier=b3_x_elevMedian_forestSpecialist_mult> b3_x_elevMedian_forestSpecialist;
  
  // pasture-x-trait interactions
  real<offset=b4_eastOnly_off, multiplier=b4_eastOnly_mult> b4_eastOnly; 
  real<offset=b4_westOnly_off, multiplier=b4_westOnly_mult> b4_westOnly;
  real<offset=b4_snsmOnly_off, multiplier=b4_snsmOnly_mult> b4_snsmOnly;
  real<offset=b4_notWandes_off, multiplier=b4_notWandes_mult> b4_notWandes;
  real<offset=b4_notEandes_off, multiplier=b4_notEandes_mult> b4_notEandes;
  
  real<offset=b4_elevMedian_off, multiplier=b4_elevMedian_mult> b4_elevMedian;
  real<offset=b4_elevBreadth_off, multiplier=b4_elevBreadth_mult> b4_elevBreadth;
  
  real<offset=b4_forestPresent_off, multiplier=b4_forestPresent_mult> b4_forestPresent;
  real<offset=b4_forestSpecialist_off, multiplier=b4_forestSpecialist_mult> b4_forestSpecialist;
  real<offset=b4_tfSpecialist_off, multiplier=b4_tfSpecialist_mult> b4_tfSpecialist;
  real<offset=b4_dryForestPresent_off, multiplier=b4_dryForestPresent_mult> b4_dryForestPresent;
  real<offset=b4_floodDrySpecialist_off, multiplier=b4_floodDrySpecialist_mult> b4_floodDrySpecialist;
  real<offset=b4_floodSpecialist_off, multiplier=b4_floodSpecialist_mult> b4_floodSpecialist;
  real<offset=b4_aridPresent_off, multiplier=b4_aridPresent_mult> b4_aridPresent;
  
  real<offset=b4_migratory_off, multiplier=b4_migratory_mult> b4_migratory;
  real<offset=b4_mass_off, multiplier=b4_mass_mult> b4_mass;
  
  real<offset=b4_dietInvert_off, multiplier=b4_dietInvert_mult> b4_dietInvert;
  real<offset=b4_dietCarn_off, multiplier=b4_dietCarn_mult> b4_dietCarn;
  real<offset=b4_dietFruitNect_off, multiplier=b4_dietFruitNect_mult> b4_dietFruitNect;
  real<offset=b4_dietGran_off, multiplier=b4_dietGran_mult> b4_dietGran;
  
  real<offset=b4_x_elevMedian_forestPresent_off, multiplier=b4_x_elevMedian_forestPresent_mult> b4_x_elevMedian_forestPresent;
  real<offset=b4_x_elevMedian_forestSpecialist_off, multiplier=b4_x_elevMedian_forestSpecialist_mult> b4_x_elevMedian_forestSpecialist;
  
  // Detection
  // Intercepts
  real<offset=mu_d0_off, multiplier=mu_d0_mult> mu_d0;
  
  real<offset=log_sigma_d0_sp_off, multiplier=log_sigma_d0_sp_mult> log_sigma_d0_sp;
  vector<offset=mu_d0, multiplier=exp(log_sigma_d0_sp)>[n_sp] d0_sp;
  
  real<offset=log_sigma_d0_fam_off, multiplier=log_sigma_d0_fam_mult> log_sigma_d0_fam;
  vector<offset=0, multiplier=exp(log_sigma_d0_fam)>[n_fam] d0_fam;
  
  // Slopes
  // Pasture effects
  real<offset=mu_d1_pasture_off, multiplier=mu_d1_pasture_mult> mu_d1_pasture;
  real<offset=log_sigma_d1_pasture_sp_off, multiplier=log_sigma_d1_pasture_sp_mult> log_sigma_d1_pasture_sp;
  vector<offset=mu_d1_pasture, multiplier=exp(log_sigma_d1_pasture_sp)>[n_sp] d1_pasture_sp;
  
  real<offset=log_sigma_d1_pasture_fam_off, multiplier=log_sigma_d1_pasture_fam_mult> log_sigma_d1_pasture_fam;
  vector<offset=0, multiplier=exp(log_sigma_d1_pasture_fam)>[n_fam] d1_pasture_fam;
  
  // Nuisance effects
  real<offset=mu_d2_time_off, multiplier=mu_d2_time_mult> mu_d2_time;
  real<offset=log_sigma_d2_time_sp_off, multiplier=log_sigma_d2_time_sp_mult> log_sigma_d2_time_sp;
  vector<offset=mu_d2_time, multiplier=exp(log_sigma_d2_time_sp)>[n_sp] d2_time_sp;
  
  real<offset=d2_obsSM_off, multiplier=d2_obsSM_mult> d2_obsSM;       
  real<offset=d2_obsDE_off, multiplier=d2_obsDE_mult> d2_obsDE;       
  real<offset=d2_obsJG_off, multiplier=d2_obsJG_mult> d2_obsJG;      
  
  // Trait effects
  real<offset=d3_mass_off, multiplier=d3_mass_mult> d3_mass;
  real<offset=d3_elevMedian_off, multiplier=d3_elevMedian_mult> d3_elevMedian;
  real<offset=d3_migratory_off, multiplier=d3_migratory_mult> d3_migratory;
  real<offset=d3_dietCarn_off, multiplier=d3_dietCarn_mult> d3_dietCarn;
  real<offset=d3_x_time_elevMedian_off, multiplier=d3_x_time_elevMedian_mult> d3_x_time_elevMedian;
} // End parameters block

model {
  // Priors and Jacobian adjustments
  // Occupancy
  mu_b0 ~ normal(0, 5);
  
  sigma_b0_spCl ~ normal(0, 3);
  b0_spCl ~ normal(0, sigma_b0_spCl);
  
  real sigma_b0_sp = exp(log_sigma_b0_sp);
  target += log_sigma_b0_sp;
  sigma_b0_sp ~ normal(0, 3);
  b0_sp ~ normal(mu_b0, sigma_b0_sp);
  
  real sigma_b0_fam = exp(log_sigma_b0_fam);
  target += log_sigma_b0_fam;
  sigma_b0_fam ~ normal(0, 3);
  b0_fam ~ normal(0, sigma_b0_fam);
  
  mu_b1_relev ~ normal(0, 4);
  
  real sigma_b1_relev_sp = exp(log_sigma_b1_relev_sp);
  target += log_sigma_b1_relev_sp;
  sigma_b1_relev_sp ~ normal(0, 3);
  b1_relev_sp ~ normal(mu_b1_relev, sigma_b1_relev_sp);
  
  mu_b1_relev2 ~ normal(0, 4);
  
  real sigma_b1_relev2_sp = exp(log_sigma_b1_relev2_sp);
  target += log_sigma_b1_relev2_sp;
  sigma_b1_relev2_sp ~ normal(0, 3);
  b1_relev2_sp ~ normal(mu_b1_relev2, sigma_b1_relev2_sp);
  
  
  b1_lowland ~ normal(0, 4);
  b1_x_lowland_relev ~ normal(0, 4);
  b1_x_lowland_relev2 ~ normal(0, 4);
  
  
  mu_b2_pasture ~ normal(0, 4);
  
  real sigma_b2_pasture_sp = exp(log_sigma_b2_pasture_sp);
  target += log_sigma_b2_pasture_sp;
  sigma_b2_pasture_sp ~ normal(0, 3);
  b2_pasture_sp ~ normal(mu_b2_pasture, sigma_b2_pasture_sp);
  
  real sigma_b2_pasture_fam = exp(log_sigma_b2_pasture_fam);
  target += log_sigma_b2_pasture_fam;
  sigma_b2_pasture_fam ~ normal(0, 3);
  b2_pasture_fam ~ normal(0, sigma_b2_pasture_fam);
  
  b3_eastOnly ~ normal(0, 4);
  b3_westOnly ~ normal(0, 4);
  b3_snsmOnly ~ normal(0, 4);
  b3_notWandes ~ normal(0, 4);
  b3_notEandes ~ normal(0, 4);
  
  b3_elevMedian ~ normal(0, 4);
  b3_elevBreadth ~ normal(0, 4);
  b3_forestPresent ~ normal(0, 4);
  b3_forestSpecialist ~ normal(0, 4);
  b3_tfSpecialist ~ normal(0, 4);
  b3_dryForestPresent ~ normal(0, 4);
  b3_floodDrySpecialist ~ normal(0, 4);
  b3_floodSpecialist ~ normal(0, 4);
  b3_aridPresent ~ normal(0, 4);
  
  b3_migratory ~ normal(0, 4);
  b3_mass ~ normal(0, 4);
  
  b3_dietInvert ~ normal(0, 4);
  b3_dietCarn ~ normal(0, 4);
  b3_dietFruitNect ~ normal(0, 4);
  b3_dietGran ~ normal(0, 4);
  
  b3_x_elevMedian_forestPresent ~ normal(0, 4);
  b3_x_elevMedian_forestSpecialist ~ normal(0, 4);
  
  
  b4_eastOnly ~ normal(0, 4);
  b4_westOnly ~ normal(0, 4);
  b4_snsmOnly ~ normal(0, 4);
  b4_notWandes ~ normal(0, 4);
  b4_notEandes ~ normal(0, 4);
  
  b4_elevMedian ~ normal(0, 4);
  b4_elevBreadth ~ normal(0, 4);
  b4_forestPresent ~ normal(0, 4);
  b4_forestSpecialist ~ normal(0, 4);
  b4_tfSpecialist ~ normal(0, 4);
  b4_dryForestPresent ~ normal(0, 4);
  b4_floodDrySpecialist ~ normal(0, 4);
  b4_floodSpecialist ~ normal(0, 4);
  b4_aridPresent ~ normal(0, 4);
  
  b4_migratory ~ normal(0, 4);
  b4_mass ~ normal(0, 4);
  
  b4_dietInvert ~ normal(0, 4);
  b4_dietCarn ~ normal(0, 4);
  b4_dietFruitNect ~ normal(0, 4);
  b4_dietGran ~ normal(0, 4);
  
  b4_x_elevMedian_forestPresent ~ normal(0, 4);
  b4_x_elevMedian_forestSpecialist ~ normal(0, 4);
  
  // Detection
  mu_d0 ~ student_t(7.763, 0, 1.566);  // This is Dorazio's suggested prior, which is approximately uniform on the probability scale between 0.01 and 0.99. See also Northrup & Gerber 2018
  
        real sigma_d0_sp = exp(log_sigma_d0_sp);
        target += log_sigma_d0_sp;
        sigma_d0_sp ~ normal(0, 2);
        d0_sp ~ normal(mu_d0, sigma_d0_sp);
  
        real sigma_d0_fam = exp(log_sigma_d0_fam);
        target += log_sigma_d0_fam;
        sigma_d0_fam ~ normal(0, 2);
        d0_fam ~ normal(0, sigma_d0_fam);
  
        mu_d1_pasture ~ normal(0, 2);
  
        real sigma_d1_pasture_sp = exp(log_sigma_d1_pasture_sp);
        target += log_sigma_d1_pasture_sp;
        sigma_d1_pasture_sp ~ normal(0, 2);
        d1_pasture_sp ~ normal(mu_d1_pasture, sigma_d1_pasture_sp);
  
        real sigma_d1_pasture_fam = exp(log_sigma_d1_pasture_fam);
        target += log_sigma_d1_pasture_fam;
        sigma_d1_pasture_fam ~ normal(0, 2);
        d1_pasture_fam ~ normal(0, sigma_d1_pasture_fam);
  
        mu_d2_time ~ normal(0, 2);
  
        real sigma_d2_time_sp = exp(log_sigma_d2_time_sp);
        target += log_sigma_d2_time_sp;
        sigma_d2_time_sp ~ normal(0, 2);
        d2_time_sp ~ normal(mu_d2_time, sigma_d2_time_sp);
  
        d2_obsSM ~ normal(0, 1); // intentionally somewhat informative
        d2_obsDE ~ normal(0, 1);
        d2_obsJG ~ normal(0, 1);        
  
        d3_mass ~ normal(0, 2);
        d3_elevMedian ~ normal(0, 2);
        d3_migratory ~ normal(0, 2);
        d3_dietCarn ~ normal(0, 2);
        d3_x_time_elevMedian ~ normal(0,2);
  
  // Likelihood computed via reduce_sum
    target += reduce_sum(
      // partial_sum function
          partial_sum, 
    
      // surveys
          det_data,
    
      // grainsize
          grainsize,
    
      // variable sizes
          n_spCl, n_sp, n_fam, 
  
      // parameters: 
        // Occupancy intercepts
          b0_spCl, b0_sp, b0_fam, 
        
        // Occupancy elevation  
          b1_relev_sp, b1_relev2_sp, 
          b1_lowland, b1_x_lowland_relev, b1_x_lowland_relev2,
        
        // Occupancy pasture
          b2_pasture_sp, b2_pasture_fam,
    
        // Occupancy traits
          b3_eastOnly, b3_westOnly, b3_snsmOnly, b3_notWandes, b3_notEandes, b3_elevMedian, b3_elevBreadth,
          b3_forestPresent, b3_forestSpecialist, b3_tfSpecialist, b3_dryForestPresent, b3_floodDrySpecialist,
          b3_floodSpecialist, b3_aridPresent, b3_migratory, b3_mass, b3_dietInvert, b3_dietCarn, b3_dietFruitNect,
          b3_dietGran, b3_x_elevMedian_forestPresent, b3_x_elevMedian_forestSpecialist, b4_eastOnly, b4_westOnly,
          b4_snsmOnly, b4_notWandes, b4_notEandes,b4_elevMedian, b4_elevBreadth, b4_forestPresent, b4_forestSpecialist,
          b4_tfSpecialist, b4_dryForestPresent, b4_floodDrySpecialist, b4_floodSpecialist, b4_aridPresent, b4_migratory,
          b4_mass, b4_dietInvert, b4_dietCarn, b4_dietFruitNect, b4_dietGran, b4_x_elevMedian_forestPresent,
          b4_x_elevMedian_forestSpecialist,
          
        // Detection intercept  
          d0_sp, d0_fam, 
          
        // Detection pasture
          d1_pasture_sp, d1_pasture_fam, 
          
        // Detection nuisance
          d2_time_sp, d2_obsSM, d2_obsDE, d2_obsJG, 
          
        // Detection traits  
          d3_mass, d3_elevMedian, d3_migratory, d3_dietCarn, d3_x_time_elevMedian,
    
      // Data
        // random effect levels
          id_spCl, id_sp, id_fam, 

        // Q and nv
          Q, nv, 
    
        // covariates
          relev, relev2, lowland, pasture, eastOnly, westOnly, snsmOnly, notWandes, notEandes,
          elevMedian, elevBreadth, forestPresent, forestSpecialist, tfSpecialist, dryForestPresent, floodDrySpecialist,
          floodSpecialist, aridPresent, migratory, mass, dietInvert, dietCarn, dietFruitNect, dietGran, time, obsSM,
          obsDE, obsJG
  ); // end reduce_sum call
} // end model block
