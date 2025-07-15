// dynamic multi-species occupancy model 

// jcu, started may, 2024.

// throughout I denote dimensions as
// species == i
// site == j
// year == k
// visit == l

data {
  
  int<lower=0> n_species; // number of species
  int<lower=1> species[n_species]; // vector of species identities
  int<lower=0> n_sites; // number of sites
  int<lower=1> sites[n_sites]; // vector of site identities
  int<lower=0> n_years; // total years
  int<lower=0> n_years_minus1;
  int<lower=0> n_surveys; // surveys per year
  int<lower=0,upper=n_surveys> y[n_species, n_sites, n_years]; // number of detections per year
  
} // end data

parameters {

  // initial state
  real psi1_0;
  vector[n_species] psi1_species_raw;
  real<lower=0> sigma_psi1_species;

  // colonization
  real gamma0;
  vector[n_species] gamma_species_raw;
  real<lower=0> sigma_gamma_species;

  // persistence
  real phi0;
  vector[n_species] phi_species_raw;
  real<lower=0> sigma_phi_species;
  
  // detection
  real p0; // intercept
  vector[n_species] p_species_raw;
  real<lower=0> sigma_p_species;
  
} // end parameters

transformed parameters {

  // logit scale psi1, gamma, phi
  real psi1[n_species, n_sites]; // odds of occurrence year 1
  real gamma[n_species, n_sites, n_years_minus1]; // odds of colonization
  real phi[n_species, n_sites, n_years_minus1]; // odds of persistence
  real p[n_species, n_sites, n_years];  // odds of detection
  
  vector[n_species] psi1_species;
  vector[n_species] gamma_species;
  vector[n_species] phi_species;
  vector[n_species] p_species;
  
  // implies: xprocess_species ~ normal(mu_xprocess_species, sigma_xprocess_species)
  psi1_species = psi1_0 + sigma_psi1_species * psi1_species_raw;
  gamma_species = gamma0 + sigma_gamma_species * gamma_species_raw;
  phi_species = phi0 + sigma_phi_species * phi_species_raw;
  p_species = p0 + sigma_p_species * p_species_raw;
  
  for(i in 1:n_species){
    for(j in 1:n_sites){    // loop across all sites
      for(k in 1:n_years_minus1){ // loop across all years
  
        psi1[i,j] = inv_logit( // probability (0-1) of occurrence in year 1 is equal to..
          psi1_species[species[i]] // a species specific intercept
          ); // end psi1[j,k]
        
        gamma[i,j,k] = inv_logit( // probability (0-1) of colonization is equal to..
          gamma_species[species[i]] // a species specific intercept
          ); // end gamma[i,j,k]
                
        phi[i,j,k] = inv_logit( // probability (0-1) of persistence is equal to..
          phi_species[species[i]] // a species specific intercept
          ); // end phi[i,j,k]

      } // end loop across all years
    } // end loop across all sites
  } // end loop across all species 
  
  for(i in 1:n_species){
    for(j in 1:n_sites){    // loop across all sites
      for(k in 1:n_years){ // loop across all years
        
        p[i,j,k] = inv_logit( // probability (0-1) of detection is equal to..
          p_species[species[i]] // a species specific intercept); // end p[j,k,l]
          ); // end p[i,j,k]  
        
      } // end loop across all years
    } // end loop across all sites
  } // end loop across all species 
   
  // construct an occurrence array
  real psi[n_species, n_sites, n_years];
  
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_years){
        
        if(k < 2){ // define initial state
          psi[i,j,k] = psi1[i,j]; 
        } else { // describe temporally autocorrelated system dynamics
          // As psi approaches 1, there's a weighted switch on phi (survival)
          // As psi approaches 0, there's a weighted switch on gamma (colonization)
          // reduce 1 from k for phi and gamma because there are n_years - 1 transitions 
          // and so there are only n_years - 1 speciesXsite "stacks" of phi and gamma
          // but phi[,,k-1] for k = 2 will actually consider the effects of 
          // e.g. flower abundacnce in year 2 (since year 2 is the first year we estimate phi)
          psi[i,j,k] = psi[i,j,k-1] * phi[i,j,k-1] + (1 - psi[i,j,k-1]) * gamma[i,j,k-1]; 
        } // end if/else
        
      } // end loop across all years
    } // end loop across all sites
  } // end loop across all species
   
} // end transformed parameters

model {
  
  // PRIORS
  
  // occupancy
  // initial state
  psi1_0 ~ normal(0, 2); 
  psi1_species_raw ~ std_normal();
  sigma_psi1_species ~ normal(0, 1);

  // colonization
  gamma0 ~ normal(0, 2); // persistence intercept
  gamma_species_raw ~ std_normal();
  sigma_gamma_species ~ normal(0, 1);

  // persistence
  phi0 ~ normal(0, 2); // global intercept
  phi_species_raw ~ std_normal();
  sigma_phi_species ~ normal(0, 1);

  // detection
  p0 ~ normal(0, 2); // global intercept
  p_species_raw ~ std_normal();
  sigma_p_species ~ normal(0, 1);

  // LIKELIHOOD
  for(i in 1:n_species){
    for (j in 1:n_sites){
      for (k in 1:n_years){
        
          if (y[i,j,k] > 0){ // lp observed 
            // detection on each visit given detection rate on each visit
            // V_NA == 0 indicates that no survey occurred. Multiplying 0 by the 
            // lpmf() statement should remove it from the target
            target += (log(psi[i,j,k]) + // present
                      binomial_lpmf(y[i,j,k] | n_surveys, p[i,j,k]));
          } else { // lp unobserved (set up for 6 annual surveys)
            // marginal likelihood of:
            // occurrence, yet...
            // non-detection on each visit given detection rate on each visit given occurrence
            target += (log_sum_exp(log(psi[i,j,k]) + 
                       binomial_lpmf(0 | n_surveys, p[i,j,k]),
                                    // or just simple no occurrence
                                    log1m(psi[i,j,k])));
          } // end if/else
          
      } // end loop across all years
    } // end loop across all sites   
  } // end loop across all species

} // end model

generated quantities{
  
  //
  // posterior predictive check (number of detections, binned by species)
  //
  int<lower=0> W_species_rep[n_species]; // sum of simulated detections
  
  int z_simmed[n_species, n_sites, n_years]; // simulate occurrence

  for(i in 1:n_species){
   for(j in 1:n_sites){
     for(k in 1:n_years){
          z_simmed[i,j,k] = bernoulli_rng(psi[i,j,k]); 
      }    
    }
  }
  
  // initialize repped detections at 0
  for(i in 1:n_species){
    W_species_rep[i] = 0;
  }
      
  // generating posterior predictive distribution
  // Predict Z at sites
  for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_years){ // loop across all years

          // detections in replicated data (us z_simmed from above)
          W_species_rep[i] = W_species_rep[i] + 
            // multiply by the NA indicator - if we didn't survey in real life
            // we don't survey in this simulation.
            (z_simmed[i,j,k] * binomial_rng(n_surveys, p[i,j,k]));
           
      } // end loop across years
    } // end loop across sites
  } // end loop across species
  
} // end generated quantities
