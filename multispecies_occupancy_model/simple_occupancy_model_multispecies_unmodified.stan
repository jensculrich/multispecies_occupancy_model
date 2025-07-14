// multi-species static occupancy model
// Jens Ulrich, Feb 2025

data {
  int<lower=1> n_sites;  // number of study sites
  int site[n_sites];
  int<lower=1> n_species; // number of species
  int species[n_species];
  int<lower=1> n_surveys;  // number of visits to each site
  int<lower=0, upper=n_surveys> detections[n_sites, n_species];  
  // 'detections' is a vector containing number of successful detections 
  // out of n_surveys maximum detections for each site
} // end data block

parameters { 
  vector[n_species] alpha; // site-level occupancy rate
  real alpha0;
  real<lower=0> sigma_alpha;
  vector[n_species] beta;  // site-level detection rate
  real beta0;
  real<lower=0> sigma_beta;
} // end parameters block

transformed parameters {
  // adding a transformed params block will make it easier to add more predictors later
  // to add predictors, for example, an site occurrence covariate, you could do something like:
  // logit_psi[i] = alpha + alpha1 * covariate_value[i];
  
  // define expected occurrence and detection rates given site predictors
  // right now there are no predictors so expected rates are equal to the intercept terms
  real logit_psi[n_sites, n_species];  // log odds  of occurrence
  real logit_p[n_sites, n_species];  // log odds of detection
  
  // occurrence predictors on a logit scale
  for (i in 1:n_sites){     // loop across all species
    for (j in 1:n_species){
      logit_psi[i,j] = alpha[species[j]] // grand mean intercept
        // add predictors of interest here if you want to do more than estimate an intercept!
        ; // end linear predictor for occurrence
        
      logit_p[i,j] = beta[species[j]] // grand mean intercept
        // add predictors of interest here if you want to do more than estimate an intercept!
        ; // end linear predictor for detection
    }
  } // end for loop across sites

} // end transformed parameters block

model {
  // PRIORS
  alpha ~ normal(alpha0, sigma_alpha); // weakly informative prior
  alpha0 ~ normal(0, 2);
  sigma_alpha ~ normal(0, 2);
  beta ~ normal(beta0, sigma_beta); // weakly informative prior
  beta0 ~ normal(0, 2);
  sigma_beta ~ normal(0, 2);
  // LIKELIHOOD
  for (i in 1:n_sites) { // loop across all sites
    for (j in 1:n_species){
      
    // if species is detected at the specific site at least once
    // then the species occurs there. lp_observed calculates
    // the probability density that species occurs given psi, plus the 
    // probability density that we successfully observed it on detections[i]/n_surveys
    if(detections[i,j] > 0) {
       // lp_observed:
       target += log_inv_logit(logit_psi[i,j]) +
                binomial_logit_lpmf(detections[i,j] | n_surveys, logit_p[i,j]);
    // else the species was never detected at the site
    } else {
       // Stan can sample the mean and sd of parameters by summing out the
      // parameter (marginalizing) across likelihood statements
      // lp_unobserved sums the probability density of:
      // 1) species occupies the site but was not detected on any survey, and
      // 2) the species does not occupy the site
      target += 
              // present but never detected
              log_sum_exp(log_inv_logit(logit_psi[i,j]) +
              binomial_logit_lpmf(0 | n_surveys, logit_p[i,j]),
              // not present
              log1m_inv_logit(logit_psi[i,j])); 
    } // end if/else ever observed
    
    }
  } // end for loop across sites
} // end model block
