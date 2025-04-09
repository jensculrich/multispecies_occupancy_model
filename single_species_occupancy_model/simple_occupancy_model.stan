// single-species static occupancy model
// Jens Ulrich, Feb 2025

data {
  int<lower=1> n_sites;  // number of study sites
  int<lower=1> n_surveys;  // number of visits to each site
  int<lower=0, upper=n_surveys> detections[n_sites];  
  // 'detections' is a vector containing number of successful detections 
  // out of n_surveys maximum detections for each site
} // end data block

parameters { 
  real alpha; // site-level occupancy rate
  real beta;  // site-level detection rate
} // end parameters block

transformed parameters {
  // adding a transformed params block will make it easier to add more predictors later
  // to add predictors, for example, an site occurrence covariate, you could do something like:
  // logit_psi[i] = alpha + alpha1 * covariate_value[i];
  
  // define expected occurrence and detection rates given site predictors
  // right now there are no predictors so expected rates are equal to the intercept terms
  real logit_psi[n_sites];  // log odds  of occurrence
  real logit_p[n_sites];  // log odds of detection
  
  // occurrence predictors on a logit scale
  for (i in 1:n_sites){     // loop across all species
    logit_psi[i] = alpha // grand mean intercept
    // add predictors of interest here if you want to do more than estimate an intercept!
    ; // end linear predictor for occurrence
  } // end for loop across sites
  // detection predictors on a logit scale
  for (i in 1:n_sites) {   // loop across all species
    logit_p[i] = beta // grand mean intercept
    // add predictors of interest here if you want to do more than estimate an intercept!
    ; // end linear predictor for detection
  } // end for loop across sites
} // end transformed parameters block

model {
  // PRIORS
  alpha ~ normal(0, 2); // weakly informative prior
  beta ~ normal(0, 2); // weakly informative prior
  // LIKELIHOOD
  for (i in 1:n_sites) { // loop across all sites
    // if species is detected at the specific site at least once
    // then the species occurs there. lp_observed calculates
    // the probability density that species occurs given psi, plus the 
    // probability density that we successfully observed it on detections[i]/n_surveys
    if(detections[i] > 0) {
       // lp_observed:
       target += log_inv_logit(logit_psi[i]) +
                binomial_logit_lpmf(detections[i] | n_surveys, logit_p[i]);
    // else the species was never detected at the site
    } else {
       // Stan can sample the mean and sd of parameters by summing out the
      // parameter (marginalizing) across likelihood statements
      // lp_unobserved sums the probability density of:
      // 1) species occupies the site but was not detected on any survey, and
      // 2) the species does not occupy the site
      target += 
              // present but never detected
              log_sum_exp(log_inv_logit(logit_psi[i]) +
              binomial_logit_lpmf(0 | n_surveys, logit_p[i]),
              // not present
              log1m_inv_logit(logit_psi[i])); 
    } // end if/else ever observed
  } // end for loop across sites
} // end model block
