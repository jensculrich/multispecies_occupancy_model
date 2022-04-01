// multi-species occupancy model
// This Stan program defines an occupancy model for 
// a community of species using detection data from
// J sites visited on K surveys
// The model has a hierarchical sturcture with species occupancy and detection
// rates linked through shared distributions and defined by a covariance matrix.
// The model also fits parameters for the effect of species type and site type on 
// occupancy rates (but currently no covariates for detection rates)

// This model draws structure from the butterfly occupancy model written by Bob Carpenter, here:
// https://mc-stan.org/users/documentation/case-studies/dorazio-royle-occupancy.html

functions {
  
  matrix cov_matrix_2d(vector sigma, real rho) {
    matrix[2,2] Sigma;
    Sigma[1,1] = square(sigma[1]);
    Sigma[2,2] = square(sigma[2]);
    Sigma[1,2] = sigma[1] * sigma[2] * rho;
    Sigma[2,1] = Sigma[1,2];
    return Sigma;
  }

  real lp_observed(int X, int K, real logit_psi, real logit_theta) {
    // if species i is observed at site j:
    return log_inv_logit(logit_psi)
      // where logit_psi[i] = (uv[i, 1] + a1_species_occ * species_cov1[i]) + 
                                // (a1_site_occ * site_cov1[j])
      + binomial_logit_lpmf(X | K, logit_theta); 
      // The log binomial probability mass of x successes (observed)
      // in K trials (survey revisits) given the
      // logit-scaled chance of success logit_theta (detection prob)
      // where logit_theta[i] = uv[i, 2]
  }

  real lp_unobserved(int K, real logit_psi, real logit_theta) {
    // if species is observed elsewhere but not at site j
    return log_sum_exp(lp_observed(0, K, logit_psi, logit_theta),
                       log1m_inv_logit(logit_psi));
        // present but undetected (with probability psi_i(1-theta_i)K) or
        // log1m_inv_logit means species i is either absent with probability 1-psi_i) 
  }

}

data {
  
  int<lower=1> J;  // sites within region
  int<lower=1> K;  // visits to sites
  int<lower=1> N;  // observed species
  int<lower=0, upper=K> X[N, J];  // visits when species i was detected at site j
  real species_cov1[N];    // unit covariate 1
  real site_cov1[J]; // unit covariate 1

}

parameters {
  
  real a1_species_occ;  // species-level occupancy effect (specialization - factor)
  real a1_site_occ;     // site-level occupancy effect (site treatment)
  
  vector[2] uv[N];                // species-level (occupancy, detection)
  real<lower=-1,upper=1> rho_uv;  // correlation of (occupancy, detection)
  vector<lower=0>[2] sigma_uv;    // sd of (occupancy, detection)
  
  real alpha; // site-level (occupancy)
  real beta;  // site-level (detection)

}

transformed parameters {
  
  real logit_psi[N, J];  // log odds  of occurrence
  real logit_theta[N, J];  // log odds of detection
  
  for (i in 1:N){     // loop across all species
      for (j in 1:J){    // loop across all sites
        logit_psi[i, j] = alpha + // grand mean intercept
        uv[i, 1] + // species-level effect on the intercept
        a1_species_occ * species_cov1[i] + // effect of species category on occupancy
        a1_site_occ * site_cov1[j]; // effect of site category on occupancy
      }
  }
  
  for (i in 1:N) {   // loop across all species
    for (j in 1:J) {     // loop across all sites
          logit_theta[i, j] = beta + // grand mean intercept
          uv[i, 2]; // species-level effect on the intercept
    }
  }
  
}

model {
  // priors
  a1_species_occ ~ cauchy(0, 2.5);
  a1_site_occ ~ cauchy(0, 2.5);
  
  sigma_uv ~ cauchy(0, 2.5);
  (rho_uv + 1) / 2 ~ beta(2, 2);
  uv ~ multi_normal(rep_vector(0, 2), cov_matrix_2d(sigma_uv, rho_uv));
  
  alpha ~ cauchy(0, 2.5);
  beta ~ cauchy(0, 2.5);
  
  // likelihood
  // Stan can sample mean and sd of parameters by summing out the
  // parameter (marginalizing) across likelihood statements
  for (i in 1:N) { // loop across all species
    for (j in 1:J) { // loop across all sites
      if (X[i, j] > 0) // if species is detected at the specific site 1 or more times
        target += lp_observed(X[i, j], K, logit_psi[i, j], logit_theta[i, j]);
        
      else // species i is detected at another site, but not site j
        target += lp_unobserved(K, logit_psi[i, j], logit_theta[i, j]); 
    }
  }
  
}

generated quantities {
  
  //  a set of simulated (u,v) pairs for a new species, 
  // along with the logit-scaled psi and theta variables. 
  // These are used for posterior predictive checking to see if the 
  // actual species are distributed as characterized by their prior.
  vector[2] sim_uv;
  real logit_psi_sim;
  real logit_theta_sim;
  real sim_species_category;
  real sim_site_category;

  sim_uv = multi_normal_rng(rep_vector(0,2),
                             cov_matrix_2d(sigma_uv, rho_uv));
                             
  sim_species_category = bernoulli_rng(0.5);
  sim_site_category = bernoulli_rng(0.5);
                     
  logit_psi_sim = alpha + sim_uv[1] + // grand mean intercept
        a1_species_occ * sim_species_category + // effect of species category on occupancy
        a1_site_occ * sim_site_category;
        
  logit_theta_sim = beta + sim_uv[2];
  
}
