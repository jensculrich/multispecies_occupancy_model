// multi-species occupancy model
// This Stan program defines a model for occupancy
// a community of species using detection data from
// J sites visited on K surveys

// adapted from model for butterfly occupancy written by Bob Carpenter, here:
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
                                // (ab[j, 1] + a1_site_occ * site_cov1[j])
      + binomial_logit_lpmf(X | K, logit_theta); 
      // The log binomial probability mass of x successes (observed)
      // in K trials (survey revisits) given the
      // logit-scaled chance of success logit_theta (detection prob)
      // where logit_theta[i] = uv[i, 2] + beta
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
  
  vector[2] ab[J];                // site-level (occupancy, detection)
  real<lower=-1,upper=1> rho_ab;  // correlation of (occupancy, detection)
  vector<lower=0>[2] sigma_ab;    // sd of (occupancy, detection)

}

transformed parameters {
  real logit_psi[N, J];  // log odds  of occurrence
  real logit_theta[N, J];  // log odds of detection
  
  for (i in 1:N){     // loop across all species
      for (j in 1:J){    // loop across all sites
        logit_psi[i, j] = (uv[i, 1] + a1_species_occ * species_cov1[i]) + (ab[j, 1] + a1_site_occ * site_cov1[j]);
      }
  }
  
  for (i in 1:N) {   // loop across all species
    for (j in 1:J) {     // loop across all sites
          logit_theta[i, j] = (uv[i, 2]) + (ab[j, 2]); // no detection covariates
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
  
  sigma_ab ~ cauchy(0, 2.5);
  (rho_ab + 1) / 2 ~ beta(2, 2);
  ab ~ multi_normal(rep_vector(0, 2), cov_matrix_2d(sigma_ab, rho_ab));
  
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
