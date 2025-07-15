### single species static occupancy model
## Jens Ulrich, Feb 2025

library(tidyverse)

# build simulation function

simulate_data <- function(n_sites,
                          site,
                          n_species,
                          species,
                          n_surveys,
                          n_years,
                          n_years_minus1,
                          years,
                          years_minus1,
                          psi1_0,
                          sigma_psi1_species,
                          gamma0,
                          sigma_gamma_species,
                          phi0,
                          sigma_phi_species,
                          p0,
                          sigma_p_species
){
  
  ## ilogit and logit functions
  ilogit <- function(x) exp(x)/(1+exp(x))
  logit <- function(x) log(x/(1-x))
  
  ## predictor center scaling function
  center_scale <- function(x) {
    (x - mean(x)) / sd(x)
  }
  
  # prepare arrays for z and y
  z <- array(NA, dim = c(n_species, n_sites, n_years)) # latent presence/absence
  y <- array(NA, dim = c(n_species, n_sites, n_years, n_surveys)) # observed data
  
  ## --------------------------------------------------
  ## Create random effects
  
  ## species-specific random intercepts

  psi1_species <- rnorm(n=n_species, mean=0, sd=sigma_psi1_species)
  
  gamma_species <- rnorm(n=n_species, mean=0, sd=sigma_gamma_species)
  
  phi_species <- rnorm(n=n_species, mean=0, sd=sigma_phi_species)
  
  p_species <- rnorm(n=n_species, mean=0, sd=sigma_p_species)
  
  ## --------------------------------------------------
  ## Create expected values
  
  # generate p with heterogeneity
  logit_p <- array(NA, dim = c(n_species, n_sites, n_years, n_surveys)) 
  
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_years){
        for(l in 1:n_surveys){
          
          logit_p[i,j,k,l] = 
            p0 +
            p_species[i] 
        }
      }
    }  
  }
  
  
  # generate ecological expected values with heterogeneity
  logit_psi1 <- array(NA, dim = c(n_species, n_sites)) 
  logit_gamma <- array(NA, dim = c(n_species, n_sites, n_years_minus1)) 
  logit_phi <- array(NA, dim = c(n_species, n_sites, n_years_minus1)) 
  
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_years_minus1){
        
        logit_psi1[i,j] = 
          psi1_0 +
          psi1_species[i] 
        
        logit_gamma[i,j,k] = # gamma for transition (starting for between years 1 and 2)
          gamma0 +
          gamma_species[i] 
        
        logit_phi[i,j,k] = 
          phi0 +
          phi_species[i]
        
      }
    }
  }
  
  
  
  # generate initial presence/absence states
  for(i in 1:n_species){
    for(j in 1:n_sites){
      z[i,j,1] <- rbinom(n=1, size=1, prob=ilogit(logit_psi1[i,j])) 
    }
  }
  
  true_occupancy <- apply(z[,,],c(1,3),sum ) / n_sites # true occupancy proportion in year 1     
  
  # generate presence/absence in subsequent years
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 2:n_years){
        
        # use z as a switch so we are estimating 
        exp_z <- z[i,j,k-1] * ilogit(logit_phi[i,j,k-1]) + # survival if z=1
          (1 - z[i,j,k-1]) * ilogit(logit_gamma[i,j,k-1]) # or colonization if z=0
        
        # and then assign z stochastically
        # some sites may transition if they are colonized or local extinction occurs
        # but might otherwise retain their state across years
        # look at year starting with 2, for the first transition (phi or gamma[,,1])
        z[i,j,k] <- rbinom(n = 1, size = 1, prob = exp_z) 
        
      }
    }    
  }
  
  # detection / non-detection data
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_years){
        for(l in 1:n_surveys){
          y[i,j,k,l] <- rbinom(n = 1, size = 1, prob = z[i,j,k]*ilogit(logit_p[i,j,k,l]))
        }
      }
    }
  }
  
  y ; str(y)
  
  sum((y/n_surveys) / sum(z)) # proportion of times detection given presence
  
  # which species never occurred
  species_not_occurring <- length(which(apply(true_occupancy, 1, sum) == 0))
  print(paste0("Number of species never occurring: ", species_not_occurring))
  
  # which species were never detected (even though they occurred)
  species_not_observed <- length(which(apply(y, 1, sum, na.rm=TRUE) == 0))
  print(paste0("Number of species occurring but never observed: ", species_not_observed))

  ## --------------------------------------------------
  # Return stuff
  return(list(
    y = y
  ))
  
} # end simulation function
