###-----------------------------------------------------------------------------
### simulate some test data

source("./multispecies_occupancy_model/sim_dynocc.R")

## Define study dimensions 
n_sites <- 50 # number of sites
sites <- seq(1:n_sites)
n_species <- 30
species <- seq(1:n_species)
n_surveys <- 5 # number of repeated visits to each site
n_years <- 5
n_years_minus1 <- n_years - 1
years <- seq(1:n_years)
years_minus1 <- seq(1:n_years_minus1)

## Define simulation params
psi1_0 <- 0 #  site-level initial occupancy rate (logit scaled)
sigma_psi1_species <- 1.5 # variance among species
gamma0 <- -2 #  site-level colonization rate (logit scaled)
sigma_gamma_species <- 1 # variance among species
phi0 <- 2 #  site-level persistence rate (logit scaled)
sigma_phi_species <- 1 # variance among species
p0 <- -1 #  site-level detection rate (logit scaled)
sigma_p_species <- 1 # variance among species

my_data <- simulate_data(n_sites,
  sites,
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
)

y <- my_data$y

# for the binomial model, sum y across all surveys:
# For row, the sum or mean is over dimensions dims+1".
y <- rowSums(y, dims = 3)

## --------------------------------------------------
### Prep data and tweak model options

stan_data <- c("y", 
               "species", "sites", 
               "n_species", "n_sites", "n_years", "n_years_minus1", "n_surveys"
) 

## Parameters monitored 
params <- c("psi1_0", 
            "sigma_psi1_species",

            "gamma0", 
            "sigma_gamma_species",
            
            "phi0", 
            "sigma_phi_species",
            
            "p0", 
            "sigma_p_species"
)

# MCMC settings
n_iterations <- 300
n_thin <- 1
n_burnin <- 150
n_chains <- 4
n_cores <- n_chains
delta = 0.95

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(psi1_0 = runif(1, -1, 1),
       sigma_psi1_species = runif(1, 0, 1),
       gamma0 = runif(1, -2, 0),
       sigma_gamma_species = runif(1, 0, 1),
       phi0 = runif(1, 2, 3),
       sigma_phi_species = runif(1, 0, 1),
       p0 = runif(1, -1, 1),
       sigma_p_species = runif(1, 0, 1)
  )
)

# targets
parameter_values <-  c(
  psi1_0, 
  sigma_psi1_species,

  gamma0, 
  sigma_gamma_species,
  
  phi0, 
  sigma_phi_species,
  
  p0, 
  sigma_p_species
  
)

targets <- as.data.frame(cbind(params, parameter_values))

## --------------------------------------------------
### Run model

stan_model <- "./multispecies_occupancy_model/multispecies_dynocc.stan"

## Call Stan from R
stan_out <- stan(stan_model,
                 data = stan_data, 
                 init = inits, 
                 pars = params,
                 chains = n_chains, iter = n_iterations, 
                 warmup = n_burnin, thin = n_thin,
                 seed = 1,
                 control=list(adapt_delta=delta),
                 open_progress = FALSE,
                 cores = n_cores)

print(stan_out, digits = 3, 
      pars = c("psi1_0", 
               "sigma_psi1_species",
               
               "gamma0", 
               "sigma_gamma_species",
               
               "phi0", 
               "sigma_phi_species",
               
               "p0", 
               "sigma_p_species"
      ))
