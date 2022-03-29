### Multi-species occupancy model
## Jens Ulrich
## Project for FRST 509 course
## Constructing an occupancy model for multiple species

library(MASS)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(rstan)

## drawing from example here, written by Bob Carpenter: 
## https://mc-stan.org/users/documentation/case-studies/dorazio-royle-occupancy.html


### simulate test data

### Define study dimensions some data
n_species <- 80 # total number of species in the supercommunity
n_sites <- 40 # number of sites
n_surveys <- 12 # number of visits

## specify site-level occupancy and detection params
# here the sites are considered exchangeable and considered to share a single parameter
# for site-level occupancy and site-level detection
alpha <- -1 #  site-level occupancy (logit scaled)
beta <- -1.5 #  site-level detection (logit scaled)

# inv logit transformation function
ilogit <- function(u) { return(1 / (1 + exp(-u))); }

# inv logit transformation function to view transformed values
transformed_alpha <- ilogit(alpha) # site-level occupancy 
transformed_beta <- ilogit(beta) # site-level detection 


## generate covariate values
# site_cov1 is treatment group (0 = restored habitat; 1 = not restored habitat)
site_unit_covars <- data.frame(site_cov1 = rep(0:1, each = n_sites/2)) # could try updating with categorical covariate
# just keep the one cov1 for now
site_cov1 <- site_unit_covars$site_cov1
a1_site_occ <- 1.5 # slope of site covariate on site-level occupancy effect
# b1_site_occ <- 0.5 # slope of site covariate on site-level detection effect
# continuous version
# unit_covars <- data.frame(cov1 = rnorm(n_sites)) # this is a continuous covariate (like flower abundance e.g.)

# (here the sites are considered exchangeable and considered to share a single parameter) that determine species occupancy and detectability.
# if wanted to, could make sites non-exchangeable
#alpha_sigma <- .25 # small / negligible variation in site intercepts (occupancy)
#beta_sigma <- .25 # small / negligible variation in site intercepts (detection)
#rho_ab <- 0 # uncorrelated

## Simulate site-level occupancy and detection if you want them to 
# be non-equivalent
# vector to fill with occurrence and detection rates (site-level)
# ab <- matrix(NA, nrow = n_sites, ncol = 2)

mu_ab <- c(alpha, beta)

sigma_ab <- matrix(NA, nrow = 2, ncol = 2)
sigma_ab[1,1] = (alpha_sigma)^2
sigma_ab[2,2] = (beta_sigma)^2
sigma_ab[1,2] = alpha_sigma * beta_sigma * rho_ab
sigma_ab[2,1] = alpha_sigma * beta_sigma * rho_ab

ab <- MASS::mvrnorm(n = n_sites, mu = mu_ab, Sigma = sigma_ab, empirical = FALSE)
ab
(rho_ab_simmed <- cor(ab[,1], ab[,2])) # correlation of site occupancy and detection


# visualize site parameters
hist(ab[,1], xlim = c(-3,3))
hist(ab[,2], xlim = c(-3,0))
# don't want to create any strong correlation between 
# site occupancy and detection effects
plot(ab[,2] ~ ab[,1])
lm1 <- lm(ab[,2] ~ ab[,1])
abline(lm1) 
rho_ab_simmed <- cor(ab[,1], ab[,2]) # correlation of species occupancy and detection


## simulate species-level parameters
# cov1 is treatment group (0 = restored habitat; 1 = not restored habitat)
num_high_occurrence_species <- 30
species_unit_covars <- data.frame(species_cov1 = rep(c(0, 1), 
                                                     times = c(n_species-num_high_occurrence_species, num_high_occurrence_species)))

# just keep the one cov1 for now
species_cov1 <- species_unit_covars$species_cov1
a1_species_occ <- 3 # slope on occupancy
# b1_species_occ <- .5 # slope on detection

# for each species i, there are a pair of parameters for species-specific 
# occupancy and detection
# now simulate some species level occupancy (u) and detection (v) values

# means for sigma u and sigma v will be centered at 0; these are not 
# estimated by the model written by BC, but in theory could be by tweaking
# the STAN model code if this ended up to be a variable of interest?
u <- 0 # species-level occupancy (logit scaled)
v <- 0 # species-level detection (logit scaled)
sigma_u <- 0.5 # std deviation of species occupancy
sigma_v <- 0.5 # std deviation of species detection
rho_uv <- 0.8

mu_uv <- c(u, v)

sigma_uv <- matrix(NA, nrow = 2, ncol = 2)

sigma_uv[1,1] = (sigma_u)^2
sigma_uv[2,2] = (sigma_v)^2
sigma_uv[1,2] = sigma_u * sigma_v * rho_uv
sigma_uv[2,1] = sigma_u * sigma_v * rho_uv

uv <- MASS::mvrnorm(n = n_species, mu = mu_uv, Sigma = sigma_uv, empirical = FALSE)
uv
(rho_uv_simmed <- cor(uv[,1], uv[,2])) # correlation of site occupancy and detection


hist(uv[,1], xlim = c(-6,6)) # most species at about 50% of sites, some at 10 or 90
hist(uv[,2], xlim = c(-6,6))
plot(uv[,2] ~ uv[,1])
lm2 <- lm(uv[,2] ~ uv[,1])
abline(lm2)

# the output correlation between occupancy and detection is ~.613
rho_uv_simmed <- cor(uv[,1], uv[,2]) # correlation of species occupancy and detection

# and then each survey will need to be simulated for each 
# species across each site

# First simulation step is to simulate if the sites are occupied,
Z <- matrix(NA, nrow = n_species, ncol = n_sites)
# set.seed(666)
for(i in 1:n_species){ # loop across species
  for(j in 1:n_sites){ # loop across sites
    Z[i, j] <- rbinom(1, 1, ilogit((uv[i, 1] + a1_species_occ * species_cov1[i]) + 
                                     (alpha + a1_site_occ * site_cov1[j]))) # 
    # Z[i, j] <- rbinom(1, 1, ilogit(uv_filtered[i,1] + (ab[j,1]))) 
    # simulate occupancy states
    # for species i at site j, draw occupancy from a binomial dist with prob 
    # a[i] + u[i]
  }
}

# and then simulation of the surveys at each site (given simulated occupancy)
det_data_3D <- array(NA, c(n_species, n_sites, n_surveys))
# set.seed(1)
for(i in 1:n_species){ # loop across available species
  for(j in 1:n_sites){ # loop across sites
    for(k in 1:n_surveys){ 
      det_data_3D[i, j, k] <- Z[i, j] * rbinom(1, 1, ilogit(uv[i, 2] + beta)) 
      # det_data_3D[i, j, k] <- Z[i, j] * rbinom(1, 1, ilogit(uv_filtered[i,2] + ab[j,2])) 
    } # multiply by the occupancy state so that you can never detect species
    # if it is not present at the site
  }
}

# Now add the values in the 3D array to get a single 2D matrix of summed detections
det_data_2D <- matrix(NA, nrow = n_species, ncol = n_sites)
for(i in 1:n_species){ # loop across available species
  for(j in 1:n_sites){ # loop across sites
    det_data_2D[i, j] <- sum(det_data_3D[i, j, 1:n_surveys])
  }
}

# now visualize the simulated detection data that is dependent on 
# simulated occupancy states
# the first thing the authors do is reshape and plot the data itself
df_x2 <- melt(det_data_2D)
colnames(df_x2)[1] <- "species"
colnames(df_x2)[2] <- "site"
colnames(df_x2)[3] <- "detections"

head(df_x2)

# plot the data as in BC's example
# I chose more sites and more species to have a higher sample size
# notice 4 rows without any detections
# so only 61 species are detected from my simulated surveys
# but I would like my model to estimate 61 + 4 = 65 species as the total diversity
# that is, E_N & E_N_2 should ~65 and not 61.
# so my model should estimate (with these param values and seed values)
detections_heatmap_plot2 <-
  ggplot(df_x2, aes(as.factor(site), as.factor(species))) +
  geom_tile(aes(fill = detections), colour = "white") +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "site number", y = "species number") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Detections of Species at Sites over Visits")

plot(detections_heatmap_plot2)
print(row_sums <- rowSums(det_data_2D))

# need to remove rows with no detections from simulated data in det_data_2D 
# for model to run, (those 4 species that are not recovered)
det_data_2D_df <- as.data.frame(det_data_2D) 
# det_data_2D_df <- det_data_2D_df[apply(det_data_2D_df[,-1], 1, function(x) !all(x==0)),]

det_data_2D_matrix <- as.matrix(det_data_2D_df)

# data to feed to the model
X <- det_data_2D_matrix # detection data
N <- n_species # _ species detected (but there are _ available)
J <- n_sites # number of sites - set at 20
K <- n_surveys # number of surveys - set at 15
species_cov1 <- species_cov1
site_cov1 <- site_cov1

my_dat <- c("X", "N", "J", "K", "species_cov1", "site_cov1")
param_names <- c("a1_species_occ",
                "a1_site_occ",
                "alpha", "beta",
                "var sp occ", "var sp det", "rho_uv")
parameters <- c(a1_species_occ,
                a1_site_occ,
                alpha, beta,
                sigma_u, sigma_v, rho_uv_simmed)
targets <- as.data.frame(cbind(param_names, parameters))

stan_model <- stan_model("./multispecies_occupancy_with_covariate_species_and_site.stan")

# keep iter and chains small just for purposes of making sure the model runs
n_iter <- 10000
n_chains <- 4
n_cores <- 4
sim_fit <- sampling(stan_model, data = my_dat, 
                    iter = n_iter, warmup = n_iter*0.6, 
                    chains = n_chains, cores = n_cores, 
                    control=list(adapt_delta=0.975))

options(width="120")
print(sim_fit, c("a1_species_occ",
                  "a1_site_occ",
                  "alpha", "beta",
                  "sigma_uv", "rho_uv", 
                  "lp__"))
#"uv[1,1]", "uv[1,2]", "uv[25,1]", "uv[25,2]"))
# alpha should be 1.25
# beta should be -1.75
# Omega should be 0.6
# sigma_uv[1] (species variance in occupancy) should be 1.0
# sigma_uv[2] (species variance in detection) should be 1.25
# rho_uv - this is where I'm getting thrown off and having a hard time in 
# simulating the data (causing differences in my simulated dataset versus butterflies ex)
# It should be the case that more abundant species (by occupancy)
# should correspondingly be easier to detect
# E_N_2 should be 40 - the number of species simulated as "available"
launch_shinystan(sim_fit)

traceplot(sim_fit,
          c("a1_species_occ",
            "a1_site_occ",
            "alpha", "beta",
            "sigma_uv", "rho_uv", 
            "uv[20,1]", "uv[20,2]", "uv[53,1]", "uv[53,2]",
            "lp__"),
          inc_warmup=FALSE) +
  coord_cartesian(xlim = c(9000, 10000)) +
  scale_x_continuous(breaks=c(9000, 9500, 10000))

traceplot(sim_fit, "uv", inc_warmup=TRUE) +
  coord_cartesian(xlim = c(1, 50))

pairs(sim_fit, pars = c("a1_species_occ",
                        "a1_site_occ",
                        "alpha", "beta",
                        "sigma_uv", "rho_uv"))




### Real data
## They apply their model to estimating total numbers of butterfly species 
x_csv <- read.csv("./clean_real_data_2.csv", row.names = 1)
n <- dim(x_csv)[1] #  distinct butterfly species detected (28)
J <- dim(x_csv)[2] # distinct sites visited (20) # each row is a site with number of detections
x <- matrix(NA, n, J) # _ matrix where _ is the number of visits in which species i was observed at site j
for (i in 1:n) {
  for (j in 1:J) {
    x[i,j] <- x_csv[i,j]
  }
}
K <- 6  # visits to each site (18) (read from paper figure 5)
S <- 30 # size of supercommunity out of which community is drawn (50) # more than adequate judging from plot in figure 6
# unit_covars <- data.frame(cov1 = rep(0:1, each = J/2)) # could try updating with categorical covariate
head(x_csv)

x_csv <- x_csv %>%
  relocate(balaclava, bobolink, gordon, kensington, moberly, quilchena,
           rupert, slocan, winona, locarno, falaise, china.creek, killarney,
           memorial.south, oak.meadows, prince.of.wales, queen.elizabeth,
           west.memorial)
# relocate so that all the control sites are on the left and restored on the right
# (not part of the model just for visualization)

# visualize the real detection matrix
df_x <- melt(x)
colnames(df_x)[1] <- "species"
colnames(df_x)[2] <- "site"
colnames(df_x)[3] <- "detections"

head(df_x)

detections_heatmap_plot <-
  ggplot(df_x, aes(as.factor(site), as.factor(species))) +
  geom_tile(aes(fill = detections), colour = "white") +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "site number", y = "species number") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Detections of Species at Sites over Visits")

plot(detections_heatmap_plot)

site_unit_covars <- data.frame(site_cov1 = rep(0:1, each = J/2)) # could try updating with categorical covariate
# just keep the one cov1 for now
site_cov1 <- site_unit_covars$site_cov1

num_high_occurrence_species <- 6
species_unit_covars <- data.frame(species_cov1 = rep(c(1, 0), 
                                                     times = c(num_high_occurrence_species, n-num_high_occurrence_species)))

# just keep the one cov1 for now
species_cov1 <- species_unit_covars$species_cov1

X <- x_csv
N <- nrow(x) # 15 species detected (but there are maybe more available)
J <- ncol(x)
K <- 6
S <- 30
species_cov1 <- species_cov1
site_cov1 <- site_cov1

my_dat <- c("X", "N", "J", "K", "S", "species_cov1", "site_cov1")

stan_model <- stan_model("./multispecies_occupancy_with_covariate_species_and_site.stan")

# keep iter and chains small just for purposes of making sure the model runs
n_iter <- 12000
n_chains <- 4
n_cores <- 4
sim_fit <- sampling(stan_model, data = my_dat, 
                    iter = n_iter, warmup = n_iter*0.6, 
                    chains = n_chains, cores = n_cores, 
                    control=list(adapt_delta=0.975))

options(width="120")
print(sim_fit, c("a1_species_occ",
                 "a1_site_occ",
                 "alpha", "beta",
                 "sigma_uv", "rho_uv", 
                 "lp__"))
#"uv[1,1]", "uv[1,2]", "uv[25,1]", "uv[25,2]"))
# alpha should be 1.25
# beta should be -1.75
# Omega should be 0.6
# sigma_uv[1] (species variance in occupancy) should be 1.0
# sigma_uv[2] (species variance in detection) should be 1.25
# rho_uv - this is where I'm getting thrown off and having a hard time in 
# simulating the data (causing differences in my simulated dataset versus butterflies ex)
# It should be the case that more abundant species (by occupancy)
# should correspondingly be easier to detect
# E_N_2 should be 40 - the number of species simulated as "available"
launch_shinystan(sim_fit)

traceplot(sim_fit,
          c("a1_species_occ",
            "a1_site_occ",
            "alpha", "beta",
            "sigma_uv", "rho_uv", 
            "uv[1,1]", "uv[1,2]", "uv[10,1]", "uv[10,2]",
            "lp__"),
          inc_warmup=FALSE) +
  coord_cartesian(xlim = c(9000, 10000)) +
  scale_x_continuous(breaks=c(9000, 9500, 10000))

traceplot(sim_fit, "uv", inc_warmup=TRUE) +
  coord_cartesian(xlim = c(1, 50))

pairs(sim_fit, pars = c("a1_species_occ",
                        "a1_site_occ",
                        "alpha", "beta",
                        "sigma_uv", "rho_uv"))
