### single species static occupancy model
## Jens Ulrich, Feb 2025

library(rstan)

###-----------------------------------------------------------------------------
### simulate some test data

## Define study dimensions 
n_sites <- 100 # number of sites
n_surveys <- 5 # number of repeated visits to each site

# inv logit transformation function
ilogit <- function(u) { return(1 / (1 + exp(-u))); }

## specify site-level occupancy and detection params
# here the sites are considered exchangeable and considered to share a single parameter
# for site-level occupancy and site-level detection
alpha <- 0 #  site-level occupancy rate (logit scaled)
beta <- -1 #  site-level detection rate (logit scaled)

# inv logit transformation function to view transformed values (percent odds of occurence or detection)
transformed_alpha <- ilogit(alpha) # site-level occupancy 
transformed_beta <- ilogit(beta) # site-level detection 

# (here the sites are considered exchangeable and considered to share a single parameter) 
# that determines species occupancy and a single param that determines species detectability.
# one could change these things by adding predictor terms to psi and p 

# occupancy linear predictor
psi <- alpha # add params here if you'd like more than just an intercept!
# e.g. 
# psi <- vector(length = n_sites)
# for(i in 1:n_sites){psi[i] <- alpha + alpha1*covariate1[i]} 

# detection linear predictor
p <- beta # add params here if you'd like more than just an intercept!
# e.g. 
# p <- vector(length = n_sites)
# for(i in 1:n_sites){p[i] <- beta + beta1*covariate1[i]} 

## First simulation step is to simulate if the sites are occupied,
Z <- matrix(NA, nrow = n_sites)

for(i in 1:n_sites){ # loop across sites
    Z[i] <- rbinom(1, 1, ilogit(psi)) # 
    # simulate occupancy states
    # for each site j, draw occupancy from a binomial dist with prob 
}

## and then simulate detections on surveys of each site (given simulated occupancy)
det_data_array <- array(NA, c(n_sites, n_surveys))

for(i in 1:n_sites){ # loop across sites
    for(j in 1:n_surveys){ # for each site loop across surveys
      det_data_array[i, j] <- Z[i] * rbinom(1, 1, ilogit(p)) 
      # multiply by the occupancy state so that you can never detect species
      # if it is not present at the site
      # only take 1 rbinom draw per survey
    } 
}

# Now add the values in the  array to get a single vector of summed detections
# i.e., across all surveys to a site, how many successfully detected the species?
det_data_vector <- vector(length = n_sites)

for(i in 1:n_sites){ # loop across all sites
  det_data_vector[i] <- sum(det_data_array[i, 1:n_surveys]) # and sum the detections
}

###-----------------------------------------------------------------------------
### Run the model

# prepare the data to feed to the model
detections <- det_data_vector # detection data
n_sites <- n_sites # number of sites 
n_surveys <- n_surveys # number of surveys
my_data <- c("detections", "n_sites", "n_surveys") # data
# params to monitor
param_names <- c("alpha", "beta") 

# 'targets' allows you to later check agreement between model estimates and the true param values
parameters <- c(alpha, beta) # true values used to simulate the data
targets <- as.data.frame(cbind(param_names, parameters)) 

# load the stan model
stan_model <- stan_model("./simple_occupancy_model.stan")

# specify model run time and number of chains
params = param_names
n_iter <- 4000
n_chains <- 4
n_cores <- 4

# fit the model
sim_fit <- sampling(stan_model, data = my_data, 
                    pars = params,
                    iter = n_iter, warmup = n_iter*0.5, 
                    chains = n_chains, cores = n_cores)

# view outputs and some simple model fitting diagnostics
options(width="120")

print(sim_fit, c("alpha", "beta"))

traceplot(sim_fit, c("alpha", "beta")) 

pairs(sim_fit, pars = c("alpha", "beta"))

###-----------------------------------------------------------------------------
### Compare parameter estimates to parameter targets

library(ggplot2)
library(tidyverse)

fit_summary <- rstan::summary(sim_fit)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

X <- as.factor(seq(1:nrow(targets))) # 4 ecological params of interest

estimates_lower_95 <- c(
  fit_summary$summary[1,4], # alpha
  fit_summary$summary[2,4] # beta
)
estimates_upper_95 <- c(
  fit_summary$summary[1,8], # alpha
  fit_summary$summary[2,8] # beta
)
estimates_lower_50 <- c(
  fit_summary$summary[1,5], # alpha
  fit_summary$summary[2,5] # beta
)
estimates_upper_50 <- c(
  fit_summary$summary[1,7], # alpha
  fit_summary$summary[2,7] # beta
)

# bind the true values and the BCI's from the model into a df
df_estimates <- as.data.frame(cbind(X, targets, 
                                    estimates_lower_95, estimates_upper_95,
                                    estimates_lower_50, estimates_upper_50))
df_estimates$parameters <- as.numeric(df_estimates$parameters)

# now plot the true values and BCI's
(p <- ggplot(df_estimates) +
    theme_bw() +
    scale_x_discrete(name="", breaks = c(1, 2),
                     labels=c(bquote(alpha), 
                              bquote(beta)
                     )
    ) +
    scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                       limits = c(-3, 3)) +
    guides(color = guide_legend(title = "")) +
    geom_hline(yintercept = 0, lty = "dashed") +
    theme(legend.text=element_text(size=10),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 20, angle=0, vjust=0),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    coord_flip()
)

# add error bars
p <- p +
  geom_errorbar(aes(x=X, ymin=estimates_lower_95, ymax=estimates_upper_95),
                color="black",width=0.1,size=1,alpha=1) +
  geom_errorbar(aes(x=X, ymin=estimates_lower_50, ymax=estimates_upper_50),
                color="black",width=0,size=3,alpha=1) +
  geom_point(aes(x=X, y=parameters),
             size = 10, alpha = 1, shape = 10, colour = "firebrick2") 

# draw the plot
p

# in the plot (p) red targets indicate the true values of the parameters
# that we used for the simulation
# ideally the BCI's will overlap with these true values!

# the plotting scheme for p isn't set up to automatically update if you
# add more covariates. You'd have to edit the plotting code accordingly.
