##### Description #####
# Code to explore using stan models for population growth state-space analysis

##### Libraries #####
library(tidyverse)
library(rstan)

##### Options #####
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

##### Load data #####
load("4_Results/MARK_analyses/rmark_best_models/CA_B_derived.RData")
load("4_Results/MARK_analyses/rmark_best_models/TY_B_derived.RData")
load("4_Results/MARK_analyses/rmark_best_models/UL_B_derived.RData")
load("4_Results/MARK_analyses/rmark_best_models/LL_B_derived.RData")
load("2_Data_manipulation/data_files/climate_data.RData")

##### Deal with climate data #####
#CA begins in March 2008 - month 13
n_occ <- dim(CA_B_derived)[1] / 2
dry_wet = c(1,1,1,0,0,0,0,0,0,0,1,1) #1 = dry season, begins at March
mar_08_col <- which(colnames(clim_data) == "Mar_2008")
jags_clim_data <- clim_data[,mar_08_col:(mar_08_col + n_occ-1)] %>% 
  as_tibble() %>%
  rbind(drywet = rep(dry_wet, times = n_occ %/% 12 + 1)[1:n_occ]) %>%
  as.matrix() %>% t() %>% as.data.frame() %>%
  mutate(V1 = V1 / max(V1), V2 = V2 / max(V2), V3 = V3 / max(V3))
colnames(jags_clim_data) <- c("prec", "tmin", "tmax", "dry")

##### Function to run stan

ssm_stan_simple <- function (data, clim, sex, update = FALSE, model = NULL) {
  
  ni <- 2000    #iterations
  nb <- 1000    #burn-in
  nc <- 3         #chains
  nt <- 1         #thinning
  
  data <- filter(data, Sex == sex)
  n_occ <- dim(data)[1]
  
  stan.data <- list(logy = log(data$N), # log population size
                    Time = n_occ, # number of occasions
                    sigma_obs = max(data$N_se / data$N), # observation standard error for each time period
                    prec = clim$prec, # average precipitation
                    tmin = clim$tmin, # average min temp
                    tmax = clim$tmax, # average max temp
                    dry = clim$dry) # 1 if dry season, 0 if wet
  
  inits <- lapply(1:nc, function(i)
    list(sigma_proc = runif(1, 0, 1),
         #sigma_obs = runif(1, 0, 1),
         b = rnorm(1, 1, 0.1),
         logN_est = c(rnorm(1, stan.data$logy[1], 0.01),
                      rep(NA, length(data$Occasion) - 1))
         )
    )
  
  parameters <- c("b",
                  #"c",
                  "sigma2_proc",
                  #"sigma2_obs",
                  "N_est")
  
  ssm <- stan("3_Analyses/state_space_models/SSM_simple_stan/ssm_simple.stan",
              data = stan.data, init = inits, pars = parameters,
              chains = nc, iter = ni, warmup = nb, thin = nt,
              seed = 1, cores = 3)

  return(ssm)
}

ssm_stan_dd <- function (data, clim, sex, update = FALSE, model = NULL) {
  
  ni <- 2000    #iterations
  nb <- 1000    #burn-in
  nc <- 3         #chains
  nt <- 1         #thinning
  
  data <- filter(data, Sex == sex)
  n_occ <- dim(data)[1]
  
  stan.data <- list(logy = log(data$N), # log population size
                    Time = n_occ, # number of occasions
                    sigma_obs = max(data$N_se / data$N), # observation standard error for each time period
                    prec = clim$prec, # average precipitation
                    tmin = clim$tmin, # average min temp
                    tmax = clim$tmax, # average max temp
                    dry = clim$dry) # 1 if dry season, 0 if wet
  
  inits <- lapply(1:nc, function(i)
    list(sigma_proc = runif(1, 0, 1),
         #sigma_obs = runif(1, 0, 1),
         b = rnorm(1, 1, 0.1),
         dd = rgamma(1, 0.25, 0.25),
         logN_est = c(rnorm(1, stan.data$logy[1], 0.01),
                      rep(NA, length(data$Occasion) - 1))
    )
  )
  
  parameters <- c("b",
                  "dd",
                  "sigma2_proc",
                  #"sigma2_obs",
                  "N_est")
  
  ssm <- stan("3_Analyses/state_space_models/SSM_simple_stan/ssm_dd.stan",
              data = stan.data, init = inits, pars = parameters,
              chains = nc, iter = ni, warmup = nb, thin = nt,
              seed = 1, cores = 3)
 
  return(ssm)
}

# Test model on CA and M - does not run in time though
CA_M_ssm_stan_simple <- ssm_stan_simple(CA_B_derived, jags_clim_data, "M")
CA_M_ssm_stan_dd <- ssm_stan_dd(CA_B_derived, jags_clim_data, "M")

  
  
  
  





