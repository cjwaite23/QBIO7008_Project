##### Description #####
# Script to fit the Jolly-Seber restricted dynamic occupancy
# state space formulation to the data with Stan

# Stan script code obtained from Ito, 2022

##### Libraries #####
library(tidyverse)
library(rstan)
library(gdata)

##### Load data #####
load("3_Analyses/state_space_models/capture_histories.RData")

CA_CH_M <- CA_CH[CA_sex == 0,]
CA_CH_F <- CA_CH[CA_sex == 1,]
TY_CH_M <- TY_CH[TY_sex == 0,]
TY_CH_F <- TY_CH[TY_sex == 1,]
LL_CH_M <- LL_CH[LL_sex == 0,]
LL_CH_F <- LL_CH[LL_sex == 1,]
UL_CH_M <- UL_CH[UL_sex == 0,]
UL_CH_F <- UL_CH[UL_sex == 1,]

##### Stan model 1 #####
# define script separately
##### Fit model #####
JS_ssm_stan <- function(CH, ni, nt, nb, nc) {
  n_occ <- dim(CH)[2]
  aug_N <- round(dim(CH)[1] / 6)
  CH.aug <- rbind(CH, matrix(0, ncol = dim(CH)[2], nrow = aug_N))
  
  dat <- list(y = CH.aug, 
              M = dim(CH.aug)[1], 
              n_occasions = n_occ)
  inits <- lapply(1:nc, function(i)
    list(beta = runif(5, 0, 1),
         p = runif(n_occ, 0, 1),
         phi = runif(n_occ, 0, 1),
         gamma = runif(n_occ, 0, 1)))
  params <- c("psi", "p", "b", "Nsuper", "N", "B", "gamma")
  
  js_occ <- stan("3_Analyses/SSMs/jolly-seber/stan_models/js_ssm.stan",
                 data = dat, init = inits, pars = params,
                 chains = nc, iter = ni, warmup = nb, thin = nt)
  return(js_occ)
}

start <- Sys.time()
TY_F_stan_ssm <- JS_ssm_stan(TY_CH_F, 30, 1, 10, 1)
end <- Sys.time()
end - start

# Model was never able to be run to completion and convergence so no
# saving lines are shown here
