##### Description #####
# Script to fit the Jolly-Seber restricted dynamic occupancy
# state space formulation to the data with JAGS

# JAGS script code obtained from Kery & Schaub 2011

##### Libraries #####
library(tidyverse)
library(jagsUI)
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

##### JAGS model  ######
##### Defining jags normal model #####
sink(file = "3_Analyses/SSMs/jolly-seber/jags_models/js_ssm.jags")
cat("
  model {of 
    # priors and constraints
    # survival probability coefficients
    for (t in 1:(Time - 1)) {
      phi[t] ~ dunif(0, 1)
    }
    # capture probabilities
    for (t in 1:Time) {
      p[t] ~ dunif(0, 1)
      gamma[t] ~ dunif(0, 1)
    }
    
    # Likelihood
    for (i in 1:M) {
      # First occasion
      # State process
      z[i,1] ~ dbern(gamma[1])
      mu1[i] <- z[i, 1] * p[1]
      # Observation process
      y[i, 1] ~ dbern(mu1[i])
      # Subsequent occasions
      for (t in 2:Time) {
        # Stateprocess
        q[i, t-1] <- 1 - z[i, t-1] #Availability for recruitment
        mu2[i, t] <- phi[t-1] * z[i, t-1] + gamma[t] * prod(q[i, 1:(t-1)])
        z[i, t] ~ dbern(mu2[i, t])
        # Observation process
        mu3[i, t] <- z[i, t] * p[t]
        y[i, t] ~ dbern(mu3[i, t])
      } #t
    } #i
    
    # Calculate derived population parameters
    for (t in 1:Time) {
      qgamma[t] <- 1 - gamma[t]
    }
    cprob[1] <- gamma[1]
    for (t in 2:Time) {
      cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
    } #t
    psi <- sum(cprob[]) # Inclusion probability
    for (t in 1:Time) {
      b[t] <- cprob[t] / psi # Entry probability
    } #t

    for (i in 1:M) {
      recruit[i, 1] <- z[i, 1]
      for (t in 2:Time) {
        recruit[i, t] <- (1 - z[i, t-1]) * z[i, t]
      } #t
    } #i

    for (t in 1:Time) {
      N[t] <- sum(z[1:M, t])  # Actual population size
      B[t] <- sum(recruit[1:M, t])  #Number of entries
    } #t
    
    for (i in 1:M) {
      Nind[i] <- sum(z[i, 1:Time])
      Nalive[i] <- 1 - equals(Nind[i], 0)
    } #i

    Nsuper <- sum(Nalive[]) # Superpopulation size
  }
",fill=TRUE)
sink()

##### Fit model #####
JS_ssm <- function(CH, ni, nt, nb, nc) {
  
  n_occ <- dim(CH)[2]
  #Create augmented capture-history matrix with aug_N all-zero rows
  aug_N <- round(dim(CH)[1] / 6)
  CH.aug <- rbind(CH, matrix(0, ncol = dim(CH)[2], nrow = aug_N))
  
  z_prior <- CH.aug
  for (row in 1:dim(z_prior)[1]) {
    if (length(which(z_prior[row,] == 1)) > 0) {
      first <- min(which(z_prior[row,] == 1))
      last <- max(which(z_prior[row,] == 1))
      z_prior[row, first:last] <- 1
    }
  }
  
  # Bundle data
  jags.data <- list(y = CH.aug,
                    Time = dim(CH.aug)[2],
                    M = dim(CH.aug)[1]
  )
  
  # Initial values
  inits <- function() {list(phi = runif(dim(CH.aug)[2] - 1, 0, 1),
                            p = runif(dim(CH.aug)[2], 0, 1),
                            z = z_prior
  )
  }
  
  # Parameters monitored
  parameters <- c("psi", "NSuper", "p", "phi", "b", "N", "B", "gamma")
  
  mod_file = "3_Analyses/SSMs/jolly-seber/jags_models/js_ssm.jags"
  
  js_ssm <- jags(data = jags.data, inits = inits, parameters.to.save = parameters,
                 model.file = mod_file,
                 n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                 parallel = TRUE)
  return(js_ssm)
}

##### Run the data
start <- Sys.time()
TY_F_jags_ssm <- JS_ssm(TY_CH_F, 30, 1, 10, 1)
end <- Sys.time()
end - start

# Model was never able to be run to completion and convergence so no
# saving lines are shown here