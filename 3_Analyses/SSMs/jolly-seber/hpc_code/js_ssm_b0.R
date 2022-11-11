##### Description #####
# Overall script to run density-independent (b0) jags models on the HPC

##### Libraries #####
library(tidyverse)
library(jagsUI)

##### Import Data #####
load("/home/s4532001/7008_guppies/data/capture_histories.RData")
river_data <- list(CA_CH, LL_CH, TY_CH, UL_CH)
rivers <- c("CA", "LL", "TY", "UL")
# River 1: Caigual (CA) - 6843
# River 2: Lower La Laja (LL) - 12737
# River 3: Taylor (TY) - 12774
# River 4: Upper La Laja (UL) - 23556

##### Get river choice from array input #####
args <- commandArgs(trailingOnly = FALSE)
riv <- as.integer(args[7])

##### Define jags model for constant phi #####
sink(file = "/home/s4532001/7008_guppies/jags_models/js_ssm_b0.jags")
cat("
  model {
    # priors and constraints
    for (t in 1:(Time-1)) {
      phi_exp[t] <- exp(beta0)
      phi_logit[t] <- phi_exp[t] / (1 + phi_exp[t])
    } #t
  
    # survival probability coefficients
    beta0 ~ dnorm(0, 0.1)
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
        mu2[i, t] <- phi_logit[t-1] * z[i, t-1] + gamma[t] * prod(q[i, 1:(t-1)])
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
##########
##### Define jags function #####

JS_ssm <- function(CH, ni, nt, nb, nc) {
  
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
                    M = dim(CH.aug)[1])
  
  # Initial values
  inits <- function() {list(beta0 = runif(1, -1, 1),
                            p = runif(dim(CH.aug)[2], 0, 1),
                            z = z_prior
  )
  }
  
  # Parameters monitored
  parameters <- c("psi", "p", "beta0", "b", "Nsuper", "N", "B", "gamma")
  
  # # MCMC settings
  # ni <- 120
  # nt <- 1
  # nb <- 20
  # nc <- 3
  
  js_ssm <- jags(data = jags.data, inits = inits, parameters.to.save = parameters,
                  model.file = "/home/s4532001/7008_guppies/jags_models/js_ssm_b0.jags",
                  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                  parallel = TRUE)
  # js_ssm <- autojags(data = jags.data, inits = inits, parameters.to.save = parameters,
  #                    model.file = "/home/s4532001/7008_guppies/jags_models/js_ssm_b0.jags",
  #                    n.chains = nc, n.burnin = nb, n.thin = nt,
  #                    parallel = TRUE, Rhat.limit = 1.001)
  
  return(js_ssm)
}

start <- Sys.time()
js_ssm <- JS_ssm(river_data[[riv]], ni = 60, nt = 1, nb = 40, nc = 3)
end <- Sys.time()
end - start

save(js_ssm, file = 
       paste("/home/s4532001/7008_guppies/output/models/", rivers[riv], "_js_ssm_b0.RData", sep = ""))

