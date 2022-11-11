##### Description #####
# Overall script to run jags models on the HPC

##### Libraries #####
library(tidyverse)
library(jagsUI)

##### Import Data #####
load("/home/s4532001/7008_guppies/data/capture_histories.RData")
rivers <- c("CA", "LL", "TY", "UL")
river_data <- list(CA_CH, LL_CH, TY_CH, UL_CH)
start_month <- c("Mar_2008", "Mar_2007", "Mar_2007", "Mar_2008")
# River 1: Caigual (CA) - 6843
# River 2: Lower La Laja (LL) - 12737
# River 3: Taylor (TY) - 12774
# River 4: Upper La Laja (UL) - 23556

#load("1_Data/worldclim_data/climate_data.RData")

##### Get river choice from array input #####
args <- commandArgs(trailingOnly = FALSE)
riv <- as.integer(args[7])

##### Define jags model for full specd phi #####
sink(file = "/home/s4532001/7008_guppies/jags_models/js_ssm_phit.jags")
cat("
  model {
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

##########
##### Define jags model for full speed phi NOT REALLY THAT FAST #####
sink(file = "/home/s4532001/7008_guppies/jags_models/js_ssm_phit2.jags")
cat("
  model {
    # priors and constraints
    # survival probability coefficients
    for (t in 1:T_less_one) {
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
      q[1:T_less_one] <- 1 - z[i, 1:T_less_one] #Availability for recruitment
      r[1:T_less_one] <- phi[1:T_less_one] * z[i, 1:T_less_one]
      for (t in 2:Time) {
        # Stateprocess
        z[i, t] ~ dbern(r[t-1] + gamma[t] * prod(q[1:(t-1)]))
        # Observation process
        y[i, t] ~ dbern(z[i, t] * p[t])
      } #t
    } #i
    
    # Calculate derived population parameters
    qgamma <- 1 - gamma
    
    cprob[1] <- gamma[1]
    for (t in 2:Time) {
      cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
    } #t
    psi <- sum(cprob[]) # Inclusion probability
    b <- cprob / psi # Entry probability

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
  }
",fill=TRUE)
sink()
##########
##### Define jags function #####

JS_ssm <- function(riv, ni, nt, nb, nc, update = FALSE, obj = NULL) {
  
  river <- rivers[riv]
  CH <- river_data[[riv]][1:3206, 1:48]
  n_occ <- dim(CH)[1]
  # #Deal with climate data and start months
  # dry_wet = c(1,1,1,0,0,0,0,0,0,0,1,1) #1 = dry season, begins at March
  # start_col <- which(colnames(clim_data) == start_month[riv])
  # clim <- clim_data[,start_col:(start_col + n_occ-1)] %>% 
  #   as_tibble() %>%
  #   rbind(drywet = rep(dry_wet, times = n_occ %/% 12 + 1)[1:n_occ]) %>%
  #   as.matrix() %>% t() %>% as.data.frame()
  # colnames(clim) <- c("prec", "tmin", "tmax", "dry")
  
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
                    M = dim(CH.aug)[1],
                    T_less_one = dim(CH.aug)[2] - 1
  )
  
  # Initial values
  inits <- function() {list(phi = runif(dim(CH.aug)[2] - 1, 0, 1),
                            p = runif(dim(CH.aug)[2], 0, 1),
                            z = z_prior
  )
  }
  
  # Parameters monitored
  parameters <- c("psi", "p", "phi", "b", "N", "B", "gamma")
  
  mod_file = "/home/s4532001/7008_guppies/jags_models/js_ssm_phit.jags"
  
  if (!update) {
    js_ssm <- jags(data = jags.data, inits = inits, parameters.to.save = parameters,
                   model.file = mod_file,
                   n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                   parallel = TRUE)
  } else {
      js_ssm <- jagsUI::update(object = obj, parameters.to.save = parameters,
                       n.iter = ni, n.thin = nt,
                       parallel = TRUE)
  }
  return(js_ssm)
}

start <- Sys.time()
js_ssm <- JS_ssm(riv, ni = 6000, nt = 1, nb = 5000, nc = 3)
save(js_ssm, file = 
       paste("/home/s4532001/7008_guppies/output/models/", rivers[riv], "_js_ssm_phit_48_base.RData", sep = ""))
while (all(unlist(js_ssm$Rhat) < 1.001)) {
  js_ssm <- JS_ssm(riv, ni = 1000, nt = 1, nb = 5000, nc = 3, 
                   update = TRUE, obj = js_ssm)
  save(js_ssm, file = 
         paste("/home/s4532001/7008_guppies/output/models/", rivers[riv], "_js_ssm_phit_48_upd.RData", sep = ""))
}
end <- Sys.time()
end - start

