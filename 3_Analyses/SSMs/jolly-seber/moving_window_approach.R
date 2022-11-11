##### Description #####
# Script to try implementing my 'moving window' approach in JAGS for the 
# state-space formulation of Jolly-Seber
# This code is limited by language constraints of JAGS - it is unable to consider
# the same variable as data and parameters

##### Libraries #####
library(tidyverse)
library(jagsUI)
library(stringr)

##### Import Data #####
data <- read_csv("~/QBIO7008 Guppy Project/7008_large_files/1_Data/Individual_differences_determine_the_strength_of_ecological_interactions.csv")
names(data)[1] <- "obs"

##### Create capture histories #####
mr_data <- data %>%
  select(obs, stream, sampling, ID, individual_id, sex) %>%
  #create binary "captured" (1) or "not_captured" (0) column
  mutate(capture = 1 - is.na(ID)) %>%
  dplyr::select(-ID)
rm(data)
enc_his <- mr_data %>%
  dplyr::select(-obs) %>%
  pivot_wider(id_cols = c(individual_id, stream, sex),
              names_from = sampling,
              values_from = capture)
rm(mr_data)
rivers <- unique(enc_his$stream)
sexes <- unique(enc_his$sex)

##### Function to run take Sex and Stream and produce data-frame of fixed values #####
fix_zs <- function(Stream, Sex, enc_his) {
  # filter out stream and sex combos
  if (Stream == "CA" | Stream == "TY") {
    CH <- enc_his %>%
      filter(sex == Sex, stream == Stream) %>%
      dplyr::select(-(1:15)) %>%
      as.matrix()
  } else {
    CH <- enc_his %>%
      filter(sex == Sex, stream == Stream) %>%
      dplyr::select(-(1:3)) %>%
      as.matrix()
  }
  n_occ <- dim(CH)[2]
  
  # Produce a matrix with 1's and 0's for fixed z and NA for unfixed
  Fixed <- CH
  for (row in 1:dim(Fixed)[1]) {
      first <- min(which(Fixed[row,] == 1))
      last <- max(which(Fixed[row,] == 1))
      Fixed[row, first:last] <- 1
      if (first == 1) {
        Fixed[row, (last+1):(last+3)] <- NA
      } else if (first == 2) {
        Fixed[row, 1] <- NA
        Fixed[row, (last+1):(last+3)] <- NA
      } else if (first == 3) {
        Fixed[row, 1:2] <- NA
        Fixed[row, (last+1):(last+3)] <- NA
      } else if (last == n_occ - 2) {
        Fixed[row, (first-3):(first-1)] <- NA
        Fixed[row, (n_occ - 1):n_occ] <- NA
      } else if (last == n_occ - 1) {
        Fixed[row, (first-3):(first-1)] <- NA
        Fixed[row, n_occ] <- NA
      } else if (last == n_occ) {
        Fixed[row, (first-3):(first-1)] <- NA
      } else {
        Fixed[row, (first-3):(first-1)] <- NA
        Fixed[row, (last+1):(last+3)] <- NA
      }
  }
  # Turn this matrix into a dataframe of each individual
  Fixed_df <- Fixed %>% as.data.frame() %>%
    mutate(ID = 1:dim(Fixed)[1]) %>%
    pivot_longer(cols = 1:n_occ, names_to = "time", values_to = "fixed") %>% 
    select(ID, time, fixed) %>%
    filter(!is.na(fixed)) %>%
    mutate(time = as.numeric(time) - 12*(Stream == "CA" | Stream == "TY"))
  
  # Augment CH matrix
  aug_N <- round(dim(CH)[1] / 6)
  CH.aug <- rbind(CH, matrix(0, ncol = n_occ, nrow = aug_N))
  
  # Create matrix of initial conditions
  z_inits <- ifelse(is.na(Fixed), 0, NA) %>%
    rbind(matrix(0, ncol = n_occ, nrow = aug_N))

  return(list(CH.aug = CH.aug, Fixed_df = Fixed_df, z_inits = z_inits))
}

##### Write the jags model #####
##### Define jags model for full specd phi #####
sink(file = "3_Analyses/SSMs/jolly-seber/jags_models/js_ssm_window.jags")
cat("
  model {
    # Fix all known z's
    for (x in 1:X) {
      z[fixed_ID[x], fixed_time[x]] <- fixed_value[x]
    }
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

##### Define function to run jags model
JS_ssm <- function(riv, sex, ni, nt, nb, nc) {
  
  data <- fix_zs(riv, sex, enc_his)
  
  # Bundle data
  jags.data <- list(y = data$CH.aug,
                    Time = dim(data$CH.aug)[2],
                    M = dim(data$CH.aug)[1],
                    T_less_one = dim(data$CH.aug)[2] - 1,
                    fixed_ID = data$Fixed_df$ID,
                    fixed_time = data$Fixed_df$time,
                    fixed_value = data$Fixed_df$fixed,
                    X = dim(data$Fixed_df)[1],
                    data1par0 = ifelse(is.na(data$z_inits), 1, 0)
  )
  
  # Initial values
  inits <- function() {list(phi = runif(dim(data$CH.aug)[2] - 1, 0, 1),
                            p = runif(dim(data$CH.aug)[2], 0, 1),
                            z = data$z_inits
  )
  }
  
  # Parameters monitored
  parameters <- c("NSuper", "psi", "p", "phi", "b", "N", "B", "gamma")
  
  mod_file = "3_Analyses/state_space_models/js_jags_window/js_ssm_window.jags"
  
  js_ssm <- jags(data = jags.data, inits = inits, parameters.to.save = parameters,
                   model.file = mod_file,
                   n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                   parallel = TRUE)
  return(js_ssm)
}

start <- Sys.time()
CA_M_ssm <- JS_ssm("CA", "M", 30, 1, 10, 3)
end <- Sys.time()
end - start





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

CA_M_jags_data <- fix_zs("CA", "M", enc_his)

##### Look at prevalence of 10001 or 100001 etc #####
enc_his_str <- mr_data %>%
  dplyr::select(-obs) %>%
  pivot_wider(id_cols = c(individual_id, stream, sex),
              names_from = sampling,
              values_from = capture) %>%
  unite(col = history, 4:100,
        sep = "", remove = TRUE) %>%
  mutate(history = ifelse(stream == "CA" | stream == "TY",
                          substr(history, 25, 109),
                          history),
         n10001 = str_extract(history, "10001"),
         n100001 = str_extract(history, "100001"),
         n1000001 = str_extract(history, "1000001"),
         n10000001 = str_extract(history, "10000001"),
         greater3 = !is.na(n10001) | !is.na(n100001) | !is.na(n1000001),
         greater4 = !is.na(n100001) | !is.na(n1000001))
