##### Description #####
# exploring the application of the Kalman Filter to the Data

##### Libraries #####
library(tidyverse)
library(dlm)

#### Load Data #####
load("4_Results/MARK_analyses/rmark_best_models/CA_B_derived.RData")
load("4_Results/MARK_analyses/rmark_best_models/TY_B_derived.RData")
load("4_Results/MARK_analyses/rmark_best_models/UL_B_derived.RData")
load("4_Results/MARK_analyses/rmark_best_models/LL_B_derived.RData")
n_occ85 <- dim(TY_B_derived)[1] / 2
n_occ97 <- dim(UL_B_derived)[1] / 2

##### Define function to modify the data and fit the Kalman filter #####

Kf <- function(data, sex, n_occ) {
  data <- filter(data, Sex == sex) %>%
    mutate(logN = log(N),
           logN_se = N_se/N + 0.0001 * (N_se/N == 0)) 
  
  #Step 1: Create a function that defines the model as a function of the parameter values to estimate.
  # m0 = initial value of the state, w[0]
  # C0 = variance of the initial value
  # GG = matrix multiplier in the State Eq: w[t] = GG*w[t-1] + epsilon[t]
  # W = state process variance
  # FF = matrix multiplier in the Obs. Eq: g[t] = FF*w[t] + eta[t]
  # V = observation error variance
  m1_dlm <- function(theta) {
    dlm(m0 = data$N[1], 
        C0 = 0, 
        GG = theta[1], 
        W = exp(theta[2]), 
        FF = 1,
        V = exp(theta[3])
    )
  }

  #Fit the MLE model
  m1_dlm_fit <- dlmMLE(y = data$N[-1],
                       parm = c(0,0,0),
                       build = m1_dlm)
  
  # Extract the parameters
  m1_dlm_par <- m1_dlm(m1_dlm_fit$par)
  # Use Kalman filter to estimate true values
  m1_dlm_filter <- dlmFilter(y = data$N[-1], mod = m1_dlm_par)
  # Obtain predicted pop sizes
  data$N_pred <- (m1_dlm_filter$m)
  # Plot data
  plot <- ggplot(data = data) +
    geom_point(aes(x = Occasion, y = N), alpha = 0.7) +
    geom_line(aes(x = Occasion, y = N), alpha = 0.7) +
    geom_point(aes(x = Occasion, y = N_pred), col = "blue", alpha = 0.7, shape = 2) +
    geom_line(aes(x = Occasion, y = N_pred), col = "blue", alpha = 0.7, linetype = 3)
  return(list(pars = m1_dlm_par, plot = plot))
}

# Specify stream, sex and occasions
data <- CA_B_derived
sex = "M"
n_occ <- n_occ85
# Actually run the model
Kf(UL_B_derived, "F", n_occ85)
