##### Description #####
# Running State-Space Population Models on data from MARK using JAGS

##### Libraries #####
library(tidyverse)
library(jagsUI)
library(car)
library(bestglm)
library(visreg)
library(gridExtra)
library(patchwork)

##### Load data for each stream and sex ##### 
# dfs are grouped by stream here
load("4_Results/MARK_analyses/rmark_best_models/CA_B_derived.RData")
load("4_Results/MARK_analyses/rmark_best_models/TY_B_derived.RData")
load("4_Results/MARK_analyses/rmark_best_models/UL_B_derived.RData")
load("4_Results/MARK_analyses/rmark_best_models/LL_B_derived.RData")

##### Load climatic data #####
load("2_Data_manipulation/data_files/climate_data.RData")

##### Manipulate climate data #####
# Climate data needs to be clipped to be the same size as the sampling periods
# for each of the 85 month sites and the 97 month sites
n_occ85 <- dim(TY_B_derived)[1] / 2
n_occ97 <- dim(UL_B_derived)[1] / 2
#CA begins in March 2008 - month 13
dry_wet = c(1,1,1,0,0,0,0,0,0,0,1,1) #1 = dry season, begins at March
mar_08_col <- which(colnames(clim_data) == "Mar_2008")
mar_07_col <- which(colnames(clim_data) == "Mar_2007")

clim85 <- clim_data[,mar_08_col:(mar_08_col + n_occ85-1)] %>% 
  as_tibble() %>%
  rbind(drywet = rep(dry_wet, times = n_occ85 %/% 12 + 1)[1:n_occ85]) %>%
  as.matrix() %>% t() %>% as.data.frame()
colnames(clim85) <- c("prec", "tmin", "tmax", "dry")

clim97 <- clim_data[,mar_07_col:(mar_07_col + n_occ97-1)] %>% 
  as_tibble() %>%
  rbind(drywet = rep(dry_wet, times = n_occ97 %/% 12 + 1)[1:n_occ97]) %>%
  as.matrix() %>% t() %>% as.data.frame()
colnames(clim97) <- c("prec", "tmin", "tmax", "dry")

##### Create JAGS SSM model #####
sink(file = "3_Analyses/SSMs/pop_growth/jags_models/ssm_final.jags")
cat("
  model {
    # Priors and constraints
    
    logN.est[1] ~ dnorm(init_logN, 1)
    b ~ dnorm(1, 0.01)
  
    sigma.proc ~ dunif(0, 2)
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
    
    # sigma.obs ~ dunif(0, 2)
    # sigma2.obs <- pow(sigma.obs, 2)
    tau.obs <- pow(sigma.obs, -2)
    
    # Likelihood
    # State process
    for (t in 1:(Time - 1)) {
      eps[t] ~ dnorm(0, tau.proc)
      logN.est[t+1] <- b + logN.est[t] + eps[t]
    }
    # Observation process
    for (t in 1:Time) {
      logy[t] ~ dnorm(logN.est[t], tau.obs[t])
    }
    
    #Population size on real scale
    for (t in 1:Time) {
      N.est[t] <- exp(logN.est[t])
    }
  }  
  ", fill =TRUE)
sink()

##### Define function to run the jags model #####
ssm_jags <- function (data, sex, n_occ) {
  
  # Filter data by sex
  data <- filter(data, Sex == sex)
  
  # Provide data, time period, initial pop size and observation error as data
  jags.data <- list(logy = log(data$N), # log population size
                    Time = n_occ, # number of occasions
                    init_logN = log(data$N[1]), # initial population size
                    sigma.obs = data$N_se / data$N + 0.0001 * (data$N_se / data$N == 0) # observation standard error for each time period
  )
  
  #function that defines initial conditions.
  #sigma.obs is commented out because we are now fixing it in the data
  inits <- function() {list(sigma.proc = runif(1, 0, 1),
                            #sigma.obs = runif(1, 0, 1),
                            b = rnorm(1, 1, 0.1),
                            logN.est = c(rnorm(1, jags.data$init_logN, 0.01),
                                         rep(NA, length(data$Occasion) - 1))
  )
  }
  #which parameters do we wish to track the values of
  parameters <- c("b", 
                  "sigma2.proc",
                  #"sigma2.obs",
                  "N.est")
  
  # set run hyperparameters
  ni <- 60000    #iterations
  nb <- 30000    #burn-in
  nc <- 3         #chains
  nt <- 1         #thinning
  
  # call the model
  ssm <- jags(data = jags.data, inits = NULL, parameters.to.save = parameters,
              model.file = "3_Analyses/SSMs/pop_growth/jags_models/ssm_final.jags",
              n.iter = ni, n.chains = nc, n.burnin = nb, n.thin = nt, parallel = TRUE, store.data = TRUE)
  # # Optional code here that calls JAGS and retries if we get an error message
  # # errors may be due to unsuitable initial conditions
  # boolFalse <- FALSE
  # while (boolFalse == FALSE) {
  #   tryCatch({
  #     ssm <- jags(data = jags.data, inits = NULL, parameters.to.save = parameters,
  #                 model.file = "3_Analyses/SSMs/pop_growth/jags_models/ssm_final.jags",
  #                 n.iter = ni, n.chains = nc, n.burnin = nb, n.thin = nt, parallel = TRUE, store.data = TRUE)
  #     boolFalse <- TRUE
  #   }, error=function(e){
  #   }, finally={})
  # }
  
  return(ssm)
}


##### Run the jags model on the data until convergence #####
# we define convergence as all RHat values being less than 1.001
# simply need to do a find and replace on the below stream labels and sex
LL_F_ssm <- ssm_jags(LL_B_derived, sex = "F", n_occ = n_occ97)
LL_M_ssm <- ssm_jags(LL_B_derived, sex = "M", n_occ = n_occ97)

save(LL_F_ssm, file = "4_Results/SSMs/pop_growth/final_riv_sex_models/LL_F_ssm.RData")
save(LL_M_ssm, file = "4_Results/SSMs/pop_growth/final_riv_sex_models/LL_M_ssm.RData")

##### Optionally load previous results #####
# load(file = "4_Results/SSMs/pop_growth/final_riv_sex_models/LL_F_ssm.RData")
# load(file = "4_Results/SSMs/pop_growth/final_riv_sex_models/LL_M_ssm.RData")

##### Generate predictions of true pop. size #####
N_preds_F <- tibble(N_est = LL_F_ssm$mean$N.est, N_LCL_est = LL_F_ssm$q2.5$N.est, N_UCL_est = LL_F_ssm$q97.5$N.est)
N_preds_M <- tibble(N_est = LL_M_ssm$mean$N.est, N_LCL_est = LL_M_ssm$q2.5$N.est, N_UCL_est = LL_M_ssm$q97.5$N.est)

ssms_covs <- LL_B_derived %>% cbind(rbind(N_preds_F, N_preds_M))

##### Plot the results for just a single stream #####
ggplot(ssms_covs) +
  geom_ribbon(aes(x = Occasion, ymin = N_LCL, ymax = N_UCL, fill = Sex), alpha = 0.2) +
  geom_point(aes(x = Occasion, y = N, col = Sex), alpha = 0.5, shape = 1) +
  scale_fill_manual(values = c("plum1", "turquoise")) +
  geom_line(aes(x = Occasion, y = N_est, col = Sex)) +
  geom_point(aes(x = Occasion, y = N_est, col = Sex)) +
  geom_errorbar(aes(x = Occasion, ymin = N_LCL_est, ymax = N_UCL_est, col = Sex)) +
  scale_colour_manual(values = c("maroon3", "dodgerblue")) +
  theme_classic()
ggsave(filename = "5_Figures/SSMs/pop_growth/final_riv_sex_figure/LL_plot.pdf",
       height = 8, width = 12)


##### Grouping the data #####
load(file = "4_Results/SSMs/pop_growth/final_riv_sex_models/LL_F_ssm.RData")
load(file = "4_Results/SSMs/pop_growth/final_riv_sex_models/LL_M_ssm.RData")
load(file = "4_Results/SSMs/pop_growth/final_riv_sex_models/UL_F_ssm.RData")
load(file = "4_Results/SSMs/pop_growth/final_riv_sex_models/UL_M_ssm.RData")
load(file = "4_Results/SSMs/pop_growth/final_riv_sex_models/TY_F_ssm.RData")
load(file = "4_Results/SSMs/pop_growth/final_riv_sex_models/TY_M_ssm.RData")
load(file = "4_Results/SSMs/pop_growth/final_riv_sex_models/CA_F_ssm.RData")
load(file = "4_Results/SSMs/pop_growth/final_riv_sex_models/CA_M_ssm.RData")

LL_preds_F <- tibble(N_est = LL_F_ssm$mean$N.est, N_LCL_est = LL_F_ssm$q2.5$N.est, N_UCL_est = LL_F_ssm$q97.5$N.est)
LL_preds_M <- tibble(N_est = LL_M_ssm$mean$N.est, N_LCL_est = LL_M_ssm$q2.5$N.est, N_UCL_est = LL_M_ssm$q97.5$N.est)
LL_ssm_preds <- LL_B_derived %>% cbind(rbind(LL_preds_F, LL_preds_M))
UL_preds_F <- tibble(N_est = UL_F_ssm$mean$N.est, N_LCL_est = UL_F_ssm$q2.5$N.est, N_UCL_est = UL_F_ssm$q97.5$N.est)
UL_preds_M <- tibble(N_est = UL_M_ssm$mean$N.est, N_LCL_est = UL_M_ssm$q2.5$N.est, N_UCL_est = UL_M_ssm$q97.5$N.est)
UL_ssm_preds <- UL_B_derived %>% cbind(rbind(UL_preds_F, UL_preds_M))
TY_preds_F <- tibble(N_est = TY_F_ssm$mean$N.est, N_LCL_est = TY_F_ssm$q2.5$N.est, N_UCL_est = TY_F_ssm$q97.5$N.est)
TY_preds_M <- tibble(N_est = TY_M_ssm$mean$N.est, N_LCL_est = TY_M_ssm$q2.5$N.est, N_UCL_est = TY_M_ssm$q97.5$N.est)
TY_ssm_preds <- TY_B_derived %>% cbind(rbind(TY_preds_F, TY_preds_M)) %>% mutate(Occasion = Occasion + 12)
CA_preds_F <- tibble(N_est = CA_F_ssm$mean$N.est, N_LCL_est = CA_F_ssm$q2.5$N.est, N_UCL_est = CA_F_ssm$q97.5$N.est)
CA_preds_M <- tibble(N_est = CA_M_ssm$mean$N.est, N_LCL_est = CA_M_ssm$q2.5$N.est, N_UCL_est = CA_M_ssm$q97.5$N.est)
CA_ssm_preds <- CA_B_derived %>% cbind(rbind(CA_preds_F, CA_preds_M)) %>% mutate(Occasion = Occasion + 12)

ssm_preds <- rbind(LL_ssm_preds, UL_ssm_preds, TY_ssm_preds, CA_ssm_preds)
save(ssm_preds, file = "4_Results/SSMs/pop_growth/final_riv_sex_models/ssm_preds.RData")

##### Plotting the data together #####
load("4_Results/SSMs/pop_growth/final_riv_sex_models/ssm_preds.RData")

ggplot(ssm_preds) +
  geom_ribbon(aes(x = Occasion, ymin = N_LCL_est, ymax = N_UCL_est, fill = Sex), alpha = 0.3) +
  geom_line(aes(x = Occasion, y = N_est, col = Sex)) +
  #geom_point(aes(x = Occasion, y = N_est, col = Sex), alpha = 0.8) +
  geom_point(aes(x = Occasion, y = N, col = Sex), shape = 1) +
  facet_wrap(vars(Stream), nrow = 2, scales = "free",
             labeller = as_labeller(c("CA" = "Caigual","LL" = "Lower La Laja", "UL" = "Upper La Laja", "TY" = "Taylor"))) +
  coord_cartesian(ylim = c(0, NA), xlim = c(0,100)) +
  scale_fill_manual(values = c("maroon3", "dodgerblue")) +
  scale_colour_manual(values = c("maroon3", "dodgerblue")) +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) +
  ylab("Population Size") + xlab("Months from March 2008")
ggsave(filename = "5_Figures/SSMs/pop_growth/final_riv_sex_figures/all_plot.pdf",
       height = 8, width = 12)
ggsave(filename = "5_Figures/SSMs/pop_growth/final_riv_sex_figures/all_plot.jpg",
       height = 8, width = 12, units = "in", dpi = 1000)

##### Create climate data #####
load("2_Data_manipulation/data_files/climate_data.RData")

#CA begins in March 2008 - month 13
dry_wet = c(1,1,1,0,0,0,0,0,0,0,1,1) #1 = dry season, begins at March
mar_08_col <- which(colnames(clim_data) == "Mar_2008")
clim <- clim_data[,mar_08_col:(mar_08_col + 97 - 1)] %>% 
  as_tibble() %>%
  rbind(drywet = rep(dry_wet, times = 97 %/% 12 + 1)[1:97]) %>%
  as.matrix() %>% t() %>% as.data.frame() %>%
  mutate(Occasion = 1:97) %>%
  rename("prec" = "V1", "tmin" = "V2", "tmax" = "V3", "dry" = "V4")

##### Performing analyses of covariates explanation #####
# We are creating linear models of growth rate as a function of covariates and N
# Firstly we need to get the right data:
#   we need to generate growth rates for all time steps not the first steps
growth_rate <- rep(NA, times = dim(ssm_preds)[1])
for (i in 2:dim(ssm_preds)[1]) {
  if (ssm_preds$Occasion[i-1] != 97) {
    growth_rate[i-1] <- ssm_preds$N_est[i] / ssm_preds$N_est[i-1]
  }
}

# create dataset of all response and predictors
reg_data <- ssm_preds %>%
  dplyr::select(Occasion, Stream, Sex, N_est, B, phi) %>%
  mutate(N_rate = growth_rate) %>%
  left_join(clim, by = "Occasion") %>%
  mutate(tav = (tmin + tmax) / 2,
         Canopy = as.factor(ifelse(Stream == "LL" | Stream == "CA", "Unthinned", "Thinned")))
save(reg_data, file = "4_Results/SSMs/pop_growth/final_riv_sex_models/reg_data.RData")

##### Analysis of growth rate #####
# N_rate will be a total growth rate, so we need to collate population sizes by sex
total_Ns <- ssm_preds %>%
  dplyr::select(Occasion, Stream, N_est) %>%
  group_by(Occasion, Stream) %>%
  summarise(N_total = sum(N_est)) %>%
  arrange(Stream)
# Calculate total growth rates
total_Ns$N_rate <- NA
for (i in 2:dim(total_Ns)[1]) {
  if (total_Ns$Occasion[i-1] != 97) {
    total_Ns$N_rate[i-1] <- total_Ns$N_total[i] / total_Ns$N_total[i-1]
  }
}
# Augment N_total dataset
reg_data_tot <- total_Ns %>%
  left_join(clim, by = "Occasion") %>%
  mutate(tav = (tmin + tmax) / 2,
         Canopy = as.factor(ifelse(Stream == "LL" | Stream == "CA", "Unthinned", "Thinned")))

# Model selection
N_rate_df <- reg_data_tot %>% ungroup() %>%
  dplyr::select(N_total, Canopy, prec, tav, N_rate) %>% na.omit() %>%
  mutate("Canopy:tav" = (as.integer(Canopy) - 1) * tav,
         "Canopy:N_total" = (as.integer(Canopy) - 1) * N_total,
         "Canopy:prec" = (as.integer(Canopy) - 1) * prec) %>%
  dplyr::select(N_total, Canopy, prec, tav, "Canopy:tav", "Canopy:prec", "Canopy:N_total", N_rate) %>%
  as.data.frame()

# inspect best models
N_rate_best_AIC <- bestglm(N_rate_df, IC = "AIC", TopModels = 10)
N_rate_best_AIC$BestModels
summary(N_rate_best_AIC$BestModel)

N_rate_best <- lm(N_rate ~ N_total + Canopy + prec + tav + Canopy*prec + Canopy*tav,
                  data = N_rate_df, na.action = na.omit)
summary(N_rate_best)
vif(N_rate_best, type = "predictor")


###### Plotting regression results #####
##### Growth Rate added variable plots #####
avPlots(N_rate_best)
##### Plot partial residuals with visreg #####
#N_rate_best is our focus
par(mfrow=c(1,3))
N_plot <- visreg(N_rate_best, "N_total", by = "Canopy", overlay = TRUE,
                 xlab = "Population Size", ylab = "Growth Rate", legend = FALSE, gg = TRUE) +
  labs(title = "P < 0.001")
prec_plot <- visreg(N_rate_best, "prec", by = "Canopy", overlay = TRUE,
                    xlab = "Monthly Rainfall (mm)", ylab = "Growth Rate", legend = FALSE, gg = TRUE) +
  labs(title = "P < 0.01")
tav_plot <- visreg(N_rate_best, "tav", by = "Canopy", overlay = TRUE,
                   xlab = "Monthly Av. Temp (°C)", ylab = "Growth Rate", legend = FALSE, gg = TRUE) +
  labs(title = "P < 0.001")

cp <- c("#7570b3", "#66a61e")

N_plot + prec_plot + tav_plot +
  plot_layout(guides = "collect") &
  scale_fill_manual(values = alpha(cp, 0.4)) &
  scale_colour_manual(values = cp) &
  theme_classic() +
  theme(legend.position = "bottom") &
  plot_annotation(tag_levels = "A")
ggsave(filename = "5_Figures/SSMs/pop_growth/final_riv_sex_figures/visreg_plot.pdf",
       height = 5, width = 10, units = "in")
ggsave(filename = "5_Figures/SSMs/pop_growth/final_riv_sex_figures/visreg_plot.jpg",
       height = 5, width = 10, units = "in", dpi = 1000)
