##### Description
# R script to perform JS Popan models of guppy data for each stream and sex
#  using RMark

# Choose your choice of stream after loading data and then do a find and replace
# with that stream's initials replacing the initials of the stream in row 59

# Caigual - CA
# Taylor - TY
# Upper Lalaja - UL
# Lower Lalaja - LL

##### Libraries #####
library(tidyverse)
library(RMark)
library(gdata)

##### Load rmark histories #####
load("2_Data_manipulation/data_files/rmark_histories.RData")
##### Keep data for CA or TY #####
keep(TY_F_history, TY_M_history, sure = TRUE)
begin_time <- 13
n_occ <- 85
##### Keep data for UL or LL #####
keep(LL_F_history, LL_M_history, sure = TRUE)
begin_time <- 1
n_occ <- 97

##### Fit JS Popan Formulation to a capture history #####
### Using mark.wrapper.parallel
# 15mins for 8 models of 3631 ind - don't run
popan_wrap <- function(history) {
  #Phi formulae
  Phi.dot <- list(formula = ~ 1)
  Phi.time <- list(formula = ~ time)
  #p formulae
  p.dot <- list(formula = ~ 1)
  p.time <- list(formula = ~ time)
  #pent formulae
  pent.dot <- list(formula = ~ 1)
  pent.time <- list(formula = ~ time)
  #N formulae
  N.dot <- list(formula = ~ 1)
  
  cml <- create.model.list("POPAN")
  rmark_processed <- process.data(history, model = "POPAN")
  ddl <- make.design.data(rmark_processed)
  #null <- mark(rmark_processed, ddl, model = "POPAN", output = FALSE)
  
  return(
    mark.wrapper.parallel(cml, data = rmark_processed, ddl = ddl,
                          #initial = null,
                          cpus = 4, delete = TRUE)
  )
}

##### Model selection for chosen stream #####
start <- Sys.time()
TY_F_model_selection <- popan_wrap(TY_F_history)
TY_M_model_selection <- popan_wrap(TY_M_history)
end <- Sys.time()
end-start

TY_F_model_selection$model.table
### Best model is always model no. 8:
# Phi ~ time, p ~ time, pent ~ time, N ~ 1
TY_F_best_model <- TY_F_model_selection$model.8
save(TY_F_best_model, file = "4_Results/MARK_analyses/rmark_best_models/TY_F_best_model.RData")

TY_M_model_selection$model.table
### Best model is always model no. 8:
# Phi ~ time, p ~ time, pent ~ time, N ~ 1
TY_M_best_model <- TY_M_model_selection$model.8
save(TY_M_best_model, file = "4_Results/MARK_analyses/rmark_best_models/TY_M_best_model.RData")


##### Store parameters and derived values #####
# Define function to store these parameters and derived values
derived_df <- function(rmark_proc, best_model, sex, stream, n_occ) {
  
  derived_N <- best_model$results$derived$`N Population Size` %>%
    as_tibble() %>%
    rename(N = estimate, N_se = se, N_LCL = lcl, N_UCL = ucl)
  derived_B <- rbind(NA, best_model$results$derived$`B* Gross Births`) %>%
    as_tibble() %>%
    rename(B = estimate, B_se = se, B_LCL = lcl, B_UCL = ucl)
  derived_phi <- rbind(best_model$results$real[1:(n_occ-1),1:4], NA) %>%
    as_tibble() %>%
    rename(phi = estimate, phi_se = se, phi_LCL = lcl, phi_UCL = ucl)
  derived_p <- best_model$results$real[n_occ:(2*n_occ - 1), 1:4] %>%
    as_tibble() %>%
    rename(p = estimate, p_se = se, p_LCL = lcl, p_UCL = ucl)
  
  return(cbind(tibble(Occasion = 1:n_occ,
                      Stream = stream,
                      Sex = sex),
               derived_N, derived_B, derived_phi, derived_p))
}

TY_F_processed <- process.data(TY_F_history, model = "POPAN")
TY_M_processed <- process.data(TY_M_history, model = "POPAN")

TY_F_derived <- derived_df(TY_F_processed, TY_F_best_model, "F", "TY", n_occ)
TY_M_derived <- derived_df(TY_M_processed, TY_M_best_model, "M", "TY", n_occ)
TY_B_derived <- rbind(TY_F_derived, TY_M_derived)
save(TY_B_derived, file = "4_Results/MARK_analyses/rmark_best_models/TY_B_derived.RData")

##### Visualise Parameters #####
ggplot(TY_B_derived) +
  geom_ribbon(aes(x = Occasion, ymin = N_LCL, ymax = N_UCL, fill = Sex), alpha = 0.3) +
  geom_line(aes(x = Occasion, y = N, col = Sex)) +
  geom_point(aes(x = Occasion, y = N, col = Sex)) +
  scale_colour_manual(values = c("magenta", "turquoise")) +
  scale_fill_manual(values = c("magenta", "turquoise")) +
  theme_bw()

