##### Description #####
# R Script to convert mark-recapture data of guppies into 
# .inp files for analysis with MARK
# and into R data frames

##### Libraries #####
library(tidyverse)
library(rmark)

##### Import Data #####
data <- read_csv("C:/Users/callu/Documents/QBIO7008 Guppy Project/7008_large_files/1_Data/Individual_differences_determine_the_strength_of_ecological_interactions.csv")
names(data)[1] <- "obs"

##### Prepare mark-recapture data in encounter history format
mr_data <- data %>%
  select(obs, stream, sampling, ID, individual_id, sex) %>%
  #create binary "captured" (1) or "not_captured" (0) column
  mutate(capture = 1 - is.na(ID)) %>%
  select(-ID)

#create data-frame of the stream and sex of each individual ID
ind_data <- mr_data %>%
  select(individual_id, stream, sex) %>%
  unique()

#create data-frame of all encounter histories for each river and sex combination
rivers <- unique(mr_data$stream)
sexes <- unique(mr_data$sex)

for (r in rivers) {
  for (s in sexes) {
    enc_his <- mr_data %>%
      filter(stream == r & sex == s)
    n_months <- max(enc_his$sampling) - min(enc_his$sampling) + 1
    enc_his <- enc_his %>%
      select(individual_id, sampling, capture) %>%
      pivot_wider(id_cols = c(individual_id),
                  names_from = sampling,
                  values_from = capture) %>%
      unite(col = history, 2:(1+n_months),
            sep = "", remove = TRUE) %>%
      mutate(individual_id = paste("/*", individual_id, "*/", sep=" "),
             end = "1;")
    write.table(enc_his, file = paste("2_Data_manipulation/MARK_encounter_histories/", paste(r, s, "history.inp", sep = "_"), sep = ""),
                sep = "   ", quote = FALSE, col.names = FALSE, row.names = FALSE)
  }
}


##### Creating files just for each river (sex as covariate) #####
for (r in rivers) {
  enc_his <- mr_data %>%
    filter(stream == r)
  n_months <- max(enc_his$sampling) - min(enc_his$sampling) + 1
  enc_his <- enc_his %>%
    select(individual_id, sex, sampling, capture) %>%
    pivot_wider(id_cols = c(individual_id, sex),
                names_from = sampling,
                values_from = capture) %>%
    unite(col = history, 3:(2+n_months),
          sep = "", remove = TRUE) %>%
    mutate(individual_id = paste("/*", individual_id, "*/", sep=" "),
           male = as.integer(sex == "M"), female = as.integer(sex == "F"), end = ";") %>%
    select(-sex)
  write.table(enc_his, file = paste("2_Data_manipulation/MARK_encounter_histories/", paste(r, "history.inp", sep = "_"), sep = ""),
              sep = "   ", quote = FALSE, col.names = FALSE, row.names = FALSE)
}

##### Creating one for file with river and sex as covariates #####
n_months <- max(enc_his$sampling) - min(enc_his$sampling) + 1
enc_his <- enc_his %>%
  select(individual_id, sex, sampling, capture) %>%
  pivot_wider(id_cols = c(individual_id, sex),
              names_from = sampling,
              values_from = capture) %>%
  unite(col = history, 3:(2+n_months),
        sep = "", remove = TRUE) %>%
  mutate(individual_id = paste("/*", individual_id, "*/", sep=" "),
         male = as.integer(sex == "M"), female = as.integer(sex == "F"), end = ";") %>%
  select(-sex)
write.table(enc_his, file = paste("2_Data_manipulation/MARK_encounter_histories/", paste(r, s, "history.inp", sep = "_"), sep = ""),
            sep = "   ", quote = FALSE, col.names = FALSE, row.names = FALSE)

##### Save .inp objects as R data.frames #####
##### Convert INP files into rmark objects
UL_F_history <- convert.inp("2_Data_manipulation/data_files/MARK_encounter_histories/UL_F_history.inp")
UL_M_history <- convert.inp("2_Data_manipulation/data_files/MARK_encounter_histories/UL_M_history.inp")
LL_F_history <- convert.inp("2_Data_manipulation/data_files/MARK_encounter_histories/LL_F_history.inp")
LL_M_history <- convert.inp("2_Data_manipulation/data_files/MARK_encounter_histories/LL_M_history.inp")
CA_F_history <- convert.inp("2_Data_manipulation/data_files/MARK_encounter_histories/CA_F_history.inp")
CA_M_history <- convert.inp("2_Data_manipulation/data_files/MARK_encounter_histories/CA_M_history.inp")
TY_F_history <- convert.inp("2_Data_manipulation/data_files/MARK_encounter_histories/TY_F_history.inp")
TY_M_history <- convert.inp("2_Data_manipulation/data_files/MARK_encounter_histories/TY_M_history.inp")
UL_history <- convert.inp("2_Data_manipulation/data_files/MARK_encounter_histories/UL_history.inp", group.df = data.frame(sex = c("M", "F")), use.comments = TRUE)
LL_history <- convert.inp("2_Data_manipulation/data_files/MARK_encounter_histories/LL_history.inp", group.df = data.frame(sex = c("M", "F")), use.comments = TRUE)
CA_history <- convert.inp("2_Data_manipulation/data_files/MARK_encounter_histories/CA_history.inp", group.df = data.frame(sex = c("M", "F")), use.comments = TRUE)
TY_history <- convert.inp("2_Data_manipulation/data_files/MARK_encounter_histories/TY_history.inp", group.df = data.frame(sex = c("M", "F")), use.comments = TRUE)

save.image("2_Data_manipulation/data_files/rmark_histories.RData")


