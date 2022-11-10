##### Description #####
# This code generates capture histories for each river to be analysed by 
# state-space formulations

##### Libraries #####
library(tidyverse)

##### Data Manipulation #####
data <- read_csv("C:/Users/callu/Documents/QBIO7008 Guppy Project/7008_large_files/1_Data/Individual_differences_determine_the_strength_of_ecological_interactions.csv")
names(data)[1] <- "obs"
# Get data columns of interest
mr_data <- data %>%
  select(obs, stream, sampling, ID, individual_id, sex) %>%
  #create binary "captured" (1) or "not_captured" (0) column
  mutate(capture = 1 - is.na(ID)) %>%
  select(-ID)

# Set of four rivers:
rivers <- unique(mr_data$stream)

# Filter data into one dataframe per river with pivoting of data to wide format
for (r in rivers) {
  enc_his <- mr_data %>%
    filter(stream == r)
  n_months <- max(enc_his$sampling) - min(enc_his$sampling) + 1
  enc_his <- enc_his %>%
    select(individual_id, sex, sampling, capture) %>%
    pivot_wider(id_cols = c(individual_id, sex),
                names_from = sampling,
                values_from = capture)
  assign(paste(r, "sex", sep = "_"), as.integer(enc_his$sex == "M"))
  enc_his <- enc_his %>%
    select(-individual_id, -sex) %>%
    as.matrix()
  assign(paste(r, "CH", sep = "_"), enc_his)
}

# Save capture history (CH) and sex data
save(CA_CH, LL_CH, UL_CH, TY_CH,
     CA_sex, LL_sex, UL_sex, TY_sex, 
     file = "2_Data_manipulation/data_files/capture_histories_df.RData")
