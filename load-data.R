# This script loads the data for the charity matching data set

library(dplyr)

# Load and preprocess data
df <- read.csv("data/charitable_withdummyvariables.csv")
outcome_variable_name <- "out_gavedum"
treatment_variable_name <- c("treatment", "treat_ratio2", "treat_ratio3")
cont_covariates <- c("pwhite", "pblack", "ave_hh_sz", "median_hhincome", "pop_propurban")
bin_covariates <- c("female", "couple", "red0", "redcty")
covariate_names <- c(cont_covariates, bin_covariates)

cont_cov_descps <- c("Proportion of white", 
                     "Proportion of black", 
                     "Average household size", 
                     "Median household income", 
                     "Proportion of urban population")

bin_cov_descps <- c("Female", "Couple", "Red state", "Red county")

covariate_descps <- c(cont_cov_descps, bin_cov_descps)

df[df == -999] <- NA
df <- df %>% 
  dplyr::select(c(outcome_variable_name, 
                  treatment_variable_name, 
                  covariate_names)) %>% 
  dplyr::rename(Y = out_gavedum) %>%
  na.omit()
df$median_hhincome <- df$median_hhincome / 10000

# Collect the treatment into one variable
df <- df %>% mutate(treatment_lvl = case_when(treatment == 0 ~ 0, 
                                              treat_ratio2 == 1 ~ 2, 
                                              treat_ratio3 == 1 ~ 3, 
                                              TRUE ~ 1))
