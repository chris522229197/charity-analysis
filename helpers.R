# This script contains helper functions

# Split entire data set into interested treatment group

convert_df <- function(dataset, ctrl, trt) {
  dataset <- dataset %>%
    filter(treatment_lvl == ctrl | treatment_lvl == trt) %>%
    mutate(W = case_when(treatment_lvl == ctrl ~ 0,
                         treatment_lvl == trt ~ 1)) %>%
    select(-c(treatment_lvl,treatment,treat_ratio2,treat_ratio3))
  return (dataset)
}

# Difference of the mean outcome between treated and control group
# This is the classical way to find ATE in a randomized trial

difference_in_means <- function(dataset) {
  # NOTE: Code from Exploring Causal Inference in Experimental and Observational
  # Studies - Part 1 (Susan Athey, Stefan Wager, and Nicolaj Norgaard Muhlbach)
  # Filter treatment / control observations, pulls outcome variable as a vector
  y1 <- dataset %>% dplyr::filter(W == 1) %>% dplyr::pull(Y) # Outcome in treatment grp
  y0 <- dataset %>% dplyr::filter(W == 0) %>% dplyr::pull(Y) # Outcome in control group
  n1 <- sum(dataset[,"W"]) # Number of obs in treatment
  n0 <- sum(1 - dataset[,"W"]) # Number of obs in control
  # Difference in means is ATE
  tauhat <- mean(y1) - mean(y0)
  # 95% Confidence intervals
  se_hat <- sqrt( var(y0)/(n0-1) + var(y1)/(n1-1) )
  lower_ci <- tauhat - 1.96 * se_hat
  upper_ci <- tauhat + 1.96 * se_hat
  return(c(ATE = tauhat, lower_ci = lower_ci, upper_ci = upper_ci))
}

# ATE with OLS for unconfoundedness assumption

ate_condmean_ols <- function(dataset) {
  # NOTE: Code from Exploring Causal Inference in Experimental and Observational
  # Studies - Part 1 (Susan Athey, Stefan Wager, and Nicolaj Norgaard Muhlbach)
  df_mod_centered = data.frame(scale(dataset, center = TRUE, scale = FALSE))
  # Running OLS with full interactions is like running OLS separately on
  # the treated and controls. If the design matrix has been pre-centered,
  # then the W-coefficient corresponds to the ATE.
  lm.interact = lm(Y ~ . * W, data = df_mod_centered)
  tau.hat = as.numeric(coef(lm.interact)["W"])
  se.hat = as.numeric(sqrt(vcovHC(lm.interact)["W", "W"]))
  c(ATE=tau.hat, lower_ci = tau.hat - 1.96 * se.hat, upper_ci = tau.hat + 1.96 * se.hat)
}

# Calculate propensity score

prop_score <- function(dataset) {
  X <- dataset %>% select(-c(Y,W))
  p_logistic.fit <- glm(dataset$W ~ as.matrix(X), family = "binomial")
  p_logistic <- predict(p_logistic.fit, type = "response")
  return(p_logistic)
}

# IPW to account for propensity
# Only requires propensity regression

ipw <- function(dataset, p) {
  # NOTE: Code from Exploring Causal Inference in Experimental and Observational
  # Studies - Part 1 (Susan Athey, Stefan Wager, and Nicolaj Norgaard Muhlbach)
  W <- dataset$W
  Y <- dataset$Y
  G <- ((W - p) * Y) / (p * (1 - p))
  tau.hat <- mean(G)
  se.hat <- sqrt(var(G) / (length(G) - 1))
  c(ATE=tau.hat, lower_ci = tau.hat - 1.96 * se.hat, upper_ci = tau.hat + 1.96 * se.hat)
}

# Plot the calibration of propensity score

plot_calibration <- function(p, W) {
  {plot(smooth.spline(p, W, df = 4))
  abline(0,1)}
}

# AIPW for doubly robust
# Requires outcome regression and propensity regression
# Some freedom for what functions/models we use

aipw_ols <- function(dataset, p) {
  # NOTE: Code from Exploring Causal Inference in Experimental and Observational
  # Studies - Part 1 (Susan Athey, Stefan Wager, and Nicolaj Norgaard Muhlbach)
  ols.fit = lm(Y ~ W * ., data = dataset)
  dataset.treatall = dataset
  dataset.treatall$W = 1
  treated_pred = predict(ols.fit, dataset.treatall)
  dataset.treatnone = dataset
  dataset.treatnone$W = 0
  control_pred = predict(ols.fit, dataset.treatnone)
  actual_pred = predict(ols.fit, dataset)
  G <- treated_pred - control_pred +
    ((dataset$W - p) * (dataset$Y - actual_pred)) / (p * (1 - p))
  tau.hat <- mean(G)
  se.hat <- sqrt(var(G) / (length(G) - 1))
  c(ATE=tau.hat, lower_ci = tau.hat - 1.96 * se.hat, upper_ci = tau.hat + 1.96 * se.hat)
}

# For ATE, show the results of different methods. Talk about how they compare to each other.
# Sensitivity analysis with data set reduction.