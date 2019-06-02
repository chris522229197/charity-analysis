# This script contains helper functions

# Split entire data set into interested treatment group
split_control_treated <- function(control_lab, treated_lab, df, 
                                  covariates, outcome_var = "Y", 
                                  treatment_var = "treatment_lvl") {
  colnames(df)[colnames(df) == outcome_var] <- "Y"
  colnames(df)[colnames(df) == treatment_var] <- "W"
  subdf <- df %>% filter(W %in% c(control_lab, treated_lab))
  subdf$W[subdf$W == control_lab] <- 0
  subdf$W[subdf$W == treated_lab] <- 1
  subdf <- subdf %>% select_(.dots = c("Y", "W", covariates))
  Y <- subdf[ , "Y"]
  W <- subdf[ , "W"]
  X <- as.matrix(subdf[ , covariates, drop = FALSE])
  return(list("df" = subdf, "Y" = Y, "W" = W, "X" = X))
}

# Difference of the mean outcome between treated and control group
# This is the classical way to find ATE in a randomized trial

# ATE with OLS for unconfoundedness assumption

# Calculate propensity score

# IPW to account for propensity
# Only requires propensity regression

# Plot the calibration of propensity score

# APIW for doubly robust
# Requires outcome regression and propensity regression
# Some freedom for what functions/models we use

# For ATE, show the results of different methods. Talk about how they compare to each other.
# Sensitivity analysis with data set reduction.


# S-learner OOB estimates with regression forest
fit_sf <- function(X, W, Y, num_trees = 50) {
  # Code from Lecture 5 of Machine Learning and Causal Inference (Stefan Wager)
  sf <- regression_forest(cbind(X, W), Y, num.trees = num_trees)
  pred.sf.0 = predict(sf, cbind(X, 0))$predictions
  pred.sf.1 = predict(sf, cbind(X, 1))$predictions
  preds.sf.oob = predict(sf)$predictions
  pred.sf.0[W==0] = preds.sf.oob[W==0]
  pred.sf.1[W==1] = preds.sf.oob[W==1]
  preds.sf = pred.sf.1 - pred.sf.0
  return(preds.sf)
}

# T-learner OOB estimates with regression forest
fit_tf <- function(X, W, Y, num_trees = 50) {
  # Code from Lecture 5 of Machine Learning and Causal Inference (Stefan Wager)
  tf0 = regression_forest(X[W==0,], Y[W==0], num.trees = num_trees)
  tf1 = regression_forest(X[W==1,], Y[W==1], num.trees = num_trees)
  tf.preds.0 = predict(tf0, X)$predictions
  tf.preds.1 = predict(tf1, X)$predictions
  tf.preds.0[W==0] = predict(tf0)$predictions #OOB
  tf.preds.1[W==1] = predict(tf1)$predictions #OOB
  preds.tf = tf.preds.1 - tf.preds.0
  return(preds.tf)
}

# X-learner OOB estimates with regression forest
fit_xf <- function(X, W, Y, num_trees = 50) {
  # Code from Lecture 5 of Machine Learning and Causal Inference (Stefan Wager)
  tf0 = regression_forest(X[W==0,], Y[W==0], num.trees = num_trees)
  yhat0 = predict(tf0, X[W==1,])$predictions
  xf1 = regression_forest(X[W==1,], Y[W==1]-yhat0, num.trees = num_trees)
  xf.preds.1 = predict(xf1, X)$predictions
  xf.preds.1[W==1] = predict(xf1)$predictions
  
  tf1 = regression_forest(X[W==1,], Y[W==1], num.trees = num_trees)
  yhat1 = predict(tf1, X[W==0,])$predictions
  xf0 = regression_forest(X[W==0,], yhat1-Y[W==0], num.trees = num_trees)
  xf.preds.0 = predict(xf0, X)$predictions
  xf.preds.0[W==0] = predict(xf0)$predictions
  
  propf = regression_forest(X, W, tune.parameters = TRUE, num.trees = num_trees)
  ehat = predict(propf)$predictions
  preds.xf = (1 - ehat) * xf.preds.1 + ehat * xf.preds.0
  return(preds.xf)
}

# Causal forest OOB estimates
fit_cf <- function(X, W, Y, num_trees = 50) {
  cf <- causal_forest(X, Y, W, num.trees = num_trees)
  oob_pred <- predict(cf, estimate.variance = TRUE)
  return(oob_pred$predictions)
}

# Find Y tilde using regression forest for the R loss computation
find_Y_tilde <- function(X, Y, num_trees = 50) {
  Y_forest <- regression_forest(X, Y, num.trees = num_trees)
  Yhat <- predict(Y_forest)$predictions
  return(Y - Yhat)
}

# Find W tilde using regression forest for the R loss computation
find_W_tilde <- function(X, W, num_trees = 50) {
  W_forest <- regression_forest(X, W, num.trees = num_trees)
  What <- predict(W_forest)$predictions
  return(W - What)
}

# Compute the R loss for each sample
rloss <- function(Y_tilde, W_tilde, estimated_tau) {
  # Robinson's method as baseline
  robinson_ols <- lm(Y_tilde ~ W_tilde)
  robinson_tau <- robinson_ols$coefficients["W_tilde"]
  delta <- (Y_tilde - estimated_tau*W_tilde)^2 - (Y_tilde - robinson_tau*W_tilde)^2
  return(delta)
}
