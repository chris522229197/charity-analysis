# This script contains helper functions

library(kableExtra)
library(ggplot2)
library(dplyr)
library(grf)
library(purrr)
library(scales)
library(caret)
library(causalTree)
library(causalToolbox)
library(rpart)
library(rpart.plot)
library(treeClust)
library(car)
library(glmnet)
library(sandwich)

# Split entire data set into interested treatment group
split_control_treated <- function(control_lab, treated_lab, df, 
                                  covariates, outcome_var = "Y", 
                                  treatment_var = "treatment_lvl", 
                                  id_var = "ID") {
  colnames(df)[colnames(df) == outcome_var] <- "Y"
  colnames(df)[colnames(df) == treatment_var] <- "W"
  colnames(df)[colnames(df) == id_var] <- "ID"
  
  subdf <- df %>% filter(W %in% c(control_lab, treated_lab))
  subdf$W[subdf$W == control_lab] <- 0
  subdf$W[subdf$W == treated_lab] <- 1
  subdf <- subdf %>% select_(.dots = c("Y", "W", "ID", covariates))
  Y <- subdf[ , "Y"]
  W <- subdf[ , "W"]
  ID <- subdf[ , "ID"]
  X <- as.matrix(subdf[ , covariates, drop = FALSE])
  return(list("df" = subdf, "Y" = Y, "W" = W, "X" = X, "ID" = ID))
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

p_logistic <- function(dataset) {
  X <- dataset %>% select(-c(Y,W))
  p_logistic.fit <- glm(dataset$W ~ as.matrix(X), family = "binomial")
  p_logistic <- predict(p_logistic.fit, type = "response")
  return(p_logistic)
}

p_rf <- function(dataset, num_trees = NUM_TREES) {
  X <- dataset %>% select(-c(Y,W))
  cf <- causal_forest(X, dataset$Y, dataset$W, num.trees = num_trees)
  return(cf$W.hat)
}

p_glmnet <- function(dataset, k = 10, alpha = 0) {
  X <- as.matrix(dataset %>% select(-c(Y,W)))
  cv <- cv.glmnet(X, dataset$W, nfolds = k, alpha = alpha)
  fitted_glmnet <- glmnet(X, dataset$W, alpha = alpha, lambda = cv$lambda.min)
  p_glmnet <- predict(fitted_glmnet, X, s = cv$lambda.min, type = "response")
  return(p_glmnet)
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

plot_calibration <- function(p, W, ...) {
  {plot(smooth.spline(p, W, df = 4), ...)
  abline(0,1)}
}

# Plot propensity score against covariates

plot_yx <- function(y, X, ylab, xlab, ...) {
  if (is.factor(X)) {
    boxplot(y ~ X, xlab = xlab, ylab = ylab,...)
  } else {
    plot(X, y, xlab = xlab, ylab = ylab, ...)
  }
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

aipw_rf <- function(dataset, num_trees = NUM_TREES) {
  X <- dataset %>% select(-c(Y,W))
  cf <- causal_forest(X, dataset$Y, dataset$W, num.trees = num_trees)
  p <- cf$W.hat
  ate_cf_aipw = average_treatment_effect(cf)
  c(ATE=ate_cf_aipw["estimate"],
    lower_ci=ate_cf_aipw["estimate"] - 1.96 * ate_cf_aipw["std.err"],
    upper_ci=ate_cf_aipw["estimate"] + 1.96 * ate_cf_aipw["std.err"])
}

aipw_glmnet <- function(dataset, p, k = 10, alpha = 0) {
  Y <- dataset$Y
  W <- dataset$W
  X <- data.matrix(dataset %>% select(-c(Y,W)))
  XW <- cbind(X,W)
  cv <- cv.glmnet(XW, Y, nfolds = k, alpha = alpha)
  fitted_glmnet <- glmnet(XW, Y, alpha = alpha, lambda = cv$lambda.min)
  dataset.treatall = cbind(X, 1)
  treated_pred = as.vector(predict(fitted_glmnet, dataset.treatall))
  dataset.treatnone = cbind(X, 0)
  control_pred = as.vector(predict(fitted_glmnet, dataset.treatnone))
  actual_pred = predict(fitted_glmnet, XW)
  G <- treated_pred - control_pred +
    ((dataset$W - p) * (dataset$Y - actual_pred)) / (p * (1 - p))
  tau.hat <- mean(G)
  se.hat <- sqrt(var(G) / (length(G) - 1))
  c(ATE=tau.hat, lower_ci = tau.hat - 1.96 * se.hat, upper_ci = tau.hat + 1.96 * se.hat)
}

# For ATE, show the results of different methods. Talk about how they compare to each other.
# Sensitivity analysis with data set reduction.

calculate_ATE <- function(df1, df2, df3, ratio = 1) {
  df1 <- df1 %>% sample_n(nrow(df1)*ratio)
  df2 <- df2 %>% sample_n(nrow(df2)*ratio)
  df3 <- df3 %>% sample_n(nrow(df3)*ratio)
  
  # Calculate RCT baseline
  tauhat_rct1 <- difference_in_means(df1)
  tauhat_rct2 <- difference_in_means(df2)
  tauhat_rct3 <- difference_in_means(df3)
  
  # ATE using Direct Conditional Regression with OLS
  tauhat_ols1 <- ate_condmean_ols(df1)
  tauhat_ols2 <- ate_condmean_ols(df2)
  tauhat_ols3 <- ate_condmean_ols(df3)
  
  # Calculate propensity scores
  p_logistic1 = p_logistic(df1)
  p_logistic2 = p_logistic(df2)
  p_logistic3 = p_logistic(df3)
  
  p_rf1 = p_rf(df1)
  p_rf2 = p_rf(df2)
  p_rf3 = p_rf(df3)
  
  p_glmnet1 = p_glmnet(df1)
  p_glmnet2 = p_glmnet(df2)
  p_glmnet3 = p_glmnet(df3)
  
  # ATE using IPW with logistic propensity score
  tauhat_logistic_ipw1 <- ipw(df1, p_logistic1)
  tauhat_logistic_ipw2 <- ipw(df2, p_logistic2)
  tauhat_logistic_ipw3 <- ipw(df3, p_logistic3)
  
  # ATE using IPW with random forest propensity score
  tauhat_rf_ipw1 <- ipw(df1, p_rf1)
  tauhat_rf_ipw2 <- ipw(df2, p_rf2)
  tauhat_rf_ipw3 <- ipw(df3, p_rf3)
  
  # ATE using IPW with glmnet propensity score
  tauhat_glmnet_ipw1 <- ipw(df1, p_glmnet1)
  tauhat_glmnet_ipw2 <- ipw(df2, p_glmnet2)
  tauhat_glmnet_ipw3 <- ipw(df3, p_glmnet3)
  
  # ATE using AIPW with logistic propensity score and OLS
  tauhat_logistic_ols_aipw1 <- aipw_ols(df1, p_logistic1)
  tauhat_logistic_ols_aipw2 <- aipw_ols(df2, p_logistic2)
  tauhat_logistic_ols_aipw3 <- aipw_ols(df3, p_logistic3)
  
  # ATE using AIPW with random forest propensity score and regression
  tauhat_rf_aipw1 <- aipw_rf(df1)
  tauhat_rf_aipw2 <- aipw_rf(df2)
  tauhat_rf_aipw3 <- aipw_rf(df3)
  
  # ATE using AIPW with glmnet propensity score and regression
  tauhat_glmnet_aipw1 <- aipw_glmnet(df1, p_glmnet1)
  tauhat_glmnet_aipw2 <- aipw_glmnet(df2, p_glmnet2)
  tauhat_glmnet_aipw3 <- aipw_glmnet(df3, p_glmnet3)

  all_estimators = rbind(
    RCT_gold_standard_1 = tauhat_rct1,
    RCT_gold_standard_2 = tauhat_rct2,
    RCT_gold_standard_3 = tauhat_rct3,
    linear_regression_1 = tauhat_ols1,
    linear_regression_2 = tauhat_ols2,
    linear_regression_3 = tauhat_ols3,
    IPW_logistic_1 = tauhat_logistic_ipw1,
    IPW_logistic_2 = tauhat_logistic_ipw2,
    IPW_logistic_3 = tauhat_logistic_ipw3,
    IPW_glmnet_1 = tauhat_glmnet_ipw1,
    IPW_glmnet_2 = tauhat_glmnet_ipw2,
    IPW_glmnet_3 = tauhat_glmnet_ipw3,
    AIPW_logistic_ols_1 = tauhat_logistic_ols_aipw1,
    AIPW_logistic_ols_2 = tauhat_logistic_ols_aipw2,
    AIPW_logistic_ols_3 = tauhat_logistic_ols_aipw3,
    AIPW_rf_1 = tauhat_rf_aipw1,
    AIPW_rf_2 = tauhat_rf_aipw2,
    AIPW_rf_3 = tauhat_rf_aipw3,
    AIPW_glmnet_1 = tauhat_glmnet_aipw1,
    AIPW_glmnet_2 = tauhat_glmnet_aipw2,
    AIPW_glmnet_3 = tauhat_glmnet_aipw3
  )
  all_estimators <- data.frame(all_estimators)
  all_estimators <- add_rownames(all_estimators, "Estimator")
  all_estimators$treatment_lvl <- substr(all_estimators$Estimator, nchar(all_estimators$Estimator), nchar(all_estimators$Estimator))
  return(all_estimators)
}

# Barplot of ATE Estimates
plot_ATE_estimates <- function(estimators) {
  f <- ggplot(estimators, aes(x = Estimator, y = ATE, color = treatment_lvl)) + 
    geom_line() + 
    geom_point() + 
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.5)
  
  
  f + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major.x = element_blank()) + 
    scale_x_discrete(labels=c(rep("Difference in means", 3), 
                              rep("OLS", 3), 
                              rep("IPW (logistic)", 3), 
                              rep("IPW (ridge)",3), 
                              rep("AIPW (logistic-OLS)", 3), 
                              rep("AIPW (forest-forest)", 3), 
                              rep("AIPW (ridge-ridge)", 3),
                              rep("AIPW ridge ridge", 3))) +
    labs(color = "Matching ratios") +
    ylab("Estimated ATE") + 
    scale_color_discrete(labels = c("1:1", "2:1", "3:1"))
}

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

# S-learner OOB estimates with glmnet
# Default is ridge regression
# Lasso might be too sparse and ignore the coefficient for W 
fit_sglmnet <- function(X, W, Y, ID, k = 10, alpha = 0) {
  XW <- cbind(X, W)
  # k-fold cross validation
  test_folds <- createFolds(Y, k)
  preds_sglmnet <- list()
  for (i in 1:k) {
    # Split data
    test_fold <- test_folds[[i]]
    XW_train <- XW[-test_fold, ]
    Y_train <- Y[-test_fold]
    X_test <- X[test_fold, ]
    # CV for hyperparameter tuning
    cv <- cv.glmnet(XW_train, Y_train, nfolds = k, alpha = alpha)
    # Fit the model with the best hyperparameter
    fitted_glmnet <- glmnet(XW_train, Y_train, lambda = cv$lambda.min, alpha = alpha)
    # Prediction on held-out set
    preds0 <- as.vector(predict(fitted_glmnet, newx = cbind(X_test, 0)))
    preds1 <- as.vector(predict(fitted_glmnet, newx = cbind(X_test, 1)))
    preds_df <- data.frame("pred" = preds1 - preds0, "ID" = ID[test_fold])
    preds_sglmnet[[i]] <- preds_df
  }
  preds_sglmnet <- Reduce(rbind, preds_sglmnet)
  return(preds_sglmnet)
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

# T-learner OOB estimates with glmnet
fit_tglmnet <- function(X, W, Y, ID, k = 10, alpha = 0) {
  # k-fold cross validation
  test_folds <- createFolds(Y, k)
  preds_tglmnet <- list()
  for (i in 1:k) {
    # Split data
    test_fold <- test_folds[[i]]
    X_train <- X[-test_fold, ]
    W_train <- W[-test_fold]
    Y_train <- Y[-test_fold]
    X_test <- X[test_fold, ]
    # CVs for hyperparameter tuning
    cv0 <- cv.glmnet(X_train[W_train == 0, ], Y_train[W_train == 0], nfolds = k, 
                     alpha = alpha)
    cv1 <- cv.glmnet(X_train[W_train == 1, ], Y_train[W_train == 1], nfolds = k, 
                     alpha = alpha)
    # Fit the models with the best hyperparameters
    fitted_glmnet0 <- glmnet(X_train[W_train == 0, ], Y_train[W_train == 0], alpha = alpha, 
                             lambda = cv0$lambda.min)
    fitted_glmnet1 <- glmnet(X_train[W_train == 1, ], Y_train[W_train == 1], alpha = alpha, 
                             lambda = cv1$lambda.min)
    # Prediction on held-out set
    preds0 <- as.vector(predict(fitted_glmnet0, newx = X_test))
    preds1 <- as.vector(predict(fitted_glmnet1, newx = X_test))
    preds_df <- data.frame("pred" = preds1 - preds0, "ID" = ID[test_fold])
    preds_tglmnet[[i]] <- preds_df
  }
  preds_tglmnet <- Reduce(rbind, preds_tglmnet)
  return(preds_tglmnet)
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

# X-learner OOB estimates with glmnet
fit_xglmnet <- function(X, W, Y, ID, k = 10, alpha = 0) {
  # k-fold cross validation
  test_folds <- createFolds(Y, k)
  preds_xglmnet <- list()
  for (i in 1:k) {
    # Split data
    test_fold <- test_folds[[i]]
    X_train <- X[-test_fold, ]
    W_train <- W[-test_fold]
    Y_train <- Y[-test_fold]
    X_test <- X[test_fold, ]
    
    # Treated prediction
    t_learner0_cv <- cv.glmnet(X_train[W_train == 0, ], Y_train[W_train == 0], 
                               alpha = alpha, nfolds = k)
    t_learner0 <- glmnet(X_train[W_train == 0, ], Y_train[W_train == 0], alpha = alpha, 
                         lambda = t_learner0_cv$lambda.min)
    yhat0 <- predict(t_learner0, newx = X_train[W_train == 1, ])
    x_learner1_cv <- cv.glmnet(X_train[W_train == 1, ], Y_train[W_train == 1] - yhat0, 
                               alpha = alpha, nfolds = k)
    x_learner1 <- glmnet(X_train[W_train == 1, ], Y_train[W_train == 1] - yhat0, 
                         alpha = alpha, lambda = x_learner1_cv$lambda.min)
    x_learner_preds1 <- as.vector(predict(x_learner1, newx = X_test))
    
    # Control prediction
    t_learner1_cv <- cv.glmnet(X_train[W_train == 1, ], Y_train[W_train == 1], 
                               alpha = alpha, nfolds = k)
    t_learner1 <- glmnet(X_train[W_train == 1, ], Y_train[W_train == 1], alpha = alpha, 
                         lambda = t_learner1_cv$lambda.min)
    yhat1 <- predict(t_learner1, newx = X_train[W_train == 0, ])
    x_learner0_cv <- cv.glmnet(X_train[W_train == 0, ], yhat1 - Y_train[W_train == 0], 
                               alpha = alpha, nfolds = k)
    x_learner0 <- glmnet(X_train[W_train == 0, ], yhat1 - Y_train[W_train == 0], 
                         alpha = alpha, lambda = x_learner0_cv$lambda.min)
    x_learner_preds0 <- as.vector(predict(x_learner0, newx = X_test))
    
    # Propensity score regression
    prop_glmnet_cv <- cv.glmnet(X_train, W_train, alpha = alpha, nfolds = k)
    prop_glmnet <- glmnet(X_train, W_train, alpha = alpha, 
                          lambda = prop_glmnet_cv$lambda.min)
    ehat <- as.vector(predict(prop_glmnet, X_test))
    preds_x_learner <- (1 - ehat) * x_learner_preds1 + ehat * x_learner_preds0
    preds_xglmnet[[i]] <- data.frame("pred" = preds_x_learner, "ID" = ID[test_fold])
  }
  preds_xglmnet <- Reduce(rbind, preds_xglmnet)
  return(preds_xglmnet)
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

# Compute the OOB estimates and R losses for different forest-based models
compare_forests <- function(X, W, Y, ID, num_trees = 50) {
  # HTE estimates with the models
  message("Fitting the S-learner forest...")
  sf <- fit_sf(X, W, Y, num_trees)
  
  message("Fitting the T-learner forest...")
  tf <- fit_tf(X, W, Y, num_trees)
  
  message("Fitting the X-learner forest...")
  xf <- fit_xf(X, W, Y, num_trees)
  
  message("Fitting the causal forest...")
  cf <- fit_cf(X, W, Y, num_trees)
  
  # R losses
  Y_tilde <- find_Y_tilde(X, Y, num_trees)
  W_tilde <- find_W_tilde(X, W, num_trees)
  
  rloss_fn <- purrr::partial(rloss, Y_tilde = Y_tilde, W_tilde = W_tilde)
  sf_rl <- rloss_fn(sf)
  tf_rl <- rloss_fn(tf)
  xf_rl <- rloss_fn(xf)
  cf_rl <- rloss_fn(cf)
  
  # Pack the results together
  estimates <- list(sf, tf, xf, cf)
  rls <- list(sf_rl, tf_rl, xf_rl, cf_rl)
  learners <- list("S-forest", "T-forest", "X-forest", "Causal forest")
  
  # Organize into a data frame output
  output_lst <- list()
  for (i in 1:length(estimates)) {
    output_df <- data.frame("estimate" = estimates[[i]], 
                            "rloss" = rls[[i]], 
                            "learner" = learners[[i]], 
                            "Y" = Y, "W" = W, "ID" = ID)
    output_df <- cbind(output_df, as.data.frame(X))
    output_lst[[i]] <- output_df
  }
  output <- Reduce(rbind, output_lst)
  output$learner <- factor(output$learner, levels = unlist(learners))
  return(output)
}

# Compute the OOB estimates and R losses for different glmnet-based models
compare_glmnets <- function(X, W, Y, ID, k = 10, alpha = 0, num_trees = 50) {
  # HTE estimates with the models
  message("Fitting the S-learner glmnet...")
  sglm <- fit_sglmnet(X, W, Y, ID, k, alpha)
  
  message("Fitting the T-learner glmnet...")
  tglm <- fit_tglmnet(X, W, Y, ID, k, alpha)
  
  message("Fitting the X-learner glmnet...")
  xglm <- fit_xglmnet(X, W, Y, ID, k, alpha)
  
  # R losses
  Y_tilde <- find_Y_tilde(X, Y, num_trees)
  W_tilde <- find_W_tilde(X, W, num_trees)
  
  rloss_fn <- purrr::partial(rloss, Y_tilde = Y_tilde, W_tilde = W_tilde)
  sglm_rl <- rloss_fn(sglm$pred)
  tglm_rl <- rloss_fn(tglm$pred)
  xglm_rl <- rloss_fn(xglm$pred)
  
  # Pack the results together
  estimates <- list(sglm, tglm, xglm)
  rls <- list(sglm_rl, tglm_rl, xglm_rl)
  learners <- list("S-glmnet", "T-glmnet", "X-glmnet")
  
  # Organize into a data frame output
  output_lst <- list()
  for (i in 1:length(estimates)) {
    output_df <- data.frame("estimate" = estimates[[i]]$pred, 
                            "rloss" = rls[[i]], 
                            "learner" = learners[[i]], 
                            "Y" = Y, "W" = W, "ID" = estimates[[i]]$ID)
    output_df <- cbind(output_df, as.data.frame(X))
    output_lst[[i]] <- output_df
  }
  output <- Reduce(rbind, output_lst)
  output$learner <- factor(output$learner, levels = unlist(learners))
  return(output)
}

# Omnibus test for heterogeneity
# Return an lm object
test_hetero <- function(X, W, Y, tauhat, num_trees = 50) {
  Y_forest <- regression_forest(X, Y, num.trees = num_trees)
  W_forest <- regression_forest(X, W, num.trees = num_trees)
  Yhat <- predict(Y_forest)$predictions
  What <- predict(W_forest)$predictions
  # Get the response and predictors
  mean_tauhat <- mean(tauhat)
  testing_df <- data.frame("Y_tilde" = Y - Yhat, 
                           "avg" = mean_tauhat * (W - What), 
                           "hetero" = (tauhat - mean_tauhat) * (W - What))
  # Perform the test by fitting an OLS
  testing_lm <- lm(Y_tilde ~ avg + hetero + 0, data = testing_df)
  return(testing_lm)
}

# Extract and essential info for testing heterogeneity
summarize_test_hetero <- function(testing_lm) {
  lm_summary <- summary(testing_lm)$coefficients
  beta <- lm_summary["hetero", "Estimate"]
  beta_se <- lm_summary["hetero", "Std. Error"]
  onesided_pval <- lm_summary["hetero", "Pr(>|t|)"] / 2
  return(list("beta" = beta, "beta_se" = beta_se, "onesided_pval" = onesided_pval))
}

# Summarize HTE learner results
summarize_hte_results <- function(hte_results, comparison) {
  hte_summary <- hte_results %>% 
    filter(compare == comparison) %>%
    group_by(learner) %>%
    summarise(MRL = mean(rloss), SDRL = sd(rloss))
  return(hte_summary)
}

# Test whether there is a meaningful difference in doubly robust scores across high and low covariate groups
test_covariate_hetero <- function(df, x_name, drs_name = "drs", x_type = "continuous") {
  drs <- df[ , drs_name]
  x <- df[ , x_name]
  if (x_type == "continuous") {
    high_idx <- x > median(x)
    return(t.test(drs[!high_idx], drs[high_idx]))
  } else {
    pos_idx <- x == 1
    return(t.test(drs[!pos_idx], drs[pos_idx]))
  }
}

# Test whether there is a meaningful difference in covariate across high and low CATE groups
test_CATE <- function(df, cate_name, x, variable = "continuous") {
  cate <- df %>% select(cate_name)
  x <- df %>% select(x)
  med_cate <- median(cate[[1]])
  idx <- cate %>% transmute(case_when(cate >= med_cate ~ 1,
                                      cate < med_cate ~ 0))
  x0 <- x[idx == 0, ]
  x1 <- x[idx == 1, ]
  
  if (variable == "continuous") {
    t.test(x0, x1)  
  } else {
    num_of_successes_in_group_0 <- sum(x0)
    num_of_successes_in_group_1 <- sum(x1)
    count_of_group_0 <- length(x0)
    count_of_group_1 <- length(x1)
    prop.test(c(num_of_successes_in_group_0, num_of_successes_in_group_1), 
              c(count_of_group_0, count_of_group_1), 
              alternative = "two.sided",
              correct = FALSE)
  }
}

# Compute doubly robust scores
compute_drs <- function(X, W, Y, cate, What, Yhat) {
  # Estimate E[Y|X, W=0] and E[Y|X, W=1]
  muhat0 <- Yhat - What * cate
  muhat1 <- Yhat + (1 - What) * cate
  
  # Compute doubly robust scores for the training set
  resid <- Y - W * muhat1 - (1 - W) * muhat0
  weights <- W - What / (What * (1 - What))
  drs <- cate + weights * resid
  return(drs)
}

# Estimate the doubly robust score for binary assignment with training
# and new data with forest-based methods.
# This does not account for the treatment cost.
estimate_binary_drs <- function(X_train, W_train, Y_train, cate_train, 
                                X_new, W_new, Y_new, cate_new, 
                                num_trees = 50, compute_new = TRUE) {
  # Fit the outcome and propensity forests
  Y_forest <- regression_forest(X_train, Y_train, num.trees = num_trees)
  W_forest <- regression_forest(X_train, W_train, num.trees = num_trees)
  
  # Find Y and W estimates for the training set
  Yhat_train <- predict(Y_forest)$predictions # OOB
  What_train <- predict(W_forest)$predictions # OOB
  # Doubly robust scores for the training set
  drs_train <- compute_drs(X_train, W_train, Y_train, cate_train, 
                           What_train, Yhat_train)
  
  # Find Y and W estimates for the new data
  # Also doubly robust scores
  if (compute_new) {
    Yhat_new <- predict(Y_forest, newdata = X_new)$predictions
    What_new <- predict(W_forest, newdata = X_new)$predictions
    drs_new <- compute_drs(X_new, W_new, Y_new, cate_new, 
                           What_new, Yhat_new)
  } else {
    drs_new <- NULL
  }
  return(list("train" = drs_train, "new" = drs_new))
}

# Compute estimated benefit of plug-in policy over random policy
benefit_plugin <- function(cate, cost, net_drs) {
  plugin_assign <- as.numeric(cate >= cost) # assignment
  benefit <- (2 * plugin_assign - 1) * net_drs
  return(benefit)
}

# Summarize the policy learning results
summarize_opt_policy <- function(Ahats_kfold) {
  k_summary_dfs <- list()
  for (k in unique(Ahats_kfold$k_fold)) {
    k_df <- Ahats_kfold %>% filter(k_fold == k)
    k_summary_df <- data.frame("Ahat_policy" = mean(k_df$policy_Ahat), 
                               "Ahat_plugin" = mean(k_df$plugin_Ahat), 
                               "Ahat_everyone" = mean(k_df$everyone_Ahat))
    k_summary_dfs[[k]] <- k_summary_df
  }
  k_summary <- Reduce(rbind, k_summary_dfs)
  summary_df <- data.frame("Ahat_policy_mean" = mean(k_summary$Ahat_policy), 
                           "Ahat_policy_sd" = sd(k_summary$Ahat_policy), 
                           "Ahat_plugin_mean" = mean(k_summary$Ahat_plugin), 
                           "Ahat_plugin_sd" = sd(k_summary$Ahat_plugin), 
                           "Ahat_everyone_mean" = mean(k_summary$Ahat_everyone), 
                           "Ahat_everyone_sd" = sd(k_summary$Ahat_everyone))
  return(summary_df)
}

# Optimal policy learning with tree-based methods accounting for treatment cost
# with cross validation
opt_policy_forest <- function(X, W, Y, cate, ID, cost = 0, k = 10, num_trees = 50, 
                              detailed = TRUE) {
  # Initialize
  k_trees <- list()
  Ahats_test <- list()
  # k-fold cross validation
  test_folds <- createFolds(Y, k)
  for (i in 1:k) {
    # Split data
    test_fold <- test_folds[[i]]
    
    X_train <- X[-test_fold, ]
    W_train <- W[-test_fold]
    Y_train <- Y[-test_fold]
    cate_train <- cate[-test_fold]
    ID_train <- ID[-test_fold]
    
    X_test <- X[test_fold, ]
    W_test <- W[test_fold]
    Y_test <- Y[test_fold]
    cate_test <- cate[test_fold]
    ID_test <- ID[test_fold]
    
    # Estimate the net doubly robust scores (accounting for treatment cost)
    drs <- estimate_binary_drs(X_train, W_train, Y_train, cate_train, 
                               X_test, W_test, Y_test, cate_test, 
                               num_trees = num_trees)
    net_drs_train <- drs$train - cost
    net_drs_test <- drs$new - cost
    
    # Prepare the data frame for policy learning
    covaraites <- colnames(X_train)
    df_train <- as.data.frame(X_train)
    df_train$label <- factor(sign(net_drs_train))
    df_train$weight <- abs(net_drs_train)
    
    # Fit the optimal policy tree
    fmla <- as.formula(paste0("label ~ ", paste0(covaraites, collapse = " + ")))
    policy_tree <- rpart(formula = fmla, data = df_train, weights = df_train$weight)
    k_trees[[i]] <- policy_tree # store the tree
    
    # Predic optimal treatment assignment on the test fold
    df_test <- as.data.frame(X_test)
    df_test$label <- factor(sign(net_drs_test))
    df_test$weight <- abs(net_drs_test)
    policy_assign_test <- as.numeric(predict(policy_tree, newdata = df_test)[ , "1"] > 0.5)
    
    # Calculate optimal policy value over random policy
    policy_val_test <- (2 * policy_assign_test - 1) * net_drs_test
    
    # Calculate plug-in policy value over random policy
    plugin_val_test <- benefit_plugin(cate_test, cost, net_drs_test)
    
    # Combine the results
    Ahats_test[[i]] <- data.frame("ID" = ID_test, 
                                  "policy_Ahat" = policy_val_test, 
                                  "plugin_Ahat" = plugin_val_test, 
                                  "everyone_Ahat" = net_drs_test, 
                                  "k_fold" = i)
    message(paste0("Finished ", i , "/", k, " folds"))
  }
  Ahats_test <- Reduce(rbind, Ahats_test)
  Ahats_test$k_fold <- factor(Ahats_test$k_fold)
  if (!detailed) {
    Ahats_test <- summarize_opt_policy(Ahats_test)
  }
  return(list("Ahats_test" = Ahats_test, "k_trees" = k_trees))
}
