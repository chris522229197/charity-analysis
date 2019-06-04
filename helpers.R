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

p_rf <- function(dataset, num_trees = 50) {
  X <- dataset %>% select(-c(Y,W))
  cf <- causal_forest(X, dataset$Y, dataset$W, num.trees = num_trees)
  return(cf$W.hat)
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

aipw_rf <- function(dataset, num_trees = 50) {
  X <- dataset %>% select(-c(Y,W))
  cf <- causal_forest(X, dataset$Y, dataset$W, num.trees = num_trees)
  p <- cf$W.hat
  ate_cf_aipw = average_treatment_effect(cf)
  c(ATE=ate_cf_aipw["estimate"],
    lower_ci=ate_cf_aipw["estimate"] - 1.96 * ate_cf_aipw["std.err"],
    upper_ci=ate_cf_aipw["estimate"] + 1.96 * ate_cf_aipw["std.err"])
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
  
  # ATE using IPW with logistic propensity score
  tauhat_logistic_ipw1 <- ipw(df1, p_logistic1)
  tauhat_logistic_ipw2 <- ipw(df2, p_logistic2)
  tauhat_logistic_ipw3 <- ipw(df3, p_logistic3)
  
  # ATE using IPW with random forest propensity score
  tauhat_rf_ipw1 <- ipw(df1, p_rf1)
  tauhat_rf_ipw2 <- ipw(df2, p_rf2)
  tauhat_rf_ipw3 <- ipw(df3, p_rf3)
  
  # ATE using AIPW with logistic propensity score and OLS
  tauhat_logistic_ols_aipw1 <- aipw_ols(df1, p_logistic1)
  tauhat_logistic_ols_aipw2 <- aipw_ols(df2, p_logistic2)
  tauhat_logistic_ols_aipw3 <- aipw_ols(df3, p_logistic3)
  
  # ATE using AIPW with random forest propensity score and regression
  tauhat_rf_aipw1 <- aipw_rf(df1)
  tauhat_rf_aipw2 <- aipw_rf(df2)
  tauhat_rf_aipw3 <- aipw_rf(df3)

  plot_calibration(p_logistic1, df1$W)
  plot_calibration(p_logistic2, df2$W)
  plot_calibration(p_logistic3, df3$W)

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
    AIPW_logistic_ols_1 = tauhat_logistic_ols_aipw1,
    AIPW_logistic_ols_2 = tauhat_logistic_ols_aipw2,
    AIPW_logistic_ols_3 = tauhat_logistic_ols_aipw3,
    AIPW_rf_1 = tauhat_rf_aipw1,
    AIPW_rf_2 = tauhat_rf_aipw2,
    AIPW_rf_3 = tauhat_rf_aipw3
  )
  all_estimators <- data.frame(all_estimators)
  all_estimators <- add_rownames(all_estimators, "Estimator")
  all_estimators$treatment_lvl <- substr(all_estimators$Estimator, nchar(all_estimators$Estimator), nchar(all_estimators$Estimator))
  print(all_estimators)

  f <- ggplot(all_estimators, aes(x = Estimator, y = ATE, ymin = lower_ci, ymax = upper_ci))
  f + geom_crossbar(aes(color = treatment_lvl),
                    position = position_dodge(1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
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
