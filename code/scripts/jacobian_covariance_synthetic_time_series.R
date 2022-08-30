# For each point in time for a given synthetic time series compute:
# (1) Jacobian matrix, (2) covariance and correlation matrix of
# perturbations, (3) determinant of covariance matrix

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/lotka_volterra.R")
source("code/functions/predator_prey.R")
source("code/functions/food_chain.R")
source("code/functions/smap_jacobian.R")
if(!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
if(!require(rootSolve)) {install.packages("rootSolve"); library(rootSolve)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(plotly)) {install.packages("plotly"); library(plotly)}
if(!require(expm)) {install.packages("expm"); library(expm)}
if(!require(MASS)) {install.packages("MASS"); library(MASS)}
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(quadprog)) {install.packages("quadprog"); library(quadprog)}
if(!require(ggsci)) {install.packages("ggsci"); library(ggsci)}
if(!require(doBy)) {install.packages("doBy"); library(doBy)}
if(!require(plyr)) {install.packages("plyr"); library(plyr)}
if(!require(rEDM)) {install.packages("rEDM"); library(rEDM)}

# load the synthetic time series ------------------------------
# number of species (2, 3 or 4)
n_sp <- 2
# model to use (predator_prey, food_chain or lotka_volterra)
func_name <- "predator_prey"
# final time series length
ts_length <- 500
# sampling frequency
sampling_freq <- 20
# time step to evolve perturbations (other values in additional SI analyses)
if (func_name == "predator_prey")
  k <- 3
if (func_name == "food_chain")
  k <- 0.5
if (func_name == "lotka_volterra") {
  if (n_sp == 2)
    k <- 3
  if (n_sp == 4)
    k <- 3
}
# whether k is fixed or variable
k_type <- "fixed"
# initial covariance matrix
cov_initial <- diag(1, nrow = n_sp)
# whether to use analytical or smap Jacobian
jacobian <- "analytical"
# whether to use whole data set or infer sequentially using a moving training set
inference <- "whole"
# whether to normalize data for s-map
normalize <- TRUE
# fraction of time series to train s-map
frac_train <- 0.5
# kernel parameter for s-map
theta_seq <- c(0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 
               1, 2, 3, 4, 5, 6, 7, 8)
# amount of observational noise to add to time series (0, 0.1 or 0.2)
noise <- 0
# amount of noise in k and initial covariance matrix 
noise_k_sigma <- 0
# load model settings
source("code/scripts/model_settings.R")
# load result files
load(file = paste("data/synthetic_time_series/time_series_", func_name, "_", n_sp, "_sp_", 
                  ts_length, "_points_", sampling_freq, "_sampling_freq_", noise, "_noise",
                  ".RData", sep = ""))

# compute time-varying Jacobian matrix ------------------------------ 
# list to store jacobian matrices
J <- list()
# time points to train s-map
train_points <- 1:ceiling(nrow(ts) * frac_train)
# time points to use for further analyses
final_points <- ceiling(nrow(ts) * frac_train):nrow(ts)
# time series to train s-map
ts_train <- ts[train_points, ]
# time series to infer Jacobians
ts <- ts[final_points, ]
# update time series length
ts_length <- nrow(ts)
# analytical Jacobian matrix
if (jacobian == "analytical") {
  J <- dlply(ts[1:(ts_length-1), ], "time", function(x) jacobian.full(y = unlist(c(x[2:(n_sp + 1)])), 
                                                                      func = func,
                                                                      parms = parms))
}
# compute Jacobian with the s-map 
if (jacobian == "smap") {
  # infer Jacobian using whole time series
  if (inference == "whole") {
    ts_smap <- ts
    # normalize time series
    if (normalize) {
      ts_smap[ , -1] <- scale(ts_smap[ , -1])
    } 
    # select theta with cross-validation
    rmse_seq <- sapply(theta_seq, function(x) mean(smap_jacobian(x, ts = ts_smap)[[2]]))
    theta <- theta_seq[which.min(rmse_seq)]
    # perform s-map with optimal theta
    smap_results <- smap_jacobian(ts_smap, theta = theta)
    # s-map jacobians
    J <- smap_results[[1]]
    J <- J[-length(J)]
  }
  # infer Jacobian sequentially
  if (inference == "sequential") {
    for (i in 1:(ts_length-1)) {
      print(i)
      ts_train_curr <- ts_train
      # normalize time series
      if (normalize) {
        ts_train_curr[ , -1] <- scale(ts_train[ , -1])
      } 
      # select theta using cross validation
      rmse_seq <- sapply(theta_seq, function(x) mean(smap_jacobian(x, ts = ts_train_curr)[[2]]))
      theta <- theta_seq[which.min(rmse_seq)]
      # fit jacobian to past points
      J_train <- smap_jacobian(ts_train_curr, theta = theta)[[1]]
      # save last jacobian
      J[[i]] <- J_train[[length(J_train)]]
      # add a point to training data
      ts_train <- rbind(ts_train, ts[i + 1, ])
      # remove first point of training data
      ts_train <- ts_train[-1, ]
    }
  }
}

# compute measures from covariance matrix of perturbations for each point in time ------------------------------ 
# to store results
mean_rate_change <- c()
k_vec <- c()
cov_final_list <- list()
cor_final_list <- list()
cov_det <- c()
VCR <- c()
lambda_cov_sd <- c()
for (i in 1:(ts_length-1)) { 
  # current mean rate of change of the system
  mean_rate_change[i] <- mean(abs((as.numeric(ts[i+1, -1]) - as.numeric(ts[i, -1])) / as.numeric(ts[i, -1])))
  # time steps to evolve perturbations
  if (k_type == "fixed") {
    k_vec[i] <- k
    times_variable <- seq(0, k_vec[i], by = time_step)
  }
  if (k_type == "variable") {
    # determine time sequence to evolve perturbations
    k_vec[i] <- round(1 / sqrt(mean_rate_change[i]), digits = 1)
    times_variable <- seq(0, k_vec[i], by = time_step)
    # if time is 0, determine a very small time
    if (k_vec[i] == 0) {
      k_vec[i] <- 2 * time_step 
      times_variable <- seq(0, k_vec[i], by = time_step)
    }
  }
  # current abundances
  state_curr <- as.numeric(ts[i, -1])
  names(state_curr) <- names(state)
  # add noise to k and initial covariance matrix
  k_vec[i] <- k_vec[i] + rnorm(1, 0, noise_k_sigma * k_vec[i])
  diag(cov_initial) <- diag(cov_initial) + rnorm(n_sp, 0, noise_k_sigma)
  # matrix exponential
  M <- expm(k_vec[i] * J[[i]])
  # compute covariance matrix using matrix exponential
  cov_final <- as.matrix(M %*% cov_initial %*% t(M))
  # determinant of covariance matrix (volume expansion)
  cov_det[i] <- det(cov_final)
  # trace of Jacobian matrix (Volume Contraction Rate)
  VCR[i] <- sum(diag(J[[i]]))
  # standard deviation of eigenvalues of covariance matrix
  lambda_cov_sd[i] <- sd(eigen(cov_final)$values)
  # compute correlation matrix from covariance matrix
  cor_final <- solve(diag(sqrt(diag(cov_final)))) %*% cov_final %*% solve(diag(sqrt(diag(cov_final))))
  cor_final_list[[i]] <- cor_final
  # store covariance matrix
  cov_final_list[[i]] <- cov_final
}

# build data frame to store results ------------------------------
# time indeces and species abundances
df <- ts[1:(ts_length-1), ]
# add Jacobian coefficients
J_df <- data.frame(matrix(unlist(J), nrow = length(J), byrow = TRUE))
names(J_df) <- apply(expand.grid(1:n_sp, 1:n_sp), 1, 
                     function(x) paste("j", paste(x, collapse = ""), sep = "_"))
# add correlations 
cor_df <- data.frame(matrix(unlist(cor_final_list), nrow = length(cor_final_list), byrow = TRUE))
names(cor_df) <- apply(expand.grid(1:n_sp, 1:n_sp), 1, 
                       function(x) paste("cor", paste(x, collapse = ""), sep = "_"))
# add covariances 
cov_df <- data.frame(matrix(unlist(cov_final_list), nrow = length(cov_final_list), byrow = TRUE))
names(cov_df) <- apply(expand.grid(1:n_sp, 1:n_sp), 1, 
                       function(x) paste("cov", paste(x, collapse = ""), sep = "_"))
# merge data frames
full_df <- cbind(df, J_df, cor_df, cov_df)
# add other variables
full_df$mean_rate_change <- mean_rate_change
full_df$k <- k_vec
full_df$cov_det <- cov_det
full_df$VCR <- VCR
full_df$lambda_cov_sd <- lambda_cov_sd

# saving objects as RData ------------------------------
if (jacobian == "analytical") {
  save(full_df, file = paste("results/synthetic_time_series/jacobian_covariance_", func_name, "_", n_sp, "_sp_",
                             k_type, "_time_points_", k, "_k_", "analytical_J", ".RData", sep = ""))
}
if (jacobian == "smap") {
  if (normalize) {
    save(full_df, file = paste("results/synthetic_time_series/jacobian_covariance_", func_name, "_", n_sp, "_sp_",
                               k_type, "_time_points_", k, "_k_", noise, "_noise_", "smap_J_normalized", ".RData", sep = ""))
  } else {
    save(full_df, file = paste("results/synthetic_time_series/jacobian_covariance_", func_name, "_", n_sp, "_sp_",
                               k_type, "_time_points_", k, "_k_", noise, "_noise_", "smap_J", ".RData", sep = ""))
  }
}
