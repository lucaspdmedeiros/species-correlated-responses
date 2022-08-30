# Infer Jacobian matrix over time with the S-map for a given empirical 
# time series and save the results

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/smap_jacobian.R")
if(!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
if(!require(rootSolve)) {install.packages("rootSolve"); library(rootSolve)}
if(!require(plyr)) {install.packages("plyr"); library(plyr)}
if(!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
if(!require(tidyr)) {install.packages("tidyr"); library(tidyr)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(plotly)) {install.packages("plotly"); library(plotly)}
if(!require(viridis)) {install.packages("viridis"); library(viridis)}
if(!require(corrplot)) {install.packages("corrplot"); library(corrplot)}
if(!require(latex2exp)) {install.packages("latex2exp"); library(latex2exp)}
if(!require(RColorBrewer)) {install.packages("RColorBrewer"); library(RColorBrewer)}
if(!require(ppcor)) {install.packages("ppcor"); library(ppcor)}
if(!require(rEDM)) {install.packages("rEDM"); library(rEDM)}
if(!require(reticulate)) {install.packages("reticulate"); library(reticulate)}

# loading time series and defining settings ------------------------------ 
# whether to use whole data set or infer sequentially using a moving training set
inference <- "whole"
# whether to normalize data for s-map
normalize <- TRUE
# proportion of data as training set
training_prop <- 0.5
# kernel parameter for s-map
theta_seq <- c(0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 
               1, 2, 3, 4, 5, 6, 7, 8)
# empirical data (blasius_2020_C1_interpolated or becks_2005_panel_g)
data <- "blasius_2020_C1_interpolated"
# load data
ts <- read.csv(paste("data/empirical_time_series/", data, ".csv", sep = ""), header = TRUE)
sp_names <- names(ts)[-1]
n_sp <- ncol(ts) - 1
# define time series section to use for s-map
if (data == "blasius_2020_C1_interpolated") {
  t1 <- 0.95
  t2 <- ts$time[nrow(ts)]
}
if (data == "becks_2005_panel_g") {
  t1 <- 8
  t2 <- ts$time[nrow(ts)]
}
ts <- ts[ts$time >= t1 & ts$time <= t2, ]
ts <- na.omit(ts)
ts_smap <- ts
# change time series names
names(ts_smap) <- c("time", paste("x", 1:n_sp, sep = ""))

# performing inference of Jacobian matrix using whole data set ------------------------------ 
if (inference == "whole") {
  # normalize time series
  if (normalize) {
    mean_ts <- apply(ts_smap[ , -1], 2, mean, na.rm = TRUE)
    sd_ts <- apply(ts_smap[ , -1], 2, sd, na.rm = TRUE)
    ts_smap[ , -1] <- t((t(ts_smap[ , -1]) - mean_ts) / sd_ts)
  }
  # select theta with cross-validation
  rmse_seq <- sapply(theta_seq, function(x) mean(smap_jacobian(x, ts = ts_smap)[[2]]))
  theta <- theta_seq[which.min(rmse_seq)]
  # perform s-map with optimal theta
  smap_results <- smap_jacobian(ts_smap, theta = theta)
  # s-map jacobians
  smap_jacobians <- smap_results[[1]]
  # s-map intercepts
  smap_intercept <- smap_results[[3]]
}

# performing inference of Jacobian matrix sequentially using only past data ------------------------------ 
if (inference == "sequential") {
  # last training point to start with
  start_last_train <- floor(nrow(ts) * training_prop)
  # number of time points
  n <- length((start_last_train + 1):nrow(ts))
  # to store results
  smap_jacobians <- list()
  # for loop moving forward one point at a time
  for (i in 1:n) {
    print(i)
    # define current training and test sets
    training_set_curr <- i:(start_last_train - 1 + i)
    # prepare time series for s-map
    ts_smap <- ts[training_set_curr, ]
    names(ts_smap) <- c("time", paste("x", 1:n_sp, sep = ""))
    ts_smap_curr <- ts_smap
    # normalize training set
    if (normalize) {
      mean_train <- apply(ts_smap_curr[ , -1], 2, mean)
      sd_train <- apply(ts_smap_curr[ , -1], 2, sd)
      ts_smap_curr[ , -1] <- t((t(ts_smap_curr[ , -1]) - mean_train) / sd_train)
    }
    # select theta with cross-validation
    rmse_seq <- sapply(theta_seq, function(x) mean(smap_jacobian(x, ts = ts_smap_curr)[[2]]))
    theta <- theta_seq[which.min(rmse_seq)]
    # perform s-map with optimal theta
    smap_results <- smap_jacobian(ts_smap_curr, theta = theta)
    # s-map jacobians
    smap_J <- smap_results[[1]]
    # s-map intercepts
    smap_intercept <- smap_results[[3]]
    # jacobian at last time step
    smap_jacobians[[i]] <- smap_J[[length(smap_J)]]
  }
}

# save results ------------------------------ 
if (normalize) {
  save(smap_jacobians, file = paste("results/empirical_time_series/smap_jacobians_", data,
                                    "_times_", t1, "_", t2, "_normalized.RData", sep = ""))
} else {
  save(smap_jacobians, file = paste("results/empirical_time_series/smap_jacobians_", data, 
                                    "_times_", t1, "_", t2, ".RData", sep = ""))
}
