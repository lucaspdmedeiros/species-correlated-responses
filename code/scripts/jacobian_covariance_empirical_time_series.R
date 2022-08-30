# For each point in time for a given empirical time series compute:
# (1) Jacobian matrix, (2) covariance and correlation matrix of
# perturbations, (3) determinant of covariance matrix

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
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

# loading time series and s-map results ------------------------------ 
# whether to use whole data set or infer sequentially using a moving training set
inference <- "whole"
# proportion of data as training set
training_prop <- 0.5
# whether to normalize data for s-map
normalize <- TRUE
# empirical data (blasius_2020_C1_interpolated or becks_2005_panel_g)
data <- "blasius_2020_C1_interpolated"
# time step to evolve perturbations
k <- 1
# whether k is fixed or variable
k_type <- "fixed"
# whether to plot Jacobian matrix coefficients
plot_J <- FALSE
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
# normalize time series
if (normalize) {
  mean_ts <- apply(ts[ , -1], 2, mean, na.rm = TRUE)
  sd_ts <- apply(ts[ , -1], 2, sd, na.rm = TRUE)
  ts[ , -1] <- t((t(ts[ , -1]) - mean_ts) / sd_ts)
}
# load results
if (normalize) {
  load(file = paste("results/empirical_time_series/smap_jacobians_", data,
                    "_times_", t1, "_", t2, "_normalized.RData", sep = ""))
} else {
  load(file = paste("results/empirical_time_series/smap_jacobians_", data, 
                    "_times_", t1, "_", t2, ".RData", sep = ""))
}
# initial covariance matrix
n_sp <- ncol(ts) - 1
cov_initial <- diag(1, nrow = n_sp)
# list of Jacobian matrices
J <- smap_jacobians
# number of Jacobian matrices
n <- length(J)
if (inference == "sequential") {
  # starting time point
  start_last_train <- floor(nrow(ts) * training_prop)
  # time series containing just test set
  ts <- ts[(start_last_train + 1):nrow(ts), ]
}

# plotting time series of Jacobian matrix coefficients ------------------------------ 
if (plot_J) {
  # data frame with Jacobian coefficients
  J_df <- data.frame(matrix(unlist(J), nrow = length(J), byrow = TRUE))
  names(J_df) <- apply(expand.grid(1:n_sp, 1:n_sp), 1, 
                       function(x) paste("j", paste(x, collapse = ""), sep = "_"))
  J_df$time <- ts$time
  # plot
  plot_df <- gather(J_df, "coefficient", "value", -time)
  ggplot(data = plot_df, aes(x = time, y = value, group = coefficient, color = coefficient)) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
    geom_line(size = 1.5) +
    xlab("Time") +
    ylab("Jacobian coefficient") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1),
          axis.text.y = element_text(size = 16),
          axis.title = element_text(size = 18),
          axis.text.x = element_text(size = 16),
          legend.title = element_blank(),
          legend.text = element_text(size = 16),
          legend.key.size = unit(0.4, "cm"))
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
for (i in 1:(nrow(ts) - 1)) { 
  # current time point
  curr_t <- ts$time[i]
  # current mean rate of change of the system
  mean_rate_change[i] <- mean(abs((as.numeric(ts[i+1, -1]) - as.numeric(ts[i, -1])) / as.numeric(ts[i, -1])))
  # time steps to evolve perturbations
  if (k_type == "fixed") {
    k_vec[i] <- k
  }
  if (k_type == "variable") {
    # determine time sequence to evolve perturbations
    k_vec[i] <- round(1 / sqrt(mean_rate_change[i]), digits = 1)
    # if time is 0, determine a very small time
    if (k_vec[i] == 0) {
      k_vec[i] <- 0.1 
    }
  }
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
df <- ts[-nrow(ts), ]
# add Jacobian coefficients
J_df <- data.frame(matrix(unlist(J), nrow = length(J), byrow = TRUE))
names(J_df) <- apply(expand.grid(1:n_sp, 1:n_sp), 1, 
                     function(x) paste("j", paste(x, collapse = ""), sep = "_"))
J_df <- J_df[-nrow(J_df), ]
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
full_df$k_vec <- k_vec
full_df$cov_det <- cov_det
full_df$VCR <- VCR
full_df$lambda_cov_sd <- lambda_cov_sd

# saving objects as RData ------------------------------
if (normalize) {
  save(full_df, file = paste("results/empirical_time_series/jacobian_covariance_", data, 
                             "_times_", t1, "_", t2, "_", k_type, 
                             "_time_points_normalized", ".RData", sep = ""))
} else {
  save(full_df, file = paste("results/empirical_time_series/jacobian_covariance_", data,
                             "_times_", t1, "_", t2, "_", k_type, 
                             "_time_points", ".RData", sep = ""))
}
