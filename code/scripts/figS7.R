# Code for Fig S7: illustration of how to use the S-map to perform our community
# sensitivity decomposition for the 2-species predator-prey model

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/predator_prey.R")
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
if(!require(RColorBrewer)) {install.packages("RColorBrewer"); library(RColorBrewer)}
if(!require(latex2exp)) {install.packages("latex2exp"); library(latex2exp)}
if(!require(scales)) {install.packages("scales"); library(scales)}
if(!require(viridis)) {install.packages("viridis"); library(viridis)}
if(!require(ggrepel)) {install.packages("ggrepel"); library(ggrepel)}
if(!require(MetBrewer)) {install.packages("MetBrewer"); library(MetBrewer)}
if(!require(rEDM)) {install.packages("rEDM"); library(rEDM)}

# settings ------------------------------
# whether to save plots
save_plots <- FALSE
# number of species
n_sp <- 2
# model to use
func_name <- "predator_prey"
# final time series length
ts_length <- 500
# sampling frequency
sampling_freq <- 20
# time steps to evolve perturbed points
k <- 3
# amount of observational noise to add to time series
noise <- 0.1
# load model settings
source("code/scripts/model_settings.R")
# fraction of time series to use
frac_train <- 0.5
# whether k is fixed or variable
k_type <- "fixed"
# initial covariance matrix
cov_initial <- diag(1, nrow = n_sp)
# whether to use analytical or smap Jacobian
jacobian <- "smap"
# whether to use whole data set or infer sequentially using a moving training set
inference <- "whole"
# whether to normalize data for s-map
normalize <- TRUE
# kernel parameter for s-map
theta_seq <- c(0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 
               1, 2, 3, 4, 5, 6, 7, 8)
# time point to use as example
t_star <- 76
# load result files
load(file = paste("data/synthetic_time_series/time_series_", func_name, "_", n_sp, "_sp_", 
                  ts_length, "_points_", sampling_freq, "_sampling_freq_", noise, "_noise",
                  ".RData", sep = ""))
# time points to use for further analyses
final_points <- ceiling(nrow(ts) * frac_train):nrow(ts)
# time series to compute Jacobians
ts <- ts[final_points, ]
# update time series length
ts_length <- nrow(ts)

# predator-prey time series ------------------------------
# build new data frame for plotting
ts_plot <- ts
ts_plot$time <- 1:nrow(ts_plot)
plot_df <- gather(ts_plot, "variable", "value", -time)
# change variable names and plot
plot_df$variable[plot_df$variable == "x1"] <- "Predator"
plot_df$variable[plot_df$variable == "x2"] <- "Prey"
fig <- ggplot(data = subset(plot_df, variable == "Predator" | variable == "Prey"), 
              aes(x = time, y = value, group = variable, linetype = variable)) +
  geom_line(size = 1) +
  geom_vline(size = 2, xintercept = t_star, color = "gray50", alpha = 0.5) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  xlab("Time") +
  ylab("Abundance") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1, "cm"),
        legend.position = "top")
if (save_plots) {
  ggsave(paste("figs/figS7/figS7_predator_prey_ts.pdf", sep = ""), 
         fig, width = 22, height = 10, units = "cm")
}

# compute time-varying Jacobian matrix ------------------------------ 
# list to store jacobian matrices
J <- list()
# update time series length
ts_length <- nrow(ts)
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
}

# plot Jacobian matrix elements ------------------------------ 
# data frame with Jacobian elements
J_df <- data.frame(matrix(unlist(J), nrow = length(J), byrow = TRUE))
names(J_df) <- apply(expand.grid(1:n_sp, 1:n_sp), 1, 
                     function(x) paste("j", paste(x, collapse = ""), sep = "_"))
J_df$time <- 1:length(J)
# plot
plot_df <- gather(J_df, "element", "value", -time)
fig <- ggplot(data = plot_df, aes(x = time, y = value, group = element, color = element)) +
  geom_line(size = 1) +
  geom_vline(size = 2, xintercept = t_star, color = "gray50", alpha = 0.5) +
  xlab("Time") +
  ylab("Jacobian matrix element") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1, "cm"),
        legend.position = "top")
if (save_plots) {
  ggsave(paste("figs/figS7/figS7_jacobian_coefficients_ts.pdf", sep = ""), 
         fig, width = 22, height = 10, units = "cm")
}

# plot state space with S-map weights for selected point ------------------------------ 
# compute weights
ts_no_time <- ts[ , -1]
curr_state <- as.numeric(ts_no_time[t_star, ])
distances <- apply(ts_no_time, 1, function(x) sqrt(sum((x - curr_state)^2)))
mean_distance <- mean(distances)
weights <- exp(-theta * (distances / mean_distance))
ts_plot <- ts
ts_plot$weights <- weights
# plot
fig <- ggplot(data = ts_plot, aes(x = x1, y = x2, fill = weights)) +
  geom_point(size = 2.5, shape = 21) +
  xlab(TeX("Predator ($N_1$)")) +
  ylab(TeX("Prey ($N_2$)")) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1), breaks = c(0, 0.5, 1),
                      name = "S-map\nweight") +
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.size = unit(0.45, "cm"))
if (save_plots) {
  ggsave(paste("figs/figS7/figS7_s-map_weights.pdf", sep = ""), 
         fig, width = 13, height = 13, units = "cm")
}

# plot contribution of species correlations ------------------------------ 
# to store results
mean_rate_change <- c()
k_vec <- c()
cov_final_list <- list()
cor_final_list <- list()
cov_det <- c()
VCR <- c()
lambda_cov_sd <- c()
for (i in 1:(ts_length-1)) { 
  # current abundances
  state_curr <- as.numeric(ts[i, -1])
  names(state_curr) <- names(state)
  # initial covariance matrix
  cov_initial <- diag(1, nrow = n_sp)
  # matrix exponential
  M <- expm(k * J[[i]])
  # compute covariance matrix using matrix exponential
  cov_final <- as.matrix(M %*% cov_initial %*% t(M))
  # determinant of covariance matrix (volume expansion)
  cov_det[i] <- det(cov_final)
  # trace of Jacobian matrix (Volume Contraction Rate)
  VCR[i] <- sum(diag(J[[i]]))
  # compute correlation matrix from covariance matrix
  cor_final <- solve(diag(sqrt(diag(cov_final)))) %*% cov_final %*% solve(diag(sqrt(diag(cov_final))))
  cor_final_list[[i]] <- cor_final
  # store covariance matrix
  cov_final_list[[i]] <- cov_final
}
# time indeces and species abundances
df <- ts[1:(ts_length-1), ]
# community sensitivity (determinant of covariance matrix)
df$community_sensitivity <- cov_det
# contribution of individual species (product of variances)
cov_df <- data.frame(matrix(unlist(cov_final_list), nrow = length(cov_final_list), byrow = TRUE))
names(cov_df) <- apply(expand.grid(1:n_sp, 1:n_sp), 1, 
                       function(x) paste("cov", paste(x, collapse = ""), sep = "_"))
sub_df <- cov_df[ , paste("cov", paste(1:n_sp, 1:n_sp, sep = ""), sep = "_")]
df$species_contrib <- apply(sub_df, 1, function(x) prod(x))
# contribution of species correlations
df$correlation_contrib <- df$community_sensitivity / df$species_contrib
# build new data frame for plotting
df$time <- 1:nrow(df)
plot_df <- gather(df[ , c("time", "community_sensitivity", "species_contrib", "correlation_contrib")], 
                  "variable", "value", -time)
# plot
fig <- ggplot(data = df, aes(x = time, y = correlation_contrib)) +
  geom_line(size = 1, color = "#83C282") +
  geom_vline(size = 2, xintercept = t_star, color = "gray50", alpha = 0.5) +
  xlab("Time") +
  ylab(expression(atop("Contribution of species",
                       "correlations (|"~bold(P)~"|) in log"))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position = "none")
if (save_plots) {
  ggsave(paste("figs/figS7/figS7_contrib_species_corr_ts.pdf", sep = ""), 
         fig, width = 22, height = 9, units = "cm")
}
