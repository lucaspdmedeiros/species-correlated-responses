# Code for Fig S1: accuracy of covariance matrix in describing the distribution of 
# perturbed abundances for a given population dynamics model

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/lotka_volterra.R")
source("code/functions/predator_prey.R")
source("code/functions/food_chain.R")
if(!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
if(!require(rootSolve)) {install.packages("rootSolve"); library(rootSolve)}
if(!require(plyr)) {install.packages("plyr"); library(plyr)}
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
if(!require(expm)) {install.packages("expm"); library(expm)}
if(!require(MetBrewer)) {install.packages("MetBrewer"); library(MetBrewer)}

# load results ------------------------------
# whether to save plots
save_plots <- FALSE
# number of species (2, 3 or 4)
n_sp <- 2
# model to use (predator_prey, food_chain or lotka_volterra)
func_name <- "predator_prey"
# standard deviation of perturbation distribution (i.e., 15% of mean sd of abundances)
if (func_name == "predator_prey") {
  pert_sd <- 2.902995
  k <- 3
}
if (func_name == "food_chain") {
  pert_sd <- 3.749998
  k <- 0.5
}
if (func_name == "lotka_volterra") {
  pert_sd <- 0.01279386
  k <- 3
}
# final time series length
ts_length <- 500
# sampling frequency
sampling_freq <- 20
# amount of observational noise added to time series
noise <- 0
# fraction of time series to train s-map
frac_train <- 0.5
# load model settings
source("code/scripts/model_settings.R")
# perturbation magnitude (sd multiplier)
pert_magnitude <- 0.15
# initial covariance matrix
cov_initial <- diag(pert_sd^2, nrow = n_sp) 
# time steps to evolve perturbed points
k <- 3
# whether number of time points to evolve perturbations is fixed or variable
k_type <- "fixed"
# load result files
load(file = paste("data/synthetic_time_series/time_series_", func_name, "_", n_sp, "_sp_", 
                  ts_length, "_points_", sampling_freq, "_sampling_freq_", noise, "_noise",
                  ".RData", sep = ""))
load(file = paste("results/perturbation_analyses/unperturbed_final_points_", func_name, "_", n_sp, "_sp_",
                  pert_magnitude, "_pert_magn_", k_type, "_time_step_k_", k, ".RData", sep = ""))
load(file = paste("results/perturbation_analyses/perturbed_initial_final_points_", func_name, "_", n_sp, "_sp_", 
                  pert_magnitude, "_pert_magn_", k_type, "_time_step_k_", k, ".RData", sep = ""))
if (k_type == "variable") {
  load(file = paste("results/perturbation_analyses/time_steps_k_", func_name, "_", n_sp, "_sp_",
                    k_type, "_time_step_k_", k, ".RData", sep = ""))
}
ts <- ts[250:(nrow(ts)-1), ]

# compute analytical Jacobian matrix ------------------------------ 
J <- dlply(ts, "time", function(x) jacobian.full(y = unlist(c(x[2:(n_sp + 1)])), 
                                                 func = func,
                                                 parms = parms))

# accuracy between inferred covariance matrix and covariance matrix computed from perturbations ------------------------------ 
cov_infer <- list()
cov_true <- list()
norm_cov_diff <- c()
for (i in 1:length(J)) {
  print(i)
  # current abundances
  state_curr <- as.numeric(ts[i, -1])
  names(state_curr) <- names(state)
  # matrix exponential
  M <- expm(k * J[[i]])
  # compute covariance matrix using matrix exponential
  cov_infer[[i]] <- as.matrix(M %*% cov_initial %*% t(M))
  # data frame with perturbed points at final time point
  curr_df <- df_full[(df_full$time == i) & (df_full$type == "final"), 1:n_sp]
  # data covariance matrix
  cov_true[[i]] <- cov(curr_df)
  # compute distance between covariance matrices
  A <- cov_true[[i]] - cov_infer[[i]]
  norm_cov_diff[i] <- sqrt(sum(diag(A %*% t(A))))
}
# convert lists to data frames
cov_infer_df <- data.frame(matrix(unlist(cov_infer), nrow = length(cov_infer), byrow = TRUE))
names(cov_infer_df) <- apply(expand.grid(1:n_sp, 1:n_sp), 1, 
                             function(x) paste("cov", paste(x, collapse = ""), sep = "_"))
cov_true_df <- data.frame(matrix(unlist(cov_true), nrow = length(cov_true), byrow = TRUE))
names(cov_true_df) <- apply(expand.grid(1:n_sp, 1:n_sp), 1, 
                            function(x) paste("cov", paste(x, collapse = ""), sep = "_"))
# correlation between inferred and true covariance matrix elements
cor_cov <- c()
for(i in 1:ncol(cov_infer_df)) {
  cor_cov[i] <- cor(cov_infer_df[ , i], cov_true_df[ , i])
}
mean(cor_cov)
# prepare data frame for plotting
cov_infer_df$type <- "Theoretical covariance matrix"
cov_infer_df$time <- 1:nrow(cov_infer_df)
cov_true_df$type <- "Covariance matrix of perturbed abundances"
cov_true_df$time <- 1:nrow(cov_true_df)
cov_df <- rbind(cov_infer_df, cov_true_df)
plot_df <- gather(cov_df, "covariance", "value", -c("time", "type"))
# plot
fig <- ggplot() +
  geom_line(data = plot_df, aes(x = time, y = value, group = type, color = type),
            size = 1.5) +
  scale_color_manual(values = c(met.brewer("Hiroshige", n = 100)[1], 
                                met.brewer("Hiroshige", n = 100)[80])) +
  facet_wrap(~ covariance, nrow = n_sp, scales = "free") +
  xlab("Time") +
  ylab("Pairwise covariance") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 14),
        strip.text = element_text(size = 20),
        strip.background = element_rect(fill = "white", size = 1.5),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.2, "cm"),
        legend.position = "top")
# save plot
if (save_plots) {
  ggsave(paste("figs/figSI_", func_name, "_theoretical_vs_perturbation_covariance", 
               ".pdf", sep = ""), 
         fig, width = 28, height = 20, units = "cm")
}
