# Code for Fig S9: contribution of species correlations across state space under 
# a given population dynamics model

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/lotka_volterra.R")
source("code/functions/predator_prey.R")
source("code/functions/food_chain.R")
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

# settings ------------------------------
# whether to save plots
save_plots <- FALSE
# number of species (2, 3 or 4)
n_sp <- 2
# model to use (predator_prey, food_chain or lotka_volterra)
func_name <- "predator_prey"
# final time series length
ts_length <- 500
# sampling frequency
sampling_freq <- 20
# time step to evolve perturbations
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
# whether number of time points to evolve perturbations is fixed or variable
k_type <- "fixed"
# amount of observational noise to add to time series
noise <- 0
# load model settings
source("code/scripts/model_settings.R")
# fraction of time series to use
frac_train <- 0.5
# load result files
load(file = paste("data/synthetic_time_series/time_series_", func_name, "_", n_sp, "_sp_", 
                  ts_length, "_points_", sampling_freq, "_sampling_freq_", noise, "_noise",
                  ".RData", sep = ""))
# results for analytical jacobian
load(file = paste("results/synthetic_time_series/jacobian_covariance_", func_name, "_", n_sp, "_sp_",
                  k_type, "_time_points_", k, "_k_", "analytical_J", ".RData", sep = ""))

# compute community sensitivity decomposition ------------------------------
# community sensitivity (determinant of covariance matrix)
full_df$community_sensitivity <- full_df$cov_det
# contribution of individual species (product of variances)
sub_df <- full_df[ , paste("cov", paste(1:n_sp, 1:n_sp, sep = ""), sep = "_")]
full_df$species_contrib <- apply(sub_df, 1, function(x) prod(x))
# contribution of species correlations
full_df$correlation_contrib <- full_df$community_sensitivity / full_df$species_contrib

# plotting contribution of species correlations across state space ------------------------------
# time points to use
train_points <- 1:ceiling(nrow(ts) * frac_train)
# time points to use for further analyses
final_points <- ceiling(nrow(ts) * frac_train):nrow(ts)
# time series to use
ts <- ts[final_points, ]
# update time series length
ts_length <- nrow(ts)
# build new data frame for plotting
ts_plot <- ts[-nrow(ts), ]
ts_plot$time <- 1:nrow(ts_plot)
ts_plot$correlation_contrib <- log(full_df$correlation_contrib)
# plot for 2-species predator-prey
if (func_name == "predator_prey") {
  fig <- ggplot() +
    geom_point(data = ts_plot, aes(x = x1, y = x2, 
                                   fill = correlation_contrib), 
               size = 3.5, shape = 21) +
    scale_fill_gradient2(low = "#567F55", 
                         mid = "#83C282", 
                         midpoint = mean(ts_plot$correlation_contrib), 
                         high = "#FFFFFF",
                         limits = c(min(ts_plot$correlation_contrib), 
                                    max(ts_plot$correlation_contrib)),
                         name = "Contribution\nof species\ncorrelations") +
    xlab(TeX("Predator ($N_1$)")) +
    ylab(TeX("Prey ($N_2$)")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.5),
          axis.text.y = element_text(size = 16),
          axis.title = element_text(size = 20),
          axis.text.x = element_text(size = 16),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.key.size = unit(0.6, "cm"))
}
# plot for 3-species food chain
if (func_name == "food_chain") {
  fig <- ggplot() +
    geom_point(data = ts_plot, aes(x = x1, y = x3, 
                                        fill = correlation_contrib), 
               size = 3.5, shape = 21) +
    scale_fill_gradient2(low = "#567F55", 
                         mid = "#83C282", 
                         midpoint = mean(ts_plot$correlation_contrib), 
                         high = "#FFFFFF",
                         limits = c(min(ts_plot$correlation_contrib), 
                                    max(ts_plot$correlation_contrib)),
                         name = "Contribution\nof species\ncorrelations") +
    xlab(TeX("Producer ($N_1$)")) +
    ylab(TeX("Secondary consumer ($N_3$)")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.5),
          axis.text.y = element_text(size = 16),
          axis.title = element_text(size = 20),
          axis.text.x = element_text(size = 16),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.key.size = unit(0.6, "cm"))
}
# plot for 2- and 4-species lotka-volterra
if (func_name == "lotka_volterra") {
  if (n_sp == 2) {
    fig <- ggplot() +
      geom_point(data = ts_plot, aes(x = x2, y = x1, 
                                     fill = correlation_contrib),
                 size = 3.5, shape = 21) +
      scale_fill_gradient2(low = "#567F55", 
                           mid = "#83C282", 
                           midpoint = mean(ts_plot$correlation_contrib), 
                           high = "#FFFFFF",
                           limits = c(min(ts_plot$correlation_contrib), 
                                      max(ts_plot$correlation_contrib)),
                           name = "Contribution\nof species\ncorrelations") +
      xlab(TeX("Predator ($N_2$)")) +
      ylab(TeX("Prey ($N_1$)")) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(size = 1.5),
            axis.text.y = element_text(size = 16),
            axis.title = element_text(size = 20),
            axis.text.x = element_text(size = 16),
            plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 14),
            legend.key.size = unit(0.6, "cm"))
  }
  if (n_sp == 4) {
    fig <- ggplot() +
      geom_point(data = ts_plot, aes(x = x4, y = x3, 
                                     fill = correlation_contrib),
                 size = 3.5, shape = 21) +
      scale_fill_gradient2(low = "#567F55", 
                           mid = "#83C282", 
                           midpoint = mean(ts_plot$correlation_contrib), 
                           high = "#FFFFFF",
                           limits = c(min(ts_plot$correlation_contrib), 
                                      max(ts_plot$correlation_contrib)),
                           name = "Contribution\nof species\ncorrelations") +
      xlab(TeX("Competitor 4 ($N_4$)")) +
      ylab(TeX("Competitor 3 ($N_3$)")) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(size = 1.5),
            axis.text.y = element_text(size = 16),
            axis.title = element_text(size = 20),
            axis.text.x = element_text(size = 16),
            plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 14),
            legend.key.size = unit(0.6, "cm"))
  }
}
# save plot
if (save_plots) {
  ggsave(paste("figs/figSI_", func_name, 
               "_contrib_sp_pairs_analytical_state_space_.pdf", sep = ""), 
         fig, width = 18, height = 14, units = "cm")
}
