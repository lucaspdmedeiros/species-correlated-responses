# Code for Figs 4, S6, and S7: analytical and s-map contribution of species correlations to
# community sensitivity over time and across state space for a given population dynamics model

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
if(!require(reshape2)) {install.packages("reshape2"); library(reshape2)}
if(!require(MetBrewer)) {install.packages("MetBrewer"); library(MetBrewer)}
if(!require(scales)) {install.packages("scales"); library(scales)}

# load results ------------------------------
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
# whether k is fixed or variable
k_type <- "fixed"
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
# whether to normalize data for s-map
normalize <- TRUE
# amount of observational noise to add to time series (0.1 or 0.2)
noise <- 0.1
# load model settings
source("code/scripts/model_settings.R")
# results for analytical jacobian
load(file = paste("results/synthetic_time_series/jacobian_covariance_", func_name, "_", n_sp, "_sp_",
                  k_type, "_time_points_", k, "_k_", "analytical_J", ".RData", sep = ""))
full_df_analytical <- full_df

# compute community sensitivity decomposition ------------------------------
# community sensitivity (determinant of covariance matrix)
full_df_analytical$community_sensitivity <- full_df_analytical$cov_det
# contribution of individual species (product of variances)
sub_df <- full_df_analytical[ , paste("cov", paste(1:n_sp, 1:n_sp, sep = ""), sep = "_")]
full_df_analytical$species_contrib <- apply(sub_df, 1, function(x) prod(x))
# contribution of species correlations
full_df_analytical$correlation_contrib <- full_df_analytical$community_sensitivity / full_df_analytical$species_contrib
# creating discrete variable for contribution of species correlations
full_df_analytical$correlation_contrib_disc <- "Intermediate"
full_df_analytical$correlation_contrib_disc[log(full_df_analytical$correlation_contrib) > quantile(log(full_df_analytical$correlation_contrib), probs = 0.75)] <- "Low"
full_df_analytical$correlation_contrib_disc[log(full_df_analytical$correlation_contrib) < quantile(log(full_df_analytical$correlation_contrib), probs = 0.25)] <- "High"
# results for s-map jacobian
if (normalize) {
  load(file = paste("results/synthetic_time_series/jacobian_covariance_", func_name, "_", n_sp, "_sp_",
                    k_type, "_time_points_", k, "_k_", noise, "_noise_", "smap_J_normalized", ".RData", sep = ""))
} else {
  load(file = paste("results/synthetic_time_series/jacobian_covariance_", func_name, "_", n_sp, "_sp_",
                    k_type, "_time_points_", k, "_k_", noise, "_noise_", "smap_J", ".RData", sep = ""))
}
full_df_smap <- full_df
# community sensitivity
full_df_smap$community_sensitivity <- full_df_smap$cov_det
# contribution of individual species
sub_df <- full_df_smap[ , paste("cov", paste(1:n_sp, 1:n_sp, sep = ""), sep = "_")]
full_df_smap$species_contrib <- apply(sub_df, 1, function(x) prod(x))
# contribution of species correlations
full_df_smap$correlation_contrib <- full_df_smap$community_sensitivity / full_df_smap$species_contrib
# creating discrete variable for contribution of species correlations
full_df_smap$correlation_contrib_disc <- "Intermediate"
full_df_smap$correlation_contrib_disc[log(full_df_smap$correlation_contrib) > quantile(log(full_df_smap$correlation_contrib), probs = 0.75)] <- "Low"
full_df_smap$correlation_contrib_disc[log(full_df_smap$correlation_contrib) < quantile(log(full_df_smap$correlation_contrib), probs = 0.25)] <- "High"
# computing correlation between analytical and inferred contribution of species correlations
cor(log(full_df_analytical$correlation_contrib), log(full_df_smap$correlation_contrib))

# plots of contribution of species correlations across state space (right panels) ------------------------------
# plot for 2-species predator-prey model
if (func_name == "predator_prey") {
  plot_df_analytical <- data.frame(time = 1:nrow(full_df),
                                   x1 = scale(full_df$x1),
                                   x2 = scale(full_df$x2),
                                   correlation_contrib = full_df_analytical$correlation_contrib,
                                   correlation_contrib_disc = full_df_analytical$correlation_contrib_disc,
                                   type = "Analytical (computed from model)")
  plot_df_smap <- data.frame(time = 1:nrow(full_df),
                             x1 = scale(full_df$x1),
                             x2 = scale(full_df$x2),
                             correlation_contrib = full_df_smap$correlation_contrib,
                             correlation_contrib_disc = full_df_smap$correlation_contrib_disc,
                             type = "S-map (inferred from noisy time series)")
  plot_df_smap_sub <- subset(plot_df_smap, correlation_contrib_disc != "Intermediate")
  fig <- ggplot() +
    geom_point(data = plot_df_smap_sub, aes(x = x1, y = x2, 
                                            fill = correlation_contrib_disc,
                                            shape = correlation_contrib_disc), size = 4) +
    scale_fill_manual(values = c("#567F55", "#83C282")) +
    scale_shape_manual(values = c(21, 24)) +
    xlab(TeX("Predator ($N_1$)")) +
    ylab(TeX("Prey ($N_2$)")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.5),
          axis.text.y = element_text(size = 17),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 17),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          legend.position = "none")
  if (save_plots) {
    ggsave(paste("figs/fig4/fig4_", func_name, "_state_dependency_smap_", 
                 noise, "_noise.pdf", sep = ""), 
           fig, width = 16, height = 16, units = "cm")
  }
}
# plot for 3-species food chain model
if (func_name == "food_chain") {
  plot_df_analytical <- data.frame(time = 1:nrow(full_df),
                                   x1 = scale(full_df$x1),
                                   x2 = scale(full_df$x2),
                                   x3 = scale(full_df$x3),
                                   correlation_contrib = full_df_analytical$correlation_contrib,
                                   correlation_contrib_disc = full_df_analytical$correlation_contrib_disc,
                                   type = "Analytical (computed from model)")
  plot_df_smap <- data.frame(time = 1:nrow(full_df),
                             x1 = scale(full_df$x1),
                             x2 = scale(full_df$x2),
                             x3 = scale(full_df$x3),
                             correlation_contrib = full_df_smap$correlation_contrib,
                             correlation_contrib_disc = full_df_smap$correlation_contrib_disc,
                             type = "S-map (inferred from noisy time series)")
  plot_df_smap_sub <- subset(plot_df_smap, correlation_contrib_disc != "Intermediate")
  fig <- ggplot() +
    geom_point(data = plot_df_smap_sub, aes(x = x1, y = x3, 
                                            fill = correlation_contrib_disc,
                                            shape = correlation_contrib_disc), size = 4) +
    scale_fill_manual(values = c("#567F55", "#83C282")) +
    scale_shape_manual(values = c(21, 24)) +
    xlab(TeX("Producer ($N_1$)")) +
    ylab(TeX("Secondary consumer ($N_3$)")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.5),
          axis.text.y = element_text(size = 17),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 17),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          legend.position = "none")
  if (save_plots) {
    ggsave(paste("figs/fig4/fig4_", func_name, "_state_dependency_smap_", 
                 noise, "_noise.pdf", sep = ""), 
           fig, width = 16, height = 16, units = "cm")
  }
}
# plot for 4-species lotka-volterra model
if (func_name == "lotka_volterra") {
  plot_df_analytical <- data.frame(time = 1:nrow(full_df),
                                   x1 = scale(full_df$x1),
                                   x2 = scale(full_df$x2),
                                   x3 = scale(full_df$x3),
                                   x4 = scale(full_df$x4),
                                   correlation_contrib = full_df_analytical$correlation_contrib,
                                   correlation_contrib_disc = full_df_analytical$correlation_contrib_disc,
                                   type = "Analytical (computed from model)")
  plot_df_smap <- data.frame(time = 1:nrow(full_df),
                             x1 = scale(full_df$x1),
                             x2 = scale(full_df$x2),
                             x3 = scale(full_df$x3),
                             x4 = scale(full_df$x4),
                             correlation_contrib = full_df_smap$correlation_contrib,
                             correlation_contrib_disc = full_df_smap$correlation_contrib_disc,
                             type = "S-map (inferred from noisy time series)")
  plot_df_smap_sub <- subset(plot_df_smap, correlation_contrib_disc != "Intermediate")
  fig <- ggplot() +
    geom_point(data = plot_df_smap_sub, aes(x = x4, y = x3, 
                                            fill = correlation_contrib_disc,
                                            shape = correlation_contrib_disc), size = 4) +
    scale_fill_manual(values = c("#567F55", "#83C282")) +
    scale_shape_manual(values = c(21, 24)) +
    xlab(TeX("Competitor 4 ($N_4$)")) +
    ylab(TeX("Competitor 3 ($N_3$)")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.5),
          axis.text.y = element_text(size = 17),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 17),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          legend.position = "none")
  if (save_plots) {
    ggsave(paste("figs/fig4/fig4_", func_name, "_state_dependency_smap_", 
                 noise, "_noise.pdf", sep = ""), 
           fig, width = 16, height = 16, units = "cm")
  }
}

# plots of contribution of species correlations across time (left panels) ------------------------------
# data frame for plotting
plot_df <- rbind(plot_df_analytical, plot_df_smap)
plot_df$correlation_contrib_disc[plot_df$correlation_contrib_disc == "Intermediate"] <- NA
# plot
fig <- ggplot() +
  geom_line(data = plot_df, aes(x = time, y = correlation_contrib),
            size = 0.4) +
  geom_point(data = plot_df, aes(x = time, y = correlation_contrib, 
                                 fill = correlation_contrib_disc,
                                 shape = correlation_contrib_disc), size = 3) +
  facet_wrap(~type, nrow = 2, scales = "free") +
  scale_fill_manual(values = c("#567F55", "#83C282")) +
  scale_shape_manual(values = c(21, 24)) +
  xlab("Time") +
  ylab(expression(atop("Contribution of species",
                       "correlations (|"~bold(P)~"|) in log"))) +
  guides(fill = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4))) +    
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        strip.text = element_text(size = 17),
        strip.background = element_rect(fill = "white", size = 1.5),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.size = unit(0.8, "cm"),
        legend.position = "top")
# save plot
if (save_plots) {
  ggsave(paste("figs/fig4/fig4_", func_name, "_contrib_sp_pairs_analytical_vs_smap_", 
               noise, "_noise.pdf", sep = ""), 
         fig, width = 24, height = 14, units = "cm")
}
# compute accuracy of contribution of species correlations inferred with s-map
plot_df_analytical$correlation_contrib_disc[plot_df_analytical$correlation_contrib_disc == "Intermediate"] <- NA
plot_df_smap$correlation_contrib_disc[plot_df_smap$correlation_contrib_disc == "Intermediate"] <- NA
(accuracy <- sum(plot_df_analytical$correlation_contrib_disc == plot_df_smap$correlation_contrib_disc, na.rm = TRUE) / 
    sum(!is.na(plot_df_analytical$correlation_contrib_disc)))
# randomization test
accuracy_random <- c()
for (i in 1:1000) {
  shuffled <- plot_df_analytical$correlation_contrib_disc == sample(plot_df_smap$correlation_contrib_disc)
  accuracy_random[i] <- sum(shuffled, na.rm = TRUE) / 
    sum(!is.na(plot_df_analytical$correlation_contrib_disc))
}
mean(accuracy_random > accuracy)
# plot results of randomization test
accuracy_df <- data.frame(accuracy = accuracy_random)
fig <- ggplot() +
  geom_histogram(data = accuracy_df, aes(x = accuracy)) +
  geom_vline(xintercept = accuracy, color = "red", linetype = "dashed", size = 1.5) +
  xlab("S-map accuracy") +
  ylab("Count") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        strip.text = element_text(size = 17),
        strip.background = element_rect(fill = "white", size = 1.5),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
