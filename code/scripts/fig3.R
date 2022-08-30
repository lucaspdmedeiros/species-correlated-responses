# Code for Fig 3: impact of species correlated responses on 
# community sensitivity over time for a given population dynamics model

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

# load synthetic time series ------------------------------
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
# whether to use analytical or smap Jacobian
jacobian <- "analytical"
# amount of observational noise to add to time series
noise <- 0
# load model settings
source("code/scripts/model_settings.R")
# load result files
if (jacobian == "analytical") {
  load(file = paste("results/synthetic_time_series/jacobian_covariance_", func_name, "_", n_sp, "_sp_",
                    k_type, "_time_points_", k, "_k_", "analytical_J", ".RData", sep = ""))
}
if (jacobian == "smap") {
  if (normalize) {
    load(file = paste("results/synthetic_time_series/jacobian_covariance_", func_name, "_", n_sp, "_sp_",
                      k_type, "_time_points_", k, "_k_", noise, "_noise_", "smap_J_normalized", ".RData", sep = ""))
  } else {
    load(file = paste("results/synthetic_time_series/jacobian_covariance_", func_name, "_", n_sp, "_sp_",
                      k_type, "_time_points_", k, "_k_", noise, "_noise_", "smap_J", ".RData", sep = ""))
  }
}

# compute community sensitivity decomposition ------------------------------
# community sensitivity (determinant of covariance matrix)
full_df$community_sensitivity <- full_df$cov_det
# contribution of individual species (product of variances)
sub_df <- full_df[ , paste("cov", paste(1:n_sp, 1:n_sp, sep = ""), sep = "_")]
full_df$species_contrib <- apply(sub_df, 1, function(x) prod(x))
# contribution of species correlations
full_df$correlation_contrib <- full_df$community_sensitivity / full_df$species_contrib

# plot Jacobian and correlation matrices at different points in time ------------------------------
# points in time to use
if (func_name == "predator_prey") {
  t <- c(45, 60, 75, 90)
}
if (func_name == "food_chain") {
  t <- c(75, 77, 79, 81)
}
if (func_name == "lotka_volterra") {
  t <- c(114, 136, 158, 180)
}
# min and max of Jacobian coefficients and correlations (used for plotting)
min_j <- min(full_df[t, apply(expand.grid(1:n_sp, 1:n_sp), 1, 
                              function(x) paste("j", paste(x, collapse = ""), sep = "_"))])
max_j <- max(full_df[t, apply(expand.grid(1:n_sp, 1:n_sp), 1, 
                              function(x) paste("j", paste(x, collapse = ""), sep = "_"))])
min_cor <- min(full_df[t, apply(t(combn(1:n_sp, 2)), 1, 
                                function(x) paste("cor", paste(x, collapse = ""), sep = "_"))])
max_cor <- max(full_df[t, apply(t(combn(1:n_sp, 2)), 1, 
                                function(x) paste("cor", paste(x, collapse = ""), sep = "_"))])
# obtain list of Jacobian matrices
J <- list()
names_J <- apply(expand.grid(1:n_sp, 1:n_sp), 1, 
                 function(x) paste("j", paste(x, collapse = ""), sep = "_"))
for (i in 1:nrow(full_df)) {
  curr_df <- full_df[i, names_J]
  J[[i]] <- matrix(as.numeric(curr_df), nrow = n_sp, ncol = n_sp)
}
# obtain list of correlation matrices
cor_final <- list()
names_cor <- apply(expand.grid(1:n_sp, 1:n_sp), 1, 
                   function(x) paste("cor", paste(x, collapse = ""), sep = "_"))
for (i in 1:nrow(full_df)) {
  curr_df <- full_df[i, names_cor]
  cor_final[[i]] <- matrix(as.numeric(curr_df), nrow = n_sp, ncol = n_sp)
}
# generate matrix heatmap for each point in time
for (i in 1:length(t)) {
  # current Jacobian matrix
  curr_J <- J[[t[i]]]
  rownames(curr_J) <- paste("sp", 1:n_sp, sep = "")
  colnames(curr_J) <- paste("sp", 1:n_sp, sep = "")
  # building data frame by melting matrix
  molten_J <- melt(curr_J)
  names(molten_J) <- c("row", "column", "value")
  # plot (uncomment legend code to generate plot with legend)
  fig <- ggplot(data = molten_J, aes(x = column, y = row, fill = value)) + 
    geom_tile() +
    coord_equal() +
    scale_fill_gradient2(low = met.brewer("Hiroshige", n = 100)[1], mid = "#FFFFFF", 
                         high = met.brewer("Hiroshige", n = 100)[80], midpoint = 0, 
                         limits = c(min_j, max_j)) +
    scale_y_discrete(limits = rev(levels(molten_J$row)), expand = c(0, 0)) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    ylab("") +
    xlab("") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 4),
          legend.position = "none",
          #legend.position = "bottom",
          #legend.title = element_blank(),
          #legend.text = element_text(size = 14),
          #legend.key.height = unit(0.6, 'cm'),
          #legend.key.width = unit(1, 'cm'),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
  # save plot
  if (save_plots) {
    ggsave(paste("figs/fig3/fig3_", func_name, "_J_t_", t[i], ".pdf", sep = ""), 
           fig, width = 14, height = 14, units = "cm")
  }
  # current correlation matrix
  curr_cor <- cor_final[[t[i]]]
  diag(curr_cor) <- NA
  rownames(curr_cor) <- paste("sp", 1:n_sp, sep = "")
  colnames(curr_cor) <- paste("sp", 1:n_sp, sep = "")
  # building data frame by melting matrix
  molten_cor <- melt(curr_cor)
  names(molten_cor) <- c("row", "column", "value")
  # plot (uncomment legend code to generate plot with legend)
  fig <- ggplot(data = molten_cor, aes(x = column, y = row, fill = value)) + 
    geom_tile() +
    coord_equal() +
    scale_fill_gradient2(low = met.brewer("Hiroshige", n = 100)[1], mid = "#FFFFFF", 
                         high = met.brewer("Hiroshige", n = 100)[80], midpoint = 0, 
                         limits = c(min_cor-0.0000001, max_cor+0.0000001)) +
    scale_y_discrete(limits = rev(levels(molten_J$row)), expand = c(0, 0)) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    ylab("") +
    xlab("") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 4),
          legend.position = "none",
          #legend.position = "bottom",
          #legend.title = element_blank(),
          #legend.text = element_text(size = 14),
          #legend.key.height = unit(0.6, 'cm'),
          #legend.key.width = unit(1, 'cm'),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
  # save plot
  if (save_plots) {
    ggsave(paste("figs/fig3/fig3_", func_name, "_cor_t_", t[i], ".pdf", sep = ""), 
           fig, width = 14, height = 14, units = "cm")
  }
}

# plot of contribution of species correlations over time (central panels) ------------------------------
# use only a subset of time series if model is food_chain
if (func_name == "food_chain") {
  full_df <- subset(full_df, time < 360 & time > 310)
}
# redefine time points
full_df$time <- full_df$time - 249
# plot
fig <- ggplot(data = full_df, aes(x = time, y = correlation_contrib)) +
  geom_line(size = 1.5, color = "#83C282") +
  geom_vline(xintercept = t, color = "#000000", linetype = "dashed", size = 1.2) +
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
        axis.title = element_text(size = 19),
        axis.text.x = element_text(size = 16))
# save plot
if (save_plots) {
  ggsave(paste("figs/fig3/fig3_", func_name, "_correlation_contrib_analytical.pdf", sep = ""), 
         fig, width = 24, height = 9, units = "cm")
}
