# Code for Fig 1: illustration of how species correlated responses change
# over a 2-species predator-prey cycle

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/predator_prey.R")
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
# set seed to replicate figure
set.seed(1)
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
# number of perturbations
n_pert <- 200
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
# time points to use for further analyses
final_points <- ceiling(nrow(ts) * frac_train):nrow(ts)
# time series to compute Jacobians
ts <- ts[final_points, ]
# update time series length
ts_length <- nrow(ts)

# inset: predator-prey time series ------------------------------
# build new data frame for plotting
ts_plot <- ts
ts_plot$time <- 1:nrow(ts_plot)
plot_df <- gather(ts_plot, "variable", "value", -time)
# change variable names and plot
plot_df$variable[plot_df$variable == "x1"] <- "Predator"
plot_df$variable[plot_df$variable == "x2"] <- "Prey"
fig <- ggplot(data = subset(plot_df, variable == "Predator" | variable == "Prey"), 
              aes(x = time, y = value, group = variable, linetype = variable)) +
  geom_line(size = 1.2) +
  geom_vline(size = 1.2, xintercept = c(45, 76), linetype = "dashed") +
  scale_linetype_manual(values = c("solid", "dotted")) +
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
        legend.key.size = unit(0.55, "cm"),
        legend.position = c(0.90, 0.17))
if (save_plots) {
  ggsave(paste("figs/fig1/fig1_predator_prey_ts.pdf", sep = ""), 
         fig, width = 20, height = 8, units = "cm")
}

# panel a (species uncorrelated responses): two opposite perturbations ------------------------------
# time point to use
t <- 294
# current state
curr_state <- as.numeric(ts[ts$time == t, -1])
names(curr_state) <- paste("x", 1:n_sp, sep = "")
# performing two opposite perturbations
curr_state_high <- curr_state + c(7, 0)
curr_state_low <- curr_state + c(-7, 0)
# times to evolve perturbations
times <- seq(0, k, by = time_step)
# generate time series
ts_pert_high <- as.data.frame(ode(y = curr_state_high, times = times, func = func, parms = parms, method = "ode45"))
ts_pert_low <- as.data.frame(ode(y = curr_state_low, times = times, func = func, parms = parms, method = "ode45"))
# sample points according to sampling frequency
ts_pert_high <- ts_pert_high[seq(1, nrow(ts_pert_high), by = 10), ]
ts_pert_low <- ts_pert_low[seq(1, nrow(ts_pert_low), by = 10), ]
# time series for plot
ts_plot <- ts
ts_plot$type <- "Unperturbed"
ts_plot$time <- NA
ts_pert_high$type <- "Predator increase"
ts_pert_low$type <- "Predator decrease"
ts_full <- rbind(ts_plot, ts_pert_high, ts_pert_low)
# plot
fig <- ggplot(data = ts_full, aes(x = x1, y = x2, color = type, alpha = time)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c(met.brewer("Hiroshige", n = 100)[1], met.brewer("Hiroshige", n = 100)[80], "#000000")) + 
  scale_alpha(guide = "none") +
  xlab(TeX("Predator ($N_1$)")) +
  ylab(TeX("Prey ($N_2$)")) +
  scale_x_continuous(limits = c(20, 80)) + 
  scale_y_continuous(limits = c(5, 82)) + 
  coord_equal() +
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
        legend.key.size = unit(0.4, "cm"))
if (save_plots) {
  ggsave(paste("figs/fig1/fig1_predator_prey_single_time_", t, ".pdf", sep = ""), 
         fig, width = 16, height = 16, units = "cm")
}

# panel c (species uncorrelated responses): many random perturbations ------------------------------
# list to store perturbations
pert <- list()
# perform and evolve n_pert perturbations 
df_initial <- data.frame()
df_final <- data.frame()
for (i in 1:n_pert) {
  pert[[i]] <- rnorm(n_sp, mean = 0, sd = 3)
  state_initial <- curr_state + pert[[i]]
  df_initial <- rbind(df_initial, state_initial)
  ts_pert <- as.data.frame(ode(y = state_initial, times = times, func = func, parms = parms, method = "ode45"))
  state_final <- as.numeric(ts_pert[nrow(ts_pert), -1])
  df_final <- rbind(df_final, state_final)
}
names(df_initial) <- paste("x", 1:n_sp, sep = "")
names(df_final) <- paste("x", 1:n_sp, sep = "")
# data frame with initial and final perturbed points
df_pert <- rbind(df_initial, df_final)
df_pert$type <- "Randomly perturbed"
df_pert$time <- c(rep(1, n_pert), rep(3, n_pert))
# computing species correlation using random perturbations
cor(subset(df_pert, time == 3)$x1, subset(df_pert, time == 3)$x2)
# time series for plot
ts_plot <- ts[ , -1]
ts_plot$type <- "Unperturbed"
ts_plot$time <- NA
ts_full <- rbind(ts_plot, df_pert)
ts_full$size <- 0
ts_full$size[ts_full$type == "Unperturbed"] <- 1
# plot
fig <- ggplot(data = ts_full, aes(x = x1, y = x2, color = type, alpha = time, size = size)) +
  geom_point() +
  scale_color_manual(values = c("#BF85BA", "#000000")) + 
  scale_alpha_continuous(range = c(0.2, 0.7), guide = "none") +
  scale_size_continuous(range = c(1.5, 2.5), guide = "none") +
  xlab(TeX("Predator ($N_1$)")) +
  ylab(TeX("Prey ($N_2$)")) +
  scale_x_continuous(limits = c(20, 80)) + 
  scale_y_continuous(limits = c(5, 82)) + 
  coord_equal() +
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
        legend.key.size = unit(0.4, "cm"))
if (save_plots) {
  ggsave(paste("figs/fig1/fig1_predator_prey_multiple_time_", t, ".pdf", sep = ""), 
         fig, width = 16, height = 16, units = "cm")
}

# panel b (species correlated responses): two opposite perturbations ------------------------------
# time point to use
t <- 325
# current state
curr_state <- as.numeric(ts[ts$time == t, -1])
names(curr_state) <- paste("x", 1:n_sp, sep = "")
# performing two opposite perturbations
curr_state_high <- curr_state + c(7, 0)
curr_state_low <- curr_state + c(-7, 0)
# times to evolve perturbations
times <- seq(0, k, by = time_step)
# generate time series
ts_pert_high <- as.data.frame(ode(y = curr_state_high, times = times, func = func, parms = parms, method = "ode45"))
ts_pert_low <- as.data.frame(ode(y = curr_state_low, times = times, func = func, parms = parms, method = "ode45"))
# sample points according to sampling frequency
ts_pert_high <- ts_pert_high[seq(1, nrow(ts_pert_high), by = 10), ]
ts_pert_low <- ts_pert_low[seq(1, nrow(ts_pert_low), by = 10), ]
# time series for plot
ts_plot <- ts
ts_plot$type <- "Unperturbed"
ts_plot$time <- NA
ts_pert_high$type <- "Predator increase"
ts_pert_low$type <- "Predator decrease"
ts_full <- rbind(ts_plot, ts_pert_high, ts_pert_low)
# plot
fig <- ggplot(data = ts_full, aes(x = x1, y = x2, color = type, alpha = time)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c(met.brewer("Hiroshige", n = 100)[1], met.brewer("Hiroshige", n = 100)[80], "#000000")) + 
  scale_alpha(guide = "none") +
  xlab(TeX("Predator ($N_1$)")) +
  ylab(TeX("Prey ($N_2$)")) +
  scale_x_continuous(limits = c(20, 80)) + 
  scale_y_continuous(limits = c(5, 82)) + 
  coord_equal() +
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
        legend.key.size = unit(0.4, "cm"))
if (save_plots) {
  ggsave(paste("figs/fig1/fig1_predator_prey_single_time_", t, ".pdf", sep = ""), 
         fig, width = 16, height = 16, units = "cm")
}

# panel d (species correlated responses): many random perturbations ------------------------------
# perform and evolve n_pert perturbations 
df_initial <- data.frame()
df_final <- data.frame()
for (i in 1:n_pert) {
  state_initial <- curr_state + pert[[i]]
  df_initial <- rbind(df_initial, state_initial)
  ts_pert <- as.data.frame(ode(y = state_initial, times = times, func = func, parms = parms, method = "ode45"))
  state_final <- as.numeric(ts_pert[nrow(ts_pert), -1])
  df_final <- rbind(df_final, state_final)
}
names(df_initial) <- paste("x", 1:n_sp, sep = "")
names(df_final) <- paste("x", 1:n_sp, sep = "")
# data frame with initial and final perturbed points
df_pert <- rbind(df_initial, df_final)
df_pert$type <- "Randomly perturbed"
df_pert$time <- c(rep(1, n_pert), rep(3, n_pert))
# computing species correlation using random perturbations
cor(subset(df_pert, time == 3)$x1, subset(df_pert, time == 3)$x2)
# time series for plot
ts_plot <- ts[ , -1]
ts_plot$type <- "Unperturbed"
ts_plot$time <- NA
ts_full <- rbind(ts_plot, df_pert)
ts_full$size <- 0
ts_full$size[ts_full$type == "Unperturbed"] <- 1
# plot
fig <- ggplot(data = ts_full, aes(x = x1, y = x2, color = type, alpha = time, size = size)) +
  geom_point() +
  scale_color_manual(values = c("#BF85BA", "#000000")) + 
  scale_alpha_continuous(range = c(0.2, 0.7), guide = "none") +
  scale_size_continuous(range = c(1.5, 2.5), guide = "none") +
  xlab(TeX("Predator ($N_1$)")) +
  ylab(TeX("Prey ($N_2$)")) +
  scale_x_continuous(limits = c(20, 80)) + 
  scale_y_continuous(limits = c(5, 82)) + 
  coord_equal() +
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
        legend.key.size = unit(0.4, "cm"))
if (save_plots) {
  ggsave(paste("figs/fig1/fig1_predator_prey_multiple_time_", t, ".pdf", sep = ""), 
         fig, width = 16, height = 16, units = "cm")
}
