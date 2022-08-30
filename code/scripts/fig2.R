# Code for Fig 2: illustration of our decomposition of community sensitivity
# for the 2-species predator-prey cycle

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
if(!require(DiagrammeR)) {install.packages("DiagrammeR"); library(DiagrammeR)}
if(!require(reshape2)) {install.packages("reshape2"); library(reshape2)}
if(!require(DiagrammeRsvg)) {install.packages("DiagrammeRsvg"); library(DiagrammeRsvg)}
if(!require(rsvg)) {install.packages("rsvg"); library(rsvg)}

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
# number of perturbations
n_pert <- 200
# time step to evolve perturbations
k <- 3
# whether k is fixed or variable
k_type <- "fixed"
# whether to use analytical or smap Jacobian
jacobian <- "analytical"
# perturbation magnitude
pert_sd <- 3
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
load(file = paste("results/synthetic_time_series/jacobian_covariance_", func_name, "_", n_sp, "_sp_",
                  k_type, "_time_points_", k, "_k_", jacobian, "_J", ".RData", sep = ""))
# time points to use for further analyses
final_points <- ceiling(nrow(ts) * frac_train):nrow(ts)
# time series to infer Jacobians
ts <- ts[final_points, ]
# update time series length
ts_length <- nrow(ts)

# panel b, left: evolution of random perturbations at two points along the cycle ------------------------------
# first time point (species uncorrelated responses)
t1 <- 294
# current state
curr_state <- as.numeric(ts[ts$time == t1, -1])
names(curr_state) <- paste("x", 1:n_sp, sep = "")
# times to evolve perturbations
times <- seq(0, k, by = time_step)
# evolve unperturbed state
sol <- as.data.frame(ode(y = curr_state, times = times, func = func, parms = parms, method = "ode45"))
future_state <- as.numeric(sol[nrow(sol), -1])
names(future_state) <- paste("x", 1:n_sp, sep = "")
# jacobian and covariance matrices at current state
J <- jacobian.full(y = curr_state, func = func, parms = parms)
cov_initial <- diag(rep(pert_sd, n_sp)^2)
M <- expm(k * J)
cov_final <- M %*% cov_initial %*% t(M) 
# sample points from multivariate normal to plot ellipse
df_cov_uncorrelated <- as.data.frame(mvrnorm(n = 1000, mu = future_state, Sigma = cov_final))
# perform and evolve n_pert perturbations 
pert <- list()
df_initial <- data.frame()
df_final <- data.frame()
for (i in 1:n_pert) {
  pert[[i]] <- rnorm(n_sp, mean = 0, sd = pert_sd)
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
df_pert$time <- c(rep(1, n_pert), rep(3, n_pert))
df_uncorrelated <- df_pert
# second time point (species correlated responses)
t2 <- 325
# current state
curr_state <- as.numeric(ts[ts$time == t2, -1])
names(curr_state) <- paste("x", 1:n_sp, sep = "")
# times to evolve perturbations
times <- seq(0, k, by = time_step)
# evolve unperturbed state
sol <- as.data.frame(ode(y = curr_state, times = times, func = func, parms = parms, method = "ode45"))
future_state <- as.numeric(sol[nrow(sol), -1])
names(future_state) <- paste("x", 1:n_sp, sep = "")
# jacobian and covariance matrices at current state
J <- jacobian.full(y = curr_state, func = func, parms = parms)
cov_initial <- diag(rep(pert_sd, n_sp)^2)
M <- expm(k * J)
cov_final <- M %*% cov_initial %*% t(M) 
# sample points from multivariate normal to plot ellipse
df_cov_correlated <- as.data.frame(mvrnorm(n = 1000, mu = future_state, Sigma = cov_final))
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
df_pert$time <- c(rep(1, n_pert), rep(3, n_pert))
df_correlated <- df_pert
# merge both data frames with perturbations
df <- rbind(df_uncorrelated, df_correlated)
# merge both data frames with sampled points (to plot ellipse)
df_cov_uncorrelated$type <- "uncorrelated"
df_cov_correlated$type <- "correlated"
df_cov <- rbind(df_cov_uncorrelated, df_cov_correlated)
# plot
fig <- ggplot() +
  geom_point(data = ts, aes(x = x1, y = x2), color = "#000000", size = 2.5) +
  geom_point(data = df, aes(x = x1, y = x2, alpha = time),
             color = "#BF85BA", size = 1.5) +
  stat_ellipse(data = subset(df_cov, type == "uncorrelated"), aes(x = x1, y = x2), 
               color = "#000000", type = "norm", size = 1.5) +
  stat_ellipse(data = subset(df_cov, type == "correlated"), aes(x = x1, y = x2), 
               color = "#000000", type = "norm", size = 1.5) +
  scale_alpha_continuous(range = c(0.2, 0.7), guide = "none") +
  xlab(TeX("Predator ($N_1$)")) +
  ylab(TeX("Prey ($N_2$)")) +
  scale_x_continuous(limits = c(20, 80)) + 
  scale_y_continuous(limits = c(0, 85)) + 
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
if (save_plots) {
  ggsave(paste("figs/fig2/fig2_predator_prey_multiple_time_", t1, "_", t2, "_ellipse.pdf", sep = ""), 
         fig, width = 16, height = 16, units = "cm")
}

# panel b, right: community sensitivity decomposition over time ------------------------------
# community sensitivity (determinant of covariance matrix)
full_df$community_sensitivity <- full_df$cov_det
# contribution of individual species (product of variances)
sub_df <- full_df[ , paste("cov", paste(1:n_sp, 1:n_sp, sep = ""), sep = "_")]
full_df$species_contrib <- apply(sub_df, 1, function(x) prod(x))
# contribution of species correlations
full_df$correlation_contrib <- full_df$community_sensitivity / full_df$species_contrib
# build new data frame for plotting
full_df$time <- 1:nrow(full_df)
plot_df <- gather(full_df[ , c("time", "community_sensitivity", "species_contrib", "correlation_contrib")], 
                  "variable", "value", -time)
# plot
fig <- ggplot(data = plot_df, 
              aes(x = time, y = value, group = variable, color = variable)) +
  scale_color_manual(values = c("#BF85BA", "#83C282", "#E3A15A")) + 
  geom_line(size = 1.2) +
  geom_vline(size = 1.2, xintercept = c(45, 76), linetype = "dashed") +
  xlab("Time") +
  ylab("Community sensitivity\ndecomposition\n(log scale)") +
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
  ggsave(paste("figs/fig2/fig2_predator_prey_comm_decomp.pdf", sep = ""), 
         fig, width = 22, height = 9, units = "cm")
}
