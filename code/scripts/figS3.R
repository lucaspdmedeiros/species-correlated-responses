# Code for Fig S3: contribution of species correlations as a function of the maximum
# squared species correlation for a given population dynamics model

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
if(!require(plyr)) {install.packages("plyr"); library(plyr)}
if(!require(scales)) {install.packages("scales"); library(scales)}

# load the synthetic time series ------------------------------
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
# whether number of time points to evolve perturbations is fixed or variable
k_type <- "fixed"
# whether to use analytical or smap Jacobian
jacobian <- "analytical"
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
# amount of observational noise to add to time series
noise <- 0
# load model settings
source("code/scripts/model_settings.R")
# load result files
load(file = paste("results/synthetic_time_series/jacobian_covariance_", func_name, "_", n_sp, "_sp_",
                  k_type, "_time_points_", k, "_k_", jacobian, "_J", ".RData", sep = ""))
# contribution of species correlations
sub_df <- full_df[ , paste("cov", paste(1:n_sp, 1:n_sp, sep = ""), sep = "_")]
full_df$correlation_contrib <- full_df$cov_det / apply(sub_df, 1, function(x) prod(x)) 

# plot contribution of species correlations as a function of maximum squared correlation ------------------------------
if (func_name == "predator_prey") {
  # maximum squared correlation
  full_df$max_cor <- full_df[ , c("cor_12")]^2
  # plot
  fig <- ggplot(data = full_df, aes(x = max_cor, y = correlation_contrib)) +
    geom_point(size = 3, fill = "#929292", shape = 21) +
    xlab(expression(atop("Maximum squared species",
                         "correlation (max("~rho[ij]^2~"))"))) +
    ylab(expression(atop("Contribution of species",
                         "correlations (|"~bold(P)~"|) in log"))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "l") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.5),
          axis.text.y = element_text(size = 18),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 18))
}
if (func_name == "food_chain") {
  # maximum squared correlation
  full_df$max_cor <- apply(full_df[ , c("cor_12", "cor_13", "cor_23")]^2, 1, function(x) max(x))
  # plot
  fig <- ggplot(data = full_df, aes(x = max_cor, y = correlation_contrib)) +
    geom_point(size = 3, fill = "#929292", shape = 21) +
    xlab(expression(atop("Maximum squared species",
                         "correlation (max("~rho[ij]^2~"))"))) +
    ylab(expression(atop("Contribution of species",
                         "correlations (|"~bold(P)~"|) in log"))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "l") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.5),
          axis.text.y = element_text(size = 18),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 18))
}
if (func_name == "lotka_volterra") {
  # maximum squared correlation
  full_df$max_cor <- apply(full_df[ , c("cor_12", "cor_13", "cor_14", 
                                        "cor_23", "cor_24", "cor_34")]^2, 1, function(x) max(x))
  # plot
  fig <- ggplot(data = full_df, aes(x = max_cor, y = correlation_contrib)) +
    geom_point(size = 3, fill = "#929292", shape = 21) +
    xlab(expression(atop("Maximum squared species",
                         "correlation (max("~rho[ij]^2~"))"))) +
    ylab(expression(atop("Contribution of species",
                         "correlations (|"~bold(P)~"|) in log"))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "l") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.5),
          axis.text.y = element_text(size = 18),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 18))
}
# save plot
if (save_plots) {
  ggsave(paste("figs/figSI_", func_name, "_pairs_contrib_vs_max_cor.pdf", sep = ""), 
         fig, width = 16, height = 14, units = "cm")
}
# computing correlation between contribution of species correlations and maximum squared correlation
cor(full_df$max_cor, log(full_df$correlation_contrib))
