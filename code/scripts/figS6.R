# Code for Fig S6: contribution of species correlations over time for 
# different values of k for a given population dynamics model

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
if(!require(viridis)) {install.packages("viridis"); library(viridis)}

# load results for synthetic time series ------------------------------
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
# time step to evolve perturbations
k_vec <- c(0.5, 1.5, 3)
# whether to use analytical or smap Jacobian
jacobian <- "analytical"
# amount of observational noise to add to time series
noise <- 0
# load model settings
source("code/scripts/model_settings.R")
# load result files
plot_df <- data.frame()
for (i in 1:length(k_vec)) {
  load(file = paste("results/synthetic_time_series/jacobian_covariance_", func_name, "_", n_sp, "_sp_",
                    k_type, "_time_points_", k_vec[i], "_k_", "analytical_J", ".RData", sep = ""))
  plot_df <- rbind(plot_df, full_df)
}

# compute community sensitivity decomposition ------------------------------
# community sensitivity (determinant of covariance matrix)
plot_df$community_sensitivity <- plot_df$cov_det
# contribution of individual species (product of variances)
sub_df <- plot_df[ , paste("cov", paste(1:n_sp, 1:n_sp, sep = ""), sep = "_")]
plot_df$species_contrib <- apply(sub_df, 1, function(x) prod(x))
# contribution of species correlations
plot_df$correlation_contrib <- plot_df$community_sensitivity / plot_df$species_contrib

# plot of contribution of species correlations for analytical results ------------------------------
# subset time series if model is food chain
if (func_name == "food_chain") {
  plot_df <- subset(plot_df, time < 360 & time > 310)
}
# redefine time points
plot_df$time <- plot_df$time - 249
plot_df$k <- as.factor(plot_df$k)
# plot
fig <- ggplot(data = plot_df, aes(x = time, y = correlation_contrib)) +
  geom_line(size = 1.5, color = "#83C282") +
  facet_wrap( ~ k, nrow = 3, labeller = "label_both", scales = "free") +
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
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = "white", size = 1.5),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position = "top",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.2, "cm"))
if (save_plots) {
  ggsave(paste("figs/figSI_", func_name, "_correlation_contrib_analytical_multiple_k.pdf", sep = ""), 
         fig, width = 24, height = 16, units = "cm")
}
