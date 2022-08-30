# Code for Figs 5a and S9a: contribution of species correlations to community sensitivity 
# over time and across state space for the 2-species predator-prey time series
# from Blasius et al (2020) Nature

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
if(!require(scales)) {install.packages("scales"); library(scales)}
if(!require(viridis)) {install.packages("viridis"); library(viridis)}
if(!require(ggrepel)) {install.packages("ggrepel"); library(ggrepel)}

# loading time series and s-map results ------------------------------ 
# whether to save plots
save_plots <- FALSE
# whether to normalize data for s-map
normalize <- TRUE
# empirical data
data <- "blasius_2020_C1_interpolated"
# whether k is fixed or variable
k_type <- "fixed"
# load data
ts <- read.csv(paste("data/empirical_time_series/", data, ".csv", sep = ""), header = TRUE)
sp_names <- names(ts)[-1]
n_sp <- ncol(ts) - 1
# define time series section to use for s-map
t1 <- 0.95
t2 <- ts$time[nrow(ts)]
ts <- ts[ts$time >= t1 & ts$time <= t2, ]
# normalize time series
if (normalize) {
  mean_ts <- apply(ts[ , -1], 2, mean, na.rm = TRUE)
  sd_ts <- apply(ts[ , -1], 2, sd, na.rm = TRUE)
  ts[ , -1] <- t((t(ts[ , -1]) - mean_ts) / sd_ts)
}
if (normalize) {
  load(file = paste("results/empirical_time_series/smap_jacobians_", data,
                    "_times_", t1, "_", t2, "_normalized.RData", sep = ""))
} else {
  load(file = paste("results/empirical_time_series/smap_jacobians_", data, 
                    "_times_", t1, "_", t2, ".RData", sep = ""))
}
J_list <- smap_jacobians
# number of Jacobian matrices
n <- length(J_list)
if (normalize) {
  load(file = paste("results/empirical_time_series/jacobian_covariance_", data, 
                    "_times_", t1, "_", t2, "_", k_type, 
                    "_time_points_normalized", ".RData", sep = ""))
} else {
  load(file = paste("results/empirical_time_series/jacobian_covariance_", data,
                    "_times_", t1, "_", t2, "_", k_type, 
                    "_time_points", ".RData", sep = ""))
}

# compute community sensitivity decomposition ------------------------------
# community sensitivity (determinant of covariance matrix)
full_df$community_sensitivity <- full_df$cov_det
# contribution of individual species (product of variances)
sub_df <- full_df[ , paste("cov", paste(1:n_sp, 1:n_sp, sep = ""), sep = "_")]
full_df$species_contrib <- apply(sub_df, 1, function(x) prod(x))
# contribution of species correlations
full_df$correlation_contrib <- full_df$community_sensitivity / full_df$species_contrib

# top left panel: abundance time series ------------------------------
# taking only a subset of time series
ts_plot <- ts[ts$time > 20 & ts$time < 135, ]
plot_df <- gather(ts_plot, "variable", "value", -time)
# change variable names and plot
plot_df$variable[plot_df$variable == "algae"] <- "Prey (algae)"
plot_df$variable[plot_df$variable == "rotifers"] <- "Predator (rotifers)"
fig <- ggplot(data = plot_df, 
              aes(x = time, y = value, group = variable, linetype = variable)) +
  geom_line(size = 1) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  xlab("Time (days)") +
  ylab("Abundance (normalized)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1.2, "cm"),
        legend.position = "top")
if (save_plots) {
  if (normalize) {
    ggsave(paste("figs/fig5/fig5_", data, "_times_", t1, "_", t2, 
                 "_ts_normalized.pdf", sep = ""), 
           fig, width = 26, height = 9, units = "cm")
  } else {
    ggsave(paste("figs/fig5/fig5_", data, "_times_", t1, "_", t2, 
                 "_ts.pdf", sep = ""), 
           fig, width = 26, height = 9, units = "cm")
  }
}

# center panel: contribution of species correlations across state space ------------------------------
# creating discrete variable for contribution of species correlations
full_df$correlation_contrib_disc <- "Intermediate"
full_df$correlation_contrib_disc[log(full_df$correlation_contrib) > quantile(log(full_df$correlation_contrib), probs = 0.75)] <- "Low"
full_df$correlation_contrib_disc[log(full_df$correlation_contrib) < quantile(log(full_df$correlation_contrib), probs = 0.25)] <- "High"
# plot
fig <- ggplot(data = full_df, aes(x = rotifers, y = algae, 
                                  shape = correlation_contrib_disc, fill = correlation_contrib_disc, 
                                  alpha = correlation_contrib_disc)) +
  geom_point(size = 3.2) +
  scale_fill_manual(values = c("#567F55", "#D5D5D5", "#83C282"),
                    name = "Contribution of\nspecies correlations") +
  scale_shape_manual(values = c(21, 22, 24),
                     name = "Contribution of\nspecies correlations") +
  scale_alpha_manual(values = c(1, 0.5, 1),
                     name = "Contribution of\nspecies correlations") +
  coord_equal() +
  scale_x_continuous(limits = c(-2, 4)) +
  scale_y_continuous(limits = c(-2, 4)) +
  guides(fill = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4)),
         alpha = guide_legend(override.aes = list(size = 4))) +   
  xlab("Predator (rotifers)") +
  ylab("Prey (algae)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.8, "cm"),
        legend.position = "top")
if (save_plots) {
  if (normalize) {
    ggsave(paste("figs/fig5/fig5_", data, "_times_", t1, "_", t2, 
                 "_state_space_pairs_contrib_smap_normalized.pdf", sep = ""), 
           fig, width = 18, height = 18, units = "cm")
  } else {
    ggsave(paste("figs/fig5/fig5_", data, "_times_", t1, "_", t2, 
                 "_state_space_pairs_contrib_smap.pdf", sep = ""), 
           fig, width = 18, height = 18, units = "cm")
  }
}

# right panel: prey abundance as a function of contribution of species correlations ------------------------------
# t-test of prey abundance as a function of contribution of species correlations
model <- lm(algae ~ correlation_contrib_disc, 
            data = subset(full_df, correlation_contrib_disc == "High" | correlation_contrib_disc == "Low"))
summary(model)
t.test(full_df$algae[full_df$correlation_contrib_disc == "High"], full_df$algae[full_df$correlation_contrib_disc == "Low"])
# plot
fig <- ggplot(data = subset(full_df, correlation_contrib_disc == "High" | correlation_contrib_disc == "Low"),
              aes(x = correlation_contrib_disc, y = algae, color = correlation_contrib_disc)) +
  geom_boxplot(size = 1.5, outlier.size = 2.5) +
  scale_color_manual(values = c("#567F55", "#83C282")) +
  scale_y_continuous(limits = c(-2, 4)) +
  xlab("") +
  ylab("Prey (algae)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position = "none")
if (save_plots) {
  if (normalize) {
    ggsave(paste("figs/fig5/fig5_", data, "_times_", t1, "_", t2, 
                 "_boxplot_pairs_contrib_smap_normalized.pdf", sep = ""), 
           fig, width = 10, height = 18, units = "cm")
  } else {
    ggsave(paste("figs/fig5/fig5_", data, "_times_", t1, "_", t2, 
                 "_boxplot_pairs_contrib_smap.pdf", sep = ""), 
           fig, width = 10, height = 18, units = "cm")
  }
}

# bottom left panel: time series of contribution of species correlations ------------------------------
# taking only a subset of time series
full_df <- full_df[full_df$time > 20 & full_df$time < 135, ]
# plot
fig <- ggplot() +
  geom_line(data = full_df, aes(x = time, y = correlation_contrib),
            size = 0.7) +
  geom_point(data = full_df, aes(x = time, y = correlation_contrib, 
                                 fill = correlation_contrib_disc,
                                 shape = correlation_contrib_disc), size = 3.2) +
  scale_fill_manual(values = c("#567F55", "#D5D5D5", "#83C282")) +
  scale_shape_manual(values = c(21, 22, 24)) +
  xlab("Time (days)") +
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
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.size = unit(0.8, "cm"),
        legend.position = "top")
if (save_plots) {
  if (normalize) {
    ggsave(paste("figs/fig5/fig5_", data, "_times_", t1, "_", t2, 
                 "_ts_pairs_contrib_smap_normalized.pdf", sep = ""), 
           fig, width = 26, height = 9, units = "cm")
  } else {
    ggsave(paste("figs/fig5/fig5_", data, "_times_", t1, "_", t2, 
                 "_ts_pairs_contrib_smap.pdf", sep = ""), 
           fig, width = 26, height = 9, units = "cm")
  }
}

# plot of full decomposition of community sensitivity over time ------------------------------
# create data frame
plot_df <- gather(full_df[ , c("time", "community_sensitivity", "species_contrib", "correlation_contrib")], 
                  "variable", "value", -time)
# plot
fig <- ggplot(data = plot_df, 
       aes(x = time, y = value, group = variable, color = variable)) +
  scale_color_manual(values = c("#BF85BA", "#83C282", "#E3A15A")) + 
  geom_line(size = 1.2) +
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
