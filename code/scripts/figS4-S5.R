# Code for Figs S4 and S5: relationship between log of determinant of a correlation matrix 
# and the maximum squared correlation as well as matrix dimension for randomly generated 
# correlation matrices

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
if(!require(clusterGeneration)) {install.packages("clusterGeneration"); library(clusterGeneration)}

# settings ------------------------------
# set seed to replicate figure
set.seed(2)
# whether to save plots
save_plots <- FALSE
# number of random correlation matrices to generate
n <- 500
# maximum dimension
S_max <- 10
# to store results
full_df <- data.frame()

# generate matrices, compute correlation between the two variables, and plot ------------------------------
for (i in 2:S_max) {
  max_cor_A <- c()
  det_A <- c()
  for (j in 1:n) {
    # generate matrix
    A <- rcorrmatrix(i, alphad = 1)
    # compute maximum squared correlation
    max_cor_A[j] <- max(A[row(A)!=col(A)]^2)
    # compute determinant
    det_A[j] <- det(A)
  }
  # print correlation between the two variables
  print(cor(max_cor_A, log(det_A)))
  # store results
  df <- data.frame(S = i,
                   matrix = 1:n,
                   max_cor = max_cor_A, 
                   det = det_A)
  full_df <- rbind(full_df, df)
}
# plot of determinant vs maximum correlation
fig <- ggplot(data = full_df, aes(x = max_cor, y = det)) +
  geom_point(size = 2, fill = "#929292", shape = 21) +
  facet_wrap( ~ S, nrow = 3, labeller = "label_both") +
  xlab(expression("Maximum squared correlation (max("~rho[ij]^2~"))")) +
  ylab(expression("Determinant of correlation matrix (|"~bold(P)~"|) in log")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = "white", size = 1.5),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
# save plot
if (save_plots) {
  ggsave("figs/figSI_det_cor_mat_vs_max_cor_random_mats.pdf", 
         fig, width = 28, height = 24, units = "cm")
}
# plot of determinant vs number of species
fig <- ggplot() +
  geom_jitter(data = full_df, aes(x = S, y = det),
              size = 1.5, alpha = 0.2, width = 0.25) +
  geom_boxplot(data = full_df, aes(x = S, y = det, group = S),
               size = 1, alpha = 0.4, outlier.shape = NA) +
  xlab("Number of species (S)") +
  ylab(expression("Determinant of correlation matrix (|"~bold(P)~"|) in log")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l") +
  scale_x_continuous(breaks = 2:S_max) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.5),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = "white", size = 1.5),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
# save plot
if (save_plots) {
  ggsave("figs/figSI_det_cor_mat_vs_number_species.pdf", 
         fig, width = 18, height = 16, units = "cm")
}
