# Data treatment for the predator-prey time series from Blasius et al 2020

# cleaning wd, loading packages ------------------------------ 
rm(list = ls(all = TRUE))
if(!require(stats)) {install.packages("stats"); library(stats)}

# loading time series ------------------------------ 
# empirical data
data <- "blasius_2020_C1"
# load data
ts <- read.csv(paste("data/empirical_time_series/", data, ".csv", sep = ""), header = TRUE)

# fitting cubic hermite interpolation ------------------------------ 
# specify time points for spline to go through (same number of points as original time series)
xout <- seq(ts$time[1], ts$time[length(ts$time)], length = nrow(ts))
# fit cubic hermite interpolation for algae
algae_spline <- spline(x = ts$time, y = ts$algae, method = "fmm", xout = xout)
# fit cubic hermite interpolation for rotifers
rotifers_spline <- spline(x = ts$time, y = ts$rotifers, method = "fmm", xout = xout)
# plot time series with splines
par(mfrow = c(2, 1), mgp = c(2, 0.8, 0), mar = 0.1 + c(3, 3, 3, 1))
plot(ts$time, ts$algae, xlab = "Time", ylab = "Algae")
lines(algae_spline, col = "red")
plot(ts$time, ts$rotifers, xlab = "Time", ylab = "Rotifers")
lines(spline(ts$time, ts$rotifers), col = "blue")

# save time series ------------------------------ 
# interpolated data
new_ts <- data.frame(time = xout,
                     algae = algae_spline$y,
                     rotifers = rotifers_spline$y)
# save as csv
write.csv(new_ts, paste("data/empirical_time_series/", data, "_interpolated.csv", sep = ""), row.names = FALSE)
