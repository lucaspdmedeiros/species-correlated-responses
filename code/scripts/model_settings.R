# Parameters and integration settings for different population dynamics models

# 2- and 4-species Lotka-Volterra ------------------------------
if (func_name == "lotka_volterra") {
  func <- lotka_volterra
  # predator prey neutral cycle
  if (n_sp == 2) {
    r <- c(0.2, -0.2)
    A <- matrix(c(0, -0.15,
                  0.1, 0),
                nrow = 2, ncol = 2, byrow = TRUE)
    # initial condition and time steps
    state <- c(x1 = 2, x2 = 1)
    n_points <- 100000
    time_step <- 0.05
    times <- seq(0, time_step * n_points, by = time_step)
  }
  # 4 sp competition chaos (Vano et al 2006)
  if (n_sp == 4) {
    r <- c(1, 0.72, 1.53, 1.27)
    A <- matrix(c(-1, -1.09, -1.52, 0, 
                  0, -1, -0.44, -1.36,
                  -2.33, 0, -1, -0.47,
                  -1.21, -0.51, -0.35, -1),
                nrow = 4, ncol = 4, byrow = TRUE)
    A <- A * r
    # initial condition and time steps
    state <- c(x1 = 0.5, x2 = 0.5, x3 = 0.5, x4 = 0.5)
    n_points <- 100000
    time_step <- 0.05
    times <- seq(0, time_step * n_points, by = time_step)
  }
  # model parameters
  parms <- c(r, A)
  names(parms) <- c(paste("r", 1:n_sp, sep = ""),
                    paste(paste("a", 1:n_sp, sep = ""), 
                          rep(1:n_sp, each = n_sp), sep = ""))
}

# 2-species predator-prey (Yodzis 1989 Introduction to theoretical ecology) ------------------------------
if (func_name == "predator_prey") {
  func <- predator_prey
  # model parameters
  parms <- c(k = 0.5, a = 0.002, h = 4, d = 0.1, r = 0.5, K = 100)
  # initial condition and time steps
  state <- c(x1 = 30, x2 = 50)
  n_points <- 100000
  time_step <- 0.05
  times <- seq(0, time_step * n_points, by = time_step)
}

# 3-species food chain (Upadhyay 2000 Mathematical and Computer Modelling) ------------------------------
if (func_name == "food_chain") {
  func <- food_chain
  # model parameters
  parms <- c(r = 4.3, k = 50, a1 = 0.1, b1 = 0.1, a2 = 0.1, b2 = 0.1, 
             s = 1, h = 0.05, l = 1, n = 0.03)
  # initial condition and time steps
  state <- c(x1 = 30, x2 = 30, x3 = 30)
  n_points <- 100000
  time_step <- 0.05
  times <- seq(0, time_step * n_points, by = time_step)
}
