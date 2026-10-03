# Data generators with known truth, shared by the LP / panel / DiD tests.

# Time series: y_t = 0.8 x_t + 0.4 x_{t-1} + e_t   (IRF: 0.8, 0.4, 0, 0, ...)
sim_lp_df <- function(n = 400L, seed = 1L) {
  set.seed(seed)
  x <- rnorm(n)
  y <- 0.8 * x + 0.4 * c(0, x[-n]) + rnorm(n, sd = 0.3)
  data.frame(y = y, x = x, z = rnorm(n))
}

# Endogenous treatment with a valid, strong instrument: true effect = 1.2,
# OLS is biased upward because u enters both D and y.
sim_iv_df <- function(n = 600L, seed = 2L) {
  set.seed(seed)
  z <- rnorm(n); u <- rnorm(n)
  D <- 0.9 * z + 0.5 * u + rnorm(n, sd = 0.5)
  y <- 1.2 * D + 0.8 * u + rnorm(n)
  data.frame(y = y, D = D, z = z)
}

# Balanced panel with unit and time effects: contemporaneous effect 0.7.
sim_panel_df <- function(N = 40L, Tt = 30L, seed = 3L) {
  set.seed(seed)
  d <- expand.grid(unit = seq_len(N), tt = seq_len(Tt))
  d <- d[order(d$unit, d$tt), ]
  d$shock <- rnorm(nrow(d))
  d$size  <- rnorm(nrow(d))
  d$y <- 0.7 * d$shock + rnorm(N)[d$unit] + rnorm(Tt)[d$tt] +
         rnorm(nrow(d), sd = 0.5)
  d
}

# Staggered absorbing adoption (some never treated); constant effect `tau`.
sim_did_df <- function(N = 150L, Tt = 20L, tau = 1, seed = 4L) {
  set.seed(seed)
  d <- expand.grid(id = seq_len(N), year = seq_len(Tt))
  d <- d[order(d$id, d$year), ]
  adopt <- sample(c(8:14, Inf), N, replace = TRUE, prob = c(rep(1, 7), 3))
  d$treat <- as.numeric(d$year >= adopt[d$id])
  d$y <- rnorm(N)[d$id] + 0.1 * d$year + tau * d$treat +
         rnorm(nrow(d), sd = 0.5)
  d
}

# Bivariate stable VAR(1) data as a named matrix.
sim_sr_data <- function(Tn = 250L, seed = 5L) {
  set.seed(seed)
  B <- rbind(c(0.5, 0.1), c(0.2, 0.4))
  y <- matrix(0, Tn, 2)
  for (t in 2:Tn) y[t, ] <- B %*% y[t - 1, ] + rnorm(2)
  colnames(y) <- c("a", "b")
  y
}
