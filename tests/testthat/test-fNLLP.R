# State-dependent and nonlinear DGP with known truth:
#   y_t = (0.9 m_t + 0.1 (1 - m_t)) s_t + b3 s_t^3 + 0.5 x_t + e_t
sim_nllp_df <- function(n = 400, seed = 4, b3 = 0.3) {
  set.seed(seed)
  s <- rnorm(n)
  x <- rnorm(n)
  m <- plogis(2 * rnorm(n))
  y <- (0.9 * m + 0.1 * (1 - m)) * s + b3 * s^3 + 0.5 * x + rnorm(n, sd = 0.3)
  data.frame(y = y, s = s, x = x, m = m, trend = seq_len(n))
}

test_that("fNLLP state design equals fLP on explicit interaction columns", {
  df <- sim_nllp_df()
  df$ms <- df$m * df$s; df$ls <- (1 - df$m) * df$s
  df$mx <- df$m * df$x; df$lx <- (1 - df$m) * df$x

  f <- fNLLP(y ~ s + x + trend, data = df, horizons = 3, shock = "s",
             type = "state", state = "m", common = "trend", n_threads = NT)
  # [1, m] spans the same space as [m, 1 - m], so the shock coefficients and
  # their HAC SEs are unchanged.
  for (k in 1:2) {
    sh <- c("ms", "ls")[k]
    g <- fLP(y ~ ms + ls + m + mx + lx + trend, data = df, horizons = 3,
             shock = sh, n_threads = NT)
    expect_equal(unname(f$irfs[, 1, k]), unname(g$irfs[, 1]), tolerance = 1e-10)
    expect_equal(unname(f$irfs_se[, 1, k]), unname(g$irfs_se[, 1]), tolerance = 1e-8)
  }
  # Truth recovery needs a DGP without the cubic term the state design omits.
  g <- fNLLP(y ~ s + x, data = sim_nllp_df(b3 = 0), horizons = 0, shock = "s",
             type = "state", state = "m", n_threads = NT)
  expect_equal(unname(g$irfs[1, 1, ]), c(0.9, 0.1), tolerance = 0.1)
  expect_equal(f$components, c("high", "low"))
  expect_equal(f$design_vars,
               c("m:(Intercept)", "m:s", "m:x", "(1-m):(Intercept)", "(1-m):s",
                 "(1-m):x", "trend"))
})

test_that("fNLLP sign and cubic designs equal fLP on explicit columns", {
  df <- sim_nllp_df()
  df$sp <- pmax(df$s, 0); df$sn <- pmin(df$s, 0); df$s3 <- df$s^3

  a <- fNLLP(y ~ s + x, data = df, horizons = 3, shock = "s", type = "sign",
             n_threads = NT)
  b <- fNLLP(y ~ s + x, data = df, horizons = 3, shock = "s", type = "cubic",
             n_threads = NT)
  for (k in 1:2) {
    ga <- fLP(y ~ sp + sn + x, data = df, horizons = 3,
              shock = c("sp", "sn")[k], n_threads = NT)
    gb <- fLP(y ~ s + s3 + x, data = df, horizons = 3,
              shock = c("s", "s3")[k], n_threads = NT)
    expect_equal(unname(a$irfs[, 1, k]), unname(ga$irfs[, 1]), tolerance = 1e-10)
    expect_equal(unname(a$irfs_se[, 1, k]), unname(ga$irfs_se[, 1]), tolerance = 1e-8)
    expect_equal(unname(b$irfs[, 1, k]), unname(gb$irfs[, 1]), tolerance = 1e-10)
    expect_equal(unname(b$irfs_se[, 1, k]), unname(gb$irfs_se[, 1]), tolerance = 1e-8)
  }
  expect_equal(unname(b$irfs[1, 1, 2]), 0.3, tolerance = 0.1)
  expect_true(all(is.nan(b$diff)))
})

test_that("fNLLP HAC covariance matches a full sandwich, fast and full paths", {
  skip_if_not_installed("sandwich")
  df <- sim_nllp_df()
  f0 <- fNLLP(y ~ s + x, data = df, horizons = 2, shock = "s", type = "state",
              state = "m", nw_lags = 2, n_threads = NT)
  f1 <- fNLLP(y ~ s + x, data = df, horizons = 2, shock = "s", type = "state",
              state = "m", nw_lags = 2, store_full = TRUE, n_threads = NT)
  n <- nrow(df)
  for (h in 0:2) {
    r  <- 1:(n - h)
    yy <- df$y[r + h]; m <- df$m[r]; s <- df$s[r]; x <- df$x[r]
    fit <- lm(yy ~ 0 + m + I(m * s) + I(m * x) + I(1 - m) + I((1 - m) * s) +
                I((1 - m) * x))
    V <- sandwich::NeweyWest(fit, lag = 2 + h + 1, prewhite = FALSE,
                             adjust = FALSE)[c(2, 5), c(2, 5)]
    for (f in list(f0, f1)) {
      expect_equal(unname(f$irfs_se[h + 1, 1, ]), unname(sqrt(diag(V))), tolerance = 1e-8)
      expect_equal(unname(f$irfs_cov[h + 1, 1]), V[1, 2], tolerance = 1e-8)
      expect_equal(unname(f$diff_se[h + 1, 1]),
                   sqrt(V[1, 1] + V[2, 2] - 2 * V[1, 2]), tolerance = 1e-8)
    }
  }
})

test_that("fNLLP balanced sample fixes the regressor dates", {
  df <- sim_nllp_df(n = 200)
  df$s[c(1:10, 191:200)] <- NA          # shock observed on rows 11-190
  f <- suppressMessages(fNLLP(y ~ s + x, data = df, horizons = 5, shock = "s",
                              type = "sign", balanced = TRUE, n_threads = NT))
  expect_equal(unname(f$nobs), rep(180, 6))
  h <- 5; r <- 11:190
  ref <- coef(lm(df$y[r + h] ~ pmax(df$s[r], 0) + pmin(df$s[r], 0) + df$x[r]))
  expect_equal(unname(f$irfs[h + 1, 1, ]), unname(ref[2:3]), tolerance = 1e-10)

  # Data end 10 rows after the last shock: with H = 11 the last regressor
  # date is dropped so y_{t+11} exists for every date kept.
  expect_message(
    g <- fNLLP(y ~ s + x, data = df, horizons = 11, shock = "s", type = "sign",
               balanced = TRUE, n_threads = NT),
    "dropping regressor dates after row 189")
  expect_equal(unname(g$nobs), rep(179, 12))

  # Fully observed data: the last H rows are reserved for the outcomes.
  full <- sim_nllp_df(n = 200)
  g <- suppressMessages(fNLLP(y ~ s + x, data = full, horizons = 4, shock = "s",
                              type = "sign", balanced = TRUE, n_threads = NT))
  expect_equal(unname(g$nobs), rep(196, 5))

  # A gap inside the outcome is an error, not a silent truncation.
  full$y[100] <- NA
  expect_error(suppressMessages(fNLLP(y ~ s + x, data = full, horizons = 4,
                                      shock = "s", type = "sign", balanced = TRUE)),
               "observed without gaps")
})

test_that("fNLLP and fLP reject formula syntax they would not estimate", {
  df <- sim_nllp_df()
  df$ypos <- exp(df$y)
  expect_error(fNLLP(log(ypos) ~ s + x, data = df, shock = "s", type = "sign"),
               "LHS must be a variable name")
  expect_error(fNLLP(y ~ s + x - 1, data = df, shock = "s", type = "sign"),
               "intercept is always included")
  expect_error(fNLLP(y ~ s + I(x^2), data = df, shock = "s", type = "sign"),
               "RHS terms must be columns")
  expect_error(fLP(log(ypos) ~ s + x, data = df, shock = "s"),
               "LHS must be a variable name")
  f <- fNLLP(c(y, x) ~ s + trend, data = df, horizons = 1, shock = "s",
             type = "sign", n_threads = NT)
  expect_equal(f$lhs_vars, c("y", "x"))
})

test_that("generated lag names do not reuse an unrelated column", {
  df <- sim_nllp_df()
  df$s_l1 <- rnorm(nrow(df))
  expect_error(fNLLP(y ~ s + l(s, 1), data = df, shock = "s", type = "sign"),
               "not lag 1 of 's'")
  df$s_l1 <- c(NA, head(df$s, -1))       # the genuine lag is accepted
  expect_s3_class(suppressMessages(fNLLP(y ~ s + l(s, 1), data = df, horizons = 1,
                                         shock = "s", type = "sign")), "fNLLP")
})

test_that("the QR fallback rank test does not depend on units", {
  set.seed(9)
  n  <- 300
  df <- data.frame(s = rnorm(n), x = rnorm(n), y = rnorm(n))
  df$x2 <- df$x + 1e-6 * rnorm(n)        # near-collinear but full rank
  df$xb <- df$x * 1e12
  a <- fNLLP(y ~ s + x + x2, data = df, horizons = 1, shock = "s",
             type = "sign", store_full = TRUE, n_threads = NT)
  b <- fNLLP(y ~ s + xb + x2, data = df, horizons = 1, shock = "s",
             type = "sign", store_full = TRUE, n_threads = NT)
  expect_equal(b$irfs, a$irfs, tolerance = 1e-8)
  expect_equal(unname(b$betas$h0["xb", 1]) * 1e12, unname(a$betas$h0["x", 1]),
               tolerance = 1e-6)
  df$x3 <- df$x + df$x2                   # exact collinearity still stops
  expect_error(fNLLP(y ~ s + x + x2 + x3, data = df, horizons = 1, shock = "s",
                     type = "sign"), "rank deficient")
})

test_that("a negative plot scale keeps lower below upper", {
  df <- sim_nllp_df()
  f <- fNLLP(y ~ s + x, data = df, horizons = 2, shock = "s", type = "sign",
             n_threads = NT)
  p <- fPlotNLLP(f, scale = -1, return_data = TRUE)
  expect_true(all(p$lower <= p$upper))
  expect_equal(p$point, -tidy.fNLLP(f, difference = FALSE)$estimate)
})

test_that("predict.fNLLP combines the two coefficients with their covariance", {
  df <- sim_nllp_df()
  fs <- fNLLP(y ~ s + x, data = df, horizons = 2, shock = "s", type = "state",
              state = "m", n_threads = NT)
  p <- predict(fs, size = 2, state = c(1, 0))
  expect_equal(p$estimate[p$state == 1], 2 * unname(fs$irfs[, 1, 1]))
  expect_equal(p$se[p$state == 0], 2 * unname(fs$irfs_se[, 1, 2]))

  fc <- fNLLP(y ~ s + x, data = df, horizons = 2, shock = "s", type = "cubic",
              n_threads = NT)
  p <- predict(fc, size = 2)
  v <- 4 * fc$irfs_se[, 1, 1]^2 + 64 * fc$irfs_se[, 1, 2]^2 + 32 * fc$irfs_cov[, 1]
  expect_equal(p$estimate, unname(2 * fc$irfs[, 1, 1] + 8 * fc$irfs[, 1, 2]))
  expect_equal(p$se, unname(sqrt(v)))
  expect_error(predict(fs, size = 1), "needs 'state'")
})

test_that("fNLLP validates its arguments", {
  df <- sim_nllp_df()
  expect_error(fNLLP(y ~ s + x, data = df, shock = "s", type = "state"),
               "needs 'state'")
  expect_error(fNLLP(y ~ s + x, data = df, shock = "s", type = "sign", state = "m"),
               "only used with type")
  expect_error(fNLLP(y ~ s + x, data = df, shock = "s", type = "state",
                     state = "m", common = "s"), "cannot be a common")
  expect_error(fNLLP(y ~ s + x, data = df, shock = "s", type = "state",
                     state = "m", common = "trend"), "not on the RHS")
  df$bad <- df$m * 2
  expect_error(fNLLP(y ~ s + x, data = df, shock = "s", type = "state",
                     state = "bad"), "must lie in \\[0, 1\\]")
})
