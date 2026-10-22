test_that("fLP returns an fLP object with the documented structure", {
  df <- sim_lp_df()
  f  <- fLP(y ~ x, data = df, horizons = 3L, shock = "x", n_threads = NT)

  expect_s3_class(f, "fLP")
  expect_true(all(c("irfs", "irfs_upper", "irfs_lower", "irfs_se", "lhs_vars",
                    "shock", "horizons", "conf", "nobs") %in% names(f)))
  expect_equal(dim(f$irfs), c(4L, 1L))
  expect_equal(rownames(f$irfs), as.character(0:3))
  expect_equal(f$lhs_vars, "y")
  expect_equal(f$shock, "x")
  expect_equal(f$horizons, 0:3)
  expect_equal(dim(f$nobs), c(4, 1))
  expect_equal(as.vector(f$nobs), nrow(df) - 0:3)
  expect_equal(coef(f), f$irfs)
})

test_that("fLP matches OLS of y_{t+h} on the shock at every horizon", {
  df <- sim_lp_df()
  n  <- nrow(df)
  f  <- fLP(y ~ x, data = df, horizons = 3L, shock = "x", n_threads = NT)

  for (h in 0:3) {
    ref <- unname(coef(lm(df$y[(1 + h):n] ~ df$x[1:(n - h)]))[2])
    expect_equal(unname(f$irfs[h + 1, 1]), ref, tolerance = 1e-8)
  }
})

test_that("fLP recovers a known impulse response", {
  f <- fLP(y ~ x, data = sim_lp_df(), horizons = 3L, shock = "x",
           n_threads = NT)
  expect_equal(unname(f$irfs[, 1]), c(0.8, 0.4, 0, 0), tolerance = 0.1)
})

test_that("fLP bands bracket the estimate and widen with the confidence level", {
  df <- sim_lp_df()
  f  <- fLP(y ~ x, data = df, horizons = 3L, shock = "x", n_threads = NT)
  expect_true(all(f$irfs_lower < f$irfs & f$irfs < f$irfs_upper))
  expect_true(all(f$irfs_se > 0))

  m <- fLP(y ~ x, data = df, horizons = 3L, shock = "x", conf = c(68, 95),
           n_threads = NT)
  expect_true(is.list(m$irfs_upper))
  expect_setequal(names(m$irfs_upper), c("68", "95"))
  expect_true(all(m$irfs_upper[["95"]] > m$irfs_upper[["68"]]))
  expect_true(all(m$irfs_lower[["95"]] < m$irfs_lower[["68"]]))
})

test_that("fLP lag operator and macros give identical results to explicit terms", {
  df <- sim_lp_df()
  a <- suppressMessages(
    fLP(y ~ x + l(y, 1:2), data = df, horizons = 2L, shock = "x",
        n_threads = NT))
  expect_equal(as.vector(a$nobs), nrow(df) - 2 - 0:2)

  df$y_l1 <- c(NA, df$y[-nrow(df)])
  df$y_l2 <- c(NA, NA, df$y[seq_len(nrow(df) - 2L)])
  b <- suppressMessages(
    fLP(y ~ x + y_l1 + y_l2, data = df, horizons = 2L, shock = "x",
        n_threads = NT))
  expect_equal(unname(a$irfs), unname(b$irfs), tolerance = 1e-10)

  ctrl <- c("y_l1", "y_l2")
  m <- suppressMessages(
    fLP(y ~ x + ..ctrl, data = df, horizons = 2L, shock = "x",
        n_threads = NT))
  expect_equal(unname(m$irfs), unname(b$irfs), tolerance = 1e-10)
})

test_that("fLP supports several outcome variables", {
  df <- sim_lp_df()
  df$y2 <- 0.3 * df$x + rnorm(nrow(df))
  m <- fLP(c(y, y2) ~ x, data = df, horizons = 2L, shock = "x",
           n_threads = NT)
  expect_equal(dim(m$irfs), c(3L, 2L))
  expect_equal(colnames(m$irfs), c("y", "y2"))
  expect_equal(m$lhs_vars, c("y", "y2"))
})

test_that("fLP shock_size = 'sd' rescales by the sample standard deviation", {
  df <- sim_lp_df()
  u <- fLP(y ~ x, data = df, horizons = 1L, shock = "x", n_threads = NT)
  s <- fLP(y ~ x, data = df, horizons = 1L, shock = "x", shock_size = "sd",
           n_threads = NT)
  expect_equal(unname(s$irfs[, 1] / u$irfs[, 1]), rep(sd(df$x), 2),
               tolerance = 1e-8)
})

test_that("fLP cumulative = TRUE changes the outcome to a long difference", {
  df <- sim_lp_df()
  n  <- nrow(df)
  cu <- fLP(y ~ x, data = df, horizons = 3L, shock = "x", cumulative = TRUE,
            n_threads = NT)
  # h = 0: y_t - y_{t-1} regressed on x_t
  ref <- unname(coef(lm(diff(df$y) ~ df$x[-1]))[2])
  expect_equal(unname(cu$irfs[1, 1]), ref, tolerance = 1e-8)
  expect_true(isTRUE(cu$cumulative))
})

test_that("fLP store_full keeps full coefficient matrices", {
  f <- fLP(y ~ x, data = sim_lp_df(), horizons = 2L, shock = "x",
           store_full = TRUE, n_threads = NT)
  expect_true(all(c("betas", "ses") %in% names(f)))
})

test_that("tidy.fLP returns one row per horizon and outcome", {
  df <- sim_lp_df()
  f  <- fLP(y ~ x, data = df, horizons = 3L, shock = "x", n_threads = NT)
  td <- tidy.fLP(f)
  expect_s3_class(td, "data.frame")
  expect_equal(nrow(td), 4L)
  expect_true(all(c("horizon", "lhs", "shock", "estimate", "se", "lower",
                    "upper") %in% names(td)))
  expect_equal(td$estimate, as.vector(f$irfs))
  expect_equal(td$horizon, 0:3)
})

test_that("print.fLP produces output and returns its input invisibly", {
  f <- fLP(y ~ x, data = sim_lp_df(), horizons = 2L, shock = "x",
           n_threads = NT)
  expect_output(print(f), "Local Projections")
  expect_invisible(print(f))
})

test_that("fLP validates its arguments", {
  df <- sim_lp_df()
  expect_error(fLP(y ~ x, data = df, horizons = 3L, n_threads = NT),
               "'shock' is required")
  expect_error(fLP(y ~ x, data = df, horizons = 3L, shock = "zz",
                   n_threads = NT), "not found on the RHS")
  expect_error(fLP(y ~ x, data = df, horizons = -1L, shock = "x",
                   n_threads = NT), "non-negative integer")
  expect_error(fLP(y ~ x, data = df, horizons = 3L, shock = "x", conf = 150,
                   n_threads = NT), "confidence levels")
  expect_error(fLP(y ~ x, data = df, horizons = 3L, shock = c("x", "z"),
                   n_threads = NT), "single character string")
})

test_that("fLP is invariant to the number of threads", {
  skip_on_cran()
  df <- sim_lp_df()
  a <- fLP(y ~ x, data = df, horizons = 6L, shock = "x", n_threads = 1L)
  b <- fLP(y ~ x, data = df, horizons = 6L, shock = "x", n_threads = 2L)
  expect_equal(a$irfs, b$irfs, tolerance = 1e-12)
  expect_equal(a$irfs_se, b$irfs_se, tolerance = 1e-12)
})

# Missing values: y_{t+h} is read from the full data before rows are dropped
# (Stata's F.h y), each outcome keeps its own sample, and nobs is per horizon.
sim_na_df <- function(T = 200, seed = 1) {
  set.seed(seed)
  e <- rnorm(T); y <- numeric(T)
  for (t in 2:T) y[t] <- 0.6 * y[t - 1] + 0.8 * e[t] + rnorm(1, sd = 0.5)
  data.frame(y = y, s = e)
}
lead_na <- function(x, h) if (h == 0) x else c(x[-(1:h)], rep(NA, h))
lag_na  <- function(x) c(NA, x[-length(x)])

test_that("fLP uses outcomes observed after the shock ends", {
  skip_if_not_installed("sandwich")
  df <- sim_na_df()
  df$s[191:200] <- NA
  H <- 6
  f <- suppressMessages(fLP(y ~ s + l(y, 1), data = df, horizons = H,
                            shock = "s", nw_lags = 2, n_threads = NT))
  for (h in 0:H) {
    m <- lm(lead_na(df$y, h) ~ df$s + lag_na(df$y))
    expect_equal(unname(f$irfs[h + 1, 1]), unname(coef(m)[2]), tolerance = 1e-10)
    se <- sqrt(sandwich::NeweyWest(m, lag = 2 + h + 1, prewhite = FALSE,
                                   adjust = FALSE)[2, 2])
    expect_equal(unname(f$irfs_se[h + 1, 1]), se, tolerance = 1e-8)
  }
  expect_equal(as.vector(f$nobs), rep(189, H + 1))

  cu <- suppressMessages(fLP(y ~ s + l(y, 1), data = df, horizons = H,
                             shock = "s", cumulative = TRUE, n_threads = NT))
  ref <- sapply(0:H, function(h)
    unname(coef(lm(I(lead_na(df$y, h) - lag_na(df$y)) ~ df$s + lag_na(df$y)))[2]))
  expect_equal(unname(cu$irfs[, 1]), ref, tolerance = 1e-10)
  expect_equal(as.vector(cu$nobs), rep(189, H + 1))
})

test_that("fLP gives each outcome its own sample", {
  df <- sim_na_df()
  df$y2 <- 0.3 * df$s + rnorm(nrow(df))
  df$y2[c(1:5, 150:200)] <- NA
  both <- suppressMessages(fLP(c(y, y2) ~ s + l(y, 1), data = df, horizons = 4,
                               shock = "s", n_threads = NT))
  one  <- suppressMessages(fLP(y ~ s + l(y, 1), data = df, horizons = 4,
                               shock = "s", n_threads = NT))
  expect_equal(both$irfs[, "y"], one$irfs[, "y"], tolerance = 1e-12)
  expect_equal(both$irfs_se[, "y"], one$irfs_se[, "y"], tolerance = 1e-12)
  ref <- sapply(0:4, function(h)
    unname(coef(lm(lead_na(df$y2, h) ~ df$s + lag_na(df$y)))[2]))
  expect_equal(unname(both$irfs[, "y2"]), ref, tolerance = 1e-10)
  expect_equal(colnames(both$nobs), c("y", "y2"))
  expect_equal(as.vector(both$nobs[, "y"]), 199 - 0:4)
  expect_equal(as.vector(both$nobs[, "y2"]),
               sapply(0:4, function(h) sum(complete.cases(lead_na(df$y2, h), df$s,
                                                          lag_na(df$y)))))
})

test_that("fLP keeps dates whose outcome is missing only at t", {
  df <- sim_na_df()
  df$y[1:10] <- NA
  f <- fLP(y ~ s, data = df, horizons = 3, shock = "s", n_threads = NT)
  ref <- sapply(0:3, function(h) unname(coef(lm(lead_na(df$y, h) ~ df$s))[2]))
  expect_equal(unname(f$irfs[, 1]), ref, tolerance = 1e-10)
  # Dates 11-h .. 200-h at horizon h; the old complete-case rule kept 11 .. 200-h.
  expect_equal(as.vector(f$nobs), rep(190, 4))
})

test_that("fLP rejects interior gaps and too-short horizon samples", {
  df <- sim_na_df()
  g <- df; g$s[100] <- NA
  expect_error(suppressMessages(fLP(y ~ s, data = g, horizons = 2, shock = "s")),
               "regressors have 1 interior gap")
  g <- df; g$y[100] <- NA
  expect_error(fLP(y ~ s, data = g, horizons = 2, shock = "s"),
               "outcome 'y' has 1 interior gap")
  expect_error(fLP(y ~ s, data = df[1:10, ], horizons = 8, shock = "s"),
               "outcome 'y' has 2 at horizon 8")
})
