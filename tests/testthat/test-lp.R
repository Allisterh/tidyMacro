test_that("fLP_cpp recovers a known impulse response", {
  set.seed(21)
  n     <- 600
  shock <- matrix(rnorm(n), ncol = 1)
  # y_t = 0.8 * shock_t + 0.4 * shock_{t-1} + noise  => IRF(0) = 0.8, IRF(1) = 0.4
  y <- 0.8 * shock + 0.4 * rbind(0, shock[-n, , drop = FALSE]) +
       matrix(rnorm(n, sd = 0.2), ncol = 1)

  res <- tidyMacro:::fLP_cpp(Y = y, X = shock, H = 3L, shock_col = 0L,
                             conf_level = 0.90, nw_lags_base = 0L,
                             n_threads = NT)
  irf <- as.numeric(res$irfs)

  expect_length(irf, 4L)
  expect_equal(irf[1], 0.8, tolerance = 0.1)
  expect_equal(irf[2], 0.4, tolerance = 0.1)
  expect_equal(irf[4], 0.0, tolerance = 0.1)
  expect_true(all(as.numeric(res$irfs_lower) <= irf + 1e-12))
  expect_true(all(as.numeric(res$irfs_upper) >= irf - 1e-12))
})

test_that("fLP_cpp impact coefficient equals the OLS slope (engine adds an intercept)", {
  set.seed(22)
  n <- 300
  x <- matrix(rnorm(n), ncol = 1)
  y <- matrix(1.5 * x + rnorm(n), ncol = 1)
  res <- tidyMacro:::fLP_cpp(y, x, H = 0L, shock_col = 0L, conf_level = 0.90,
                             nw_lags_base = 0L, n_threads = NT)
  ols <- unname(coef(lm(y ~ x))[2])
  expect_equal(as.numeric(res$irfs)[1], ols, tolerance = 1e-8)
})

test_that("fLP_cpp gives identical results for 1 and 2 threads", {
  skip_on_cran()
  skip_if(parallel::detectCores() < 2L)
  set.seed(23)
  n <- 400
  x <- matrix(rnorm(n), ncol = 1)
  y <- matrix(0.7 * x + rnorm(n), ncol = 1)
  a <- tidyMacro:::fLP_cpp(y, x, H = 8L, shock_col = 0L, conf_level = 0.9,
                           nw_lags_base = 0L, n_threads = 1L)
  b <- tidyMacro:::fLP_cpp(y, x, H = 8L, shock_col = 0L, conf_level = 0.9,
                           nw_lags_base = 0L, n_threads = 2L)
  expect_equal(a$irfs, b$irfs, tolerance = 1e-12)
  expect_equal(a$irfs_se, b$irfs_se, tolerance = 1e-12)
})
