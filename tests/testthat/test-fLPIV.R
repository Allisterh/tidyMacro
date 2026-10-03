iv_fit <- function(df = sim_iv_df(), H = 2L, ...) {
  fLPIV(y ~ D, instruments = ~ z, data = df, endog = "D", horizons = H,
        cumulative = FALSE, nw_lags_iv = 0L, n_threads = NT, ...)
}

test_that("fLPIV returns an object of class c('fLPIV', 'fLP')", {
  iv <- iv_fit()
  expect_s3_class(iv, "fLPIV")
  expect_s3_class(iv, "fLP")
  expect_true(all(c("irfs", "irfs_se", "irfs_lower", "irfs_upper",
                    "Fstat_fs", "rsqr_fs", "endog", "instrument_vars") %in%
                    names(iv)))
  expect_equal(dim(iv$irfs), c(3L, 1L))
  expect_equal(iv$endog, "D")
})

test_that("fLPIV impact effect equals textbook just-identified 2SLS", {
  df <- sim_iv_df()
  iv <- iv_fit(df)
  # with an intercept, 2SLS slope = cov(z, y) / cov(z, D)
  expect_equal(unname(iv$irfs[1, 1]), cov(df$z, df$y) / cov(df$z, df$D),
               tolerance = 1e-8)
})

test_that("fLPIV removes the endogeneity bias that OLS has", {
  df  <- sim_iv_df()
  iv  <- iv_fit(df)
  ols <- unname(coef(lm(y ~ D, data = df))[2])
  expect_lt(abs(iv$irfs[1, 1] - 1.2), abs(ols - 1.2))
  expect_equal(unname(iv$irfs[1, 1]), 1.2, tolerance = 0.15)
})

test_that("fLPIV first-stage diagnostics are sensible for a strong instrument", {
  iv <- iv_fit()
  expect_length(iv$Fstat_fs, 3L)
  expect_length(iv$rsqr_fs, 3L)
  expect_true(all(iv$Fstat_fs > 10))
  expect_true(all(iv$rsqr_fs > 0 & iv$rsqr_fs < 1))
})

test_that("fLPIV bands bracket the estimate", {
  iv <- iv_fit()
  expect_true(all(iv$irfs_lower < iv$irfs & iv$irfs < iv$irfs_upper))
  expect_true(all(iv$irfs_se > 0))
})

test_that("fLPIV shock_size = 'sd' scales by the standard deviation of D", {
  df <- sim_iv_df()
  u <- iv_fit(df, H = 1L)
  s <- iv_fit(df, H = 1L, shock_size = "sd")
  expect_equal(unname(s$irfs[1, 1] / u$irfs[1, 1]), sd(df$D), tolerance = 1e-8)
})

test_that("fLPIV defaults to cumulative responses", {
  df <- sim_iv_df()
  cu <- fLPIV(y ~ D, instruments = ~ z, data = df, endog = "D",
              horizons = 2L, nw_lags_iv = 0L, n_threads = NT)
  expect_true(isTRUE(cu$cumulative))
  expect_false(isTRUE(all.equal(unname(cu$irfs), unname(iv_fit(df)$irfs))))
})

test_that("fLP-family methods work on fLPIV objects", {
  iv <- iv_fit()
  td <- tidy.fLP(iv)
  expect_s3_class(td, "data.frame")
  expect_equal(nrow(td), 3L)
  expect_equal(coef(iv), iv$irfs)
  expect_output(print(iv))
})

test_that("fLPIV validates endog", {
  df <- sim_iv_df()
  expect_error(fLPIV(y ~ D, instruments = ~ z, data = df, horizons = 2L,
                     n_threads = NT), "'endog' is required")
  expect_error(fLPIV(y ~ D, instruments = ~ z, data = df, endog = "zz",
                     horizons = 2L, n_threads = NT), "not on the RHS")
})

test_that("fLPIV is invariant to the number of threads", {
  skip_on_cran()
  df <- sim_iv_df()
  a <- iv_fit(df, H = 5L)
  b <- fLPIV(y ~ D, instruments = ~ z, data = df, endog = "D", horizons = 5L,
             cumulative = FALSE, nw_lags_iv = 0L, n_threads = 2L)
  expect_equal(a$irfs, b$irfs, tolerance = 1e-12)
})
