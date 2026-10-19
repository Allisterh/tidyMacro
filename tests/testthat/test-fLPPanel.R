panel_fit <- function(d = sim_panel_df(), H = 2L, ...) {
  fLPPanel(y ~ shock | unit + tt, data = d, panel_id = c("unit", "tt"),
           shock = "shock", horizons = H, n_threads = NT, ...)
}

test_that("fLPPanel returns an object of class c('fLPPanel', 'fLP')", {
  p <- panel_fit()
  expect_s3_class(p, "fLPPanel")
  expect_s3_class(p, "fLP")
  expect_equal(dim(p$irfs), c(3L, 1L))
  expect_equal(dim(p$irfs_se), c(3L, 1L))
  expect_equal(colnames(p$irfs), "shock")
  expect_equal(coef(p), p$irfs)
})

test_that("fLPPanel impact effect equals the two-way fixed-effects OLS slope", {
  d  <- sim_panel_df()
  p  <- panel_fit(d)
  fe <- lm(y ~ shock + factor(unit) + factor(tt), data = d)
  expect_equal(unname(p$irfs[1, 1]), unname(coef(fe)["shock"]),
               tolerance = 1e-6)
})

test_that("fLPPanel recovers the true effect and a null response afterwards", {
  p <- panel_fit()
  expect_equal(unname(p$irfs[1, 1]), 0.7, tolerance = 0.05)
  expect_lt(abs(p$irfs[2, 1]), 0.15)
  expect_true(all(p$irfs_se > 0))
})

test_that("fLPPanel heterogeneous-shock interactions give one IRF per margin", {
  d <- sim_panel_df()
  d$y <- d$y + 0.5 * d$shock * d$size
  h <- fLPPanel(y ~ shock + shock:size | unit + tt, data = d,
                panel_id = c("unit", "tt"),
                shock = c("shock", "shock:size"), horizons = 1L,
                n_threads = NT)
  expect_equal(dim(h$irfs), c(2L, 2L))
  expect_equal(colnames(h$irfs), c("shock", "shock:size"))
  expect_equal(unname(h$irfs[1, ]), c(0.7, 0.5), tolerance = 0.1)
})

test_that("fLPPanel clustering choices change the standard errors", {
  d <- sim_panel_df()
  t_cl <- panel_fit(d)
  u_cl <- panel_fit(d, cluster = ~unit)
  two  <- panel_fit(d, cluster = ~unit + tt)
  expect_true(all(is.finite(c(t_cl$irfs_se, u_cl$irfs_se, two$irfs_se))))
  expect_false(isTRUE(all.equal(t_cl$irfs_se, u_cl$irfs_se)))
  expect_false(isTRUE(all.equal(u_cl$irfs_se, two$irfs_se)))
  # point estimates do not depend on how the errors are clustered
  expect_equal(t_cl$irfs, u_cl$irfs, tolerance = 1e-10)
  expect_equal(t_cl$irfs, two$irfs, tolerance = 1e-10)
})

test_that("fLPPanel ssc = 'IK' only applies to time clustering", {
  d <- sim_panel_df()
  base <- panel_fit(d)
  ss   <- panel_fit(d, ssc = "IK")
  expect_equal(base$irfs, ss$irfs, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(base$irfs_se, ss$irfs_se)))
  expect_error(panel_fit(d, cluster = ~unit, ssc = "IK"),
               "only for time clustering")
})

test_that("small_sample = TRUE is deprecated in favour of ssc = 'IK'", {
  d <- sim_panel_df()
  expect_warning(ss <- panel_fit(d, small_sample = TRUE),
                 "'small_sample = TRUE' is deprecated")
  expect_identical(ss$ssc, "IK")
  expect_equal(ss$irfs_se, panel_fit(d, ssc = "IK")$irfs_se)
  expect_error(panel_fit(d, small_sample = TRUE, ssc = "CR1"), "disagree")
})

test_that("ssc = 'CR1' matches the Stata / fixest finite-sample rule", {
  d <- sim_panel_df()
  for (cl in list(NULL, ~unit, ~unit + tt)) {
    base <- panel_fit(d, cluster = cl)
    cr1  <- panel_fit(d, cluster = cl, ssc = "CR1")

    # Only the variance changes, never the coefficient.
    expect_equal(base$irfs, cr1$irfs, tolerance = 1e-12)
    expect_true(all(cr1$irfs_se > base$irfs_se))

    # CR1 uses t degrees of freedom; the asymptotic mode uses the normal.
    # One-way clusters on its own dimension, two-way on the smaller one.
    G <- switch(cr1$cluster_type,
                time = length(unique(d$tt)),
                unit = length(unique(d$unit)),
                min(length(unique(d$unit)), length(unique(d$tt))))
    expect_true(all(is.finite(cr1$df)))
    expect_true(all(cr1$df <= G - 1))
    expect_true(all(is.infinite(base$df)))

    # Wider bands, from both the factor and the t critical value.
    expect_true(all((cr1$irfs_upper - cr1$irfs_lower) >
                    (base$irfs_upper - base$irfs_lower)))
  }
})

test_that("ssc = 'CR1' reproduces the documented one-way factor exactly", {
  d <- sim_panel_df()
  N <- nrow(d)

  # K = slopes + intercept + (L - 1) for every absorbed factor that is NOT
  # nested in the clustering dimension. Clustering on units makes the unit
  # fixed effect free, and vice versa.
  check <- function(cl, G, L_other) {
    base <- panel_fit(d, cluster = cl)
    cr1  <- panel_fit(d, cluster = cl, ssc = "CR1")
    K <- 1 + 1 + (L_other - 1)
    expect_equal(unname(cr1$irfs_se[1, 1] / base$irfs_se[1, 1]),
                 sqrt((G / (G - 1)) * ((N - 1) / (N - K))),
                 tolerance = 1e-8)
    expect_equal(unname(cr1$df[1, 1]), G - 1)
  }
  check(~tt,   length(unique(d$tt)),   length(unique(d$unit)))
  check(~unit, length(unique(d$unit)), length(unique(d$tt)))
})

test_that("fLPPanel cumulative projections differ from level projections", {
  d  <- sim_panel_df()
  lv <- panel_fit(d)
  cu <- panel_fit(d, cumulative = TRUE)
  expect_true(all(is.finite(cu$irfs)))
  expect_false(isTRUE(all.equal(lv$irfs, cu$irfs)))
})

test_that("tidy.fLPPanel returns one row per horizon", {
  p  <- panel_fit()
  td <- tidy.fLPPanel(p)
  expect_s3_class(td, "data.frame")
  expect_equal(nrow(td), 3L)
  expect_true(all(c("horizon", "shock", "estimate", "se", "pval", "lower",
                    "upper") %in% names(td)))
  expect_equal(td$estimate, as.vector(p$irfs))
})

test_that("print.fLPPanel runs", {
  expect_output(print(panel_fit()))
})

test_that("fLPPanel validates its arguments", {
  d <- sim_panel_df()
  expect_error(fLPPanel(y ~ shock | unit + tt, data = d,
                        panel_id = c("unit", "tt"), horizons = 2L,
                        n_threads = NT), "'shock' is required")
  expect_error(fLPPanel(y ~ shock | unit + tt, data = d, panel_id = "unit",
                        shock = "shock", horizons = 2L, n_threads = NT),
               "length-2 character vector")
})

test_that("fLPPanel is invariant to the number of threads", {
  skip_on_cran()
  d <- sim_panel_df()
  a <- fLPPanel(y ~ shock | unit + tt, data = d, panel_id = c("unit", "tt"),
                shock = "shock", horizons = 4L, n_threads = 1L)
  b <- fLPPanel(y ~ shock | unit + tt, data = d, panel_id = c("unit", "tt"),
                shock = "shock", horizons = 4L, n_threads = 2L)
  expect_equal(a$irfs, b$irfs, tolerance = 1e-12)
  expect_equal(a$irfs_se, b$irfs_se, tolerance = 1e-12)
})

test_that("fLPPanel rejects an existing x_l1 column that is not the lag of x", {
  d <- sim_panel_df()
  d$x <- rnorm(nrow(d))
  d$x_l1 <- 999
  expect_error(
    fLPPanel(y ~ shock + l(x, 1) | unit + tt, data = d,
             panel_id = c("unit", "tt"), shock = "shock", horizons = 1,
             n_threads = NT),
    "not lag 1"
  )
})

test_that("fLPPanel uses a separate sample at each horizon, as with pre-built leads", {
  d <- sim_panel_df()
  set.seed(11)
  d$shock[sample(nrow(d), 40)] <- NA   # rows that are only lead targets
  d$y[sample(nrow(d), 40)]     <- NA   # rows that only enter at h >= 1
  H <- 3
  p <- fLPPanel(y ~ shock | unit + tt, data = d, panel_id = c("unit", "tt"),
                shock = "shock", horizons = H, n_threads = NT)

  key <- paste(d$unit, d$tt)
  for (h in 0:H) {
    d$y_h <- d$y[match(paste(d$unit, d$tt + h), key)]
    ref <- lm(y_h ~ shock + factor(unit) + factor(tt), data = d)
    expect_equal(unname(p$irfs[h + 1, 1]), unname(coef(ref)["shock"]),
                 tolerance = 1e-8)
    expect_equal(unname(p$nobs[h + 1]), nobs(ref))
  }
})
