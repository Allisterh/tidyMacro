did_fit <- function(d = sim_did_df(), post = 4L, pre = 3L, ...) {
  fLPDID(y ~ treat, data = d, panel_id = c("id", "year"), treat = "treat",
         post = post, pre = pre, n_threads = 1L, ...)
}

test_that("fLPDID returns the documented tibble", {
  m <- did_fit()
  expect_s3_class(m, "fLPDID")
  expect_s3_class(m, "tbl_df")
  expect_named(m, c("event_time", "estimate", "se", "conf_low", "conf_high",
                    "nobs", "nclust", "ndrop"))
  expect_equal(m$event_time, -3:4)
})

test_that("fLPDID normalises event time -1 to zero", {
  m <- did_fit()
  r <- m[m$event_time == -1, ]
  expect_equal(r$estimate, 0)
  expect_equal(r$se, 0)
})

test_that("fLPDID recovers a constant treatment effect and flat pre-trends", {
  m <- did_fit()
  post <- m[m$event_time >= 0, ]
  pre  <- m[m$event_time < -1, ]
  expect_true(all(abs(post$estimate - 1) < 0.25))
  expect_true(all(abs(pre$estimate) < 0.3))
  expect_true(all(m$nclust[m$event_time != -1] == 150L))
})

test_that("fLPDID finds no effect when there is none (placebo)", {
  m <- did_fit(sim_did_df(tau = 0))
  m <- m[m$event_time != -1, ]
  expect_true(all(abs(m$estimate) < 3.5 * m$se))
})

test_that("fLPDID bands bracket the estimate and widen with conf", {
  m90 <- did_fit()
  m95 <- did_fit(conf = 95)
  k <- m90$event_time != -1
  expect_true(all(m90$conf_low[k] < m90$estimate[k] &
                  m90$estimate[k] < m90$conf_high[k]))
  expect_true(all((m95$conf_high - m95$conf_low)[k] >
                  (m90$conf_high - m90$conf_low)[k]))
})

test_that("fLPDID is deterministic and invariant to the number of threads", {
  d <- sim_did_df()
  a <- did_fit(d)
  b <- did_fit(d)
  expect_identical(a$estimate, b$estimate)
  skip_on_cran()
  c2 <- fLPDID(y ~ treat, data = d, panel_id = c("id", "year"),
               treat = "treat", post = 4L, pre = 3L, n_threads = 2L)
  expect_equal(a$estimate, c2$estimate, tolerance = 1e-12)
  expect_equal(a$se, c2$se, tolerance = 1e-12)
})

test_that("fLPDID options pmd and reweight run and change the baseline", {
  d  <- sim_did_df()
  b  <- did_fit(d, post = 3L, pre = 2L)
  pm <- did_fit(d, post = 3L, pre = 2L, pmd = TRUE)
  rw <- did_fit(d, post = 3L, pre = 2L, reweight = TRUE)
  expect_true(all(is.finite(pm$estimate)))
  expect_true(all(is.finite(rw$estimate)))
  expect_false(isTRUE(all.equal(b$estimate, pm$estimate)))
})

test_that("fLPDID accepts panel-aware lag controls", {
  d <- sim_did_df()
  m <- fLPDID(y ~ treat + l(y, 1:2), data = d, panel_id = c("id", "year"),
              treat = "treat", post = 3L, pre = 2L, n_threads = 1L)
  expect_equal(m$event_time, -2:3)
  expect_true(all(is.finite(m$estimate)))
})

test_that("fLPDID supports non-absorbing treatment with a stabilisation window", {
  d <- sim_did_df()
  set.seed(9)
  d$treat <- rbinom(nrow(d), 1, 0.2)
  d$y <- rnorm(150)[d$id] + d$treat + rnorm(nrow(d), sd = 0.5)
  m <- fLPDID(y ~ treat, data = d, panel_id = c("id", "year"),
              treat = "treat", post = 3L, pre = 2L, nonabsorbing = TRUE,
              L = 2L, n_threads = 1L)
  expect_equal(m$event_time, -2:3)
  expect_true(all(is.finite(m$estimate)))
  expect_equal(m$estimate[m$event_time == 0], 1, tolerance = 0.3)
})

test_that("fLPDID validates its arguments", {
  d <- sim_did_df()
  expect_error(fLPDID(y ~ treat, data = d, panel_id = c("id", "year"),
                      post = 3L, n_threads = 1L), "'treat' is required")
  expect_error(fLPDID(y ~ treat, data = d, treat = "treat", post = 3L,
                      n_threads = 1L), "'panel_id' is required")
  expect_error(fLPDID(y ~ treat, data = d, panel_id = c("id", "yr"),
                      treat = "treat", post = 3L, n_threads = 1L),
               "column 'yr' not found")
  expect_error(fLPDID(y ~ treat, data = d, panel_id = c("id", "year"),
                      treat = "treat", post = 3L, nonabsorbing = TRUE,
                      n_threads = 1L), "requires the stabilization window")
  d$id[1] <- NA
  expect_error(fLPDID(y ~ treat, data = d, panel_id = c("id", "year"),
                      treat = "treat", n_threads = 1L),
               "missing value\\(s\\) in unit id")
})

test_that("fPlotLPDID returns a ggplot for one or several fits", {
  d <- sim_did_df()
  a <- did_fit(d)
  b <- did_fit(d, conf = 95)
  expect_s3_class(fPlotLPDID(a), "ggplot")
  expect_s3_class(fPlotLPDID(base = a, wide = b), "ggplot")
})
