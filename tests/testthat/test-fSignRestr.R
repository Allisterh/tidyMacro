SIGN2 <- matrix(c(1, 1, -1, 1), 2, 2)   # shock 1: (+,+); shock 2: (-,+)

sr_fit <- function(y = sim_sr_data(), nsteps = 6L, ndraws = 15L, seed = 1L,
                   ...) {
  fSignRestr(y, p = 1, sign = SIGN2, nsteps = nsteps, ndraws = ndraws,
             n_threads = 1, seed = seed, ...)
}

test_that("fSignRestr returns an fSignRestr object with reshaped draws", {
  r <- sr_fit()
  expect_s3_class(r, "fSignRestr")
  expect_equal(dim(r$IRmed), c(2L, 2L, 6L))
  expect_equal(dim(r$IRall), c(2L, 2L, 6L, 15L))
  expect_equal(dim(r$VDall), c(2L, 2L, 6L, 15L))
  expect_equal(r$varnames, c("a", "b"))
  expect_equal(r$ident, "sign")
  expect_equal(r$p, 1L)
  expect_equal(r$nsteps, 6L)
  expect_true(r$accept_rate > 0 && r$accept_rate <= 1)
})

test_that("every accepted draw satisfies the impact sign restrictions", {
  r  <- sr_fit(ndraws = 30L)
  IR <- r$IRall
  expect_true(all(IR[1, 1, 1, ] >= -1e-10))
  expect_true(all(IR[2, 1, 1, ] >= -1e-10))
  expect_true(all(IR[1, 2, 1, ] <=  1e-10))
  expect_true(all(IR[2, 2, 1, ] >= -1e-10))
})

test_that("credible bands are ordered and FEVD shares are proper", {
  r <- sr_fit(ndraws = 30L)
  expect_true(all(r$IRinf <= r$IRmed + 1e-10))
  expect_true(all(r$IRmed <= r$IRsup + 1e-10))
  expect_true(all(r$VDall >= -1e-10 & r$VDall <= 1 + 1e-10))
  shares <- apply(r$VDall, c(1, 3, 4), sum)   # sum over shocks
  expect_equal(as.numeric(shares), rep(1, length(shares)), tolerance = 1e-8)
})

test_that("results are reproducible for a seed and differ across seeds", {
  expect_identical(sr_fit(seed = 1L)$IRmed, sr_fit(seed = 1L)$IRmed)
  expect_false(identical(sr_fit(seed = 1L)$IRmed, sr_fit(seed = 2L)$IRmed))
})

test_that("store_draws = FALSE drops the draws but not the summaries", {
  full <- sr_fit()
  lean <- sr_fit(store_draws = FALSE)
  expect_null(lean$IRall)
  expect_identical(full$IRmed, lean$IRmed)
})

test_that("several confidence levels are returned in bands", {
  r <- sr_fit(conf = c(68, 90))
  expect_equal(names(r$bands), c("68", "90"))
})

test_that("fSignRestr accepts data frames and custom variable names", {
  y <- sim_sr_data()
  r <- sr_fit(as.data.frame(y), varnames = c("X", "Y"))
  expect_equal(r$varnames, c("X", "Y"))
  expect_equal(dim(r$IRmed), c(2L, 2L, 6L))
})

test_that("narrative sign restrictions are accepted", {
  r <- sr_fit(ndraws = 10L,
              narrative = list(sign = list(shock = 1, period = 100, sign = 1)))
  expect_true(r$accept_rate > 0 && r$accept_rate <= 1)
  expect_equal(dim(r$IRmed), c(2L, 2L, 6L))
})

test_that("an external instrument switches to the sign+iv scheme", {
  y <- sim_sr_data()
  set.seed(6)
  Z <- matrix(rnorm(nrow(y)), ncol = 1)
  r <- fSignRestr(y, p = 1, sign = matrix(1, 2, 1), nsteps = 4L, ndraws = 10L,
                  instrument = list(Z = Z), n_threads = 1, seed = 1)
  expect_equal(r$ident, "sign+iv")
  expect_equal(dim(r$IRmed), c(2L, 2L, 4L))
})

test_that("print.fSignRestr produces output", {
  expect_output(print(sr_fit()))
})

test_that("fSignRestr validates its inputs", {
  y <- sim_sr_data()
  expect_error(fSignRestr(y, p = 1, sign = matrix(1, 3, 3), nsteps = 4L,
                          ndraws = 5L, n_threads = 1),
               "must have 2 rows")
  yb <- y; yb[3, 1] <- NA
  expect_error(fSignRestr(yb, p = 1, sign = SIGN2, n_threads = 1),
               "finite values")
  expect_error(fSignRestr(y, p = 0.5, sign = SIGN2, n_threads = 1),
               "finite integer")
})

test_that("fSignRestr is invariant to the number of threads", {
  skip_on_cran()
  y <- sim_sr_data()
  a <- fSignRestr(y, p = 1, sign = SIGN2, nsteps = 5L, ndraws = 15L,
                  n_threads = 1, seed = 3)
  b <- fSignRestr(y, p = 1, sign = SIGN2, nsteps = 5L, ndraws = 15L,
                  n_threads = 2, seed = 3)
  expect_equal(a$IRmed, b$IRmed, tolerance = 1e-12)
})
