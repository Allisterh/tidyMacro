test_that("fGenerateQ returns an orthonormal matrix", {
  set.seed(1)
  Q <- fGenerateQ(4L)
  expect_equal(crossprod(Q), diag(4), tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(Q %*% t(Q), diag(4), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("fSignRestrictions_cpp finds an impact matrix that satisfies signs", {
  sigma <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  SIGN  <- matrix(c( 1, 1,     # shock 1: +, +
                    -1, 1), 2, 2)  # shock 2: -, +
  res <- fSignRestrictions_cpp(sigma, SIGN, sr_hor = 1L, sr_rot = 5000L,
                               seed = 7L)
  expect_true(res$found)
  B <- res$B
  expect_equal(B %*% t(B), sigma, tolerance = 1e-8, ignore_attr = TRUE)
  expect_true(all(B[, 1] >= 0))
  expect_true(B[1, 2] <= 0 && B[2, 2] >= 0)
})

test_that("fSignRestrictions_cpp is reproducible for a fixed seed", {
  sigma <- sim_sigma(3L)
  SIGN  <- matrix(c(1, 0, 0,
                    0, 1, 0,
                    0, 0, 1), 3, 3)
  a <- fSignRestrictions_cpp(sigma, SIGN, sr_rot = 2000L, seed = 99L)
  b <- fSignRestrictions_cpp(sigma, SIGN, sr_rot = 2000L, seed = 99L)
  expect_identical(a, b)
})

test_that("fVARPosterior_cpp draws have the right shape and are reproducible", {
  y <- sim_var(T = 150L)
  d1 <- fVARPosterior_cpp(y, p = 1L, c = 1L, ndraws = 25L, seed = 5L)
  d2 <- fVARPosterior_cpp(y, p = 1L, c = 1L, ndraws = 25L, seed = 5L)

  expect_equal(dim(d1$beta_draws),  c(3L, 2L, 25L))
  expect_equal(dim(d1$sigma_draws), c(2L, 2L, 25L))
  expect_identical(d1, d2)

  # every covariance draw must be positive definite
  pd <- vapply(seq_len(25), function(i)
    all(eigen(d1$sigma_draws[, , i], symmetric = TRUE,
              only.values = TRUE)$values > 0), logical(1))
  expect_true(all(pd))
})

test_that("fSR_cpp: sign-restricted draws respect the impact restrictions", {
  skip_on_cran()
  y    <- sim_var(T = 200L)
  SIGN <- matrix(c(1, 1,
                  -1, 1), 2, 2)
  res <- fSR_cpp(y, p = 1L, c = 1L, SIGN = SIGN, nsteps = 6L, ndraws = 20L,
                 sr_hor = 1L, store_draws = TRUE, n_threads = NT, seed = 42L)

  expect_equal(dim(res$IRmed), c(2L, 2L, 6L))
  expect_true(res$accept_rate > 0)

  # impact responses of every accepted draw satisfy SIGN
  IR <- array(res$IRall, c(2, 2, 6, ncol(res$IRall)))
  expect_true(all(IR[1, 1, 1, ] >= -1e-10))
  expect_true(all(IR[2, 1, 1, ] >= -1e-10))
  expect_true(all(IR[1, 2, 1, ] <=  1e-10))
  expect_true(all(IR[2, 2, 1, ] >= -1e-10))

  # FEVD shares live on [0, 1]
  expect_true(all(res$VDmed >= -1e-10 & res$VDmed <= 1 + 1e-10))
})

test_that("fSR_cpp does not depend on the number of threads", {
  skip_on_cran()
  skip_if(parallel::detectCores() < 2L)
  y    <- sim_var(T = 200L)
  SIGN <- matrix(c(1, 1, -1, 1), 2, 2)
  a <- fSR_cpp(y, 1L, 1L, SIGN, nsteps = 5L, ndraws = 15L,
               store_draws = FALSE, n_threads = 1L, seed = 3L)
  b <- fSR_cpp(y, 1L, 1L, SIGN, nsteps = 5L, ndraws = 15L,
               store_draws = FALSE, n_threads = 2L, seed = 3L)
  expect_equal(a$IRmed, b$IRmed, tolerance = 1e-12)
})

test_that("fHDShock_cpp components add up to the data", {
  skip_on_cran()
  y   <- sim_var(T = 120L)
  fit <- fVAR(y, p = 1L, c = 1L)
  B   <- t(chol(fit$sigma))
  hd  <- fHDShock_cpp(y, fit$beta, B, p = 1L, c = 1L)

  rows <- 2:nrow(y)   # first p rows are NA (initial conditions)
  shocks <- apply(hd$shock, c(1, 2), sum)
  # initial condition + intercept + sum of shock contributions = y
  recon <- (hd$init + hd$const + shocks)[rows, ]
  expect_equal(recon, y[rows, ], tolerance = 1e-6, ignore_attr = TRUE)
})
