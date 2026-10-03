test_that("fCompanionMatrix has the textbook structure", {
  beta <- rbind(c(0.5, 0.1),     # lag 1 (N = 2, p = 2, c = 0)
                c(0.2, 0.4),
                c(0.1, 0.0),     # lag 2
                c(0.0, 0.05))
  cm <- fCompanionMatrix(beta, c = 0L, p = 2L)
  expect_equal(cm$N, 2L)
  comp <- cm$comp
  expect_equal(dim(comp), c(4L, 4L))
  expect_equal(unname(comp[1:2, 1:4]), unname(t(beta)))
  expect_equal(unname(comp[3:4, 1:2]), diag(2))
  expect_equal(unname(comp[3:4, 3:4]), matrix(0, 2, 2))
})

test_that("fWoldIRF: impact is identity and later horizons are powers of companion", {
  y    <- sim_var(T = 300L)
  fit  <- fVAR(y, p = 1L, c = 1L)
  wold <- fWoldIRF(fit, horizon = 5L)
  comp <- fCompanionMatrix(fit$beta, fit$c, fit$p)$comp

  expect_equal(dim(wold), c(2L, 2L, 6L))
  expect_equal(wold[, , 1], diag(2), ignore_attr = TRUE)
  expect_equal(wold[, , 2], comp[1:2, 1:2], ignore_attr = TRUE)
  expect_equal(wold[, , 3], (comp %*% comp)[1:2, 1:2], ignore_attr = TRUE)
})

test_that("stable VAR has Wold responses that decay", {
  y    <- sim_var(T = 1000L)
  fit  <- fVAR(y, p = 1L, c = 1L)
  wold <- fWoldIRF(fit, horizon = 40L)
  expect_lt(max(abs(wold[, , 41])), 1e-2)
})

test_that("fCheckStability runs without error", {
  fit <- fVAR(sim_var(T = 200L), p = 1L, c = 1L)
  expect_no_error(capture.output(fCheckStability(fit)))
})

test_that("fCholeskyIRF: impact response equals the Cholesky factor", {
  fit  <- fVAR(sim_var(T = 400L), p = 1L, c = 1L)
  S    <- t(chol(fit$sigma))
  wold <- fWoldIRF(fit, horizon = 8L)
  irf  <- fCholeskyIRF(wold, S)

  expect_equal(dim(irf), dim(wold))
  expect_equal(irf[, , 1], S, ignore_attr = TRUE)
  expect_equal(irf[, , 4], wold[, , 4] %*% S, ignore_attr = TRUE)
})

test_that("fFEVDChol shares are in [0, 1] and sum to one across shocks", {
  fit  <- fVAR(sim_var(T = 400L), p = 1L, c = 1L)
  S    <- t(chol(fit$sigma))
  irf  <- fCholeskyIRF(fWoldIRF(fit, horizon = 10L), S)
  fevd <- fFEVDChol(irf)$fevd

  expect_true(all(fevd >= -1e-12 & fevd <= 1 + 1e-12))
  totals <- apply(fevd, c(1, 3), sum)
  expect_equal(unname(totals), matrix(1, nrow(totals), ncol(totals)),
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("fBQIRF: long-run cumulative response is lower triangular", {
  fit  <- fVAR(sim_var(T = 600L), p = 1L, c = 1L)
  wold <- fWoldIRF(fit, horizon = 200L)
  C1   <- apply(wold, c(1, 2), sum)
  D1   <- t(chol(C1 %*% fit$sigma %*% t(C1)))
  K    <- solve(C1, D1)
  bq   <- fBQIRF(wold, K)

  lr <- apply(bq, c(1, 2), sum)
  expect_equal(lr[1, 2], 0, tolerance = 1e-8)   # BQ restriction
})

test_that("fPolyConvolve with a unit polynomial returns the input", {
  A  <- array(rnorm(2 * 2 * 4), c(2, 2, 4))
  I  <- array(0, c(2, 2, 1)); I[, , 1] <- diag(2)
  out <- fPolyConvolve(A, I, nlags = 4L)
  expect_equal(out, A, tolerance = 1e-12)
})

test_that("fGetBands orders lower <= median <= upper", {
  set.seed(8)
  boot  <- array(rnorm(2 * 5 * 200), c(2, 5, 200))
  bands <- fGetBands(boot, conf = 68)
  expect_true(all(bands$lower <= bands$median + 1e-12))
  expect_true(all(bands$median <= bands$upper + 1e-12))
})

test_that("fBootstrapChol returns coherent bands (small, single thread)", {
  skip_on_cran()
  y   <- sim_var(T = 200L)
  fit <- fVAR(y, p = 1L, c = 1L)
  set.seed(11)
  bs <- fBootstrapChol(y, fit, nboot = 30L, horizon = 6L,
                       conf = 90, conf2 = 68, n_threads = NT)
  expect_equal(dim(bs$upper), c(2L, 2L, 7L))
  expect_true(all(bs$lower <= bs$upper + 1e-12))
  expect_true(all(bs$upper2 <= bs$upper + 1e-12))
})
