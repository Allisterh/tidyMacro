long_run_fixture <- function() {
  fit <- fVAR(sim_var(T = 300L, seed = 41L), p = 1L, c = 1L)
  list(
    wold = fWoldIRF(fit, horizon = 8L),
    S = t(chol(fit$sigma))
  )
}

test_that("fMaxIRF returns the constrained maximum-response direction", {
  x <- long_run_fixture()
  irf <- fMaxIRF(x$wold, x$S, var_idx = 2L)
  direction <- solve(x$S, irf[, 1])

  expect_equal(dim(irf), c(2L, 9L))
  expect_equal(direction[1], 0, tolerance = 1e-12)
  expect_equal(sqrt(sum(direction^2)), 1, tolerance = 1e-10)
  expect_gte(irf[2, 9], -1e-12)
})

test_that("fUhligMaxShare and fUhligIRF use the same normalized direction", {
  x <- long_run_fixture()
  h2 <- as.numeric(fUhligMaxShare(x$wold, x$S, idx = 2L))
  impact <- as.numeric(x$S %*% h2)
  if (sum(x$wold[1, , 9] * impact) < 0) impact <- -impact
  expected <- vapply(seq_len(9L), function(h) x$wold[, , h] %*% impact,
                     numeric(2L))
  irf <- fUhligIRF(x$wold, x$S, idx = 2L)

  expect_equal(h2[1], 0, tolerance = 1e-12)
  expect_equal(sqrt(sum(h2^2)), 1, tolerance = 1e-10)
  expect_equal(dim(irf), c(2L, 9L))
  expect_equal(irf, expected, tolerance = 1e-10, ignore_attr = TRUE)
  expect_gte(irf[1, 9], -1e-12)
})

test_that("long-run identification functions validate dimensions and indices", {
  x <- long_run_fixture()
  expect_error(fMaxIRF(x$wold, x$S, 0L), "out of range")
  expect_error(fUhligMaxShare(x$wold, diag(3), 1L), "conformable")
  expect_error(fUhligIRF(x$wold, x$S, 3L), "out of range")
})
