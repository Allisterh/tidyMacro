test_that("fLagMakerMatrix returns (T-p) x (N*p) with lag-ordered columns", {
  y <- cbind(a = 1:10, b = 11:20)
  X <- fLagMakerMatrix(y, 2L)
  expect_equal(dim(X), c(8L, 4L))
  # columns: a_lag1, b_lag1, a_lag2, b_lag2; first row corresponds to t = 3
  expect_equal(as.numeric(X[1, ]), c(2, 12, 1, 11))
  expect_equal(as.numeric(X[8, ]), c(9, 19, 8, 18))
})

test_that("fVAR returns the documented structure and dimensions", {
  y   <- sim_var(T = 200L)
  fit <- fVAR(y, p = 1L, c = 1L)

  expect_named(fit, c("beta", "residuals", "sigma", "p", "c", "n_exog"),
               ignore.order = TRUE)
  expect_equal(dim(fit$beta), c(3L, 2L))            # (N*p + c) x N
  expect_equal(dim(fit$residuals), c(199L, 2L))     # (T - p) x N
  expect_equal(dim(fit$sigma), c(2L, 2L))
  expect_equal(fit$p, 1L)
  expect_equal(fit$c, 1L)
  expect_equal(fit$n_exog, 0L)
  expect_true(isSymmetric(unname(fit$sigma), tol = 1e-10))
})

test_that("fVAR recovers the true coefficients in a large sample", {
  y   <- sim_var(T = 5000L)
  fit <- fVAR(y, p = 1L, c = 1L)
  expect_equal(unname(fit$beta), true_beta_var1(), tolerance = 0.05)
})

test_that("fVAR matches equation-by-equation OLS via lm()", {
  y   <- sim_var(T = 300L)
  fit <- fVAR(y, p = 1L, c = 1L)

  Y  <- y[-1, ]
  X  <- y[-nrow(y), ]
  for (j in 1:2) {
    ref <- unname(coef(lm(Y[, j] ~ X)))   # intercept, lag y1, lag y2
    expect_equal(as.numeric(fit$beta[, j]), ref, tolerance = 1e-8)
  }
})

test_that("fVAR sigma uses the documented degrees-of-freedom correction", {
  y   <- sim_var(T = 250L)
  N <- 2L; p <- 1L; c <- 1L
  fit <- fVAR(y, p = p, c = c)
  dof <- nrow(fit$residuals) - c - N * p
  expect_equal(unname(fit$sigma), crossprod(fit$residuals) / dof,
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("fVAR works without an intercept and with exogenous regressors", {
  y <- sim_var(T = 200L, c = 0L, beta = true_beta_var1()[-1, ])
  fit0 <- fVAR(y, p = 1L, c = 0L)
  expect_equal(dim(fit0$beta), c(2L, 2L))

  set.seed(9)
  ex  <- matrix(rnorm(nrow(y)), ncol = 1)
  fitx <- fVAR(y, p = 1L, c = 1L, exog = ex)
  expect_equal(nrow(fitx$beta), 2L * 1L + 1L + 1L)
  expect_equal(fitx$n_exog, 1L)
})

test_that("fGenerateVARData keeps the first p rows and is deterministic", {
  set.seed(1)
  y0  <- matrix(rnorm(40), 20, 2)
  eps <- matrix(rnorm(38), 19, 2)
  out1 <- fGenerateVARData(y0, 1L, 1L, true_beta_var1(), eps)
  out2 <- fGenerateVARData(y0, 1L, 1L, true_beta_var1(), eps)
  expect_equal(dim(out1), dim(y0))
  expect_equal(out1[1, ], y0[1, ])
  expect_identical(out1, out2)
})

test_that("fAICBIC selects a sensible lag order", {
  y   <- sim_var(T = 1500L)          # true DGP is VAR(1)
  sel <- fAICBIC(y, pmax = 6L, c = 1L)
  expect_named(sel, c("aic", "bic", "hq"), ignore.order = TRUE)
  expect_true(all(unlist(sel) >= 1 & unlist(sel) <= 6))
  expect_equal(sel$bic, 1L)
})

test_that("fOLS agrees with lm()", {
  set.seed(5)
  X <- matrix(rnorm(300), 100, 3)
  y <- matrix(1 + X %*% c(0.5, -1, 2) + rnorm(100), ncol = 1)
  res <- fOLS(y, X, c = 1L)
  ref <- lm(y ~ X)

  expect_equal(sort(as.numeric(res$beta)), sort(unname(coef(ref))),
               tolerance = 1e-8)
  expect_equal(as.numeric(res$fitted), unname(fitted(ref)), tolerance = 1e-8)
  expect_equal(as.numeric(res$r2), summary(ref)$r.squared, tolerance = 1e-8)
  expect_equal(as.numeric(res$r2adj), summary(ref)$adj.r.squared,
               tolerance = 1e-8)
})

test_that("fMBBVAR returns residuals of the right shape", {
  set.seed(3)
  eps <- matrix(rnorm(200), 100, 2)
  out <- fMBBVAR(eps, lags = 1L, BlockSize = 5L)
  expect_equal(dim(out$eps_boot), dim(eps))
})
