# Shared strong-proxy fixture for IV inference and bootstrap tests.
proxy_fixture <- function(T = 160L, seed = 31L) {
  y <- sim_var(T = T, seed = seed)
  fit <- fVAR(y, p = 1L, c = 1L)
  set.seed(seed + 1L)
  z <- fit$residuals[, 1] + rnorm(nrow(fit$residuals), sd = 0.15)
  list(
    y = y,
    fit = fit,
    Z = matrix(z, ncol = 1L),
    instr = c(0, z),
    adjust = c(1L, nrow(fit$residuals))
  )
}

expect_nested_bands <- function(x, N, H) {
  for (nm in c("upper", "lower", "upper2", "lower2")) {
    expect_equal(dim(x[[nm]]), c(N, H))
    expect_true(all(is.finite(x[[nm]])))
  }
  expect_true(all(x$lower <= x$lower2 + 1e-12))
  expect_true(all(x$lower2 <= x$upper2 + 1e-12))
  expect_true(all(x$upper2 <= x$upper + 1e-12))
}
