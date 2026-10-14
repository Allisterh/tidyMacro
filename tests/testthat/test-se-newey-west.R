# Independent R reference for the Newey-West meat and its Andrews-Monahan
# (1992) VAR(1)-prewhitened version (singular values capped at 0.97).
nw_meat_ref <- function(g, L) {
  S <- crossprod(g); n <- nrow(g)
  for (l in seq_len(min(L, n - 1))) {
    Gl <- t(g[(l + 1):n, , drop = FALSE]) %*% g[1:(n - l), , drop = FALSE]
    S <- S + (1 - l / (L + 1)) * (Gl + t(Gl))
  }
  S
}
pw_meat_ref <- function(g, L, cap = 0.97) {
  n <- nrow(g); k <- ncol(g)
  G0 <- g[-1, , drop = FALSE]; G1 <- g[-n, , drop = FALSE]
  A <- t(qr.solve(G1, G0))
  if (max(Mod(eigen(A, only.values = TRUE)$values)) > cap) {
    s <- svd(A); A <- s$u %*% diag(pmin(s$d, cap), k) %*% t(s$v)
  }
  D <- solve(diag(k) - A)
  D %*% nw_meat_ref(G0 - G1 %*% t(A), L) %*% t(D)
}
ols_pw_se <- function(X, y, L) {
  XXi <- solve(crossprod(X))
  u <- as.numeric(y - X %*% (XXi %*% crossprod(X, y)))
  sqrt(diag(XXi %*% pw_meat_ref(X * u, L) %*% XXi))
}

test_that("fLP prewhitened SEs match the reference on fast and full paths", {
  df <- sim_lp_df(n = 300)
  n <- nrow(df)
  fast <- fLP(y ~ x, data = df, horizons = 3, shock = "x", nw_lags = 0,
              nw_prewhite = TRUE, n_threads = NT)
  full <- fLP(y ~ x, data = df, horizons = 3, shock = "x", nw_lags = 0,
              nw_prewhite = TRUE, store_full = TRUE, n_threads = NT)
  expect_equal(fast$irfs_se, full$irfs_se, tolerance = 1e-12)
  for (h in 0:3) {
    se <- ols_pw_se(cbind(1, df$x[1:(n - h)]), df$y[(h + 1):n], h + 1)
    expect_equal(unname(full$ses[[h + 1]][, "y"]), se, tolerance = 1e-10)
  }
  plain <- fLP(y ~ x, data = df, horizons = 3, shock = "x", nw_lags = 0,
               n_threads = NT)
  expect_false(isTRUE(all.equal(plain$irfs_se, fast$irfs_se)))
  expect_error(fLP(y ~ x, data = df, horizons = 1, shock = "x",
                   nw_prewhite = NA), "nw_prewhite")
})

test_that("prewhitening caps an explosive VAR(1) for the scores", {
  set.seed(3)
  n <- 300
  x <- cumsum(rnorm(n)) / 5
  y <- 0.5 * x + cumsum(rnorm(n)) / 5
  f <- fLP(y ~ x, data = data.frame(y = y, x = x), horizons = 0, shock = "x",
           nw_lags = 2, nw_prewhite = TRUE, n_threads = NT)
  X <- cbind(1, x)
  u <- as.numeric(y - X %*% solve(crossprod(X), crossprod(X, y)))
  g <- X * u
  A <- t(qr.solve(g[-n, ], g[-1, ]))
  expect_gt(max(Mod(eigen(A)$values)), 0.97)
  expect_equal(f$irfs_se[1, 1], unname(ols_pw_se(X, y, 3)[2]), tolerance = 1e-10)
})

test_that("fLPIV prewhitens the full 2SLS score vector", {
  set.seed(9)
  n <- 300
  z <- rnorm(n); w <- rnorm(n); D <- 0.7 * z + 0.2 * w + rnorm(n)
  y <- 1.2 * D + 0.3 * w + rnorm(n) + 0.5 * c(0, D[-n])
  iv <- fLPIV(y ~ D + w, instruments = ~ z, data = data.frame(y, D, z, w),
              endog = "D", horizons = 0, nw_lags_iv = 3, nw_prewhite = TRUE,
              cumulative = FALSE, n_threads = NT)
  X <- cbind(1, D, w); Zm <- cbind(1, z, w)
  Xh <- Zm %*% solve(crossprod(Zm), crossprod(Zm, X))
  e <- as.numeric(y - X %*% solve(crossprod(Xh, X), crossprod(Xh, y)))
  A <- solve(crossprod(Xh))
  V <- A %*% pw_meat_ref(Xh * e, 3) %*% A
  expect_equal(iv$irfs_se[1, 1], sqrt(V[2, 2]), tolerance = 1e-10)
})

test_that("fNLLP prewhitened SEs agree across paths and with the reference", {
  df <- sim_lp_df(n = 300)
  set.seed(4)
  df$m <- plogis(rnorm(nrow(df)))
  args <- list(y ~ x, data = df, horizons = 2, shock = "x", type = "state",
               state = "m", nw_lags = 0, nw_prewhite = TRUE, n_threads = NT)
  a <- do.call(fNLLP, args)
  b <- do.call(fNLLP, c(args, store_full = TRUE))
  expect_equal(a$irfs_se, b$irfs_se, tolerance = 1e-12)
  expect_equal(a$irfs_cov, b$irfs_cov, tolerance = 1e-12)

  n <- nrow(df) - 2
  m <- df$m[1:n]; x <- df$x[1:n]
  se <- ols_pw_se(cbind(m, m * x, 1 - m, (1 - m) * x), df$y[3:nrow(df)], 3)
  expect_equal(a$irfs_se[3, 1, ], se[c(2, 4)], tolerance = 1e-10,
               ignore_attr = TRUE)
})
