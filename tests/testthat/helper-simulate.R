# Shared helpers for the tidyMacro test suite (loaded automatically by testthat).

# A stable bivariate VAR(1) with intercept.
# fVAR() convention: beta is (N*p + c) x N, intercept in the first row.
true_beta_var1 <- function() {
  rbind(c(0.5, -0.2),          # intercept
        c(0.5,  0.1),          # lag-1 coefficients on y1
        c(0.2,  0.4))          # lag-1 coefficients on y2
}

# Simulate T x N data from a VAR(p) using the package's own generator.
sim_var <- function(T = 400L, N = 2L, p = 1L, c = 1L,
                    beta = true_beta_var1(), seed = 123L) {
  set.seed(seed)
  eps  <- matrix(rnorm((T - p) * N), T - p, N)
  y0   <- matrix(0, T, N)
  tidyMacro::fGenerateVARData(y0, p, c, beta, eps)
}

# A random symmetric positive-definite matrix.
sim_sigma <- function(N = 3L, seed = 1L) {
  set.seed(seed)
  A <- matrix(rnorm(N * N), N, N)
  crossprod(A) + diag(N)
}

# Always run C++ code single-threaded in tests (CRAN: max 2 cores).
NT <- 1L
