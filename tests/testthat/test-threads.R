# Every routine that uses OpenMP must return the same numbers whatever the
# thread count: all random draws are made on the main thread before the
# parallel region. These run with two threads, the CRAN limit.

var_fixture <- function(T = 140L) {
  y <- sim_var(T = T)
  list(y = y, fit = fVAR(y, p = 1L, c = 1L))
}

# Numeric content only: results also record the thread count and the call,
# which legitimately differ between the two runs.
numbers_of <- function(x) {
  if (!is.list(x)) return(if (is.numeric(x)) x else NULL)
  x <- x[setdiff(names(x), c("n_threads", "call", "formula", "formula_orig",
                             "verbose", "cluster"))]
  out <- lapply(x, numbers_of)
  out[!vapply(out, is.null, logical(1))]
}

# Run `f(nt)` from the same seed and compare the two thread counts.
same_across_threads <- function(f, seed = 5L) {
  run <- function(nt) {
    set.seed(seed)
    numbers_of(suppressMessages(f(nt)))
  }
  expect_identical(run(1L), run(2L))
}

test_that("identified-VAR bootstraps are thread-count invariant", {
  skip_on_cran()
  x <- var_fixture()

  same_across_threads(function(nt)
    fBootstrapChol(x$y, x$fit, nboot = 25, horizon = 3, n_threads = nt))
  same_across_threads(function(nt)
    fBootstrapCholCorrected(x$y, x$fit, nboot1 = 15, nboot2 = 25, horizon = 3,
                            n_threads = nt))
  same_across_threads(function(nt)
    fBootstrapBQ(x$y, x$fit, nboot = 25, horizon = 3, n_threads = nt))
  same_across_threads(function(nt)
    fBootstrapMax(x$y, x$fit, nboot = 25, horizon = 3, var_idx = 1,
                  n_threads = nt))
  same_across_threads(function(nt)
    fBootstrapMaxCorrected(x$y, x$fit, nboot1 = 15, nboot2 = 25, horizon = 3,
                           var_idx = 1, n_threads = nt))
  same_across_threads(function(nt)
    fBootstrapUhlig(x$y, x$fit, nboot = 25, horizon = 3, idx = 1,
                    n_threads = nt))
  same_across_threads(function(nt)
    fBootstrapUhligCorrected(x$y, x$fit, nboot1 = 15, nboot2 = 25, horizon = 3,
                             idx = 1, n_threads = nt))
})

test_that("instrument-based bootstraps are thread-count invariant", {
  skip_on_cran()
  x <- proxy_fixture(T = 140L)

  same_across_threads(function(nt)
    fBootstrapIVMBB(x$y, x$fit, x$Z, nboot = 25L, blocksize = 5L,
                    adjustZ = x$adjust, adjustu = x$adjust, policyvar = 1L,
                    horizon = 3L, n_threads = nt))
  same_across_threads(function(nt)
    fBootstrapIVInvertible(x$y, x$instr, x$fit, nboot = 25L, p = 1L, c = 1L,
                           hor = 3L, cumu = integer(), n_threads = nt))
  same_across_threads(function(nt)
    fBootstrapHetero(x$y, x$fit, x$Z,
                     indsR1 = as.integer(abs(x$Z[, 1]) > stats::median(abs(x$Z[, 1]))),
                     adjustu = x$adjust, nboot = 60L, blocksize = 5L,
                     hor = 3L, nvar = 1L, n_threads = nt))
  s_hat <- c(1, rep(0.5, ncol(x$y) - 1))
  same_across_threads(function(nt)
    fBootstrapHDIV(x$y, x$fit, x$Z, s = s_hat, nboot = 25L, blocksize = 5L,
                   adjustZ = x$adjust, adjustu = x$adjust, policyvar = 1L,
                   n_threads = nt))
})

test_that("local projections and sign restrictions are thread-count invariant", {
  skip_on_cran()
  df <- sim_lp_df()
  same_across_threads(function(nt)
    fLP(y ~ x + l(y, 1:2), data = df, horizons = 3L, shock = "x",
        n_threads = nt))
  same_across_threads(function(nt)
    fLPIV(y ~ D, instruments = ~ z, data = sim_iv_df(), endog = "D",
          horizons = 3L, n_threads = nt))
  same_across_threads(function(nt)
    fLPPanel(y ~ shock | unit + tt, data = sim_panel_df(),
             panel_id = c("unit", "tt"), shock = "shock", horizons = 2L,
             ssc = "CR1", n_threads = nt))
  set.seed(4)
  n <- 300
  nl <- data.frame(s = rnorm(n), x = rnorm(n), m = plogis(2 * rnorm(n)))
  nl$y <- (0.9 * nl$m + 0.1 * (1 - nl$m)) * nl$s + 0.5 * nl$x + rnorm(n, sd = 0.3)
  same_across_threads(function(nt)
    fNLLP(y ~ s + x, data = nl, horizons = 2, shock = "s",
          type = "state", state = "m", n_threads = nt))
})

test_that("a failed draw is reported rather than aborting the session", {
  x <- proxy_fixture(T = 120L)

  # An instrument that is orthogonal to the residuals by construction leaves
  # the second stage singular on the original sample.
  dead <- matrix(0, nrow = nrow(x$Z), ncol = 1L)
  expect_error(
    suppressMessages(fBootstrapIVMBB(
      x$y, x$fit, dead, nboot = 5L, blocksize = 5L,
      adjustZ = x$adjust, adjustu = x$adjust, policyvar = 1L,
      horizon = 2L, n_threads = 1L)),
    "singular|failed")
})
