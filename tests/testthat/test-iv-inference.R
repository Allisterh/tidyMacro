test_that("fRecoverBIV_cpp recovers a strong single-instrument impact column", {
  x <- proxy_fixture()
  r <- fRecoverBIV_cpp(x$fit$residuals, x$Z, x$fit$sigma,
                       ntotcoeff = nrow(x$fit$beta))

  expect_named(r, c("b1", "B", "sigma_b", "fs_beta", "fs_F", "fs_r2",
                    "shock_sd", "n_iv"))
  expect_equal(dim(r$B), c(2L, 2L))
  expect_equal(r$B %*% t(r$B), x$fit$sigma, tolerance = 1e-8,
               ignore_attr = TRUE)
  ratio <- r$B[, 1] / as.numeric(r$b1)
  expect_equal(ratio, rep(ratio[1], 2L),
               tolerance = 1e-8)
  expect_gt(r$fs_F, 10)
  expect_true(r$fs_r2 > 0 && r$fs_r2 < 1)
  expect_equal(r$n_iv, nrow(x$Z))
})

test_that("fRecoverBIVMulti_cpp returns coherent joint-IV diagnostics", {
  set.seed(32)
  u <- matrix(rnorm(600), 200, 3)
  u[, 3] <- 0.3 * u[, 1] - 0.2 * u[, 2] + rnorm(200, sd = 0.7)
  Z <- cbind(u[, 1] + rnorm(200, sd = 0.15),
             u[, 2] + rnorm(200, sd = 0.15))
  r <- fRecoverBIVMulti_cpp(u, Z, ntotcoeff = 4L)

  expect_named(r, c("B", "sigma_b", "fs_beta", "fs_F", "fs_r2",
                    "relEig", "n_iv"))
  expect_equal(dim(r$B), c(3L, 2L))
  expect_equal(dim(r$sigma_b), c(3L, 3L))
  expect_equal(dim(r$fs_beta), c(3L, 2L))
  expect_length(r$fs_F, 2L)
  expect_true(all(r$fs_F > 10))
  expect_true(all(r$fs_r2 > 0 & r$fs_r2 < 1))
  expect_true(all(is.finite(r$relEig)))
  expect_equal(r$n_iv, 200L)
})

test_that("fRecoverBIV functions reject misaligned inputs", {
  x <- proxy_fixture(T = 80L)
  expect_error(
    fRecoverBIV_cpp(x$fit$residuals, x$Z[-1, , drop = FALSE],
                    x$fit$sigma, 3L),
    "same number of rows"
  )
  expect_error(
    fRecoverBIVMulti_cpp(x$fit$residuals, matrix(numeric(), nrow(x$Z), 0), 3L),
    "at least one instrument"
  )
})

test_that("fMSW returns normalized weak-IV and delta-method inference", {
  x <- proxy_fixture()
  r <- fMSW(x$fit, x$Z, x$y, x$adjust, hor = 4L, nvar = 1L,
            scale = 2, confidence = 0.90, NWlags = 1L)

  expect_named(r, c("IRF", "IRFstderr", "Dmlbound", "Dmubound",
                    "MSWlbound", "MSWubound", "Waldstat", "Fstat"))
  for (nm in c("IRF", "IRFstderr", "Dmlbound", "Dmubound",
               "MSWlbound", "MSWubound")) {
    expect_equal(dim(r[[nm]]), c(2L, 5L))
  }
  expect_equal(r$IRF[1, 1], 2, tolerance = 1e-12)
  expect_equal(r$MSWlbound[1, 1], 2, tolerance = 1e-12)
  expect_equal(r$MSWubound[1, 1], 2, tolerance = 1e-12)
  expect_true(all(r$IRFstderr >= 0))
  expect_true(all(r$Dmlbound <= r$IRF + 1e-12))
  expect_true(all(r$IRF <= r$Dmubound + 1e-12))
  expect_gt(r$Fstat, 10)
  expect_equal(r$Fstat, r$Waldstat)
  expect_error(fMSW(x$fit, x$Z, x$y, x$adjust, NWlags = -1L),
               "NWlags must be >= 0")
})

test_that("fHeteroIRF normalizes impact and propagates it through the Wold MA", {
  x <- proxy_fixture()
  inds <- rep(1L, nrow(x$Z))
  r <- fHeteroIRF(x$fit, x$Z, x$adjust, inds, hor = 4L,
                  nvar = 1L, scale = 2)
  wold <- fWoldIRF(x$fit, horizon = 4L)
  expected <- vapply(seq_len(5L), function(h) wold[, , h] %*% r$b1,
                     numeric(2L))

  expect_equal(dim(r$IRF), c(2L, 5L))
  expect_equal(r$b1[1], 2, tolerance = 1e-12)
  expect_equal(r$IRF, expected, tolerance = 1e-10, ignore_attr = TRUE)
  expect_error(
    fHeteroIRF(x$fit, x$Z, x$adjust, integer(nrow(x$Z)), 2L, 1L),
    "no treatment observations"
  )
})

test_that("fHDIV follows the companion-form historical-decomposition recursion", {
  x <- proxy_fixture(T = 120L)
  s <- c(1, 0.35)
  hd <- fHDIV(x$fit$residuals, x$fit$sigma, s, x$fit$beta, c = 1L, p = 1L)
  shock <- fGetShock(x$fit$residuals, x$fit$sigma, s,
                     shockSize = 1, normalize = "unit")
  A <- t(x$fit$beta[-1, , drop = FALSE])
  ref <- matrix(0, nrow(x$fit$residuals), 2L)
  ref[1, ] <- s * shock[1]
  for (i in 2:nrow(ref)) ref[i, ] <- A %*% ref[i - 1, ] + s * shock[i]

  expect_named(hd, "HDshock")
  expect_equal(dim(hd$HDshock), dim(x$fit$residuals))
  expect_equal(hd$HDshock, ref, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("fBootstrapIVMBB returns nested bands and preserves impact normalization", {
  x <- proxy_fixture(T = 120L)
  set.seed(33)
  b <- suppressMessages(fBootstrapIVMBB(
    x$y, x$fit, x$Z, nboot = 15L, blocksize = 5L,
    adjustZ = x$adjust, adjustu = x$adjust, policyvar = 1L,
    horizon = 3L, conf = 90, conf2 = 68, n_threads = NT
  ))

  expect_nested_bands(b, 2L, 4L)
  expect_equal(dim(b$meanirf), c(2L, 4L))
  expect_equal(dim(b$medianirf), c(2L, 4L))
  expect_equal(b$meanirf[1, 1], 1, tolerance = 1e-12)
  expect_equal(b$medianirf[1, 1], 1, tolerance = 1e-12)
})

test_that("fBootstrapIVInvertible returns coherent wild-bootstrap bands", {
  x <- proxy_fixture(T = 120L)
  set.seed(34)
  b <- fBootstrapIVInvertible(
    x$y, x$instr, x$fit, nboot = 15L, p = 1L, c = 1L,
    hor = 3L, cumu = integer(), conf = 90, conf2 = 68
  )

  expect_nested_bands(b, 2L, 4L)
  expect_equal(dim(b$median), c(2L, 4L))
  expect_true(all(b$lower <= b$median + 1e-12))
  expect_true(all(b$median <= b$upper + 1e-12))
})

test_that("fBootstrapIVRecover returns coherent recoverable-IV bands", {
  x <- proxy_fixture(T = 120L)
  r <- 1L
  eps <- x$fit$residuals
  n_eff <- nrow(eps) - r
  delta <- c(0, 1, 0.2, 0.45, 0.1)
  set.seed(35)
  noise <- rnorm(n_eff, sd = 0.15)
  eps_f <- cbind(eps[seq_len(n_eff), , drop = FALSE],
                 eps[2:(n_eff + 1L), , drop = FALSE])
  instr <- numeric(nrow(x$y))
  instr[2:(nrow(x$y) - r)] <- cbind(1, eps_f) %*% delta + noise

  set.seed(36)
  b <- fBootstrapIVRecover(
    x$y, instr, x$fit, noise, delta, nboot = 15L,
    p = 1L, c = 1L, r = r, hor = 3L, cumu = integer(),
    conf = 90, conf2 = 68
  )

  expect_nested_bands(b, 2L, 4L)
  expect_equal(dim(b$median), c(2L, 4L))
  expect_true(all(b$lower <= b$median + 1e-12))
  expect_true(all(b$median <= b$upper + 1e-12))
})
