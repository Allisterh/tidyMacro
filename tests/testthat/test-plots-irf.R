plot_irf_fixture <- function() {
  cube <- array(seq_len(2 * 2 * 4) / 10, c(2, 2, 4))
  mat <- cube[, 1, ]
  list(
    cube = cube,
    mat = mat,
    cube_bands = list(
      upper = cube + 0.20, lower = cube - 0.20,
      upper2 = cube + 0.10, lower2 = cube - 0.10
    ),
    mat_bands = list(
      upper = mat + 0.20, lower = mat - 0.20,
      upper2 = mat + 0.10, lower2 = mat - 0.10
    )
  )
}

mock_sign_fit <- function(offset = 0) {
  med <- array(seq_len(2 * 2 * 4) / 10 + offset, c(2, 2, 4))
  structure(
    list(
      IRmed = med,
      bands = list("90" = list(IRinf = med - 0.2, IRsup = med + 0.2)),
      conf = 90,
      varnames = c("y1", "y2")
    ),
    class = "fSignRestr"
  )
}

fIRF_cols <- c("model", "shock", "variable", "horizon", "band", "point",
               "lower", "upper")

# Values of column `col` for one shock / variable / band, in horizon order.
pick <- function(d, col, shock = NULL, variable = NULL, band = NULL) {
  keep <- rep(TRUE, nrow(d))
  if (!is.null(shock))    keep <- keep & d$shock == shock
  if (!is.null(variable)) keep <- keep & d$variable == variable
  if (!is.null(band))     keep <- keep & d$band == band
  d <- d[keep, ]
  d[[col]][order(d$horizon)]
}

# ---- fIRF: the common format ----------------------------------------------

test_that("fIRF turns a point cube and bootstrap bands into the common format", {
  x <- plot_irf_fixture()
  d <- fIRF(x$cube, bands = c(x$cube_bands, list(conf = c(95, 68))),
            varnames = c("y1", "y2"), shocknames = c("s1", "s2"))
  expect_s3_class(d, "fIRF")
  expect_named(d, fIRF_cols)
  expect_equal(nrow(d), 2 * 2 * 4 * 2)        # shocks x vars x horizons x bands
  expect_setequal(unique(d$band), c("95%", "68%"))
  expect_equal(pick(d, "point", "s2", "y1", "95%"), as.numeric(x$cube[1, 2, ]))
  expect_equal(pick(d, "upper", "s2", "y1", "68%"), as.numeric(x$cube[1, 2, ]) + 0.1)
  expect_equal(pick(d, "lower", "s1", "y2", "95%"), as.numeric(x$cube[2, 1, ]) - 0.2)
  expect_equal(pick(d, "horizon", "s1", "y1", "95%"), 0:3)
})

test_that("fIRF defaults: names, conf levels and no bands", {
  x <- plot_irf_fixture()
  d <- fIRF(x$cube, bands = x$cube_bands)
  expect_setequal(unique(d$variable), c("y1", "y2"))
  expect_setequal(unique(d$shock), c("Shock 1", "Shock 2"))
  expect_setequal(unique(d$band), c("90%", "68%"))

  p <- fIRF(x$mat)                              # N x H: one shock, no bands
  expect_equal(nrow(p), 8L)
  expect_true(all(is.na(p$band)))
  expect_equal(unique(p$shock), "Shock")
})

test_that("fIRF cumulates the point estimate only", {
  x <- plot_irf_fixture()
  d <- fIRF(x$cube, bands = x$cube_bands, varnames = c("y1", "y2"),
            cumulate = "y1")
  expect_equal(pick(d, "point", "Shock 1", "y1", "90%"), cumsum(x$cube[1, 1, ]))
  expect_equal(pick(d, "point", "Shock 1", "y2", "90%"), as.numeric(x$cube[2, 1, ]))
  expect_equal(pick(d, "upper", "Shock 1", "y1", "90%"), x$cube[1, 1, ] + 0.2)
})

test_that("fIRF reads IV-MBB, heteroskedasticity and MSW result lists", {
  x <- plot_irf_fixture()
  iv  <- c(list(meanirf = x$mat), x$mat_bands)
  he  <- c(list(point = x$mat, meanirf = x$mat + 1), x$mat_bands)
  msw <- list(IRF = x$mat, Dmlbound = x$mat - 0.1, Dmubound = x$mat + 0.1,
              MSWlbound = x$mat - 0.2, MSWubound = x$mat + 0.2)

  d_iv <- fIRF(iv, varnames = c("y1", "y2"), shocknames = "IV", scale = 10)
  d_he <- fIRF(he, varnames = c("y1", "y2"), shocknames = "Hetero")
  d_ms <- fIRF(msw, varnames = c("y1", "y2"), shocknames = "IV")

  expect_equal(pick(d_iv, "point", variable = "y1", band = "90%"), 10 * x$mat[1, ])
  expect_equal(pick(d_he, "point", variable = "y2", band = "90%"), x$mat[2, ])
  expect_setequal(unique(d_ms$band), c("Anderson-Rubin", "Delta method"))
  expect_equal(pick(d_ms, "upper", variable = "y1", band = "Delta method"),
               x$mat[1, ] + 0.1)
  expect_error(fIRF(list(bootchol = 1, upper = 1, lower = 1)),
               "Pass the point")
  expect_error(fIRF(list(a = 1)), "unrecognised")
})

test_that("fIRF reads sign-restriction fits and rejects missing levels", {
  d <- fIRF(mock_sign_fit())
  expect_equal(nrow(d), 2 * 2 * 4)
  expect_equal(unique(d$band), "90%")
  expect_equal(pick(d, "point", "Shock 2", "y1"), mock_sign_fit()$IRmed[1, 2, ])
  expect_error(fIRF(mock_sign_fit(), conf = 68), "Band 68% is not available")
})

test_that("fIRF reads LP fits: outcomes for fLP, shock terms for fLPPanel", {
  df <- sim_lp_df()
  lp <- fLP(y ~ x, data = df, horizons = 3L, shock = "x", conf = c(90, 68),
            n_threads = 1L)
  d <- fIRF(lp)
  expect_equal(unique(d$shock), "x")
  expect_equal(unique(d$variable), "y")
  expect_equal(pick(d, "point", band = "68%"), as.numeric(lp$irfs[, 1]))
  expect_equal(pick(d, "upper", band = "90%"), as.numeric(lp$irfs_upper[["90"]][, 1]))

  pn <- fLPPanel(y ~ shock | unit + tt, data = sim_panel_df(),
                 panel_id = c("unit", "tt"), shock = "shock", horizons = 2L,
                 n_threads = 1L)
  dp <- fIRF(pn)
  expect_equal(unique(dp$variable), "y")
  expect_equal(unique(dp$shock), colnames(pn$irfs))
  expect_equal(pick(dp, "point"), as.numeric(pn$irfs[, 1]))
})

test_that("a negative scale keeps lower <= upper", {
  x <- plot_irf_fixture()
  d <- fIRF(x$cube, bands = x$cube_bands, scale = -2)
  expect_true(all(d$lower <= d$upper))
  expect_equal(pick(d, "point", "Shock 1", "y1", "90%"), -2 * x$cube[1, 1, ])

  d2 <- fIRF(x$cube, bands = x$cube_bands, varnames = c("a", "b"),
             scale = c(a = 1, b = -1))
  expect_true(all(d2$lower <= d2$upper))
  expect_equal(pick(d2, "point", "Shock 1", "b", "90%"), -x$cube[2, 1, ])
  expect_error(fIRF(x$cube, scale = c(1, 2, 3)), "length 2")
})

test_that("fIRF checks band dimensions", {
  x <- plot_irf_fixture()
  expect_error(fIRF(x$cube, bands = x$mat_bands), "does not match")
  expect_error(fIRF(x$cube, varnames = "a"), "length 2")
})

# ---- fPlot -----------------------------------------------------------------

test_that("fPlot overlays several results and labels them by argument name", {
  x <- plot_irf_fixture()
  a <- fIRF(x$cube, bands = x$cube_bands, varnames = c("y1", "y2"))
  b <- fIRF(x$cube * 2, bands = x$cube_bands, varnames = c("y1", "y2"))
  d <- fPlot(A = a, B = b, shock = 1, return_data = TRUE)
  expect_s3_class(d, "fIRF")
  expect_equal(unique(d$model), c("A", "B"))
  expect_equal(unique(d$shock), "Shock 1")
  expect_s3_class(fPlot(A = a, B = b, shock = 1), "ggplot")

  # Unnamed results get generic labels
  d2 <- fPlot(a, b, shock = 2, return_data = TRUE)
  expect_equal(unique(d2$model), c("Model 1", "Model 2"))
})

test_that("fPlot selects shocks and variables and validates them", {
  x <- plot_irf_fixture()
  a <- fIRF(x$cube, bands = x$cube_bands, varnames = c("y1", "y2"),
            shocknames = c("s1", "s2"))
  expect_equal(unique(fPlot(a, shock = "s2", return_data = TRUE)$shock), "s2")
  expect_equal(unique(fPlot(a, shock = 2:1, return_data = TRUE)$shock),
               c("s2", "s1"))
  expect_equal(unique(fPlot(a, variables = "y2", return_data = TRUE)$variable),
               "y2")
  expect_error(fPlot(a, shock = 3L), "shock must be between 1 and 2")
  expect_error(fPlot(a, shock = "s9"), "not found")
  expect_error(fPlot(a, variables = "zz"), "not found")
  expect_error(fPlot(), "at least one")
})

test_that("fPlot facets by variable for one shock and shock x variable otherwise", {
  x <- plot_irf_fixture()
  a <- fIRF(x$cube, bands = x$cube_bands)
  one <- fPlot(a, shock = 1)
  all <- fPlot(a)
  expect_equal(names(one$facet$params$facets), "variable")
  expect_equal(names(all$facet$params$facets), c("shock", "variable"))
})

test_that("fPlot draws results without bands and accepts fits directly", {
  x <- plot_irf_fixture()
  expect_s3_class(fPlot(fIRF(x$mat)), "ggplot")
  expect_s3_class(fPlot(mock_sign_fit(), shock = 1), "ggplot")
  expect_error(fPlot(ribbon_alpha = 1:3, fIRF(x$cube, bands = x$cube_bands)),
               "ribbon_alpha")
})

# ---- legacy wrappers ---------------------------------------------------------

test_that("legacy IRF plotters are wrappers around fIRF + fPlot", {
  x <- plot_irf_fixture()
  vn <- c("y1", "y2")

  d <- fPlotIRFChol(x$cube, x$cube_bands, shock = 2L, varnames = vn,
                    return_data = TRUE)
  expect_equal(d, fPlot(fIRF(x$cube, bands = x$cube_bands, varnames = vn),
                        shock = 2L, return_data = TRUE))
  expect_s3_class(fPlotIRFChol(x$cube, x$cube_bands, 2L, vn), "ggplot")
  expect_error(fPlotIRFChol(x$cube, x$cube_bands, 3L, vn),
               "shock must be between")

  d <- fPlotIRFBQ(x$cube, x$cube_bands, vn, c("s1", "s2"),
                  cumulate = 1L, shocks = 1:2, return_data = TRUE)
  expect_equal(pick(d, "point", "s1", "y1", "90%"), cumsum(x$cube[1, 1, ]))
  expect_s3_class(fPlotIRFBQ(x$cube, x$cube_bands, vn, c("s1", "s2"),
                             cumulate = 1L, shocks = 1:2), "ggplot")

  d <- fPlotIRFLR(x$mat, x$mat_bands, vn, return_data = TRUE)
  d0 <- fPlotIRFLR(x$mat, varnames = vn, return_data = TRUE)
  expect_equal(nrow(d), 16L)
  expect_true(all(is.na(d0$band)))
  expect_s3_class(fPlotIRFLR(x$mat, x$mat_bands, vn), "ggplot")

  iv <- c(list(meanirf = x$mat), x$mat_bands)
  he <- c(list(point = x$mat), x$mat_bands)
  expect_equal(pick(fPlotIRFIV(iv, vn, "IV", scale = 10, return_data = TRUE),
                    "point", variable = "y1", band = "90%"), 10 * x$mat[1, ])
  expect_equal(pick(fPlotIRFHetero(he, vn, "H", scale = 10, return_data = TRUE),
                    "point", variable = "y2", band = "90%"), 10 * x$mat[2, ])
  expect_s3_class(fPlotIRFIV(iv, vn, "IV"), "ggplot")
  expect_s3_class(fPlotIRFHetero(he, vn, "H"), "ggplot")

  msw <- list(IRF = x$mat, Dmlbound = x$mat - 0.1, Dmubound = x$mat + 0.1,
              MSWlbound = x$mat - 0.2, MSWubound = x$mat + 0.2)
  d <- fPlotIRFMSW(msw, vn, "IV", scale = -1, return_data = TRUE)
  expect_setequal(unique(d$band), c("Anderson-Rubin", "Delta method"))
  expect_true(all(d$lower <= d$upper))
  expect_s3_class(fPlotIRFMSW(msw, vn, "IV"), "ggplot")
})

test_that("fPlotIRFSign maps one or two fitted models", {
  a <- mock_sign_fit()
  b <- mock_sign_fit(0.1)
  d <- fPlotIRFSign(a, shock = 1L, compare = b, labels = c("A", "B"),
                    return_data = TRUE)
  expect_equal(nrow(d), 16L)
  expect_setequal(unique(d$model), c("A", "B"))
  expect_equal(pick(d[d$model == "A", ], "horizon", variable = "y1"), 0:3)
  expect_s3_class(fPlotIRFSign(a, shock = 1L), "ggplot")
  expect_error(fPlotIRFSign(list()), "must be an 'fSignRestr'")
  d <- fPlotIRFSign(a, shock = 1, scale = -1, return_data = TRUE)
  expect_true(all(d$lower <= d$upper))
})

test_that("bootstrap results report their band levels", {
  y <- sim_var(T = 120)
  fit <- fVAR(y, p = 1, c = 1)
  set.seed(1)
  b <- suppressMessages(fBootstrapChol(y, fit, nboot = 20, horizon = 3,
                                       conf = 95, conf2 = 80, n_threads = 1L))
  expect_equal(b$conf, c(95, 80))
  wold <- fWoldIRF(fit, horizon = 3)
  d <- fIRF(fCholeskyIRF(wold, t(chol(fit$sigma))), bands = b)
  expect_setequal(unique(d$band), c("95%", "80%"))
})
