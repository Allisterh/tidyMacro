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

test_that("fPlotIRFChol maps a selected shock to plot data", {
  x <- plot_irf_fixture()
  d <- fPlotIRFChol(x$cube, x$cube_bands, shock = 2L,
                    varnames = c("y1", "y2"), return_data = TRUE)
  expect_s3_class(d, "tbl_df")
  expect_equal(nrow(d), 8L)
  expect_equal(d$point[d$variable == "y1"], as.numeric(x$cube[1, 2, ]))
  expect_s3_class(fPlotIRFChol(x$cube, x$cube_bands, 2L, c("y1", "y2")),
                  "ggplot")
  expect_error(fPlotIRFChol(x$cube, x$cube_bands, 3L, c("y1", "y2")),
               "shock must be between")
})

test_that("fPlotIRFBQ handles multiple shocks and cumulative variables", {
  x <- plot_irf_fixture()
  d <- fPlotIRFBQ(x$cube, x$cube_bands, c("y1", "y2"), c("s1", "s2"),
                  cumulate = 1L, shocks = 1:2, return_data = TRUE)
  expect_equal(nrow(d), 16L)
  keep <- d$variable == "y1" & d$shock == "s1"
  expect_equal(d$point[keep], cumsum(x$cube[1, 1, ]))
  expect_s3_class(
    fPlotIRFBQ(x$cube, x$cube_bands, c("y1", "y2"), c("s1", "s2"),
               cumulate = 1L, shocks = 1:2),
    "ggplot"
  )
})

test_that("fPlotIRFIV and fPlotIRFHetero return scaled data and ggplots", {
  x <- plot_irf_fixture()
  iv <- c(list(meanirf = x$mat), x$mat_bands)
  he <- c(list(point = x$mat), x$mat_bands)
  d_iv <- fPlotIRFIV(iv, c("y1", "y2"), "IV", scale = 10,
                     return_data = TRUE)
  d_he <- fPlotIRFHetero(he, c("y1", "y2"), "Hetero", scale = 10,
                         return_data = TRUE)

  expect_equal(d_iv$point[d_iv$variable == "y1"], 10 * x$mat[1, ])
  expect_equal(d_he$point[d_he$variable == "y2"], 10 * x$mat[2, ])
  expect_s3_class(fPlotIRFIV(iv, c("y1", "y2"), "IV"), "ggplot")
  expect_s3_class(fPlotIRFHetero(he, c("y1", "y2"), "Hetero"), "ggplot")
})

test_that("fPlotIRFMSW exposes both delta-method and Anderson-Rubin bands", {
  x <- plot_irf_fixture()
  m <- list(
    IRF = x$mat,
    Dmlbound = x$mat - 0.1,
    Dmubound = x$mat + 0.1,
    MSWlbound = x$mat - 0.2,
    MSWubound = x$mat + 0.2
  )
  d <- fPlotIRFMSW(m, c("y1", "y2"), "IV", return_data = TRUE)
  expect_named(d, c("variable", "horizon", "point", "dm_lo", "dm_hi",
                    "ar_lo", "ar_hi", "shock"))
  expect_equal(nrow(d), 8L)
  expect_s3_class(fPlotIRFMSW(m, c("y1", "y2"), "IV"), "ggplot")
})

test_that("fPlotIRFLR works with and without bootstrap bands", {
  x <- plot_irf_fixture()
  d <- fPlotIRFLR(x$mat, x$mat_bands, c("y1", "y2"), return_data = TRUE)
  d0 <- fPlotIRFLR(x$mat, varnames = c("y1", "y2"), return_data = TRUE)
  expect_equal(nrow(d), 8L)
  expect_true(all(c("upper", "lower", "upper2", "lower2") %in% names(d)))
  expect_false(any(c("upper", "lower") %in% names(d0)))
  expect_s3_class(fPlotIRFLR(x$mat, x$mat_bands, c("y1", "y2")), "ggplot")
})

test_that("fPlotIRFSign maps one or two fitted models", {
  a <- mock_sign_fit()
  b <- mock_sign_fit(0.1)
  d <- fPlotIRFSign(a, shock = 1L, compare = b, labels = c("A", "B"),
                    return_data = TRUE)
  expect_equal(nrow(d), 16L)
  expect_setequal(as.character(unique(d$model)), c("A", "B"))
  expect_equal(d$horizon[d$model == "A" & d$variable == "y1"], 0:3)
  expect_s3_class(fPlotIRFSign(a, shock = 1L), "ggplot")
  expect_error(fPlotIRFSign(list()), "must be an 'fSignRestr'")
})
