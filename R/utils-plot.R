# Scale a confidence band. A negative scale flips the sign, so pmin/pmax keep
# lower <= upper.
.fScaleBand <- function(upper, lower, scale) {
  u <- upper * scale
  l <- lower * scale
  list(upper = pmax(u, l), lower = pmin(u, l))
}

# Responses of all variables to shock s, as an N x H matrix, from an
# N x N x H array.
.fShockMat <- function(a, s) {
  matrix(a[, s, ], dim(a)[1])
}

# Long IRF data with one block of rows per variable: variable, horizon, then
# `cols` in the given order. N x H matrices contribute their row i; other
# values (e.g. a shock label) are repeated. NULL entries are dropped.
.fIRFLong <- function(varnames, cols) {
  cols <- Filter(Negate(is.null), cols)
  N <- nrow(cols[[1]])
  H <- ncol(cols[[1]])
  if (length(varnames) != N)
    stop(sprintf("varnames must have length %d (number of variables)", N))

  rows <- lapply(seq_len(N), function(i) {
    vals <- lapply(cols, function(m) if (is.matrix(m)) as.numeric(m[i, ]) else m)
    tibble::as_tibble(c(list(variable = varnames[i], horizon = 0:(H - 1)), vals))
  })
  d <- do.call(rbind, rows)
  d$variable <- factor(d$variable, levels = varnames)
  d
}

# IRF panel plot shared by the fPlotIRF* functions: outer band (lower/upper),
# inner band (lower2/upper2) when present, point line, zero line, facets.
.fIRFPlot <- function(d, facet, facet_scales, facet_ncol,
                      ribbon_fill, ribbon_alpha, ribbon_alpha2,
                      line_color, zero_line_color) {
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$horizon))
  if ("upper" %in% names(d))
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                                  fill = ribbon_fill, alpha = ribbon_alpha)
  if ("upper2" %in% names(d))
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower2, ymax = .data$upper2),
                                  fill = ribbon_fill, alpha = ribbon_alpha2)
  p +
    ggplot2::geom_line(ggplot2::aes(y = .data$point),
                       color = line_color, linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, color = zero_line_color,
                        linetype = "dashed", linewidth = 0.6) +
    ggplot2::facet_wrap(facet, scales = facet_scales, ncol = facet_ncol)
}
