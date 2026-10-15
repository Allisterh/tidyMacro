#' Plot Local Projection Impulse Responses
#'
#' Plots impulse responses and confidence bands from an object returned by
#' \code{fLP()}. Handles both single-band and multi-band fits transparently:
#' if \code{fLP()} was called with \code{conf = c(68, 95)} the plot
#' automatically draws one ribbon per level (widest ribbon is drawn first with
#' lowest opacity so the tighter band is visually on top).
#'
#' @param x An object of class \code{"fLP"} returned by \code{fLP()}.
#' @param variables Optional character vector of response variables to plot.
#'   Defaults to all responses.
#' @param scale Numeric scalar applied to IRFs and confidence bands before
#'   plotting. Default is 1.
#' @param return_data Logical. If \code{TRUE}, return the plotting data
#'   (long format when multi-band) instead of a ggplot object.
#' @param ribbon_fill Color for the confidence band fill.
#' @param ribbon_alpha Transparency for the confidence band. For multi-band
#'   fits, either a single number (interpreted as the alpha of the widest band;
#'   inner bands get proportionally darker) or a numeric vector of the same
#'   length as the confidence levels, matched widest-to-tightest.
#' @param line_color Color for the IRF line.
#' @param zero_line_color Color for the horizontal zero line.
#' @param facet_scales Facet scales option: \code{"free_y"}, \code{"free"},
#'   \code{"free_x"}, or \code{"fixed"}.
#' @param facet_ncol Number of facet columns. Defaults to an automatic layout.
#'
#' @return A ggplot object, or (if \code{return_data = TRUE}) the
#'   \code{\link{fIRF}} tibble behind it.
#'
#' @details A thin wrapper around \code{\link{fIRF}} and \code{\link{fPlot}}.
#'   Use \code{fPlot()} directly to compare several results in one figure.
#'
#' @seealso \code{\link{fPlot}}, \code{\link{fIRF}}
#'
#' @export
fPlotLP <- function(x,
                    variables       = NULL,
                    scale           = 1,
                    return_data     = FALSE,
                    ribbon_fill     = "#407EC9",
                    ribbon_alpha    = 0.2,
                    line_color      = "#910048",
                    zero_line_color = "#707372",
                    facet_scales    = "free_y",
                    facet_ncol      = NULL) {
  if (!inherits(x, "fLP")) {
    stop("fPlotLP: 'x' must be an object returned by fLP().", call. = FALSE)
  }
  fPlot(fIRF(x, variables = variables, scale = scale),
        return_data = return_data,
        ribbon_fill = ribbon_fill, line_color = line_color,
        ribbon_alpha = ribbon_alpha, zero_line_color = zero_line_color,
        facet_scales = facet_scales, facet_ncol = facet_ncol)
}
