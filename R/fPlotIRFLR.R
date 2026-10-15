#' Plot Long-Run Identified Impulse Response Functions with Confidence Bands
#'
#' Creates a grid of IRF plots with shaded confidence bands using ggplot2.
#' Designed for long-run identified VARs where the point estimate is an N x H
#' matrix — works for both LR-Max (\code{fMaxIRF()}) and Uhlig max-share
#' (\code{fUhligIRF()}) identifications.
#'
#' @param point N x (horizon+1) matrix of point estimates, from
#'   \code{fMaxIRF()} or \code{fUhligIRF()}.
#' @param boot_result List returned by \code{fBootstrapMax()},
#'   \code{fBootstrapMaxCorrected()}, \code{fBootstrapUhlig()}, or
#'   \code{fBootstrapUhligCorrected()}. Must contain \code{upper} and
#'   \code{lower}. If \code{upper2} and \code{lower2} are present, a second
#'   (inner) confidence band is drawn automatically.
#' @param varnames Character vector of variable names (length N).
#' @param return_data Logical, if TRUE returns tibble instead of plot. Default FALSE.
#' @param ribbon_fill Color for confidence band fill (default: "#407EC9")
#' @param ribbon_alpha Transparency for outer confidence band (default: 0.2)
#' @param ribbon_alpha2 Transparency for inner confidence band (default: 0.35)
#' @param line_color Color for point estimate line (default: "#910048")
#' @param zero_line_color Color for horizontal zero line (default: "#707372")
#' @param facet_scales Facet scales option: "free", "free_y", "free_x", or
#'   "fixed" (default: "free")
#' @param facet_ncol Number of columns in facet grid (default: NULL, auto-calculated)
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
fPlotIRFLR <- function(point, boot_result = NULL, varnames,
                       return_data     = FALSE,
                       ribbon_fill     = "#407EC9",
                       ribbon_alpha    = 0.2,
                       ribbon_alpha2   = 0.35,
                       line_color      = "#910048",
                       zero_line_color = "#707372",
                       facet_scales    = "free",
                       facet_ncol      = NULL) {
    fPlot(fIRF(point, bands = boot_result, varnames = varnames),
          return_data = return_data,
          ribbon_fill = ribbon_fill, line_color = line_color,
          ribbon_alpha = if (is.null(boot_result$upper2)) ribbon_alpha
                         else c(ribbon_alpha, ribbon_alpha2),
          zero_line_color = zero_line_color,
          facet_scales = facet_scales, facet_ncol = facet_ncol)
}
