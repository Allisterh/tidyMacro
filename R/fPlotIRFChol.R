#' Plot Cholesky Impulse Response Functions with Confidence Bands
#'
#' Creates a grid of IRF plots showing responses of all variables to a single
#' shock, with shaded confidence bands using ggplot2. Designed for
#' Cholesky-identified VARs.
#'
#' @param point Array of point estimates (N x N x horizon), typically from
#'   \code{fCholeskyIRF()}.
#' @param boot_result List returned by \code{fBootstrapChol()} or
#'   \code{fBootstrapCholCorrected()}. Must contain \code{upper} and
#'   \code{lower}. If \code{upper2} and \code{lower2} are present, a second
#'   (inner) confidence band is drawn automatically.
#' @param shock Integer, which shock to plot (1-indexed, 1 to N).
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
#' @details \strong{Deprecated.} \code{fPlotIRFChol()} will be removed in the next
#'   release. Use \code{\link{fPlotIRF}} instead:
#'   \preformatted{fPlotIRF(fIRF(point, bands = boot_result, varnames = varnames), shock = shock)}
#'
#' @seealso \code{\link{fPlotIRF}}, \code{\link{fIRF}}
#'
#' @export
fPlotIRFChol <- function(point, boot_result, shock, varnames,
                          return_data     = FALSE,
                          ribbon_fill     = "#407EC9",
                          ribbon_alpha    = 0.2,
                          ribbon_alpha2   = 0.35,
                          line_color      = "#910048",
                          zero_line_color = "#707372",
                          facet_scales    = "free",
                          facet_ncol      = NULL) {
    .fDeprecatedPlot("fPlotIRFChol",
        "fPlotIRF(fIRF(point, bands = boot_result, varnames = varnames), shock = shock)")
    fPlotIRF(fIRF(point, bands = boot_result, varnames = varnames),
             shock = shock, return_data = return_data,
             ribbon_fill = ribbon_fill, line_color = line_color,
             ribbon_alpha = c(ribbon_alpha, ribbon_alpha2),
             zero_line_color = zero_line_color,
             facet_scales = facet_scales, facet_ncol = facet_ncol)
}
