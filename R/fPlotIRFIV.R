#' Plot IV/Proxy SVAR Impulse Response Functions with Confidence Bands
#'
#' Creates a grid of IRF plots with shaded confidence bands using ggplot2.
#' Accepts the full list returned by \code{fBootstrapIVMBB()} directly.
#'
#' @param result_mbb List returned by \code{fBootstrapIVMBB()}. Must contain
#'   \code{meanirf}, \code{upper}, and \code{lower}. If \code{upper2} and
#'   \code{lower2} are present, a second (inner) confidence band is drawn.
#' @param varnames Character vector of variable names (length N)
#' @param shockname Character string for the shock label
#' @param scale Numeric scalar applied to all IRF values before plotting.
#'   Useful for unit rescaling (e.g. \code{scale = 10} to express as percent
#'   of a 10-unit shock). Default 1.
#' @param return_data Logical, if TRUE returns tibble instead of plot. Default FALSE.
#' @param ribbon_fill Color for confidence band fill (default: "#407EC9")
#' @param ribbon_alpha Transparency for outer confidence band (default: 0.2)
#' @param ribbon_alpha2 Transparency for inner confidence band (default: 0.35)
#' @param line_color Color for point estimate line (default: "#910048")
#' @param zero_line_color Color for horizontal zero line (default: "#707372")
#' @param facet_scales Facet scales option: "free", "free_y", "free_x", or
#'   "fixed" (default: "free")
#' @param facet_ncol Number of columns in facet grid (default: 2)
#'
#' @return A ggplot object, or (if \code{return_data = TRUE}) the
#'   \code{\link{fIRF}} tibble behind it.
#'
#' @details \strong{Deprecated.} \code{fPlotIRFIV()} will be removed in the next
#'   release. Use \code{\link{fPlotIRF}} instead:
#'   \preformatted{fPlotIRF(fIRF(result_mbb, varnames = varnames, shocknames = shockname, scale = scale))}
#'
#' @seealso \code{\link{fPlotIRF}}, \code{\link{fIRF}}
#'
#' @export
fPlotIRFIV <- function(result_mbb, varnames, shockname,
                        scale           = 1,
                        return_data     = FALSE,
                        ribbon_fill     = "#407EC9",
                        ribbon_alpha    = 0.2,
                        ribbon_alpha2   = 0.35,
                        line_color      = "#910048",
                        zero_line_color = "#707372",
                        facet_scales    = "free",
                        facet_ncol      = 2) {
    .fDeprecatedPlot("fPlotIRFIV",
        "fPlotIRF(fIRF(result_mbb, varnames = varnames, shocknames = shockname, scale = scale))")
    fPlotIRF(fIRF(result_mbb, varnames = varnames, shocknames = shockname,
                  scale = scale),
             return_data = return_data,
             ribbon_fill = ribbon_fill, line_color = line_color,
             ribbon_alpha = if (is.null(result_mbb$upper2)) ribbon_alpha
                            else c(ribbon_alpha, ribbon_alpha2),
             zero_line_color = zero_line_color,
             facet_scales = facet_scales, facet_ncol = facet_ncol)
}
