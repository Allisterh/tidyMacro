#' Plot MSW Weak-IV Robust IRFs with Two Confidence Bands
#'
#' Creates a faceted IRF plot showing both the plug-in (delta method) and
#' Anderson-Rubin weak-IV robust confidence sets from \code{fMSW()}.
#' The outer (lighter) band is the Anderson-Rubin set; the inner (darker) band
#' is the delta-method interval.
#'
#' @param msw_result List returned by \code{fMSW()}. Must contain
#'   \code{IRF}, \code{Dmlbound}, \code{Dmubound}, \code{MSWlbound},
#'   \code{MSWubound}.
#' @param varnames Character vector of variable names (length N).
#' @param shockname Character string for the shock label.
#' @param scale Numeric scalar applied to all values before plotting. Default 1.
#' @param return_data Logical; if \code{TRUE} returns a tibble instead of a
#'   ggplot object. Default \code{FALSE}.
#' @param ribbon_fill_ar  Fill colour for Anderson-Rubin band. Default "#407EC9".
#' @param ribbon_fill_dm  Fill colour for delta-method band. Default "#407EC9".
#' @param ribbon_alpha_ar Transparency for AR outer band. Default 0.18.
#' @param ribbon_alpha_dm Transparency for delta-method inner band. Default 0.35.
#' @param line_color      Colour for point estimate line. Default "#910048".
#' @param zero_line_color Colour for zero reference line. Default "#707372".
#' @param facet_scales    Facet scales: "free", "free_y", "free_x", or "fixed".
#'   Default "free".
#' @param facet_ncol Number of columns in the facet grid. Default 2.
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
fPlotIRFMSW <- function(msw_result, varnames, shockname,
                         scale           = 1,
                         return_data     = FALSE,
                         ribbon_fill_ar  = "#407EC9",
                         ribbon_fill_dm  = ribbon_fill_ar,
                         ribbon_alpha_ar = 0.18,
                         ribbon_alpha_dm = 0.35,
                         line_color      = "#910048",
                         zero_line_color = "#707372",
                         facet_scales    = "free",
                         facet_ncol      = 2) {
    if (!identical(ribbon_fill_dm, ribbon_fill_ar))
        warning("fPlotIRFMSW: both bands now share `ribbon_fill_ar`; ",
                "they are told apart by transparency.", call. = FALSE)
    fPlot(fIRF(msw_result, varnames = varnames, shocknames = shockname,
               scale = scale),
          return_data = return_data,
          ribbon_fill = ribbon_fill_ar, line_color = line_color,
          ribbon_alpha = c(ribbon_alpha_ar, ribbon_alpha_dm),
          zero_line_color = zero_line_color,
          facet_scales = facet_scales, facet_ncol = facet_ncol)
}
