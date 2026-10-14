#' Plot Heteroskedasticity Bootstrap IRFs with Confidence Bands
#'
#' Creates a grid of IRF plots for results from \code{fBootstrapHetero()}.
#' Uses the plug-in point estimate (\code{result$point}) as the center line,
#' which is the correct center for Hall (1992) recentered confidence bands.
#'
#' @param result      List returned by \code{fBootstrapHetero()}. Must contain
#'   \code{point}, \code{upper}, \code{lower}. If \code{upper2} and
#'   \code{lower2} are present, a second (inner) band is drawn.
#' @param varnames    Character vector of variable names (length N).
#' @param shockname   Character string for the shock label.
#' @param scale       Numeric scalar applied to all IRF values. Default 1.
#' @param return_data Logical; if \code{TRUE} returns a tibble. Default \code{FALSE}.
#' @param ribbon_fill    Fill colour for confidence bands. Default \code{"#407EC9"}.
#' @param ribbon_alpha   Transparency for outer band. Default 0.2.
#' @param ribbon_alpha2  Transparency for inner band. Default 0.35.
#' @param line_color     Colour for point estimate line. Default \code{"#910048"}.
#' @param zero_line_color Colour for zero reference line. Default \code{"#707372"}.
#' @param facet_scales   Facet scales: \code{"free"}, \code{"free_y"},
#'   \code{"free_x"}, or \code{"fixed"}. Default \code{"free"}.
#' @param facet_ncol     Number of columns in facet grid. Default 2.
#'
#' @return A ggplot object (or a tibble if \code{return_data = TRUE}).
#'
#' @export
fPlotIRFHetero <- function(result, varnames, shockname,
                           scale           = 1,
                           return_data     = FALSE,
                           ribbon_fill     = "#407EC9",
                           ribbon_alpha    = 0.2,
                           ribbon_alpha2   = 0.35,
                           line_color      = "#910048",
                           zero_line_color = "#707372",
                           facet_scales    = "free",
                           facet_ncol      = 2) {

    two   <- !is.null(result$upper2) && !is.null(result$lower2)
    band  <- .fScaleBand(result$upper, result$lower, scale)
    band2 <- if (two) .fScaleBand(result$upper2, result$lower2, scale)

    irf_data <- .fIRFLong(varnames, list(
        point  = result$point * scale,
        upper  = band$upper,
        lower  = band$lower,
        shock  = shockname,
        upper2 = band2$upper,
        lower2 = band2$lower
    ))

    if (return_data)
        return(irf_data)

    .fIRFPlot(irf_data, ~ variable, facet_scales, facet_ncol, ribbon_fill,
              ribbon_alpha, ribbon_alpha2, line_color, zero_line_color)
}
