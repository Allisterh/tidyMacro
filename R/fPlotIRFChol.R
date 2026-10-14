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
#' @return A ggplot object (if \code{return_data = FALSE}) or a tibble.
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

    N <- dim(point)[1]
    if (shock < 1 || shock > N)
        stop(sprintf("shock must be between 1 and %d", N))
    if (is.null(facet_ncol))
        facet_ncol <- ceiling(sqrt(N))

    two <- !is.null(boot_result$upper2) && !is.null(boot_result$lower2)
    irf_data <- .fIRFLong(varnames, list(
        point  = .fShockMat(point, shock),
        upper  = .fShockMat(boot_result$upper, shock),
        lower  = .fShockMat(boot_result$lower, shock),
        upper2 = if (two) .fShockMat(boot_result$upper2, shock),
        lower2 = if (two) .fShockMat(boot_result$lower2, shock)
    ))

    if (return_data)
        return(irf_data)

    .fIRFPlot(irf_data, ~ variable, facet_scales, facet_ncol, ribbon_fill,
              ribbon_alpha, ribbon_alpha2, line_color, zero_line_color)
}
