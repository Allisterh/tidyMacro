#' Plot Blanchard-Quah Long-Run Identified Impulse Response Functions
#'
#' Creates a grid of IRF plots with shaded confidence bands using ggplot2.
#' Designed for BQ long-run identified VARs. Supports one or multiple shocks;
#' when a single shock is supplied the facet shows only response variables.
#'
#' @param point Array of point estimates (N x N x horizon+1), typically from
#'   \code{fBQIRF()}. Variables listed in \code{cumulate} will be cumulated
#'   inside this function to recover level responses.
#' @param boot_result List returned by \code{fBootstrapBQ()}. Must contain
#'   \code{upper} and \code{lower}. If \code{upper2} and \code{lower2} are
#'   present, a second (inner) confidence band is drawn automatically.
#' @param varnames Character vector of variable names (length N).
#' @param shocknames Character vector of shock names (length equal to
#'   \code{length(shocks)}).
#' @param cumulate Integer vector (1-based) of variable indices whose point
#'   IRFs should be cumulated. Must match the \code{cumulate} argument passed
#'   to \code{fBootstrapBQ()}.
#' @param shocks Integer vector (1-based) of shock indices to plot.
#' @param return_data Logical, if TRUE returns a tibble instead of a plot.
#'   Default FALSE.
#' @param ribbon_fill Color for confidence band fill (default: "#407EC9")
#' @param ribbon_alpha Transparency for outer confidence band (default: 0.2)
#' @param ribbon_alpha2 Transparency for inner confidence band (default: 0.35)
#' @param line_color Color for point estimate line (default: "#910048")
#' @param zero_line_color Color for horizontal zero line (default: "#707372")
#' @param facet_scales Scale option for \code{facet_wrap}: one of
#'   \code{"free"} (default), \code{"free_y"}, \code{"free_x"}, or
#'   \code{"fixed"}.
#' @param facet_ncol Number of columns in the facet layout. Default NULL.
#'
#' @return A ggplot object (if \code{return_data = FALSE}) or a tibble.
#'
#' @export
fPlotIRFBQ <- function(point, boot_result, varnames, shocknames,
                        cumulate, shocks,
                        return_data     = FALSE,
                        ribbon_fill     = "#407EC9",
                        ribbon_alpha    = 0.2,
                        ribbon_alpha2   = 0.35,
                        line_color      = "#910048",
                        zero_line_color = "#707372",
                        facet_scales    = c("free", "free_y", "free_x", "fixed"),
                        facet_ncol      = NULL) {

    facet_scales <- match.arg(facet_scales)

    N        <- dim(point)[1]
    n_shocks <- length(shocks)

    if (is.null(facet_ncol))
        facet_ncol <- if (n_shocks == 1L) ceiling(sqrt(N)) else n_shocks

    if (length(shocknames) != n_shocks)
        stop("shocknames must have the same length as shocks")
    if (any(shocks < 1) || any(shocks > N))
        stop(sprintf("All shocks must be between 1 and %d", N))
    if (any(cumulate < 1) || any(cumulate > N))
        stop(sprintf("All cumulate indices must be between 1 and %d", N))

    two <- !is.null(boot_result$upper2) && !is.null(boot_result$lower2)
    irf_data <- do.call(rbind, lapply(seq_along(shocks), function(j) {
        s <- shocks[j]
        # Cumulate the point IRF (bands are already cumulated in the bootstrap)
        pt <- .fShockMat(point, s)
        for (ci in cumulate) pt[ci, ] <- cumsum(pt[ci, ])
        d <- .fIRFLong(varnames, list(
            point  = pt,
            upper  = .fShockMat(boot_result$upper, s),
            lower  = .fShockMat(boot_result$lower, s),
            upper2 = if (two) .fShockMat(boot_result$upper2, s),
            lower2 = if (two) .fShockMat(boot_result$lower2, s)
        ))
        tibble::add_column(d, shock = shocknames[j], .after = "variable")
    }))
    irf_data$shock <- factor(irf_data$shock, levels = shocknames)

    if (return_data)
        return(irf_data)

    facet <- if (n_shocks == 1L) ~ variable else ~ shock + variable
    .fIRFPlot(irf_data, facet, facet_scales, facet_ncol, ribbon_fill,
              ribbon_alpha, ribbon_alpha2, line_color, zero_line_color) +
        ggplot2::labs(x = NULL, y = NULL, caption = NULL)
}
