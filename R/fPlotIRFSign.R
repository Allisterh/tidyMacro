#' Plot Impulse Responses from a Sign-Restricted SVAR
#'
#' Draws the posterior median and credible band of the responses to one
#' structural shock, for every variable in the system. Two models can be
#' overlaid, which is the standard way of showing what narrative restrictions
#' add to plain sign restrictions.
#'
#' @param sr An object of class \code{"fSignRestr"} from
#'   \code{\link{fSignRestr}}.
#' @param shock Integer (1-indexed) structural shock to plot (default 1).
#' @param compare Optional second \code{"fSignRestr"} object drawn on the same
#'   panels for comparison (default NULL).
#' @param labels Character vector of length 1 or 2 naming the models in the
#'   legend (default \code{c("Sign", "Sign + narrative")}).
#' @param varnames Optional character vector overriding the variable names.
#' @param scale Numeric vector of length 1 or N multiplying the responses, e.g.
#'   to express a variable in basis points (default 1).
#' @param scale_compare Numeric scaling applied to \code{compare}; defaults to
#'   \code{scale}.
#' @param return_data Logical; return the plotting tibble instead of the plot
#'   (default FALSE).
#' @param colors Length-2 character vector of colours for the two models.
#' @param ribbon_alpha Transparency of the credible bands (default 0.18).
#' @param conf Optional numeric vector of band levels to draw, in percent. Each
#'   must have been requested in the \code{fSignRestr} call. Defaults to the
#'   first level of the fitted object. Several levels are drawn as nested
#'   ribbons, the widest most transparent.
#' @param facet_scales Facet scale option passed to
#'   \code{ggplot2::facet_wrap} (default "free_y").
#' @param facet_ncol Number of facet columns (default NULL, auto).
#'
#' @return A ggplot object, or (if \code{return_data = TRUE}) the
#'   \code{\link{fIRF}} tibble behind it.
#'
#' @details \strong{Deprecated.} \code{fPlotIRFSign()} will be removed in the next
#'   release. Use \code{\link{fPlotIRF}} instead:
#'   \preformatted{fPlotIRF(sr, shock = shock), or fPlotIRF(Sign = sr, Narrative = compare, shock = shock)}
#'
#' @seealso \code{\link{fPlotIRF}}, \code{\link{fIRF}}, \code{\link{fSignRestr}}, \code{\link{fPlotVarDec}}
#'
#' @export
fPlotIRFSign <- function(sr, shock = 1, compare = NULL,
                         labels        = c("Sign", "Sign + narrative"),
                         varnames      = NULL,
                         scale         = 1,
                         scale_compare = NULL,
                         return_data   = FALSE,
                         colors        = c("#407EC9", "#E4002B"),
                         ribbon_alpha  = 0.18,
                         conf          = NULL,
                         facet_scales  = "free_y",
                         facet_ncol    = NULL) {
    .fDeprecatedPlot("fPlotIRFSign",
        "fPlotIRF(sr, shock = shock), or fPlotIRF(Sign = sr, Narrative = compare, shock = shock)")
    if (!inherits(sr, "fSignRestr")) stop("`sr` must be an 'fSignRestr' object.")
    if (is.null(conf)) conf <- sr$conf[1]

    fits <- list(fIRF(sr, conf = conf, varnames = varnames, scale = scale))
    names(fits) <- labels[1]
    if (!is.null(compare)) {
        if (!inherits(compare, "fSignRestr"))
            stop("`compare` must be an 'fSignRestr' object.")
        if (is.null(scale_compare)) scale_compare <- scale
        lab2 <- if (length(labels) >= 2) labels[2] else "Comparison"
        fits[[lab2]] <- fIRF(compare, conf = conf, varnames = varnames,
                             scale = scale_compare)
    }

    do.call(fPlotIRF, c(fits, list(
        shock = shock, return_data = return_data, colors = colors,
        ribbon_fill = colors[1], line_color = colors[1],
        ribbon_alpha = ribbon_alpha, facet_scales = facet_scales,
        facet_ncol = facet_ncol)))
}
