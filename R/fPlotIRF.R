# -----------------------------------------------------------------------
# fPlotIRF.R — one plot function for impulse responses
#
#   fPlotIRF(fit)                                  one result
#   fPlotIRF(Cholesky = irf_chol, IV = irf_iv)     several results, one colour each
#   fPlotIRF(lp_fit, return_data = TRUE)           the tidy data (fIRF format)
#
# Each argument is converted with fIRF(), so anything fIRF() understands can
# be plotted directly: fLP / fLPIV / fLPPanel / fSignRestr fits, IV-MBB,
# heteroskedasticity and MSW result lists, or fIRF() objects built from a
# point estimate and bootstrap bands.
#
# Author: Dr. Muhsin Ciftci
# -----------------------------------------------------------------------

#' Plot Impulse Responses
#'
#' Draws impulse responses with nested confidence bands, one panel per
#' variable. Several results can be passed at once; each gets its own
#' colour, which makes it easy to compare identification schemes, bootstrap
#' variants or specifications in one figure.
#'
#' Every argument in \code{...} is converted with \code{\link{fIRF}}, and the
#' \code{fIRF()} arguments below can be given directly, so plotting one
#' result never needs a nested call:
#' \code{fPlotIRF(point, bands = boot, varnames = vn)}. Comparing several
#' results is the one case that still needs \code{fIRF()}, because each
#' result carries its own bands.
#'
#' @param ... One or more results. Argument names become the legend labels,
#'   e.g. \code{fPlotIRF(Standard = a, Corrected = b)}.
#' @param bands Bootstrap result supplying the confidence bands for a point
#'   estimate, e.g. \code{fPlotIRF(point_irf, bands = boot_chol)}. Only
#'   meaningful with a single result: when several are compared, each one
#'   carries its own bands, so build them with \code{\link{fIRF}} first.
#' @param varnames,shocknames Labels for the variables and shocks, applied to
#'   every result given. They are passed to \code{\link{fIRF}} when a result
#'   still has to be converted, and rename a result that \code{fIRF()}
#'   already built, so one set of labels can serve several results.
#' @param cumulate,conf Passed to \code{\link{fIRF}} when a result still has
#'   to be converted. They change how the table is built, so they are refused
#'   for a result that \code{fIRF()} already built.
#' @param shock Shocks to plot, by position or name (within each result).
#'   \code{NULL} plots every shock. With one shock per result the panels are
#'   indexed by variable; with several, by shock and variable.
#' @param variables Variables to plot (names). \code{NULL} plots all.
#' @param scale Scalar or vector named by variable, multiplying the point
#'   estimate and bands.
#' @param colors Colours for several results (one per result). Defaults to
#'   the tidyMacro palette.
#' @param ribbon_fill,line_color Band and line colours for a single result.
#'   \code{ribbon_fill} may also give one colour per band, outer first or
#'   named by band, e.g.
#'   \code{c("Anderson-Rubin" = "red", "Delta method" = "blue")}.
#' @param ribbon_alpha Band transparency: one value per band, outer first. A
#'   single value sets the outer band; inner bands get more opaque. The
#'   default is lighter when several results overlap.
#' @param zero_line_color Colour of the dashed zero line.
#' @param facet_scales Passed to \code{ggplot2::facet_wrap()}.
#' @param facet_ncol Number of panel columns. Defaults to a square-ish grid
#'   for one shock and one column per variable for several shocks.
#' @param return_data If \code{TRUE}, return the \code{fIRF} tibble instead
#'   of the plot.
#'
#' @return A \code{ggplot} object, or an \code{fIRF} tibble when
#'   \code{return_data = TRUE}.
#'
#' @seealso \code{\link{fIRF}}
#'
#' @examples
#' set.seed(123)
#' y     <- matrix(rnorm(400), ncol = 2, dimnames = list(NULL, c("y1", "y2")))
#' fit   <- fVAR(y, p = 2, c = 1)
#' wold  <- fWoldIRF(fit, horizon = 12)
#' point <- fCholeskyIRF(wold, t(chol(fit$sigma)))
#' boot  <- fBootstrapChol(y, fit, nboot = 50, horizon = 12, n_threads = 2)
#'
#' # A point estimate and its bands, in one call
#' fPlotIRF(point, bands = boot, varnames = colnames(y), shock = 1)
#'
#' # Comparing results: each one carries its own bands, so build them first
#' corrected <- fBootstrapCholCorrected(y, fit, nboot1 = 50, nboot2 = 50,
#'                                      horizon = 12, n_threads = 2)
#' fPlotIRF(Standard    = fIRF(point, bands = boot),
#'          `Corrected` = fIRF(point, bands = corrected),
#'          varnames = colnames(y), shock = 1)
#' @export
fPlotIRF <- function(...,
                  bands           = NULL,
                  varnames        = NULL,
                  shocknames      = NULL,
                  cumulate        = NULL,
                  conf            = NULL,
                  shock           = NULL,
                  variables       = NULL,
                  scale           = 1,
                  colors          = NULL,
                  ribbon_fill     = "#407EC9",
                  line_color      = "#910048",
                  ribbon_alpha    = NULL,
                  zero_line_color = "#707372",
                  facet_scales    = "free",
                  facet_ncol      = NULL,
                  return_data     = FALSE) {
    objs <- list(...)
    if (!length(objs)) stop("fPlotIRF: supply at least one result.", call. = FALSE)
    labs <- names(objs)
    if (is.null(labs)) labs <- character(length(objs))
    .fPlotIRF_check_names(labs)

    # Arguments of fIRF() passed straight through, so a single result needs
    # no fIRF() call. `bands` belongs to one specific result, so it is only
    # meaningful when exactly one is supplied.
    conv <- list(varnames = varnames, shocknames = shocknames,
                 cumulate = cumulate, conf = conf)
    conv <- conv[!vapply(conv, is.null, logical(1))]
    # varnames and shocknames are labels, so they can also rename a result
    # that fIRF() already built; the others change how it is built.
    rebuild_only <- setdiff(names(conv), c("varnames", "shocknames"))
    if (!is.null(bands)) {
        if (length(objs) > 1L)
            stop("fPlotIRF: `bands` applies to one result only. With several ",
                 "results attach the bands to each one first, e.g.\n",
                 "  fPlotIRF(Standard = fIRF(point, bands = boot_a),\n",
                 "           Corrected = fIRF(point, bands = boot_b))",
                 call. = FALSE)
        conv$bands <- bands
    }

    d <- lapply(seq_along(objs), function(k) {
        x <- objs[[k]]
        if (inherits(x, "fIRF")) {
            if (length(rebuild_only))
                stop("fPlotIRF: ", paste(rebuild_only, collapse = ", "),
                     " cannot be applied to a result that fIRF() already ",
                     "built; pass them to fIRF() instead.", call. = FALSE)
            x <- .fPlotIRF_relabel(x, varnames, shocknames)
        } else {
            x <- do.call(fIRF, c(list(x), conv))
        }
        if (nzchar(labs[k])) {
            x$model <- labs[k]
        } else if (anyNA(x$model)) {
            x$model[is.na(x$model)] <-
                if (length(objs) == 1L) "" else paste("Model", k)
        }
        .fPlotIRF_pick_shock(x, shock)
    })
    d <- .fIRF_new(do.call(rbind, d))

    if (!is.null(variables)) {
        miss <- setdiff(variables, unique(d$variable))
        if (length(miss))
            stop(sprintf("fPlotIRF: variable(s) not found: %s",
                         paste(miss, collapse = ", ")), call. = FALSE)
        d <- d[d$variable %in% variables, , drop = FALSE]
    }
    d <- .fIRF_scale(d, scale)

    if (return_data) return(d)
    .fPlotIRF_draw(d, colors, ribbon_fill, line_color, ribbon_alpha,
                zero_line_color, facet_scales, facet_ncol)
}

# Keep the requested shocks (positions or names) of one result.
.fPlotIRF_pick_shock <- function(x, shock) {
    if (is.null(shock)) return(x)
    sh <- unique(x$shock)
    if (is.numeric(shock)) {
        if (any(shock < 1 | shock > length(sh) | shock != round(shock)))
            stop(sprintf("fPlotIRF: shock must be between 1 and %d.", length(sh)),
                 call. = FALSE)
        keep <- sh[shock]
    } else {
        miss <- setdiff(shock, sh)
        if (length(miss))
            stop(sprintf("fPlotIRF: shock(s) not found: %s (available: %s).",
                         paste(miss, collapse = ", "),
                         paste(sh, collapse = ", ")), call. = FALSE)
        keep <- shock
    }
    x <- x[x$shock %in% keep, , drop = FALSE]
    # Keep the requested order for the panels.
    x[order(match(x$shock, keep)), , drop = FALSE]
}

.fPlotIRF_palette <- c("#910048", "#407EC9", "#FF8200", "#00B0B9", "#001E60",
                    "#78BE20", "#8031A7", "#F2A900")

# Band levels, outer first: numeric levels ("90%") from widest to narrowest
# after labelled bands ("Anderson-Rubin", ...) in order of appearance.
.fPlotIRF_band_levels <- function(band) {
    lev <- unique(band[!is.na(band)])
    num <- suppressWarnings(as.numeric(sub("%$", "", lev)))
    c(lev[is.na(num)], lev[!is.na(num)][order(num[!is.na(num)], decreasing = TRUE)])
}

# Overlapping results get lighter bands so each stays visible.
.fPlotIRF_alphas <- function(alpha, k, multi = FALSE) {
    if (is.null(alpha)) {
        lo <- if (multi) 0.12 else 0.20
        hi <- if (multi) 0.22 else 0.35
        return(if (k == 1L) lo else seq(lo, hi, length.out = k))
    }
    if (length(alpha) == k) return(alpha)
    if (length(alpha) == 1L)
        return(alpha + (seq_len(k) - 1L) * min(alpha, (1 - alpha) / k))
    stop(sprintf("fPlotIRF: `ribbon_alpha` must have length 1 or %d.", k),
         call. = FALSE)
}

.fPlotIRF_draw <- function(d, colors, ribbon_fill, line_color, ribbon_alpha,
                        zero_line_color, facet_scales, facet_ncol) {
    d$model    <- factor(d$model,    levels = unique(d$model))
    d$variable <- factor(d$variable, levels = unique(d$variable))
    d$shock    <- factor(d$shock,    levels = unique(d$shock))

    models <- levels(d$model)
    multi  <- length(models) > 1L
    if (multi) {
        pal <- rep_len(if (is.null(colors)) .fPlotIRF_palette else colors,
                       length(models))
        fill_values <- line_values <- stats::setNames(pal, models)
    } else {
        fill_values <- stats::setNames(ribbon_fill, models)
        line_values <- stats::setNames(line_color, models)
    }

    band_levels <- .fPlotIRF_band_levels(d$band)
    bands <- d[!is.na(d$band), , drop = FALSE]

    # One result with one colour per band (e.g. Anderson-Rubin vs delta method)
    by_band <- !multi && length(ribbon_fill) > 1L && nrow(bands) > 0L
    if (by_band) {
        band_fill <- if (is.null(names(ribbon_fill))) {
            if (length(ribbon_fill) != length(band_levels))
                stop(sprintf("fPlotIRF: `ribbon_fill` needs one colour per band (%s).",
                             paste(band_levels, collapse = ", ")), call. = FALSE)
            stats::setNames(ribbon_fill, band_levels)
        } else {
            miss <- setdiff(band_levels, names(ribbon_fill))
            if (length(miss))
                stop(sprintf("fPlotIRF: `ribbon_fill` has no colour for: %s",
                             paste(miss, collapse = ", ")), call. = FALSE)
            ribbon_fill[band_levels]
        }
        fill_values <- band_fill
    }
    lines <- d[!duplicated(d[c("model", "shock", "variable", "horizon")]), ,
               drop = FALSE]

    p <- ggplot2::ggplot(mapping = ggplot2::aes(x = .data$horizon))
    if (nrow(bands)) {
        bands$band <- factor(bands$band, levels = band_levels)
        alphas <- stats::setNames(.fPlotIRF_alphas(ribbon_alpha, length(band_levels), multi),
                                  band_levels)
        band_guide <- if (length(band_levels) > 1L && !by_band)
            ggplot2::guide_legend(override.aes = list(
                fill = if (multi) "grey40" else ribbon_fill))
        else "none"
        p <- p +
            ggplot2::geom_ribbon(
                data = bands,
                ggplot2::aes(ymin = .data$lower, ymax = .data$upper,
                             fill = if (by_band) .data$band else .data$model,
                             alpha = .data$band,
                             group = interaction(.data$model, .data$band)),
                colour = NA) +
            ggplot2::scale_alpha_manual(values = alphas, name = NULL,
                                        guide = band_guide)
    }

    p <- p +
        ggplot2::geom_line(data = lines,
                           ggplot2::aes(y = .data$point, colour = .data$model),
                           linewidth = 0.8) +
        ggplot2::geom_hline(yintercept = 0, colour = zero_line_color,
                            linetype = "dashed", linewidth = 0.6) +
        ggplot2::scale_fill_manual(
            values = fill_values, name = NULL,
            guide = if (multi) "legend"
                    else if (by_band) ggplot2::guide_legend(override.aes = list(
                        alpha = unname(alphas[band_levels])))
                    else "none") +
        ggplot2::scale_colour_manual(values = line_values, name = NULL,
                                     guide = if (multi) "legend" else "none")

    # One shock per result: panels by variable. Several: shock x variable.
    one_shock <- all(tapply(as.character(d$shock), d$model,
                            function(s) length(unique(s))) == 1L)
    n_var <- nlevels(d$variable)
    if (one_shock) {
        if (is.null(facet_ncol)) facet_ncol <- ceiling(sqrt(n_var))
        p <- p + ggplot2::facet_wrap(~ variable, scales = facet_scales,
                                     ncol = facet_ncol)
    } else {
        if (is.null(facet_ncol)) facet_ncol <- n_var
        p <- p + ggplot2::facet_wrap(~ shock + variable, scales = facet_scales,
                                     ncol = facet_ncol)
    }
    p +
        ggplot2::scale_x_continuous(
            breaks = function(lim) unique(round(pretty(lim)))) +
        ggplot2::labs(x = "Horizon", y = NULL)
}

# Rename the variables or shocks of an existing fIRF table, in the order in
# which they appear. Used when one set of labels applies to several results.
.fPlotIRF_relabel <- function(x, varnames, shocknames) {
    swap <- function(col, new, what) {
        old <- unique(col)
        if (length(new) != length(old))
            stop(sprintf("fPlotIRF: `%s` must have length %d, to match the %s in the result.",
                         what, length(old), what), call. = FALSE)
        new[match(col, old)]
    }
    if (!is.null(varnames))
        x$variable <- swap(x$variable, as.character(varnames), "varnames")
    if (!is.null(shocknames))
        x$shock <- swap(x$shock, as.character(shocknames), "shocknames")
    x
}

# A misspelled argument name would otherwise be taken as a model label and
# silently drawn as an extra model, so reject anything close to a real one.
.fPlotIRF_check_names <- function(labs) {
    known <- c("bands", "varnames", "shocknames", "cumulate", "conf", "shock",
               "variables", "scale", "colors", "ribbon_fill", "line_color",
               "ribbon_alpha", "zero_line_color", "facet_scales",
               "facet_ncol", "return_data")
    for (nm in labs[nzchar(labs)]) {
        hit <- known[tolower(known) == tolower(nm)]
        if (!length(hit)) {
            d <- utils::adist(tolower(nm), tolower(known))[1, ]
            hit <- known[d <= 2 & d > 0]
        }
        if (length(hit))
            stop(sprintf(
                "fPlotIRF: '%s' is not an argument of fPlotIRF(); did you mean '%s'? %s",
                nm, hit[1],
                "Named arguments that are not options become model labels."),
                call. = FALSE)
    }
    invisible(TRUE)
}
