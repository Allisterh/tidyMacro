# -----------------------------------------------------------------------
# fIRF.R — one data format for every impulse response in tidyMacro
#
#   fIRF(point, bands = boot)          Cholesky, BQ, LR-max, Uhlig (VAR + bootstrap)
#   fIRF(result)                       IV-MBB, heteroskedasticity, MSW lists
#   fIRF(fit)                          fLP, fLPIV, fLPPanel, fSignRestr objects
#
# Every method returns the same long tibble (class "fIRF"):
#
#   model | shock | variable | horizon | band | point | lower | upper
#
# one row per (model, shock, variable, horizon, band). `band` labels a
# confidence/credible band ("90%", "68%", "Anderson-Rubin", ...); a result
# without bands has a single row with band = NA. fPlotIRF() draws this format,
# so a new identification scheme only needs an fIRF() method.
#
# Author: Dr. Muhsin Ciftci
# -----------------------------------------------------------------------

#' Impulse Responses in a Common Tidy Format
#'
#' Converts the output of any identification scheme into one long tibble
#' with the same columns, so that results can be plotted with
#' \code{\link{fPlotIRF}}, combined with \code{rbind()}, or used in a
#' hand-made \pkg{ggplot2} figure.
#'
#' Supported inputs:
#' \itemize{
#'   \item a point estimate (\code{N x N x H} array from
#'     \code{\link{fCholeskyIRF}} or \code{\link{fBQIRF}}, or \code{N x H}
#'     matrix for a single shock) with optional \code{bands} from
#'     \code{\link{fBootstrapChol}}, \code{\link{fBootstrapCholCorrected}},
#'     \code{\link{fBootstrapBQ}}, \code{\link{fBootstrapMax}} or
#'     \code{\link{fBootstrapUhlig}};
#'   \item the lists returned by \code{\link{fBootstrapIVMBB}} and
#'     \code{\link{fBootstrapHetero}} (both carry the sample estimate in
#'     \code{point}), and \code{\link{fMSW}} (Anderson-Rubin and
#'     delta-method bands);
#'   \item \code{fLP}, \code{fLPIV} and \code{fLPPanel} fits;
#'   \item \code{fSignRestr} fits (point = posterior median).
#' }
#'
#' @param x Point estimate, result list or fitted object (see Details).
#' @param bands Bootstrap result with \code{upper}, \code{lower} and
#'   optionally \code{upper2}, \code{lower2}, in the same dimensions as
#'   \code{x}.
#' @param varnames Character vector of variable names (length N).
#' @param shocknames Character vector of shock names (one per shock).
#' @param cumulate Variables (indices or names) whose \emph{point} estimate
#'   is cumulated over horizons. Bootstrap bands are taken as given:
#'   \code{fBootstrapBQ(cumulate = ...)} already returns cumulated bands.
#' @param scale Scalar, or one value per variable, multiplying the point
#'   estimate and the bands. A negative value flips the sign and keeps
#'   \code{lower <= upper}.
#' @param conf Band levels in percent, outer band first. By default taken
#'   from the result (\code{conf} element) or \code{c(90, 68)}. For
#'   \code{fSignRestr}, selects which stored levels to use.
#' @param variables For LP fits, the response variables to keep.
#' @param model Label for the \code{model} column; \code{fPlotIRF()} fills it
#'   from argument names when missing.
#' @param ... Passed to methods.
#'
#' @return A tibble of class \code{"fIRF"} with columns \code{model},
#'   \code{shock}, \code{variable}, \code{horizon}, \code{band},
#'   \code{point}, \code{lower}, \code{upper}.
#'
#' @seealso \code{\link{fPlotIRF}}
#'
#' @examples
#' set.seed(123)
#' y     <- matrix(rnorm(400), ncol = 2, dimnames = list(NULL, c("y1", "y2")))
#' fit   <- fVAR(y, p = 2, c = 1)
#' wold  <- fWoldIRF(fit, horizon = 12)
#' point <- fCholeskyIRF(wold, t(chol(fit$sigma)))
#' boot  <- fBootstrapChol(y, fit, nboot = 50, horizon = 12, n_threads = 2)
#'
#' irf <- fIRF(point, bands = boot, varnames = colnames(y))
#' head(irf)
#' fPlotIRF(irf, shock = 1)
#' @export
fIRF <- function(x, ...) UseMethod("fIRF")

#' @rdname fIRF
#' @export
fIRF.default <- function(x, bands = NULL, varnames = NULL, shocknames = NULL,
                         cumulate = NULL, scale = 1, conf = NULL,
                         model = NA_character_, ...) {
    if (!is.numeric(x))
        stop(sprintf("fIRF: no method for an object of class '%s'.",
                     class(x)[1]), call. = FALSE)
    point <- .fIRF_cube(x, "x")
    N <- dim(point)[1]
    if (is.null(varnames)) varnames <- rownames(x)

    if (length(cumulate)) {
        idx <- .fIRF_index(cumulate, .fIRF_names(varnames, N, "varnames", "y"),
                           "cumulate")
        for (i in idx) for (s in seq_len(dim(point)[2]))
            point[i, s, ] <- cumsum(point[i, s, ])
    }

    d <- .fIRF_build(point, .fIRF_boot_bands(bands, conf), varnames,
                     shocknames, model = model)
    .fIRF_scale(d, scale)
}

#' @rdname fIRF
#' @export
fIRF.list <- function(x, varnames = NULL, shocknames = NULL, scale = 1,
                      conf = NULL, model = NA_character_, ...) {
    nm <- names(x)
    if (all(c("IRF", "MSWlbound", "MSWubound", "Dmlbound", "Dmubound") %in% nm)) {
        # MSW: weak-IV robust (Anderson-Rubin) set outside the delta method
        bands <- list(
            "Anderson-Rubin" = list(lower = x$MSWlbound, upper = x$MSWubound),
            "Delta method"   = list(lower = x$Dmlbound,  upper = x$Dmubound))
        point <- x$IRF
    } else if (!is.null(x$point) && !is.null(x$upper)) {
        # fBootstrapIVMBB() and fBootstrapHetero(): `point` is the estimate on
        # the original sample, never a summary of the bootstrap draws.
        point <- x$point
        bands <- .fIRF_boot_bands(x, conf)
    } else if (any(c("bootchol", "bootbq", "bootmax", "bootuhlig") %in% nm)) {
        stop("fIRF: bootstrap results contain bands only. Pass the point ",
             "estimate too: fIRF(point, bands = <bootstrap result>).",
             call. = FALSE)
    } else {
        stop("fIRF: unrecognised result list. Expected the output of ",
             "fBootstrapIVMBB(), fBootstrapHetero() or fMSW().", call. = FALSE)
    }
    d <- .fIRF_build(.fIRF_cube(point, "point"), bands, varnames, shocknames,
                     model = model)
    .fIRF_scale(d, scale)
}

#' @rdname fIRF
#' @export
fIRF.fLP <- function(x, variables = NULL, scale = 1, model = NA_character_,
                     ...) {
    irf  <- x$irfs                                   # H x K
    cols <- colnames(irf)
    if (is.null(cols)) cols <- x$lhs_vars
    if (is.null(cols)) cols <- paste0("response_", seq_len(ncol(irf)))
    H <- nrow(irf)
    K <- ncol(irf)
    horizons <- if (is.null(x$horizons)) seq_len(H) - 1L else x$horizons

    # fLPPanel columns are shock terms of one outcome; otherwise columns are
    # outcomes responding to one shock.
    panel <- inherits(x, "fLPPanel")
    to_cube <- function(m) {
        a <- array(m, c(H, K, 1L))
        if (panel) aperm(a, c(3L, 2L, 1L)) else aperm(a, c(2L, 3L, 1L))
    }
    if (panel) {
        varnames   <- x$lhs_vars[1]
        shocknames <- cols
    } else {
        varnames   <- cols
        shocknames <- x$shock
    }

    keys <- if (is.list(x$irfs_upper)) names(x$irfs_upper) else format(x$conf)
    bands <- lapply(keys, function(k) {
        up <- if (is.list(x$irfs_upper)) x$irfs_upper[[k]] else x$irfs_upper
        lo <- if (is.list(x$irfs_lower)) x$irfs_lower[[k]] else x$irfs_lower
        list(lower = to_cube(lo), upper = to_cube(up))
    })
    names(bands) <- .fIRF_label(keys)

    d <- .fIRF_build(to_cube(irf), bands, varnames, shocknames,
                     horizons = horizons, model = model)
    if (!is.null(variables)) {
        miss <- setdiff(variables, unique(d$variable))
        if (length(miss))
            stop(sprintf("fIRF: variable(s) not found in LP result: %s",
                         paste(miss, collapse = ", ")), call. = FALSE)
        d <- d[d$variable %in% variables, , drop = FALSE]
    }
    .fIRF_scale(d, scale)
}

#' @rdname fIRF
#' @export
fIRF.fSignRestr <- function(x, conf = NULL, varnames = NULL, shocknames = NULL,
                            scale = 1, model = NA_character_, ...) {
    levels <- sort(as.numeric(if (is.null(conf)) x$conf else conf),
                   decreasing = TRUE)
    # A level the model was not asked for fails loudly instead of being
    # drawn silently at the wrong coverage.
    bands <- lapply(levels, function(level) {
        b <- x$bands[[format(level)]]
        if (is.null(b))
            stop(sprintf("Band %s%% is not available; refit with conf = c(%s).",
                         format(level),
                         paste(unique(c(x$conf, level)), collapse = ", ")),
                 call. = FALSE)
        list(lower = b$IRinf, upper = b$IRsup)
    })
    names(bands) <- .fIRF_label(levels)
    if (is.null(varnames)) varnames <- x$varnames

    d <- .fIRF_build(.fIRF_cube(x$IRmed, "IRmed"), bands, varnames,
                     shocknames, model = model)
    .fIRF_scale(d, scale)
}

#' @rdname fIRF
#' @export
fIRF.fIRF <- function(x, scale = 1, model = NULL, ...) {
    if (!is.null(model)) x$model <- model
    .fIRF_scale(x, scale)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# N x S x H array from an N x S x H array, an N x H matrix (one shock) or a
# length-H vector (one variable, one shock).
.fIRF_cube <- function(a, what) {
    if (is.array(a) && length(dim(a)) == 3L) return(a)
    if (is.matrix(a)) return(array(a, c(nrow(a), 1L, ncol(a))))
    if (is.numeric(a) && is.null(dim(a))) return(array(a, c(1L, 1L, length(a))))
    stop(sprintf("fIRF: `%s` must be an N x H matrix or an N x S x H array.",
                 what), call. = FALSE)
}

.fIRF_names <- function(nm, n, arg, prefix) {
    if (is.null(nm))
        return(if (n == 1L && prefix == "Shock") prefix
               else paste0(prefix, if (prefix == "Shock") " " else "",
                           seq_len(n)))
    nm <- as.character(nm)
    if (length(nm) != n)
        stop(sprintf("fIRF: `%s` must have length %d.", arg, n), call. = FALSE)
    nm
}

# Positions of `which` (indices or names) in `names`.
.fIRF_index <- function(which, names, arg) {
    idx <- if (is.character(which)) match(which, names) else as.integer(which)
    if (anyNA(idx) || any(idx < 1L | idx > length(names)))
        stop(sprintf("fIRF: `%s` must be indices between 1 and %d or variable names.",
                     arg, length(names)), call. = FALSE)
    idx
}

.fIRF_label <- function(level) paste0(format(as.numeric(level), trim = TRUE), "%")

# Named list of bands (outer first) from a bootstrap result.
.fIRF_boot_bands <- function(b, conf = NULL) {
    if (is.null(b)) return(list())
    if (is.null(b$upper) || is.null(b$lower))
        stop("fIRF: `bands` must contain `upper` and `lower`.", call. = FALSE)
    lev <- if (!is.null(conf)) conf else if (!is.null(b$conf)) b$conf else c(90, 68)
    inner <- !is.null(b$upper2) && !is.null(b$lower2)
    if (inner && length(lev) < 2L)
        stop("fIRF: two bands found; `conf` needs two levels.", call. = FALSE)
    out <- list(list(lower = b$lower, upper = b$upper))
    if (inner) out[[2L]] <- list(lower = b$lower2, upper = b$upper2)
    names(out) <- .fIRF_label(lev[seq_along(out)])
    out
}

# Long tibble from an N x S x H point cube and a named list of band cubes.
.fIRF_build <- function(point, bands, varnames, shocknames, horizons = NULL,
                        model = NA_character_) {
    d <- dim(point)
    N <- d[1]; S <- d[2]; H <- d[3]
    varnames   <- .fIRF_names(varnames, N, "varnames", "y")
    shocknames <- .fIRF_names(shocknames, S, "shocknames", "Shock")
    if (is.null(horizons)) horizons <- seq_len(H) - 1L
    if (length(horizons) != H)
        stop("fIRF: `horizons` does not match the number of horizons.", call. = FALSE)

    # Horizon varies fastest, then variable, then shock.
    flat <- function(a) as.vector(aperm(a, c(3L, 1L, 2L)))
    base <- tibble::tibble(
        model    = as.character(model),
        shock    = rep(shocknames, each = N * H),
        variable = rep(rep(varnames, each = H), times = S),
        horizon  = rep(as.integer(horizons), times = N * S),
        band     = NA_character_,
        point    = flat(point),
        lower    = NA_real_,
        upper    = NA_real_)

    if (length(bands)) {
        base <- do.call(rbind, lapply(names(bands), function(lab) {
            lo <- .fIRF_cube(bands[[lab]]$lower, "lower")
            up <- .fIRF_cube(bands[[lab]]$upper, "upper")
            if (!identical(dim(lo), d) || !identical(dim(up), d))
                stop(sprintf("fIRF: band '%s' does not match the point estimate (%s).",
                             lab, paste(d, collapse = " x ")), call. = FALSE)
            b <- base
            b$band  <- lab
            b$lower <- flat(lo)
            b$upper <- flat(up)
            b
        }))
    }
    .fIRF_new(base)
}

.fIRF_new <- function(d) {
    d <- tibble::as_tibble(d)[, c("model", "shock", "variable", "horizon",
                                  "band", "point", "lower", "upper")]
    tibble::new_tibble(d, class = "fIRF")
}

# Multiply point and bands by `scale`: a scalar, a vector named by variable,
# or one value per variable in order of appearance.
.fIRF_scale <- function(d, scale) {
    if (!is.numeric(scale) || !length(scale) || any(!is.finite(scale)))
        stop("fIRF: `scale` must be finite numeric.", call. = FALSE)
    if (length(scale) == 1L && scale == 1) return(d)
    vars <- unique(d$variable)
    sc <- if (length(scale) == 1L) {
        rep(scale, nrow(d))
    } else if (!is.null(names(scale))) {
        miss <- setdiff(vars, names(scale))
        if (length(miss))
            stop(sprintf("fIRF: `scale` has no value for: %s",
                         paste(miss, collapse = ", ")), call. = FALSE)
        unname(scale[d$variable])
    } else {
        if (length(scale) != length(vars))
            stop(sprintf("fIRF: `scale` must be a scalar or have length %d.",
                         length(vars)), call. = FALSE)
        scale[match(d$variable, vars)]
    }
    lo <- d$lower * sc
    up <- d$upper * sc
    d$point <- d$point * sc
    d$lower <- pmin(lo, up)
    d$upper <- pmax(lo, up)
    d
}
