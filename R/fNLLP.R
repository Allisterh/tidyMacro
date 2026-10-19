# -----------------------------------------------------------------------
# fNLLP.R — Nonlinear and state-dependent Local Projections
#
# Provides:
#   fNLLP()          — main user-facing function
#   print.fNLLP()    — print method
#   coef.fNLLP()     — extract the two shock coefficients
#   tidy.fNLLP()     — long data frame of components (and their difference)
#   predict.fNLLP()  — response to a shock of a given size (and state value)
#
# The formula is read exactly as in fLP(): the RHS is what is estimated, a
# constant is added internally, and l()/f()/..macro expand the same way.
# `type` only says how the shock (and, for "state", every non-common term)
# enters the regression:
#
#   sign  : shock -> max(shock, 0) + min(shock, 0)
#   cubic : shock -> shock + shock^3
#   state : every term x -> state * x + (1 - state) * x, the intercept is
#           split into state and (1 - state); terms named in `common` keep a
#           single coefficient.
#
# Author: Dr. Muhsin Ciftci
# -----------------------------------------------------------------------


# =======================================================================
# Internal helpers
# =======================================================================

# Fixed regressor sample (balanced = TRUE). The regressors decide the
# estimation dates t; the outcome is read at t + h from later rows of data,
# so `data` must run H rows past the last regressor date. Mirrors the
# startS:endS / startS+h:endS+h alignment of stlpm.m.
.fNLLP_balanced_sample <- function(data, x_vars, lhs_vars, H, cumulative) {
  keep <- which(complete.cases(data[, x_vars, drop = FALSE]))
  if (length(keep) == 0)
    stop("fNLLP: no complete rows for the regressors.")

  gaps <- which(diff(keep) != 1)
  if (length(gaps) > 0) {
    stop(sprintf(paste0(
      "fNLLP: the regressors have %d interior gap(s) in the time index ",
      "(after row(s) %s of `data`). Subset `data` to a contiguous span."),
      length(gaps), paste(keep[gaps][seq_len(min(5, length(gaps)))],
                          collapse = ", ")))
  }

  first <- keep[1]
  last  <- keep[length(keep)]

  # The outcome must be observed from `first` through last + H. Where it
  # stops earlier (fully observed data, or an outcome that ends before the
  # regressors do), the last regressor dates are dropped so that y_{t+H}
  # exists for every date kept. A gap inside the outcome is an error.
  y_ok  <- complete.cases(data[, lhs_vars, drop = FALSE])
  if (!y_ok[first])
    stop(sprintf("fNLLP: the outcome is missing at the first regressor date (row %d).",
                 first))
  y_end <- first - 1 + match(FALSE, y_ok[first:nrow(data)], nomatch = nrow(data) - first + 2) - 1
  if (y_end < nrow(data) && any(y_ok[(y_end + 1):nrow(data)]) &&
      y_end < last + H) {
    stop(sprintf(paste0(
      "fNLLP: balanced = TRUE needs the outcome observed without gaps from ",
      "row %d onwards; it is missing at row %d and observed again later."),
      first, y_end + 1))
  }
  last_eff <- min(last, y_end - H)
  if (last_eff - first + 1 < 1)
    stop(sprintf(paste0(
      "fNLLP: balanced = TRUE reads the outcome %d row(s) after each ",
      "regressor date, but the outcome ends at row %d."), H, y_end))
  if (last_eff < last) {
    message(sprintf(paste0(
      "fNLLP: dropping regressor dates after row %d so the outcome at ",
      "horizon %d is observed."), last_eff, H))
    keep <- keep[keep <= last_eff]
    last <- last_eff
  }

  y_rows <- first:(last + H)
  Ydf    <- data[y_rows, lhs_vars, drop = FALSE]

  pre <- NULL
  if (isTRUE(cumulative) && first > 1) {
    cand <- data[first - 1, lhs_vars, drop = FALSE]
    if (all(complete.cases(cand))) pre <- cand
  }

  message(sprintf(
    "fNLLP: regressors fixed at rows %d-%d of `data` (%d observations); outcomes read up to row %d.",
    first, last, length(keep), last + H))

  list(X = data[keep, x_vars, drop = FALSE], Y = Ydf, pre = pre)
}


# Column names of the regression design, in the order fNLLPDesign_cpp
# builds it.
.fNLLP_design_names <- function(type, rhs_vars, shock, state, common) {
  if (type %in% c("sign", "cubic")) {
    tr <- if (type == "sign") {
      c(sprintf("max(%s, 0)", shock), sprintf("min(%s, 0)", shock))
    } else {
      c(shock, sprintf("%s^3", shock))
    }
    out <- "(Intercept)"
    for (v in rhs_vars) out <- c(out, if (v == shock) tr else v)
    return(out)
  }
  nc <- c("(Intercept)", rhs_vars[!rhs_vars %in% common])
  c(paste0(state, ":", nc),
    paste0("(1-", state, "):", nc),
    rhs_vars[rhs_vars %in% common])
}


# =======================================================================
# Main function
# =======================================================================

#' Nonlinear and State-Dependent Local Projections
#'
#' Estimates local projections in which the shock has two coefficients per
#' horizon: one per regime (smooth state dependence), one per sign
#' (asymmetry), or a linear and a cubic term (size nonlinearity). The
#' formula syntax is that of \code{\link{fLP}}: the RHS is estimated as
#' written, a constant is added internally, and \code{l()}, \code{f()} and
#' \code{..macro} expand the same way. Both coefficients come from one
#' regression, together with their Newey-West covariance, so their
#' difference has a valid standard error.
#'
#' @inheritParams fLP
#' @param shock Character. \strong{Required.} Name of the shock variable; it
#'   must appear on the RHS.
#' @param type How the shock enters the regression:
#'   \describe{
#'     \item{\code{"state"}}{(default) every RHS term \eqn{x} and the
#'       intercept enter as \eqn{F x} and \eqn{(1 - F) x}, with \eqn{F} the
#'       column named by \code{state}. Components are named by
#'       \code{regimes}.}
#'     \item{\code{"sign"}}{the shock enters as \eqn{\max(s, 0)} and
#'       \eqn{\min(s, 0)}. Components \code{"positive"} and
#'       \code{"negative"}, each per unit of the shock.}
#'     \item{\code{"cubic"}}{the shock enters as \eqn{s} and \eqn{s^3}.
#'       Components \code{"linear"} and \code{"cubic"}; use
#'       \code{predict()} for the response to a shock of a given size.}
#'   }
#' @param state Character. For \code{type = "state"}: name of the column of
#'   \code{data} holding the regime weight \eqn{F}, in \eqn{[0, 1]}. It must
#'   not also appear on the RHS.
#' @param common Character vector of RHS terms that keep a single coefficient
#'   across regimes (e.g. a trend). Only for \code{type = "state"}; the shock
#'   cannot be common.
#' @param regimes Two distinct labels for the components of
#'   \code{type = "state"}: the first for weight \eqn{F}, the second for
#'   \eqn{1 - F}. Default \code{c("high", "low")}.
#' @param balanced Logical. If \code{TRUE}, the regressor dates are the same
#'   at every horizon and the outcome is read \eqn{h} rows later, as in the
#'   Tenreyro-Thwaites (2016) code; \code{data} must then run \code{horizons}
#'   rows past the last regressor date, and the regressors must have no
#'   interior gaps. Default \code{FALSE}: the sample shrinks with the horizon,
#'   as in \code{fLP()}.
#' @param store_full Logical. If \code{TRUE}, also return the full
#'   coefficient and standard-error matrices of every horizon (\code{betas},
#'   \code{ses}). Default \code{FALSE}.
#'
#' @return An object of class \code{"fNLLP"}, a list with
#'   \describe{
#'     \item{\code{irfs}, \code{irfs_se}}{Arrays \code{(H+1) x n_y x 2}
#'       [horizon, outcome, component] of the two shock coefficients and
#'       their Newey-West standard errors.}
#'     \item{\code{irfs_upper}, \code{irfs_lower}}{Normal bands of the same
#'       shape; named lists of arrays (one per level) when \code{conf} has
#'       several levels.}
#'     \item{\code{irfs_cov}}{\code{(H+1) x n_y} covariance of the two
#'       coefficients.}
#'     \item{\code{diff}, \code{diff_se}}{\code{(H+1) x n_y} first minus
#'       second component and its standard error; \code{NaN} for
#'       \code{type = "cubic"}.}
#'     \item{\code{nobs}}{Observations per horizon.}
#'     \item{\code{betas}, \code{ses}}{With \code{store_full = TRUE}: one
#'       matrix per horizon, rows named by \code{design_vars}.}
#'   }
#'   plus the settings (\code{type}, \code{components}, \code{design_vars},
#'   \code{shock}, \code{state}, \code{common}, \code{horizons}, \code{conf},
#'   \code{nw_lags}, \code{nw_offset}, \code{nw_prewhite}, \code{balanced},
#'   \code{cumulative}, \code{formula}, ...). \code{print}, \code{coef},
#'   \code{predict} and \code{\link{tidy.fNLLP}} methods are available.
#'
#' @references
#' Tenreyro, S., & Thwaites, G. (2016). Pushing on a string: US monetary
#' policy is less powerful in recessions. \emph{American Economic Journal:
#' Macroeconomics}, 8(4), 43--74.
#'
#' @seealso \code{\link{fLP}}, \code{\link{fPlotNLLP}},
#'   \code{\link{tidy.fNLLP}}, \code{\link{predict.fNLLP}},
#'   \code{\link{TenreyroThwaites2016}}
#'
#' @examples
#' data("TenreyroThwaites2016")
#' tt <- TenreyroThwaites2016
#' tt$GDP   <- 100 * tt$logGDP
#' tt$trend <- seq_len(nrow(tt))
#' # Expansion weight: logistic transform of standardised trailing
#' # four-quarter GDP growth
#' z <- as.numeric(stats::filter(c(NA, diff(tt$GDP)), rep(1 / 4, 4), sides = 1))
#' tt$Expansion <- plogis(3 * (z - mean(z, na.rm = TRUE)) / sd(z, na.rm = TRUE))
#' tt <- tt[tt$Date >= as.Date("1968-10-01"), ]
#'
#' fit <- fNLLP(GDP ~ RRShock + l(GDP, 1) + l(FFR, 1) + trend,
#'              data = tt, horizons = 12, shock = "RRShock",
#'              type = "state", state = "Expansion", common = "trend",
#'              regimes = c("Expansion", "Recession"), balanced = TRUE,
#'              n_threads = 2)
#' fit
#'
#' # Asymmetry between contractionary and expansionary shocks
#' fit_sign <- fNLLP(GDP ~ RRShock + l(GDP, 1) + l(FFR, 1) + trend,
#'                   data = tt, horizons = 12, shock = "RRShock",
#'                   type = "sign", balanced = TRUE, n_threads = 2)
#' fit_sign
#'
#' @export
fNLLP <- function(formula, data, horizons = 12,
                  shock,
                  type       = c("state", "sign", "cubic"),
                  state      = NULL,
                  common     = NULL,
                  regimes    = c("high", "low"),
                  conf       = 90,
                  nw_lags    = NULL,
                  nw_offset  = 1,
                  nw_prewhite = FALSE,
                  balanced   = FALSE,
                  cumulative = FALSE,
                  store_full = FALSE,
                  n_threads  = 0) {

  env  <- parent.frame()
  type <- match.arg(type)

  # =====================================================================
  # 0. Arguments
  # =====================================================================
  if (missing(shock)) {
    stop("fNLLP: 'shock' is required. ",
         "Pass the name of the impulse variable, e.g. shock = \"mp_shock\".")
  }
  if (!is.character(shock) || length(shock) != 1 || is.na(shock))
    stop("fNLLP: 'shock' must be a single character string.")

  H    <- .fLP_validate_horizon_max(horizons)
  conf <- .fLP_validate_conf_vector(conf)
  if (!is.data.frame(data)) data <- as.data.frame(data)

  for (nm in c("balanced", "cumulative", "store_full", "nw_prewhite")) {
    v <- get(nm)
    if (!is.logical(v) || length(v) != 1 || is.na(v))
      stop(sprintf("fNLLP: '%s' must be TRUE or FALSE.", nm))
  }

  if (type == "state") {
    if (is.null(state) || !is.character(state) || length(state) != 1 ||
        is.na(state)) {
      stop("fNLLP: type = \"state\" needs 'state', the name of a column of ",
           "`data` holding the regime weight in [0, 1].")
    }
    if (!state %in% names(data))
      stop(sprintf("fNLLP: state variable '%s' not found in data.", state))
    if (!is.character(regimes) || length(regimes) != 2 || anyNA(regimes) ||
        regimes[1] == regimes[2]) {
      stop("fNLLP: 'regimes' must be two distinct labels, ",
           "for the weight state and 1 - state, e.g. c(\"expansion\", \"recession\").")
    }
  } else {
    if (!is.null(state))
      stop("fNLLP: 'state' is only used with type = \"state\".")
    if (!is.null(common))
      stop("fNLLP: 'common' is only used with type = \"state\".")
  }

  # =====================================================================
  # 1. Formula -> LHS / RHS
  # =====================================================================
  parsed           <- .fLP_parse(formula, data, env, "fNLLP")
  lhs_vars         <- parsed$lhs_vars
  rhs_vars         <- parsed$rhs_vars
  data             <- parsed$data
  formula_expanded <- parsed$formula_expanded

  if (!shock %in% rhs_vars) {
    stop(sprintf(
      "fNLLP: shock variable '%s' not found on the RHS.\n  RHS terms: %s",
      shock, paste(rhs_vars, collapse = ", ")))
  }
  shock_col <- which(rhs_vars == shock) - 1

  if (type == "state" && state %in% rhs_vars)
    stop("fNLLP: the state variable cannot also be an RHS term; ",
         "the intercept is already split into state and 1 - state.")

  common_cols <- integer(0)
  if (!is.null(common)) {
    if (!is.character(common) || anyNA(common))
      stop("fNLLP: 'common' must be a character vector of RHS terms.")
    bad <- setdiff(common, rhs_vars)
    if (length(bad) > 0)
      stop(sprintf("fNLLP: 'common' term(s) not on the RHS: %s\n  RHS terms: %s",
                   paste(bad, collapse = ", "), paste(rhs_vars, collapse = ", ")))
    if (shock %in% common)
      stop("fNLLP: the shock cannot be a common term.")
    common      <- unique(common)
    common_cols <- match(common, rhs_vars) - 1
  }

  state_var  <- if (type == "state") state else character(0)
  all_needed <- unique(c(lhs_vars, rhs_vars, state_var))
  missing_v  <- setdiff(all_needed, names(data))
  if (length(missing_v) > 0)
    stop(sprintf("fNLLP: the following variables are not in data: %s",
                 paste(missing_v, collapse = ", ")))

  # =====================================================================
  # 2. Estimation sample
  # =====================================================================
  if (isTRUE(balanced)) {
    samp <- .fNLLP_balanced_sample(data, unique(c(rhs_vars, state_var)),
                                   lhs_vars, H, cumulative)
    Xdf  <- samp$X
    Ydf  <- samp$Y
  } else {
    samp <- .fLP_estimation_sample(data, all_needed, lhs_vars,
                                   cumulative = cumulative, fn = "fNLLP")
    Xdf  <- samp$raw
    Ydf  <- samp$raw[, lhs_vars, drop = FALSE]
  }

  non_numeric <- c(names(Xdf)[!vapply(Xdf, is.numeric, logical(1))],
                   names(Ydf)[!vapply(Ydf, is.numeric, logical(1))])
  if (length(non_numeric) > 0)
    stop(sprintf(
      "fNLLP: all model variables must be numeric. Non-numeric variable(s): %s",
      paste(unique(non_numeric), collapse = ", ")))

  Y <- as.matrix(Ydf[, lhs_vars, drop = FALSE])
  X <- as.matrix(Xdf[, rhs_vars, drop = FALSE])
  storage.mode(Y) <- "double"
  storage.mode(X) <- "double"

  state_v <- numeric(0)
  if (type == "state") {
    state_v <- as.double(Xdf[[state]])
    if (any(state_v < 0 | state_v > 1))
      stop(sprintf("fNLLP: state variable '%s' must lie in [0, 1].", state))
  }

  Y_pre <- NULL
  if (!is.null(samp$pre)) {
    Y_pre <- as.matrix(samp$pre[, lhs_vars, drop = FALSE])
    storage.mode(Y_pre) <- "double"
  }

  T_eff <- nrow(X)

  # =====================================================================
  # 3. Newey-West bandwidth and sample-size check
  # =====================================================================
  if (is.null(nw_lags)) {
    nw_lags <- as.integer(floor(T_eff^(1.0 / 3.0)))
  } else {
    nw_lags <- .fLP_validate_scalar_int(nw_lags, "nw_lags")
  }
  nw_offset <- .fLP_validate_scalar_int_signed(nw_offset, "nw_offset")
  n_threads <- .fLP_validate_scalar_int(n_threads, "n_threads")

  design_vars <- .fNLLP_design_names(type, rhs_vars, shock, state, common)
  min_obs <- length(design_vars) + 1 +
             if (isTRUE(balanced)) 0 else H
  min_obs <- min_obs + as.integer(isTRUE(cumulative) && is.null(Y_pre))
  if (T_eff < min_obs) {
    stop(sprintf(
      "fNLLP: not enough observations. The design has %d columns and needs at least %d rows; have %d.",
      length(design_vars), min_obs, T_eff))
  }

  # =====================================================================
  # 4. C++ engine
  # =====================================================================
  res <- fNLLP_cpp(
    Y             = Y,
    X             = X,
    H             = as.integer(H),
    shock_col     = as.integer(shock_col),
    specification = as.integer(switch(type, sign = 0, cubic = 1, state = 2)),
    state         = state_v,
    common_cols   = as.integer(common_cols),
    nw_lags_base  = as.integer(nw_lags),
    store_full    = isTRUE(store_full),
    cumulative    = isTRUE(cumulative),
    balanced      = isTRUE(balanced),
    n_threads     = as.integer(n_threads),
    nw_offset     = as.integer(nw_offset),
    verbose       = FALSE,
    Y_pre         = Y_pre,
    nw_prewhite   = nw_prewhite
  )

  # =====================================================================
  # 5. Bands (rebuilt in R from irfs_se, as in fLP) and labels
  # =====================================================================
  components <- switch(type,
                       state = regimes,
                       sign  = c("positive", "negative"),
                       cubic = c("linear", "cubic"))
  h_seq <- 0:H

  dimnames(res$irfs)    <- list(as.character(h_seq), lhs_vars, components)
  dimnames(res$irfs_se) <- dimnames(res$irfs)
  for (nm in c("irfs_cov", "diff", "diff_se"))
    dimnames(res[[nm]]) <- list(as.character(h_seq), lhs_vars)

  bands <- .fLP_bands(res$irfs, res$irfs_se, conf)
  res$irfs_upper <- bands$upper
  res$irfs_lower <- bands$lower

  res$nobs <- stats::setNames(as.integer(res$nobs), as.character(h_seq))

  if (isTRUE(store_full)) {
    for (i in seq_along(res$betas)) {
      dimnames(res$betas[[i]]) <- list(design_vars, lhs_vars)
      dimnames(res$ses[[i]])   <- list(design_vars, lhs_vars)
    }
    names(res$betas) <- paste0("h", h_seq)
    names(res$ses)   <- paste0("h", h_seq)
  }

  # =====================================================================
  # 6. Annotate and return
  # =====================================================================
  res$type         <- type
  res$components   <- components
  res$lhs_vars     <- lhs_vars
  res$rhs_vars     <- rhs_vars
  res$design_vars  <- design_vars
  res$shock        <- shock
  res$state        <- if (type == "state") state else NULL
  res$common       <- common
  res$horizons     <- h_seq
  res$conf         <- conf
  res$nw_lags      <- nw_lags
  res$nw_offset    <- nw_offset
  res$nw_prewhite  <- nw_prewhite
  res$balanced     <- isTRUE(balanced)
  res$cumulative   <- isTRUE(cumulative)
  res$store_full   <- isTRUE(store_full)
  res$n_threads    <- as.integer(n_threads)
  res$formula      <- formula_expanded
  res$formula_orig <- formula
  res$call         <- match.call()

  class(res) <- "fNLLP"
  res
}


# =======================================================================
# S3 methods
# =======================================================================

#' @export
print.fNLLP <- function(x, digits = 4, ...) {
  type_lbl <- switch(x$type,
                     state = "smooth state dependence",
                     sign  = "sign asymmetry (positive / negative)",
                     cubic = "size nonlinearity (linear + cubic)")
  cat("\nNonlinear Local Projections (fNLLP)\n")
  cat(strrep("-", 45), "\n")
  cat("Original formula : ", deparse(x$formula_orig), "\n")
  cat("Expanded formula : ", deparse(x$formula), "\n")
  cat("Type             : ", type_lbl, "\n")
  cat("Shock            : ", x$shock, "\n")
  if (identical(x$type, "state")) {
    cat("State            : ", x$state, sprintf("(%s = %s, %s = 1 - %s)",
        x$components[1], x$state, x$components[2], x$state), "\n")
    cat("Common terms     : ",
        if (length(x$common)) paste(x$common, collapse = ", ") else "none", "\n")
  }
  cat("Design columns   : ", paste(x$design_vars, collapse = ", "), "\n")
  cat("Horizons         : 0 to", max(x$horizons), "\n")
  cat("Sample           : ",
      if (isTRUE(x$balanced)) "fixed regressor dates (balanced)" else
        "shrinks with the horizon", "\n")
  cat("Cumulative       : ", isTRUE(x$cumulative), "\n")
  cat("Confidence       : ",
      paste0(format(x$conf, trim = TRUE), "%", collapse = ", "), "\n")
  cat("NW lags          : base =", x$nw_lags,
      sprintf(", offset = %d -> effective at horizon h: max(base + h + %d, 0)\n",
              x$nw_offset, x$nw_offset))
  if (isTRUE(x$nw_prewhite)) cat("Prewhitening     :  VAR(1) (Andrews-Monahan)\n")
  cat("Observations     : ", paste(unique(range(x$nobs)), collapse = " to "), "\n")

  for (v in x$lhs_vars) {
    tab <- cbind(x$irfs[, v, 1], x$irfs[, v, 2])
    colnames(tab) <- x$components
    if (!identical(x$type, "cubic")) {
      tab <- cbind(tab, difference = x$diff[, v],
                   t_diff = x$diff[, v] / x$diff_se[, v])
    } else {
      tab <- cbind(tab, t_cubic = x$irfs[, v, 2] / x$irfs_se[, v, 2])
    }
    cat("\nResponse of '", v, "' (shock = '", x$shock, "'):\n", sep = "")
    print(round(tab, digits))
  }
  invisible(x)
}


#' @export
coef.fNLLP <- function(object, ...) {
  object$irfs
}


#' Tidy Nonlinear Local Projection Estimates
#'
#' @param x An object returned by \code{\link{fNLLP}}.
#' @param difference Logical. Add rows for the difference between the two
#'   components (\code{component = "difference"}). Ignored for
#'   \code{type = "cubic"}. Default \code{TRUE}.
#' @param ... Reserved for additional arguments.
#' @return A data frame with \code{horizon}, \code{lhs}, \code{shock},
#'   \code{component} (a factor), \code{estimate}, \code{se}, and normal
#'   confidence bounds \code{lower}, \code{upper}. Multiple confidence levels
#'   use columns named \code{lower_<level>} and \code{upper_<level>}.
#' @seealso \code{\link{fNLLP}}, \code{\link{fPlotNLLP}}
#' @examples
#' data("TenreyroThwaites2016")
#' tt <- TenreyroThwaites2016
#' tt <- tt[tt$Date >= as.Date("1968-10-01"), ]
#' fit <- fNLLP(FFR ~ RRShock + l(FFR, 1), data = tt, horizons = 8,
#'              shock = "RRShock", type = "sign", balanced = TRUE,
#'              n_threads = 2)
#' head(tidy.fNLLP(fit))
#' @export
tidy.fNLLP <- function(x, difference = TRUE, ...) {
  h_seq <- x$horizons
  ny    <- length(x$lhs_vars)
  multi <- is.list(x$irfs_upper)
  keys  <- format(x$conf, trim = TRUE)

  rows <- list()
  for (k in seq_along(x$components)) {
    d <- data.frame(
      horizon   = rep(h_seq, times = ny),
      lhs       = rep(x$lhs_vars, each = length(h_seq)),
      shock     = x$shock,
      component = x$components[k],
      estimate  = as.vector(x$irfs[, , k]),
      se        = as.vector(x$irfs_se[, , k]),
      stringsAsFactors = FALSE
    )
    rows[[k]] <- d
  }
  if (isTRUE(difference) && !identical(x$type, "cubic")) {
    rows[[3]] <- data.frame(
      horizon   = rep(h_seq, times = ny),
      lhs       = rep(x$lhs_vars, each = length(h_seq)),
      shock     = x$shock,
      component = "difference",
      estimate  = as.vector(x$diff),
      se        = as.vector(x$diff_se),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)

  for (i in seq_along(x$conf)) {
    z <- stats::qnorm(0.5 * (1 + x$conf[i] / 100))
    lo <- out$estimate - z * out$se
    up <- out$estimate + z * out$se
    if (!multi) {
      out$lower <- lo
      out$upper <- up
    } else {
      out[[paste0("lower_", keys[i])]] <- lo
      out[[paste0("upper_", keys[i])]] <- up
    }
  }
  out$component <- factor(out$component,
                          levels = unique(c(x$components, "difference")))
  out$component <- droplevels(out$component)
  rownames(out) <- NULL
  out
}


# Response to a shock of size `size`:
#   sign  : b+ max(d, 0) + b- min(d, 0)
#   cubic : b1 d + b3 d^3
#   state : d (m bH + (1 - m) bL), for state value m
# Variance a^2 V11 + b^2 V22 + 2ab V12 from the HAC block of (c1, c2).
#' Responses to Shocks of Given Sizes from a Nonlinear Local Projection
#'
#' Combines the two shock coefficients of an \code{\link{fNLLP}} fit into
#' the response to a shock of size \eqn{d}: \eqn{\beta^+ \max(d, 0) +
#' \beta^- \min(d, 0)} for \code{type = "sign"}, \eqn{\beta_1 d + \beta_3
#' d^3} for \code{type = "cubic"}, and \eqn{d (F \beta_1 + (1 - F)
#' \beta_2)} at state value \eqn{F} for \code{type = "state"}. Standard
#' errors use the covariance of the two coefficients.
#'
#' @param object An object returned by \code{\link{fNLLP}}.
#' @param size Numeric vector of shock sizes. Default 1.
#' @param state Numeric vector of regime weights in \eqn{[0, 1]}. Required
#'   for \code{type = "state"}, not allowed otherwise; every combination
#'   with \code{size} is evaluated.
#' @param conf Confidence levels in percent; defaults to those of the fit.
#' @param ... Not used.
#' @return A data frame with \code{horizon}, \code{lhs}, \code{shock},
#'   \code{size}, \code{state} (for \code{type = "state"}), \code{estimate},
#'   \code{se}, and \code{lower}/\code{upper} (or
#'   \code{lower_<level>}/\code{upper_<level>} for several levels).
#' @seealso \code{\link{fNLLP}}
#' @examples
#' data("TenreyroThwaites2016")
#' tt <- TenreyroThwaites2016
#' tt <- tt[tt$Date >= as.Date("1968-10-01"), ]
#' fit <- fNLLP(FFR ~ RRShock + l(FFR, 1), data = tt, horizons = 8,
#'              shock = "RRShock", type = "cubic", balanced = TRUE,
#'              n_threads = 2)
#' # A small and a large shock
#' head(predict(fit, size = c(0.25, 2)))
#' @export
predict.fNLLP <- function(object, size = 1, state = NULL, conf = NULL, ...) {
  if (!is.numeric(size) || length(size) == 0 || anyNA(size) ||
      any(!is.finite(size)))
    stop("predict.fNLLP: 'size' must be finite numbers.")
  conf <- if (is.null(conf)) object$conf else .fLP_validate_conf_vector(conf)

  if (identical(object$type, "state")) {
    if (is.null(state))
      stop("predict.fNLLP: type = \"state\" needs 'state', the value(s) of ",
           "the regime weight in [0, 1] at which to evaluate the response.")
    if (!is.numeric(state) || anyNA(state) || any(state < 0 | state > 1))
      stop("predict.fNLLP: 'state' must lie in [0, 1].")
    grid <- expand.grid(size = size, state = state)
    a <- grid$size * grid$state
    b <- grid$size * (1 - grid$state)
  } else {
    if (!is.null(state))
      stop("predict.fNLLP: 'state' is only used with type = \"state\".")
    grid <- data.frame(size = size)
    if (identical(object$type, "sign")) {
      a <- pmax(size, 0)
      b <- pmin(size, 0)
    } else {
      a <- size
      b <- size^3
    }
  }

  h_seq <- object$horizons
  rows  <- vector("list", nrow(grid) * length(object$lhs_vars))
  idx   <- 1
  for (v in object$lhs_vars) {
    c1  <- object$irfs[, v, 1]
    c2  <- object$irfs[, v, 2]
    v11 <- object$irfs_se[, v, 1]^2
    v22 <- object$irfs_se[, v, 2]^2
    v12 <- object$irfs_cov[, v]
    for (g in seq_len(nrow(grid))) {
      est <- a[g] * c1 + b[g] * c2
      se  <- sqrt(pmax(a[g]^2 * v11 + b[g]^2 * v22 + 2 * a[g] * b[g] * v12, 0))
      d <- data.frame(horizon = h_seq, lhs = v, shock = object$shock,
                      size = grid$size[g], stringsAsFactors = FALSE)
      if (!is.null(grid$state)) d$state <- grid$state[g]
      d$estimate <- est
      d$se       <- se
      for (cl in conf) {
        z <- stats::qnorm(0.5 * (1 + cl / 100))
        nm <- if (length(conf) == 1) c("lower", "upper") else
          paste0(c("lower_", "upper_"), format(cl, trim = TRUE))
        d[[nm[1]]] <- est - z * se
        d[[nm[2]]] <- est + z * se
      }
      rows[[idx]] <- d
      idx <- idx + 1
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
