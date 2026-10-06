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

# Formula -> expanded LHS / RHS names, with l()/f() columns added to data.
# Same steps as fLP(): macro-lag, macro, lag/lead, then term labels.
.fNLLP_parse <- function(formula, data, env) {
  formula_chr <- paste(deparse(formula), collapse = "")
  formula_chr <- .fLP_expand_macro_lag(formula_chr, env)
  formula_chr <- .fLP_expand_macro(formula_chr, env)

  expanded    <- .fLP_expand_lag_terms(formula_chr, data)
  formula_chr <- expanded$formula_chr
  data        <- expanded$data

  formula_expanded <- tryCatch(
    as.formula(formula_chr, env = env),
    error = function(e) stop(
      "fNLLP: could not parse the expanded formula:\n  ", formula_chr,
      "\n  ", conditionMessage(e)
    )
  )

  .fLP_check_formula(formula_expanded, "fNLLP")

  lhs_vars  <- all.vars(formula_expanded[[2]])
  rhs_str   <- paste(deparse(formula_expanded[[3]]), collapse = " ")
  rhs_terms <- terms(as.formula(paste("~", rhs_str), env = env),
                     keep.order = TRUE)
  rhs_vars  <- attr(rhs_terms, "term.labels")

  if (length(lhs_vars) == 0)
    stop("fNLLP: formula must have at least one LHS variable.")
  if (length(rhs_vars) == 0)
    stop("fNLLP: formula must have at least one RHS variable.")

  list(lhs_vars = lhs_vars, rhs_vars = rhs_vars, data = data,
       formula_expanded = formula_expanded)
}


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

  for (nm in c("balanced", "cumulative", "store_full")) {
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
  parsed           <- .fNLLP_parse(formula, data, env)
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
    Y_pre         = Y_pre
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

  .band <- function(cl, sign) {
    res$irfs + sign * stats::qnorm(0.5 * (1.0 + cl / 100)) * res$irfs_se
  }
  if (length(conf) == 1) {
    res$irfs_upper <- .band(conf, 1)
    res$irfs_lower <- .band(conf, -1)
  } else {
    keys <- format(conf, trim = TRUE)
    res$irfs_upper <- stats::setNames(lapply(conf, .band, sign = 1), keys)
    res$irfs_lower <- stats::setNames(lapply(conf, .band, sign = -1), keys)
  }

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
