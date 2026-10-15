# -----------------------------------------------------------------------
# fLP.R — Classical Local Projections with fixest-like formula syntax
#
# Provides:
#   fLP()        — main user-facing function
#   print.fLP()  — print method
#   tidy.fLP()   — broom-compatible tidy method
#   coef.fLP()   — extract IRF coefficients
#
# Syntax mirrors fixest::feols():
#   fLP(y ~ shock + control1 + control2, data, horizons = 12, shock = "shock")
#
# Macro expansion (like fixest's ..vars):
#   controls <- c("ip_l1", "cpi_l1")
#   fLP(ip ~ mp_shock + ..controls, data, horizons = 12, shock = "mp_shock")
#
# Lag / lead operators:
#   l(var, k)     — lag var by k periods  → column "var_lk"
#   l(var, k1:k2) — lags k1 through k2    → columns "var_lk1", ..., "var_lk2"
#   f(var, k)     — lead var by k periods  → column "var_fk"
#
#   Works inside macros too:
#   l(..controls, 1:3) — lags all variables in 'controls' by 1, 2, 3
#
# Multi-equation:
#   fLP(c(ip, cpi, rate) ~ mp_shock + ..controls, ..., shock = "mp_shock")
#   (runs one LP per LHS variable; LHS vars are NOT added to the RHS)
#
# Author: Dr. Muhsin Ciftci
# -----------------------------------------------------------------------


# =======================================================================
# Internal helpers
# =======================================================================

.fLP_validate_horizon_max <- function(horizons) {
  if (!is.numeric(horizons) || length(horizons) == 0L ||
      anyNA(horizons) || any(!is.finite(horizons)) ||
      any(horizons < 0) || any(horizons != floor(horizons))) {
    stop("fLP: 'horizons' must be a non-negative integer or the sequence 0:H.")
  }

  h <- as.integer(horizons)
  if (length(h) == 1L) {
    return(h)
  }

  max_h <- max(h)
  if (!identical(sort(unique(h)), 0L:max_h)) {
    stop("fLP: vector 'horizons' must be the complete sequence 0:H.")
  }
  max_h
}

.fLP_validate_scalar_int <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x < 0 || x != floor(x)) {
    stop(sprintf("fLP: '%s' must be a non-negative integer.", name))
  }
  as.integer(x)
}

# Like .fLP_validate_scalar_int but allows negative values (e.g. nw_offset,
# which is legitimately -1 for the classic Jorda h-1 NW bandwidth rule).
.fLP_validate_scalar_int_signed <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x != floor(x)) {
    stop(sprintf("fLP: '%s' must be an integer.", name))
  }
  as.integer(x)
}

.fLP_validate_scalar_number <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x <= lower || x >= upper) {
    stop(sprintf("fLP: '%s' must be a scalar in (%s, %s).",
                 name, lower, upper))
  }
  as.double(x)
}

# Validates 'conf': one or more confidence levels in PERCENT, e.g. 90 or
# c(68, 95) for a double-band plot. Returned sorted descending (widest first)
# so downstream band-building and plotting order is deterministic. Percent is
# the package-wide convention; every band argument in tidyMacro uses it.
.fLP_validate_conf_vector <- function(conf) {
  if (!is.numeric(conf) || length(conf) == 0L || anyNA(conf) ||
      any(!is.finite(conf)) || any(conf <= 0) || any(conf >= 100)) {
    stop("fLP: 'conf' must be one or more confidence levels in (0, 100), e.g. 90 or c(68, 95).")
  }
  # A fraction is the likeliest mistake now that the scale is percent, and
  # 0.90 would otherwise pass as a 0.9% band. Reject it by name.
  if (is.numeric(conf) && length(conf) && any(conf < 1))
    stop("fLP: 'conf' is in percent, not a fraction. Use ",
         paste(format(conf[conf < 1] * 100, trim = TRUE), collapse = ", "),
         " instead of ", paste(format(conf[conf < 1], trim = TRUE), collapse = ", "),
         ".", call. = FALSE)
  conf <- sort(unique(as.double(conf)), decreasing = TRUE)
  conf
}

# -----------------------------------------------------------------------
# Step 1-helper: expand l(..macro, lags) / f(..macro, lags)
#   Called BEFORE general ..macro expansion so the inner macro token
#   is still intact and can be found in the parent environment.
#
#   l(..controls, 1:3)  →  l(ip, 1:3) + l(cpi, 1:3) + l(rate, 1:3)
# -----------------------------------------------------------------------
.fLP_expand_macro_lag <- function(formula_chr, env) {

  # Pattern: l or f, then (..macroname, anything)
  pattern <- "(l|f)\\(\\.\\.[A-Za-z_.][A-Za-z0-9_.]*,([^)]+)\\)"

  repeat {
    m <- regexpr(pattern, formula_chr, perl = TRUE)
    if (m == -1L) break

    full_match <- regmatches(formula_chr, m)
    parsed     <- regmatches(full_match,
                             regexec("(l|f)\\(\\.\\.([A-Za-z_.][A-Za-z0-9_.]*),([^)]+)\\)",
                                     full_match, perl = TRUE))[[1]]
    op       <- parsed[2L]            # "l" or "f"
    mac_name <- parsed[3L]            # macro name (without ..)
    lags_str <- trimws(parsed[4L])   # e.g. "1:3"

    if (!exists(mac_name, envir = env, inherits = TRUE)) {
      stop(sprintf("fLP: macro '..%s' not found in the calling environment.", mac_name))
    }
    vars <- get(mac_name, envir = env, inherits = TRUE)
    if (!is.character(vars) || length(vars) == 0L) {
      stop(sprintf("fLP: macro '..%s' must be a non-empty character vector.", mac_name))
    }

    # Build e.g. l(ip, 1:3) + l(cpi, 1:3) + l(rate, 1:3)
    expansion <- paste(
      sprintf("%s(%s, %s)", op, vars, lags_str),
      collapse = " + "
    )
    formula_chr <- sub(pattern, expansion, formula_chr, perl = TRUE)
  }

  formula_chr
}


# -----------------------------------------------------------------------
# Step 2-helper: expand plain ..macro tokens (not inside l/f calls)
# -----------------------------------------------------------------------
.fLP_expand_macro <- function(formula_chr, env) {

  pattern <- "\\.\\.[A-Za-z_.][A-Za-z0-9_.]*"
  matches <- unique(regmatches(formula_chr,
                               gregexpr(pattern, formula_chr))[[1]])

  for (macro in matches) {
    mac_name <- sub("^\\.\\.", "", macro)
    if (!exists(mac_name, envir = env, inherits = TRUE)) {
      stop(sprintf("fLP: macro '..%s' not found in the calling environment.", mac_name))
    }
    vars <- get(mac_name, envir = env, inherits = TRUE)
    if (!is.character(vars) || length(vars) == 0L) {
      stop(sprintf("fLP: macro '..%s' must be a non-empty character vector.", mac_name))
    }
    formula_chr <- gsub(macro,
                        paste(vars, collapse = " + "),
                        formula_chr, fixed = TRUE)
  }

  formula_chr
}


# -----------------------------------------------------------------------
# Step 3-helper: expand l(var, lags) / f(var, leads) term-by-term.
#   Creates new columns in (a copy of) data, returns updated formula
#   string and data.
#
#   Supported lag specs:
#     l(ip, 1)      → ip_l1
#     l(ip, 1:3)    → ip_l1 + ip_l2 + ip_l3
#     f(ip, 1:2)    → ip_f1 + ip_f2
# -----------------------------------------------------------------------
# ---------------------------------------------------------------------
# Build the estimation sample from a data frame of model variables.
#
# Two things the naive `raw[complete.cases(raw), ]` compaction gets wrong:
#
#   1. Interior gaps. The C++ engines treat retained rows as consecutive
#      time periods, so deleting an interior row silently makes t-1 and t+1
#      adjacent — leads, long differences and HAC autocovariances then use
#      the wrong spacing. Leading NAs (generated lags) and trailing NAs
#      (generated leads) are contiguous after removal and therefore fine;
#      an interior gap is not, and is rejected.
#
#   2. The cumulative base. A cumulative LP regresses y_{t+h} - y_{t-1}.
#      When the leading trim is driven by RHS lags, the outcome one period
#      before the first estimable date is still observed, and it is the
#      correct long-difference base for that date (LPmodel.m's
#      endo_lag1 = ENDO(lag_trim:end-1)). Dropping it costs one observation
#      at every horizon. It is returned separately as `pre` so the engine
#      can use it as a base without treating it as an estimable row.
#
# Returns list(raw = <estimation rows>, pre = <1-row LHS frame or NULL>).
# ---------------------------------------------------------------------
# Resolve the shock-size convention. Mirrors LPopt.impact in LPmodel.m:
#   "unit" (impact = 1): IRF is the response to a one-unit shock.
#   "sd"   (impact = 0): IRF is the response to a one-standard-deviation
#                        shock, the sd taken over the ESTIMATION sample
#                        (LPmodel.m z-scores/scales the lag-trimmed series,
#                        not the raw input column).
.fLP_shock_size <- function(shock_size, fn = "fLP") {
  if (is.null(shock_size)) return("unit")
  if (!is.character(shock_size) || length(shock_size) < 1L)
    stop(sprintf("%s: 'shock_size' must be \"unit\" or \"sd\".", fn))
  out <- match.arg(shock_size[1L], c("unit", "sd"))
  out
}


# Reject formula syntax the LP engines cannot honour, instead of estimating
# something else silently:
#   - LHS transformations, e.g. log(y): all.vars() would keep only "y".
#     The LHS must be a name or c(name, name, ...).
#   - intercept removal (- 1, + 0): the engines always add the intercept.
#   - RHS terms that are not plain columns, e.g. log(x), I(x^2), x:z.
# Transformed variables belong in `data` as their own columns.
.fLP_check_formula <- function(formula_expanded, fn = "fLP") {
  lhs <- formula_expanded[[2]]
  lhs_ok <- is.name(lhs) ||
    (is.call(lhs) && identical(lhs[[1]], as.name("c")) && length(lhs) > 1 &&
       all(vapply(as.list(lhs)[-1], is.name, logical(1))))
  if (!lhs_ok)
    stop(sprintf(paste0(
      "%s: the LHS must be a variable name or c(y1, y2, ...); got '%s'.\n",
      "  Create transformed outcomes (e.g. log(y)) as columns of `data` first."),
      fn, paste(deparse(lhs), collapse = " ")), call. = FALSE)

  tt <- terms(formula_expanded, keep.order = TRUE)
  if (attr(tt, "intercept") == 0)
    stop(sprintf(paste0(
      "%s: the intercept is always included; remove '- 1' / '+ 0' from the ",
      "formula."), fn), call. = FALSE)

  labs <- attr(tt, "term.labels")
  plain <- vapply(labs, function(l) is.name(str2lang(l)), logical(1))
  if (!all(plain))
    stop(sprintf(paste0(
      "%s: RHS terms must be columns of `data` (or l()/f() lags and leads); ",
      "got %s.\n  Create transformed regressors as columns of `data` first."),
      fn, paste(sprintf("'%s'", labs[!plain]), collapse = ", ")), call. = FALSE)
  invisible(TRUE)
}


# Expand a formula: l(..macro, k) first, then plain ..macro tokens, then
# l()/f() terms (lag/lead columns are added to `data`). Returns the parsed
# expanded formula and the augmented data.
.fLP_expand_formula <- function(formula, data, env, fn = "fLP") {
  formula_chr <- paste(deparse(formula), collapse = "")
  formula_chr <- .fLP_expand_macro_lag(formula_chr, env)
  formula_chr <- .fLP_expand_macro(formula_chr, env)
  expanded    <- .fLP_expand_lag_terms(formula_chr, data)

  formula_expanded <- tryCatch(
    as.formula(expanded$formula_chr, env = env),
    error = function(e) stop(
      fn, ": could not parse the expanded formula:\n  ", expanded$formula_chr,
      "\n  ", conditionMessage(e)
    )
  )
  list(formula = formula_expanded, data = expanded$data)
}

# Term labels of one side of a formula, in order. deparse() can split long
# expressions, so the pieces are collapsed before as.formula().
.fLP_term_labels <- function(side, env) {
  side_str <- paste(deparse(side), collapse = " ")
  attr(terms(as.formula(paste("~", side_str), env = env), keep.order = TRUE),
       "term.labels")
}

# Two-sided LP formula (fLP, fNLLP, fLPIV) -> expanded LHS / RHS names, with
# l()/f() columns added to data. The LHS may be y or c(y1, y2); the RHS
# excludes the intercept, which the engines always add.
.fLP_parse <- function(formula, data, env, fn = "fLP") {
  ex <- .fLP_expand_formula(formula, data, env, fn)
  .fLP_check_formula(ex$formula, fn)

  lhs_vars <- all.vars(ex$formula[[2L]])
  rhs_vars <- .fLP_term_labels(ex$formula[[3L]], env)
  if (length(lhs_vars) == 0L)
    stop(fn, ": formula must have at least one LHS variable.")
  if (length(rhs_vars) == 0L)
    stop(fn, ": formula must have at least one RHS variable.")

  list(lhs_vars = lhs_vars, rhs_vars = rhs_vars, data = ex$data,
       formula_expanded = ex$formula)
}


.fLP_estimation_sample <- function(data, all_needed, lhs_vars, cumulative,
                                   fn = "fLP") {
  raw  <- data[, all_needed, drop = FALSE]
  keep <- which(complete.cases(raw))

  if (length(keep) == 0L)
    stop(sprintf("%s: no complete rows after lag/lead expansion.", fn))

  # ---- reject interior gaps ------------------------------------------
  gaps <- which(diff(keep) != 1L)
  if (length(gaps) > 0L) {
    at <- keep[gaps]
    stop(sprintf(
      paste0("%s: the model variables have %d interior gap(s) in the time ",
             "index (after row(s) %s of `data`).\n",
             "  Rows with NA cannot simply be deleted from a time-series LP: ",
             "the engine treats the retained rows as consecutive periods, so ",
             "leads, long differences and Newey-West autocovariances would ",
             "use the wrong spacing.\n",
             "  Fill the gap, or subset `data` to a contiguous span before ",
             "calling %s()."),
      fn, length(gaps),
      paste(at[seq_len(min(5L, length(at)))], collapse = ", "), fn))
  }

  n_dropped <- nrow(raw) - length(keep)

  # ---- recover the cumulative base row --------------------------------
  pre <- NULL
  if (isTRUE(cumulative) && keep[1L] > 1L) {
    cand <- data[keep[1L] - 1L, lhs_vars, drop = FALSE]
    if (all(complete.cases(cand))) pre <- cand
  }

  if (n_dropped > 0L) {
    message(sprintf(
      "%s: removing %d row(s) with NA (from lag/lead or missing data).%s",
      fn, n_dropped,
      if (!is.null(pre))
        " Retained the preceding outcome as the cumulative base." else ""
    ))
  }

  list(raw = raw[keep, , drop = FALSE], pre = pre)
}


# Expand l(x, k) / f(x, k) terms into x_lk / x_fk columns. `shift(x, k, op)`
# returns x shifted by k > 0 rows (default: plain time-series shift; the
# panel version looks up (id, t -/+ k)). `paren` wraps multi-column
# replacements so ":" interactions distribute: l(x, 1:2):z -> (x_l1 + x_l2):z.
.fLP_expand_lag_terms <- function(formula_chr, data, shift = NULL,
                                  fn = "fLP", paren = FALSE) {

  # Match l(...) or f(...) anywhere in the string.
  # We iterate to handle multiple occurrences.
  pattern <- "(l|f)\\(([A-Za-z_.][A-Za-z0-9_.]*),([^)]+)\\)"
  T_dat   <- nrow(data)

  if (is.null(shift)) {
    shift <- function(x_raw, k, op) {
      if (k >= T_dat) {
        rep(NA_real_, T_dat)
      } else if (op == "l") {
        # lag: shift forward k rows, fill head with NA
        c(rep(NA_real_, k), x_raw[seq_len(T_dat - k)])
      } else {
        # lead: shift backward k rows, fill tail with NA
        c(x_raw[seq.int(k + 1L, T_dat)], rep(NA_real_, k))
      }
    }
  }

  repeat {
    m <- regexpr(pattern, formula_chr, perl = TRUE)
    if (m == -1L) break

    full_match <- regmatches(formula_chr, m)
    parsed     <- regmatches(full_match,
                             regexec(pattern, full_match, perl = TRUE))[[1]]
    op       <- parsed[2L]            # "l" or "f"
    var      <- trimws(parsed[3L])   # variable name
    lags_str <- trimws(parsed[4L])   # lag spec, e.g. "1:3"

    # Evaluate lag indices safely
    lag_idx <- tryCatch(
      eval(parse(text = lags_str)),
      error = function(e) stop(sprintf(
        "%s: cannot parse lag/lead spec '%s' in term '%s': %s",
        fn, lags_str, full_match, conditionMessage(e)
      ))
    )
    if (!is.numeric(lag_idx) || anyNA(lag_idx) || any(!is.finite(lag_idx)) ||
        any(lag_idx < 0L) || any(lag_idx != floor(lag_idx))) {
      stop(sprintf(
        "%s: lag/lead indices must be non-negative integers; got '%s'.",
        fn, lags_str
      ))
    }
    lag_idx <- as.integer(lag_idx)

    if (!var %in% names(data)) {
      stop(sprintf(
        "%s: variable '%s' used in lag/lead term '%s' not found in data.",
        fn, var, full_match
      ))
    }

    x_raw <- data[[var]]

    new_cols <- character(length(lag_idx))
    for (i in seq_along(lag_idx)) {
      k        <- lag_idx[i]
      col_name <- if (op == "l") sprintf("%s_l%d", var, k) else
                                 sprintf("%s_f%d", var, k)
      shifted  <- if (k == 0L) x_raw else shift(x_raw, k, op)

      # A column with the generated name may already exist (e.g. the same
      # l() term twice). Reuse it only if it IS this lag/lead; an unrelated
      # column that happens to be called x_l1 must not be estimated silently.
      if (col_name %in% names(data)) {
        old <- data[[col_name]]
        same <- length(old) == T_dat &&
          identical(is.na(old), is.na(shifted)) &&
          isTRUE(all(old[!is.na(old)] == shifted[!is.na(shifted)]))
        if (!same)
          stop(sprintf(paste0(
            "%s: `data` already has a column '%s' that is not %s of '%s'.\n",
            "  Rename that column, or drop it, before using %s."),
            fn, col_name,
            if (op == "l") sprintf("lag %d", k) else sprintf("lead %d", k),
            var, full_match), call. = FALSE)
      } else {
        data[[col_name]] <- shifted
      }
      new_cols[i] <- col_name
    }

    # Replace the lag/lead term with expanded column names
    replacement <- paste(new_cols, collapse = " + ")
    if (paren && length(new_cols) > 1L)
      replacement <- paste0("(", replacement, ")")
    formula_chr <- sub(pattern, replacement, formula_chr, perl = TRUE)
  }

  list(formula_chr = formula_chr, data = data)
}


# =======================================================================
# Main function
# =======================================================================

#' Classical Local Projections (fixest-like syntax)
#'
#' Estimates impulse responses via the Jordà (2005) local projection method.
#' The syntax mirrors \pkg{fixest}: use a formula for the model,
#' \code{..macroname} tokens to expand character-vector variable lists, and
#' \code{l()} / \code{f()} operators to create lags and leads on-the-fly.
#'
#' @param formula A formula of the form \code{y ~ shock + controls}.
#'   \itemize{
#'     \item \strong{LHS}: one variable, or \code{c(y1, y2)} for multi-equation.
#'       LHS variables are \emph{never} automatically placed on the RHS.
#'     \item \strong{RHS}: exactly what you want in the regression. No automatic
#'       lag augmentation. Use \code{l()} / \code{f()} for that.
#'     \item \strong{Macros}: write \code{..name} where \code{name} is a
#'       character vector in the calling environment.
#'     \item \strong{Lag operator}: \code{l(var, k)} or \code{l(var, k1:k2)}.
#'     \item \strong{Lead operator}: \code{f(var, k)} or \code{f(var, k1:k2)}.
#'     \item Macros inside lag/lead: \code{l(..controls, 1:3)}.
#'   }
#' @param data A \code{data.frame} or \code{tibble}. Must contain all variables
#'   referenced in the formula (before lag/lead expansion).
#' @param horizons Integer. Maximum horizon H, or the complete sequence
#'   \code{0:H}; projections run for \eqn{h = 0, 1, \ldots, H}. Default: 12.
#' @param shock Character. \strong{Required.} Name of the shock/impulse variable
#'   (must appear on the RHS). This is the variable whose coefficient at each
#'   horizon is returned as the IRF.
#' @param conf Numeric in (0, 100), in percent. One or more confidence levels
#'   for bands, e.g. \code{90} or \code{c(68, 95)}. When more than one is supplied,
#'   \code{irfs_upper} and \code{irfs_lower} become named lists of matrices
#'   (one per level, keyed by the level as a string) and the plot functions
#'   (\code{fPlotIRF()}, \code{fPlotNLLP()}) draw one ribbon per level. Bands are always rebuilt in R from \code{irfs_se},
#'   so multi-level fits cost the same as single-level fits. Default: 90.
#' @param nw_lags Integer or \code{NULL}. Base Newey-West bandwidth.
#'   At horizon \eqn{h} the effective bandwidth is
#'   \code{max(nw_lags + h + nw_offset, 0)} (Bartlett kernel).
#'   Defaults to \code{floor(T^(1/3))}.
#' @param nw_offset Integer. Shift applied to the NW bandwidth rule. Default
#'   \code{1L}, the Miranda-Agrippino & Ricco convention (bandwidth grows as
#'   \eqn{h+1}). Set to \code{0L} to reproduce the classic Jordà (2005) rule,
#'   bandwidth = \eqn{h} (LP residuals at horizon \eqn{h} are MA(h) by
#'   construction). \code{fLP}'s \code{h} is 0-indexed (h=0 is the impact
#'   horizon), which makes this the exact equivalent of Cesa-Bianchi's VAR
#'   Toolbox (\code{LPmodel.m}: \code{OLSmodel(Y,X,0,hh-1)} with 1-indexed
#'   \code{hh}, so \code{hh-1 == h}) — needed to match that toolbox's SEs
#'   exactly, e.g. in the Jordà-Taylor (2025) replication.
#' @param nw_prewhite Logical. If \code{TRUE}, apply Andrews-Monahan (1992)
#'   VAR(1) prewhitening to the regression scores before the Newey-West sum
#'   (same bandwidth), shrinking the VAR(1) matrix when its spectral radius
#'   exceeds 0.97. Default \code{FALSE} (plain Newey-West, as in the VAR
#'   Toolbox).
#' @param store_full Logical. If \code{TRUE}, return full coefficient and
#'   standard-error matrices for every horizon. Defaults to \code{FALSE} to
#'   keep the result small and fast.
#' @param cumulative Logical. If \code{TRUE}, the LHS at horizon \eqn{h} is
#'   \eqn{y_{t+h} - y_{t-1}} (cumulative response from \eqn{t-1} to \eqn{t+h}).
#'   Default \code{FALSE}. When generated lags trim the leading rows, the
#'   outcome immediately before the first estimable date is still observed
#'   and is used as that date's long-difference base, so no observation is
#'   lost — matching \code{endo_lag1} in Cesa-Bianchi's \code{LPmodel.m}.
#' @param shock_size Character. Units of the reported impulse response.
#'   \code{"unit"} (default) = response to a one-unit shock. \code{"sd"} =
#'   response to a one-standard-deviation shock, obtained by z-scoring the
#'   shock over the estimation sample. These correspond to
#'   \code{LPopt.impact = 1} and \code{LPopt.impact = 0} respectively in
#'   Cesa-Bianchi's \code{LPmodel.m}; note the standard deviation is taken
#'   over the lag-trimmed estimation sample, not the raw input column.
#' @param n_threads Integer. Number of OpenMP threads for the horizon loop.
#'   \code{0} = all cores minus one. Default \code{0}.
#'
#' @return An object of class \code{"fLP"}.
#'
#' @examples
#' \dontrun{
#' ## 1. Basic usage — explicit shock, explicit controls
#' fLP(ip ~ mp_shock + ip_l1 + ip_l2 + cpi_l1,
#'     data = df, horizons = 24, shock = "mp_shock")
#'
#' ## 2. Macro expansion
#' controls <- c("ip_l1", "ip_l2", "cpi_l1", "rate_l1")
#' fLP(ip ~ mp_shock + ..controls,
#'     data = df, horizons = 24, shock = "mp_shock")
#'
#' ## 3. Lag operator — auto-creates ip_l1, ip_l2, ip_l3 from data
#' fLP(ip ~ mp_shock + l(ip, 1:3) + l(cpi, 1:2),
#'     data = df, horizons = 24, shock = "mp_shock")
#'
#' ## 4. Macro + lag operator
#' vars <- c("ip", "cpi", "rate")
#' fLP(ip ~ mp_shock + l(..vars, 1:3),
#'     data = df, horizons = 24, shock = "mp_shock")
#'
#' ## 5. Multi-equation (ip, cpi, rate as LHS — separate regressions)
#' fLP(c(ip, cpi, rate) ~ mp_shock + l(..vars, 1:3),
#'     data = df, horizons = 24, shock = "mp_shock")
#'
#' ## 6. Lead operator
#' fLP(ip ~ mp_shock + l(ip, 1:12) + f(sentiment, 1),
#'     data = df, horizons = 24, shock = "mp_shock")
#' }
#'
#' @export
fLP <- function(formula, data, horizons = 12L,
                shock,
                conf         = 90,
                nw_lags      = NULL,
                nw_offset    = 1L,
                nw_prewhite  = FALSE,
                store_full   = FALSE,
                cumulative   = FALSE,
                shock_size   = c("unit", "sd"),
                n_threads    = 0L) {

  # Capture caller environment for macro lookup
  env <- parent.frame()

  shock_size <- .fLP_shock_size(shock_size, "fLP")

  # =====================================================================
  # 0. Require 'shock' argument
  # =====================================================================
  if (missing(shock)) {
    stop(
      "fLP: 'shock' is required. ",
      "Pass the name of the impulse variable, e.g. shock = \"mp_shock\"."
    )
  }
  if (!is.character(shock) || length(shock) != 1L || is.na(shock)) {
    stop("fLP: 'shock' must be a single character string.")
  }

  H <- .fLP_validate_horizon_max(horizons)
  conf <- .fLP_validate_conf_vector(conf)   # sorted descending (widest first)
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  # =====================================================================
  # 1-4. Expand macros and l()/f() terms; parse LHS / RHS variable names
  # =====================================================================
  parsed           <- .fLP_parse(formula, data, env, "fLP")
  lhs_vars         <- parsed$lhs_vars
  rhs_vars         <- parsed$rhs_vars
  data             <- parsed$data   # may have new lag/lead columns
  formula_expanded <- parsed$formula_expanded

  # =====================================================================
  # 5. Validate shock variable
  # =====================================================================
  if (!shock %in% rhs_vars) {
    stop(sprintf(
      "fLP: shock variable '%s' not found on the RHS.\n  RHS terms: %s",
      shock, paste(rhs_vars, collapse = ", ")
    ))
  }

  # 0-indexed column position in X (C++ convention)
  shock_col <- which(rhs_vars == shock) - 1L

  # =====================================================================
  # 6. Validate all variables exist in data
  # =====================================================================
  all_needed <- unique(c(lhs_vars, rhs_vars))
  missing_v  <- setdiff(all_needed, names(data))
  if (length(missing_v) > 0L) {
    stop(sprintf(
      "fLP: the following variables are not in data: %s",
      paste(missing_v, collapse = ", ")
    ))
  }

  # =====================================================================
  # 7. Build numeric matrices Y (T x ny) and X (T x k)
  #    Remove rows with NA in any used variable.
  # =====================================================================
  samp <- .fLP_estimation_sample(data, all_needed, lhs_vars,
                                 cumulative = cumulative, fn = "fLP")
  raw  <- samp$raw

  non_numeric <- names(raw)[!vapply(raw, is.numeric, logical(1))]
  if (length(non_numeric) > 0L) {
    stop(sprintf(
      "fLP: all model variables must be numeric. Non-numeric variable(s): %s",
      paste(non_numeric, collapse = ", ")
    ))
  }

  Y <- as.matrix(raw[, lhs_vars, drop = FALSE])
  X <- as.matrix(raw[, rhs_vars, drop = FALSE])

  storage.mode(Y) <- "double"
  storage.mode(X) <- "double"

  # Cumulative long-difference base: y_{t-1} for the first estimable date.
  # NULL whenever it does not exist (no earlier row, or a non-cumulative fit).
  Y_pre <- NULL
  if (!is.null(samp$pre)) {
    Y_pre <- as.matrix(samp$pre[, lhs_vars, drop = FALSE])
    storage.mode(Y_pre) <- "double"
  }

  # Shock-size convention. "sd" reproduces LPmodel.m's impact = 0, which
  # estimates on zscore(s) where s is the LAG-TRIMMED shock — i.e. the
  # estimation sample built above, not the raw column. Rescaling the
  # regressor (rather than the fitted coefficient) matches the toolbox
  # operation order exactly.
  if (identical(shock_size, "sd")) {
    sc_idx  <- shock_col + 1L
    shock_v <- X[, sc_idx]
    shock_sd <- stats::sd(shock_v)
    if (!is.finite(shock_sd) || shock_sd <= 0)
      stop("fLP: shock_size = \"sd\" requires a shock with positive variance.")
    X[, sc_idx] <- (shock_v - mean(shock_v)) / shock_sd
  }

  T_eff <- nrow(Y)

  # =====================================================================
  # 8. Default Newey-West base bandwidth
  # =====================================================================
  if (is.null(nw_lags)) {
    nw_lags <- as.integer(floor(T_eff^(1.0 / 3.0)))
  } else {
    nw_lags <- .fLP_validate_scalar_int(nw_lags, "nw_lags")
  }
  nw_offset <- .fLP_validate_scalar_int_signed(nw_offset, "nw_offset")

  # =====================================================================
  # 9. Sanity checks
  # =====================================================================
  if (!is.logical(store_full) || length(store_full) != 1L || is.na(store_full))
    stop("fLP: 'store_full' must be TRUE or FALSE.")
  if (!is.logical(cumulative) || length(cumulative) != 1L || is.na(cumulative))
    stop("fLP: 'cumulative' must be TRUE or FALSE.")
  if (!is.logical(nw_prewhite) || length(nw_prewhite) != 1L || is.na(nw_prewhite))
    stop("fLP: 'nw_prewhite' must be TRUE or FALSE.")
  n_threads <- .fLP_validate_scalar_int(n_threads, "n_threads")
  min_obs <- ncol(X) + 2L + H +
             as.integer(isTRUE(cumulative) && is.null(Y_pre))
  if (T_eff < min_obs) {
    stop(sprintf(
      "fLP: not enough observations. Need at least %d rows after NA removal; have %d.",
      min_obs, T_eff
    ))
  }

  # =====================================================================
  # 10. Call C++ engine (all horizon loops in native C++)
  # =====================================================================
  res_cpp <- fLP_cpp(
    Y            = Y,
    X            = X,
    H            = H,
    shock_col    = as.integer(shock_col),
    conf_level   = as.double(conf[1L]) / 100,  # C++ wants a fraction
    nw_lags_base = as.integer(nw_lags),
    store_full   = isTRUE(store_full),
    cumulative   = isTRUE(cumulative),
    n_threads    = as.integer(n_threads),
    nw_offset    = as.integer(nw_offset),
    verbose      = FALSE,
    Y_pre        = Y_pre,
    nw_prewhite  = nw_prewhite
  )

  # =====================================================================
  # 11. Annotate and return
  # =====================================================================
  h_seq <- 0L:H
  dimnames(res_cpp$irfs)    <- list(as.character(h_seq), lhs_vars)
  dimnames(res_cpp$irfs_se) <- dimnames(res_cpp$irfs)

  # Bands rebuilt from irfs_se so several confidence levels come from one fit
  bands <- .fLP_bands(res_cpp$irfs, res_cpp$irfs_se, conf)
  res_cpp$irfs_upper <- bands$upper
  res_cpp$irfs_lower <- bands$lower

  if (isTRUE(store_full)) {
    rhs_with_const <- c("(Intercept)", rhs_vars)
    for (i in seq_along(res_cpp$betas)) {
      rownames(res_cpp$betas[[i]]) <- rhs_with_const
      colnames(res_cpp$betas[[i]]) <- lhs_vars
      rownames(res_cpp$ses[[i]])   <- rhs_with_const
      colnames(res_cpp$ses[[i]])   <- lhs_vars
    }
    names(res_cpp$betas) <- paste0("h", h_seq)
    names(res_cpp$ses)   <- paste0("h", h_seq)
  }

  res_cpp$lhs_vars        <- lhs_vars
  res_cpp$rhs_vars        <- rhs_vars
  res_cpp$shock           <- shock
  res_cpp$horizons        <- h_seq
  res_cpp$conf            <- conf
  res_cpp$nw_lags         <- nw_lags
  res_cpp$nw_offset       <- nw_offset
  res_cpp$nw_prewhite     <- nw_prewhite
  res_cpp$store_full      <- isTRUE(store_full)
  res_cpp$cumulative      <- isTRUE(cumulative)
  res_cpp$shock_size      <- shock_size
  res_cpp$n_threads       <- as.integer(n_threads)
  res_cpp$nobs            <- T_eff
  res_cpp$formula         <- formula_expanded
  res_cpp$formula_orig    <- formula
  res_cpp$call            <- match.call()

  class(res_cpp) <- "fLP"
  return(res_cpp)
}


# =======================================================================
# S3 methods
# =======================================================================

# Bands irfs +/- crit(level) * se at each confidence level (percent). One
# level gives matrices (arrays); several give lists keyed by format(conf).
# crit() may return a scalar or an array shaped like irfs (e.g. t quantiles).
# Dimnames are taken from irfs.
.fLP_bands <- function(irfs, se, conf,
                       crit = function(cl) stats::qnorm(0.5 * (1 + cl / 100))) {
  one <- function(cl) {
    q <- crit(cl)
    list(upper = irfs + q * se, lower = irfs - q * se)
  }
  if (length(conf) == 1L) return(one(conf))
  b <- stats::setNames(lapply(conf, one), format(conf, trim = TRUE))
  list(upper = lapply(b, `[[`, "upper"), lower = lapply(b, `[[`, "lower"))
}

# Append the band columns of an fLP-type fit to a tidy data frame: lower and
# upper for one level, lower_<level> / upper_<level> for several.
.fLP_tidy_bands <- function(base, x) {
  if (!is.list(x$irfs_upper)) {
    base$lower <- as.vector(x$irfs_lower)
    base$upper <- as.vector(x$irfs_upper)
    return(base)
  }
  for (k in names(x$irfs_upper)) {
    base[[paste0("lower_", k)]] <- as.vector(x$irfs_lower[[k]])
    base[[paste0("upper_", k)]] <- as.vector(x$irfs_upper[[k]])
  }
  base
}


#' @export
print.fLP <- function(x, digits = 4L, ...) {
  cat("\nLocal Projections (fLP)\n")
  cat(strrep("-", 45L), "\n")
  cat("Original formula : ", deparse(x$formula_orig), "\n")
  cat("Expanded formula : ", deparse(x$formula),      "\n")
  cat("Shock            : ", x$shock,                 "\n")
  cat("LHS variables    : ", paste(x$lhs_vars, collapse = ", "), "\n")
  cat("Horizons         : 0 to", max(x$horizons),     "\n")
  cat("Cumulative       : ", isTRUE(x$cumulative),     "\n")
  cat("Confidence       : ",
      paste0(format(x$conf, trim = TRUE), "%", collapse = ", "), "\n")
  nw_off <- if (!is.null(x$nw_offset)) x$nw_offset else 1L
  cat("NW lags          : base =", x$nw_lags,
      sprintf(", offset = %d -> effective at horizon h: max(base + h + %d, 0)\n",
              nw_off, nw_off))
  cat("HAC variance     :  full Newey-West sandwich (X'X)^-1 X'\u03a9X (X'X)^-1",
      if (isTRUE(x$nw_prewhite)) ", VAR(1)-prewhitened" else "", "\n", sep = "")
  cat("Observations     : ", x$nobs, "\n")
  cat("\nIRF (shock = '", x$shock, "'):\n", sep = "")
  print(round(x$irfs, digits))
  invisible(x)
}


#' @export
coef.fLP <- function(object, ...) {
  object$irfs
}


#' Tidy Local Projection Estimates
#'
#' @param x An object returned by \code{\link{fLP}} or \code{\link{fLPIV}}.
#' @param ... Reserved for additional arguments.
#' @return A data frame with horizon, response variable, shock, estimate,
#'   standard error, and confidence bounds. Multiple confidence levels use
#'   columns named \code{lower_<level>} and \code{upper_<level>}.
#' @export
tidy.fLP <- function(x, ...) {
  h_seq <- 0L:max(x$horizons)
  base <- data.frame(
    horizon  = rep(h_seq, times = length(x$lhs_vars)),
    lhs      = rep(x$lhs_vars, each = length(h_seq)),
    shock    = x$shock,
    estimate = as.vector(x$irfs),
    se       = as.vector(x$irfs_se),
    stringsAsFactors = FALSE
  )
  .fLP_tidy_bands(base, x)
}
