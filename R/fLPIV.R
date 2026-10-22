# -----------------------------------------------------------------------
# fLPIV.R — Local Projections with an External Instrument (LP-IV)
#
# Syntax mirrors fLP():
#   fLPIV(
#     formula     = y ~ D + l(..controls, 1:6),
#     instruments = ~ l(z, 0:6),
#     data        = df,
#     endog       = "D",
#     horizons    = 48,
#     conf        = c(68, 95),
#     nw_lags_iv  = 6,          # fixed HAC bandwidth (matches Stata vce(hac nw 6))
#     cumulative  = TRUE
#   )
#
# 'formula' is the main outcome equation: LHS = outcome(s), RHS =
# endogenous treatment + exogenous controls (with l()/f()/..macro support).
# 'instruments' is a one-sided formula whose RHS is the full external
# instrument set (again with l()/f()/..macro support).
#
# 'endog' names the single RHS term that is the endogenous treatment. The
# remaining RHS terms become the exogenous control matrix C.
# -----------------------------------------------------------------------

#' Local Projections with an External Instrument (LP-IV)
#'
#' Estimates 2SLS local projections horizon-by-horizon using the delta-method
#' Newey-West HAC sandwich for the endogenous treatment effect. Matches
#' Cesa-Bianchi's VAR Toolbox \code{LPmodel.m} in LP-IV mode when
#' \code{nw_lags_iv} is set to the toolbox's \code{nlag_iv}.
#'
#' @param formula A formula \code{y ~ D + controls}. LHS may be
#'   \code{c(y1, y2, ...)} for multi-outcome LP-IV. RHS lists the endogenous
#'   treatment plus all exogenous controls. \code{l()} / \code{f()} lag/lead
#'   operators and \code{..macro} expansion behave exactly as in \code{fLP()}.
#' @param instruments A one-sided formula \code{~ z + ...} whose RHS is the
#'   full external instrument matrix. Use \code{l(z, 0:6)} to include the
#'   contemporaneous instrument plus six lags (over-identified 2SLS).
#' @param data A \code{data.frame} or \code{tibble}. Rows are consecutive
#'   periods. Treatment, controls and instruments must be observed on one
#'   block of rows; outcomes may have \code{NA}, handled per horizon and
#'   outcome as in \code{\link{fLP}} (both stages use the dates whose
#'   \eqn{y_{t+h}}, and \eqn{y_{t-1}} when cumulative, is observed).
#' @param endog Character. The single RHS term (of \code{formula}) that is
#'   the endogenous treatment. The remaining RHS terms are the exogenous
#'   control matrix.
#' @param horizons Integer or the sequence \code{0:H}. Default 12.
#' @param conf Numeric in (0, 100), in percent. Scalar or vector of levels.
#'   Vectors trigger multi-band output; see \code{\link{fLP}}.
#' @param nw_lags_iv Integer. Bartlett HAC bandwidth for the IV score
#'   variance. \code{> 0}: fixed bandwidth across horizons — this matches
#'   Stata's \code{vce(hac nw <k>)} and the toolbox's LP-IV bands.
#'   \code{= 0}: horizon-varying bandwidth \code{= h} (just-identified rule).
#'   Default \code{NULL} → use \code{ncol(Z) - 1}, i.e. the lag order implied
#'   by the instrument formula, which reproduces the toolbox default.
#' @param nw_prewhite Logical. If \code{TRUE}, apply Andrews-Monahan (1992)
#'   VAR(1) prewhitening to the regression scores before the Newey-West sum
#'   (same bandwidth), shrinking the VAR(1) matrix when its spectral radius
#'   exceeds 0.97. Default \code{FALSE} (plain Newey-West, as in the VAR
#'   Toolbox).
#' @param cumulative Logical. \code{TRUE} regresses
#'   \eqn{y_{t+h} - y_{t-1}} on the RHS. Default \code{TRUE} — the standard
#'   long-difference LP-IV convention. When generated lags trim the leading
#'   rows, the outcome immediately before the first regressor date is still
#'   observed and is used as that date's base, so no observation is lost.
#' @param shock_size Character. Units of the reported impulse response.
#'   \code{"unit"} (default) = response to a one-unit move in the treatment.
#'   \code{"sd"} = response to a one-standard-deviation move, scaling the
#'   2SLS coefficient and its standard error by \code{sd(D)} over the
#'   dates with every regressor observed (one value for all horizons).
#'   These correspond to \code{LPopt.impact = 1} and
#'   \code{LPopt.impact = 0} in Cesa-Bianchi's \code{LPmodel.m}.
#' @param n_threads Integer. OpenMP threads for the horizon loop. Default 0
#'   (= all cores minus one).
#'
#' @return An object of class \code{c("fLPIV", "fLP")} so all existing
#'   \code{fLP}-family methods (\code{print}, \code{tidy.fLP}, \code{fPlotIRF})
#'   work unchanged. Additional fields: \code{Fstat_fs} and \code{rsqr_fs}
#'   (both \code{(H+1) x n_y}: first-stage F-statistic and R² per horizon
#'   and outcome, since each outcome's sample has its own first stage),
#'   \code{nobs} (\code{(H+1) x n_y} observations), plus \code{endog} and
#'   \code{instrument_vars}. \code{Fstat_fs} is the
#'   textbook exclusion F for \eqn{D \sim [1, C, Z]}, with
#'   \eqn{T_h - k_c - n_z} residual degrees of freedom. Cesa-Bianchi's
#'   \code{LPmodel.m} reports \eqn{n_h - k_z} instead (it ignores the df
#'   consumed by the residualised controls), so the two F values differ by
#'   construction. Point estimates, standard errors and bands are unaffected.
#'
#' @examples
#' data("JordaTaylor2025IV")
#'
#' fLPIV(Unemployment ~ FFRates + l(Unemployment, 1:3),
#'       instruments = ~ l(RRShock, 0:3),
#'       data        = JordaTaylor2025IV,
#'       endog       = "FFRates",
#'       horizons    = 6,
#'       n_threads   = 2)
#'
#' controls <- c("Unemployment", "Inflation", "FFRates")
#' fLPIV(Unemployment ~ FFRates + l(..controls, 1:6),
#'       instruments = ~ l(RRShock, 0:6),
#'       data        = JordaTaylor2025IV,
#'       endog       = "FFRates",
#'       horizons    = 48,
#'       conf        = c(68, 95),
#'       nw_lags_iv  = 6,
#'       cumulative  = TRUE,
#'       n_threads   = 2)
#'
#' @export
fLPIV <- function(formula, instruments, data,
                  endog,
                  horizons   = 12L,
                  conf       = 90,
                  nw_lags_iv = NULL,
                  nw_prewhite = FALSE,
                  cumulative = TRUE,
                  shock_size = c("unit", "sd"),
                  n_threads  = 0L) {

  env <- parent.frame()

  shock_size <- .fLP_shock_size(shock_size, "fLPIV")

  # -------------------------------------------------------------------
  # 0. Basic checks
  # -------------------------------------------------------------------
  if (missing(endog)) {
    stop("fLPIV: 'endog' is required (name of the endogenous treatment RHS term).")
  }
  if (!is.character(endog) || length(endog) != 1L || is.na(endog)) {
    stop("fLPIV: 'endog' must be a single character string.")
  }
  if (missing(instruments) || !inherits(instruments, "formula")) {
    stop("fLPIV: 'instruments' must be a one-sided formula, e.g. ~ l(z, 0:6).")
  }

  H    <- .fLP_validate_horizon_max(horizons)
  conf <- .fLP_validate_conf_vector(conf)
  if (!is.data.frame(data)) data <- as.data.frame(data)

  # -------------------------------------------------------------------
  # 1. Parse 'formula' with the fLP expander
  # -------------------------------------------------------------------
  parsed    <- .fLP_parse(formula, data, env, "fLPIV")
  lhs_vars  <- parsed$lhs_vars
  rhs_vars  <- parsed$rhs_vars
  data      <- parsed$data
  fmla_full <- parsed$formula_expanded

  # -------------------------------------------------------------------
  # 2. Parse 'instruments' (one-sided) with the same expander
  # -------------------------------------------------------------------
  exp_iv  <- .fLP_expand_formula(instruments, data, env, "fLPIV")
  data    <- exp_iv$data
  iv_full <- exp_iv$formula
  iv_vars <- .fLP_term_labels(iv_full[[length(iv_full)]], env)
  if (length(iv_vars) == 0L)
    stop("fLPIV: 'instruments' must contribute at least one instrument.")

  # -------------------------------------------------------------------
  # 3. Validate & split endog / exog controls
  # -------------------------------------------------------------------
  if (!endog %in% rhs_vars) {
    stop(sprintf("fLPIV: endogenous treatment '%s' not on the RHS. RHS terms: %s",
                 endog, paste(rhs_vars, collapse = ", ")))
  }
  ctrl_vars <- setdiff(rhs_vars, endog)  # can be empty

  # -------------------------------------------------------------------
  # 4. Check all needed vars exist; regressor dates and outcome rows
  # -------------------------------------------------------------------
  x_vars <- unique(c(endog, ctrl_vars, iv_vars))
  miss_v <- setdiff(c(lhs_vars, x_vars), names(data))
  if (length(miss_v) > 0L)
    stop(sprintf("fLPIV: variable(s) not in data: %s", paste(miss_v, collapse = ", ")))
  if (!is.logical(cumulative) || length(cumulative) != 1 || is.na(cumulative))
    stop("fLPIV: 'cumulative' must be TRUE or FALSE.")

  samp <- .fLP_estimation_sample(data, x_vars, lhs_vars, H,
                                 cumulative = cumulative, fn = "fLPIV")
  non_num <- c(names(samp$X)[!vapply(samp$X, is.numeric, logical(1))],
               names(samp$Y)[!vapply(samp$Y, is.numeric, logical(1))])
  if (length(non_num) > 0L)
    stop(sprintf("fLPIV: all variables must be numeric. Non-numeric: %s",
                 paste(unique(non_num), collapse = ", ")))

  raw <- samp$X
  Y <- as.matrix(samp$Y); storage.mode(Y) <- "double"
  D <- as.numeric(raw[[endog]])
  Z <- as.matrix(raw[, iv_vars, drop = FALSE]); storage.mode(Z) <- "double"
  C <- if (length(ctrl_vars) > 0L)
         as.matrix(raw[, ctrl_vars, drop = FALSE])
       else matrix(0.0, nrow = nrow(raw), ncol = 0)
  storage.mode(C) <- "double"

  # Cumulative long-difference base (see .fLP_estimation_sample); NULL when
  # no earlier row exists or the fit is not cumulative, NA when missing.
  Y_pre <- NULL
  if (!is.null(samp$pre)) {
    Y_pre <- as.matrix(samp$pre)
    storage.mode(Y_pre) <- "double"
  }

  T_eff <- nrow(raw)

  # -------------------------------------------------------------------
  # 5. Sanity / defaults
  # -------------------------------------------------------------------
  if (is.null(nw_lags_iv)) {
    nw_lags_iv <- max(0L, ncol(Z) - 1L)   # toolbox default
  } else {
    nw_lags_iv <- .fLP_validate_scalar_int(nw_lags_iv, "nw_lags_iv")
  }
  if (!is.logical(nw_prewhite) || length(nw_prewhite) != 1L || is.na(nw_prewhite))
    stop("fLPIV: 'nw_prewhite' must be TRUE or FALSE.")
  n_threads <- .fLP_validate_scalar_int(n_threads, "n_threads")

  # At least one residual degree of freedom in every first stage.
  nobs <- .fLP_horizon_nobs(Y, Y_pre, T_eff, H, cumulative)
  .fLP_check_nobs(nobs, ncol(C) + ncol(Z) + 3, lhs_vars, "fLPIV")

  # -------------------------------------------------------------------
  # 6. Call C++ engine (with conf[1]; bands rebuilt in R from irfs_se)
  # -------------------------------------------------------------------
  res_cpp <- fLPIV_cpp(
    Y          = Y,
    D          = D,
    Z          = Z,
    C          = C,
    H          = H,
    conf_level = as.double(conf[1L]) / 100,
    nw_lags_iv = as.integer(nw_lags_iv),
    cumulative = isTRUE(cumulative),
    n_threads  = as.integer(n_threads),
    verbose    = FALSE,
    Y_pre      = Y_pre,
    nw_prewhite = nw_prewhite
  )

  # Shock-size convention. LPmodel.m's impact = 0 does not touch the
  # treatment; it scales the 2SLS coefficient and its SE by d_norm =
  # std(treat) over the lag-trimmed regressor sample, the same at every
  # horizon (beta_h = beta_h*d_norm;
  # se_h = se_h*d_norm). Applied here, before the bands are built, so the
  # bands inherit the scaling.
  if (identical(shock_size, "sd")) {
    d_norm <- stats::sd(D)
    if (!is.finite(d_norm) || d_norm <= 0)
      stop("fLPIV: shock_size = \"sd\" requires a treatment with positive variance.")
    res_cpp$irfs    <- res_cpp$irfs    * d_norm
    res_cpp$irfs_se <- res_cpp$irfs_se * d_norm
  }

  # -------------------------------------------------------------------
  # 7. Labels and bands (as in fLP)
  # -------------------------------------------------------------------
  h_seq <- 0L:H
  dimnames(res_cpp$irfs)    <- list(as.character(h_seq), lhs_vars)
  dimnames(res_cpp$irfs_se) <- dimnames(res_cpp$irfs)
  dimnames(res_cpp$Fstat_fs) <- dimnames(res_cpp$irfs)
  dimnames(res_cpp$rsqr_fs)  <- dimnames(res_cpp$irfs)

  bands <- .fLP_bands(res_cpp$irfs, res_cpp$irfs_se, conf)
  res_cpp$irfs_upper <- bands$upper
  res_cpp$irfs_lower <- bands$lower

  res_cpp$lhs_vars        <- lhs_vars
  res_cpp$rhs_vars        <- rhs_vars
  res_cpp$endog           <- endog
  res_cpp$ctrl_vars       <- ctrl_vars
  res_cpp$instrument_vars <- iv_vars
  res_cpp$shock           <- endog                # for fPlotIRF labelling
  res_cpp$horizons        <- h_seq
  res_cpp$conf            <- conf
  res_cpp$nw_lags_iv      <- nw_lags_iv
  res_cpp$nw_prewhite     <- nw_prewhite
  res_cpp$cumulative      <- isTRUE(cumulative)
  res_cpp$shock_size      <- shock_size
  res_cpp$n_threads       <- as.integer(n_threads)
  res_cpp$nobs            <- .fLP_label_nobs(res_cpp$nobs, lhs_vars)
  res_cpp$formula         <- fmla_full
  res_cpp$instruments     <- iv_full
  res_cpp$formula_orig    <- formula
  res_cpp$call            <- match.call()

  class(res_cpp) <- c("fLPIV", "fLP")
  res_cpp
}


# =======================================================================
# S3 methods
# =======================================================================

#' @export
print.fLPIV <- function(x, digits = 4L, ...) {
  cat("\nLocal Projections \u2014 IV (fLPIV)\n")
  cat(strrep("-", 45L), "\n")
  cat("Original formula : ", deparse(x$formula_orig), "\n")
  cat("Expanded formula : ", deparse(x$formula),      "\n")
  cat("Instruments      : ", deparse(x$instruments),  "\n")
  cat("Endogenous       : ", x$endog,                 "\n")
  cat("Controls         : ", if (length(x$ctrl_vars)) paste(x$ctrl_vars, collapse = ", ") else "(none)", "\n")
  cat("LHS variables    : ", paste(x$lhs_vars, collapse = ", "), "\n")
  cat("Horizons         : 0 to", max(x$horizons),     "\n")
  cat("Cumulative       : ", isTRUE(x$cumulative),     "\n")
  cat("Confidence       : ",
      paste0(format(x$conf, trim = TRUE), "%", collapse = ", "), "\n")
  cat("HAC bandwidth    :  ")
  if (isTRUE(x$nw_lags_iv > 0)) {
    cat(sprintf("fixed = %d (vce(hac nw %d))\n", x$nw_lags_iv, x$nw_lags_iv))
  } else {
    cat("horizon-varying = h (just-identified rule)\n")
  }
  if (isTRUE(x$nw_prewhite)) cat("Prewhitening     :  VAR(1) (Andrews-Monahan)\n")
  cat("Observations     : ", .fLP_format_nobs(x$nobs), "\n")
  cat("\nIRF (endog = '", x$endog, "'):\n", sep = "")
  print(round(x$irfs, digits))
  if (!is.null(x$Fstat_fs)) {
    # One line when every outcome has the same first stage, else one each.
    Fs <- as.matrix(x$Fstat_fs)
    shared <- all(apply(Fs, 2, function(f) isTRUE(all.equal(f, Fs[, 1]))))
    for (j in if (shared) 1 else seq_len(ncol(Fs))) {
      lbl <- if (shared) "" else sprintf(" [%s]", colnames(Fs)[j])
      cat("\nFirst-stage F", lbl, " (h = 0..H):  ",
          paste(round(Fs[, j], 2), collapse = ", "), "\n", sep = "")
    }
  }
  invisible(x)
}
