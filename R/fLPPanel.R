# -----------------------------------------------------------------------
# fLPPanel.R — Panel Local Projections with fixest-like formula syntax
#
# Provides:
#   fLPPanel()       — main user-facing function
#   print.fLPPanel() — print method
#   coef.fLPPanel()  — extract IRF matrix
#   tidy.fLPPanel()  — broom-compatible tidy method
#
# Formula syntax (fixest-style):
#   fLPPanel(y ~ shock + control1 + control2 | id + time, data,
#            id = "id", time = "time", shock = "shock", horizons = 12)
#
#   - LHS: single response variable.
#   - RHS-main (before "|"): shock + additional controls.
#   - RHS-FE (after "|"): categorical fixed effects (any number).
#     Omit "|" for no fixed effects.
#
# Macros (like fixest ..vars):
#   controls <- c("gdp_l1", "inf_l1", "rate_l1")
#   fLPPanel(y ~ mp_shock + ..controls | id + time, ..., shock = "mp_shock")
#
# Panel-aware lag / lead operators (respect unit boundaries):
#   l(x, k)     lag x by k periods within each unit  → column "x_lk"
#   l(x, k1:k2) → columns "x_lk1", ..., "x_lk2"
#   f(x, k)     lead x by k periods within each unit → column "x_fk"
#   l(..vars, 1:3) works too.
#
# Heterogeneous responses (Almuzara-Sancibrián s ⊗ X):
#   fLPPanel(y ~ mp_shock | id + time, ..., shock = "mp_shock",
#            interact = c("size", "leverage"))
#   returns one IRF column per interaction variable.
#
# Compatible with fPlotIRF() via class = c("fLPPanel", "fLP").
#
# Author: Dr. Muhsin Ciftci
# -----------------------------------------------------------------------


# =======================================================================
# Internal helpers
# =======================================================================

# Normalize the raw time column to an integer grid where consecutive
# periods differ by 1. Preserves calendar gaps: missing periods become
# missing integer values, so `l(x, 1)` at t = 5 returns NA if t = 4 is
# unobserved. `on_grid` flags rows whose time does not sit on the grid
# and should be dropped by the caller.
.fLPPanel_normalize_time <- function(t_raw) {
  if (inherits(t_raw, "Date") || inherits(t_raw, "POSIXct"))
    t_raw <- as.numeric(t_raw)
  if (!is.numeric(t_raw))
    stop("fLPPanel: 'time' column must be numeric, Date, or POSIXct.",
         call. = FALSE)
  t_uniq <- sort(unique(t_raw))
  if (length(t_uniq) < 2L)
    stop("fLPPanel: 'time' has fewer than 2 unique values.", call. = FALSE)
  step     <- min(diff(t_uniq))
  if (!is.finite(step) || step <= 0)
    stop("fLPPanel: could not determine positive time step.", call. = FALSE)
  t_scaled <- (t_raw - t_uniq[1L]) / step + 1
  on_grid  <- abs(t_scaled - round(t_scaled)) < 1e-2
  list(t_norm = as.integer(round(t_scaled)), on_grid = on_grid)
}


# Parse a one-sided formula containing the panel unit, the panel time
# variable, or both. Clustering is deliberately restricted to the declared
# panel dimensions so its meaning remains unambiguous.
.fLPPanel_parse_cluster <- function(cluster, panel_id, env) {
  if (is.null(cluster)) {
    cluster_vars <- panel_id[2L]
    cluster_formula <- stats::as.formula(
      paste("~", cluster_vars), env = env
    )
  } else {
    if (!inherits(cluster, "formula") || length(cluster) != 2L) {
      stop(
        "fLPPanel: 'cluster' must be a one-sided formula, e.g. ~unit, ~tt, or ~unit + tt.",
        call. = FALSE
      )
    }
    cluster_formula <- cluster
    cluster_terms <- tryCatch(
      stats::terms(cluster_formula, keep.order = TRUE),
      error = function(e) stop(
        "fLPPanel: could not parse 'cluster': ", conditionMessage(e),
        call. = FALSE
      )
    )
    cluster_labels <- attr(cluster_terms, "term.labels")
    cluster_vars   <- all.vars(cluster_formula)

    # Require plain variable terms. This rejects interactions and
    # transformations such as ~unit:tt or ~factor(unit).
    if (length(cluster_labels) != length(cluster_vars) ||
        !setequal(cluster_labels, cluster_vars)) {
      stop(
        "fLPPanel: 'cluster' may contain only plain panel variables joined by '+'.",
        call. = FALSE
      )
    }
  }

  if (length(cluster_vars) < 1L || length(cluster_vars) > 2L ||
      anyDuplicated(cluster_vars)) {
    stop(
      "fLPPanel: 'cluster' must contain one or both panel variables.",
      call. = FALSE
    )
  }
  bad <- setdiff(cluster_vars, panel_id)
  if (length(bad)) {
    stop(sprintf(
      paste0(
        "fLPPanel: clustering variable(s) must be declared in 'panel_id'. ",
        "Invalid: %s; available: %s."
      ),
      paste(bad, collapse = ", "), paste(panel_id, collapse = ", ")
    ), call. = FALSE)
  }

  cluster_mode <- if (length(cluster_vars) == 2L) {
    2L
  } else if (identical(cluster_vars, panel_id[1L])) {
    1L
  } else {
    0L
  }

  list(
    formula = cluster_formula,
    vars    = cluster_vars,
    mode    = cluster_mode,
    label   = c("time", "unit", "unit + time")[cluster_mode + 1L]
  )
}


# Panel-aware lag/lead expander — respects unit boundaries AND preserves
# calendar-time gaps. `id_vec` and `t_norm` are the id and normalized
# integer time aligned with data rows, so shifts look up (id, t + k)
# rather than adjacent rows in sort order.
.fLPPanel_expand_lag_terms <- function(formula_chr, data, id_vec, t_norm) {

  T_dat <- nrow(data)
  if (length(id_vec) != T_dat || length(t_norm) != T_dat) {
    stop("fLPPanel: internal error - id/time vectors misaligned with data.",
         call. = FALSE)
  }

  # Row lookup keyed by (id, normalized time). Duplicates would collide;
  # they are rejected upstream so match() returns a single row per key.
  key_row <- paste(id_vec, t_norm, sep = "|")

  panel_lag <- function(x_raw, k, op) {
    t_target <- if (op == "l") (t_norm - k) else (t_norm + k)
    key_look <- paste(id_vec, t_target, sep = "|")
    idx <- match(key_look, key_row)
    out <- rep(NA_real_, T_dat)
    ok  <- !is.na(idx)
    out[ok] <- x_raw[idx[ok]]
    out
  }

  .fLP_expand_lag_terms(formula_chr, data, shift = panel_lag,
                        fn = "fLPPanel", paren = TRUE)
}


# =======================================================================
# Main function
# =======================================================================

#' Panel Local Projections (fixest-like syntax)
#'
#' Estimates panel impulse responses via the local-projection estimator of
#' Almuzara & Sancibrián (2024, NY Fed SR 1090). Formula syntax mirrors
#' \pkg{fixest}: RHS variables before \code{|} are the shock and controls,
#' variables after \code{|} are absorbed fixed effects. Interactions
#' (e.g. \code{shock:size}) are computed explicitly from the formula — the
#' estimator uses \emph{exactly} the specification the user writes.
#' Supports \code{..macro} expansion and panel-aware \code{l()} / \code{f()}
#' operators.
#'
#' @param formula A formula of the form
#'   \code{y ~ shock + controls | fe1 + fe2}, e.g.
#'   \code{y ~ shock:size + l(y, 1:2) | unit + tt}.
#'   \itemize{
#'     \item \strong{LHS}: single response variable.
#'     \item \strong{RHS-main} (before \code{|}): every term you want in the
#'       regression. Interactions with \code{:} are computed as elementwise
#'       products (this is where the Almuzara-Sancibrián \eqn{s_{it} \cdot X_t}
#'       appears — write \code{shock:size} in the formula).
#'     \item \strong{RHS-FE} (after \code{|}): categorical fixed-effect
#'       variables (omit \code{|} for no fixed effects).
#'     \item \strong{Macros}: \code{..name} where \code{name} is a
#'       character vector in the calling environment.
#'     \item \strong{Panel-aware lag/lead}: \code{l(x, k)}, \code{l(x, k1:k2)},
#'       \code{f(x, k)}, \code{l(..vars, 1:3)}. Never crosses unit boundaries.
#'   }
#' @param data A \code{data.frame} or \code{tibble} in long panel form.
#' @param panel_id Length-2 character vector \code{c(id_col, time_col)}
#'   naming the unit-id column and the time column, respectively.
#' @param horizons Integer. Maximum horizon H (projections run for
#'   \eqn{h = 0, 1, \ldots, H}) or the full sequence \code{0:H}. Default 12.
#' @param shock Character (single or vector). The RHS-main term label(s)
#'   whose coefficient is the impulse response. Everything else in the
#'   RHS-main is treated as a control. Must match the term label exactly
#'   (after macro / lag expansion); interactions are labelled with
#'   \code{":"}, e.g. \code{"shock:size"}.
#' @param conf Numeric in (0, 100), in percent. One or more confidence levels,
#'   e.g. \code{90} or \code{c(68, 95)}. Bands are rebuilt from the raw
#'   SEs, so multi-level fits cost the same as single-level fits.
#' @param ssc Character. Small-sample correction applied to the
#'   cluster-robust variance:
#'   \describe{
#'     \item{\code{"none"}}{(default) the asymptotic sandwich with normal
#'       critical values (\code{df = Inf}). This is the convention of the
#'       Almuzara-Sancibrian panel LP theory, which is asymptotic in both
#'       the cross-section and the time dimension.}
#'     \item{\code{"CR1"}}{the correction used by Stata and \pkg{fixest}:
#'       every one-way component of the meat is multiplied by
#'       \eqn{G/(G-1) \cdot (N-1)/(N-K)} with its own number of clusters
#'       \eqn{G}, and inference uses \eqn{t} with \eqn{\min(G) - 1}
#'       degrees of freedom. \eqn{K} counts the regressors plus the
#'       absorbed fixed effects that are not nested in a clustering
#'       dimension. Standard errors are slightly wider than
#'       \code{"none"}; the difference is negligible with many clusters
#'       and matters when there are few. Available for every clustering
#'       dimension.}
#'     \item{\code{"IK"}}{the Imbens-Kolesar (2016, REStat) refinement, a
#'       bias-reduced (CR2-type) variance with Satterthwaite degrees of
#'       freedom. Defined only for time clustering.}
#'   }
#'   Point estimates are identical under all three; only the standard
#'   errors and the critical values change.
#' @param small_sample Deprecated. \code{TRUE} is equivalent to
#'   \code{ssc = "IK"}; use \code{ssc} instead.
#' @param cumulative Logical. If \code{TRUE}, project
#'   \eqn{\sum_{k=0}^h y_{i,t+k}} instead of \eqn{y_{i,t+h}}. Default
#'   \code{FALSE}.
#' @param p_max Integer. Maximum number of lags of the regressand and of
#'   the shock terms added as controls (at horizon h the effective number
#'   of lags is \eqn{\min(h, p_{\max})}). Default 0.
#' @param n_threads Integer. OpenMP threads for the horizon loop. 0 = all
#'   cores minus one. Default 0.
#' @param cluster One-sided formula selecting the clustering dimension:
#'   the unit variable, the time variable, or both variables declared in
#'   \code{panel_id}; for example \code{~unit}, \code{~tt}, or
#'   \code{~unit + tt}. Two-way clustering uses the
#'   Cameron-Gelbach-Miller inclusion-exclusion estimator. The default
#'   \code{NULL} is equivalent to clustering on the time variable. The
#'   \code{ssc = "IK"} refinement is available only with time clustering.
#'
#' @return An object of class \code{c("fLPPanel", "fLP")} that plugs into
#'   \code{\link{fPlotIRF}} and \code{\link{tidy.fLP}}.
#'
#' @examples
#' data("DemocracyGrowth")
#'
#' # Response of log GDP to a democratisation event, unit and year fixed
#' # effects, clustered on countries.
#' fLPPanel(lgdp ~ dem | id + year,
#'          data = DemocracyGrowth, panel_id = c("id", "year"),
#'          shock = "dem", horizons = 4, cluster = ~id, n_threads = 2)
#'
#' \donttest{
#' # The Stata / fixest finite-sample correction
#' fLPPanel(lgdp ~ dem | id + year,
#'          data = DemocracyGrowth, panel_id = c("id", "year"),
#'          shock = "dem", horizons = 12, cluster = ~id,
#'          ssc = "CR1", n_threads = 2)
#'
#' # Cumulative response with lags of the outcome as controls
#' fLPPanel(lgdp ~ dem + l(lgdp, 1:2) | id + year,
#'          data = DemocracyGrowth, panel_id = c("id", "year"),
#'          shock = "dem", horizons = 12, cumulative = TRUE,
#'          cluster = ~id + year, n_threads = 2)
#' }
#'
#' @export
fLPPanel <- function(formula, data,
                     panel_id,
                     horizons     = 12L,
                     shock,
                     conf         = 90,
                     ssc          = c("none", "CR1", "IK"),
                     small_sample = FALSE,
                     cumulative   = FALSE,
                     p_max        = 0L,
                     n_threads    = 0L,
                     cluster      = NULL) {

  env <- parent.frame()

  # ---- required arguments ----------------------------------------------
  if (missing(shock))    stop("fLPPanel: 'shock' is required.")
  if (missing(panel_id)) stop("fLPPanel: 'panel_id' is required (c(id_col, time_col)).")

  if (!is.character(panel_id) || length(panel_id) != 2L || anyNA(panel_id))
    stop("fLPPanel: 'panel_id' must be a length-2 character vector, e.g. c('unit', 'tt').")
  id_var   <- panel_id[1L]
  time_var <- panel_id[2L]

  if (!is.character(shock) || length(shock) == 0L || anyNA(shock))
    stop("fLPPanel: 'shock' must be a non-empty character vector.")

  H    <- .fLP_validate_horizon_max(horizons)
  conf <- .fLP_validate_conf_vector(conf)

  if (!is.logical(small_sample) || length(small_sample) != 1L || is.na(small_sample))
    stop("fLPPanel: 'small_sample' must be TRUE or FALSE.")
  ssc_given <- !missing(ssc)          # before match.arg(), which assigns ssc
  ssc <- match.arg(ssc)
  if (isTRUE(small_sample)) {
    if (ssc_given && !identical(ssc, "IK"))
      stop("fLPPanel: 'small_sample' and 'ssc' disagree; set 'ssc' only.",
           call. = FALSE)
    .Deprecated(msg = paste(
      "fLPPanel: 'small_sample = TRUE' is deprecated and will be removed in",
      "the next release.\nUse ssc = \"IK\" instead."))
    ssc <- "IK"
  }
  ssc_mode <- match(ssc, c("none", "CR1", "IK")) - 1L
  if (!is.logical(cumulative)   || length(cumulative)   != 1L)
    stop("fLPPanel: 'cumulative' must be TRUE or FALSE.")

  p_max     <- .fLP_validate_scalar_int(p_max, "p_max")
  n_threads <- .fLP_validate_scalar_int(n_threads, "n_threads")

  if (!is.data.frame(data)) data <- as.data.frame(data)

  if (!id_var   %in% names(data))
    stop(sprintf("fLPPanel: id column '%s' not found in data.", id_var))
  if (!time_var %in% names(data))
    stop(sprintf("fLPPanel: time column '%s' not found in data.", time_var))

  cluster_spec <- .fLPPanel_parse_cluster(cluster, panel_id, env)
  if (identical(ssc, "IK") && cluster_spec$mode != 0L) {
    stop(
      paste0(
        "fLPPanel: ssc = \"IK\" is defined only for time clustering; ",
        "use cluster = ~", time_var, ", or ssc = \"CR1\" for the ",
        "Stata / fixest correction with this clustering dimension."
      ),
      call. = FALSE
    )
  }

  # ---- normalize time first, so gap-aware lag/lead can look up
  #      (id, t + k) rather than (id, adjacent row). Drops off-grid rows
  #      up-front so downstream helpers, dedup, and the C++ engine all see
  #      the same integer time index.
  tn      <- .fLPPanel_normalize_time(data[[time_var]])
  if (any(!tn$on_grid)) {
    message(sprintf(
      "fLPPanel: dropping %d row(s) with off-grid time index.",
      sum(!tn$on_grid)))
    data    <- data[tn$on_grid, , drop = FALSE]
    tn$t_norm <- tn$t_norm[tn$on_grid]
  }
  id_vec  <- data[[id_var]]
  t_norm  <- tn$t_norm

  # Duplicate (id, time) rows would silently corrupt gap-aware shifts and
  # the C++ engine's per-unit time->row map — reject them explicitly.
  dup <- duplicated(paste(id_vec, t_norm, sep = "|"))
  if (any(dup))
    stop(sprintf(
      "fLPPanel: found %d duplicate (%s, %s) row(s); each unit-time pair must be unique.",
      sum(dup), id_var, time_var), call. = FALSE)

  # ---- formula → string -------------------------------------------------
  formula_chr <- paste(deparse(formula), collapse = " ")

  # Macro helpers from fLP (pure string operations)
  formula_chr <- .fLP_expand_macro_lag(formula_chr, env)
  formula_chr <- .fLP_expand_macro(formula_chr, env)

  # ---- split on "|" ----------------------------------------------------
  parts    <- strsplit(formula_chr, "\\|", fixed = FALSE)[[1]]
  main_str <- trimws(parts[1])
  fe_str   <- if (length(parts) > 1L) trimws(paste(parts[-1L], collapse = "|")) else ""

  # ---- panel-aware lag/lead expansion (main part only) ------------------
  main_ex  <- .fLPPanel_expand_lag_terms(main_str, data, id_vec, t_norm)
  main_str <- main_ex$formula_chr
  data     <- main_ex$data

  # ---- parse main formula ----------------------------------------------
  main_formula <- tryCatch(
    as.formula(main_str, env = env),
    error = function(e) stop(
      "fLPPanel: could not parse the main formula:\n  ", main_str,
      "\n  ", conditionMessage(e)
    )
  )
  lhs_expr <- main_formula[[2L]]
  rhs_expr <- main_formula[[3L]]

  lhs_vars <- all.vars(lhs_expr)
  if (length(lhs_vars) != 1L)
    stop("fLPPanel: exactly one LHS variable is required.")

  rhs_str   <- paste(deparse(rhs_expr), collapse = " ")
  rhs_terms <- terms(as.formula(paste("~", rhs_str), env = env),
                     keep.order = TRUE)
  rhs_labels <- attr(rhs_terms, "term.labels")
  # Raw base columns (e.g. 'shock', 'size' for term 'shock:size')
  rhs_base_vars <- all.vars(rhs_terms)

  # ---- validate shock terms -------------------------------------------
  bad_shocks <- setdiff(shock, rhs_labels)
  if (length(bad_shocks))
    stop(sprintf(
      "fLPPanel: shock term(s) not found in RHS-main: %s. Available: %s",
      paste(bad_shocks, collapse = ", "),
      paste(rhs_labels, collapse = ", ")))

  control_labels <- setdiff(rhs_labels, shock)

  # ---- parse FE variables ----------------------------------------------
  if (nzchar(fe_str)) {
    fe_expr <- tryCatch(
      as.formula(paste("~", fe_str), env = env),
      error = function(e) stop(
        "fLPPanel: could not parse the fixed-effects part '", fe_str,
        "': ", conditionMessage(e)
      )
    )
    fe_vars <- attr(terms(fe_expr, keep.order = TRUE), "term.labels")
  } else {
    fe_vars <- character(0L)
  }

  # ---- collect all needed base columns and drop NA rows ---------------
  all_needed <- unique(c(lhs_vars, rhs_base_vars, fe_vars, id_var, time_var))
  miss <- setdiff(all_needed, names(data))
  if (length(miss))
    stop(sprintf("fLPPanel: columns not found in data: %s",
                 paste(miss, collapse = ", ")))

  raw <- data[, all_needed, drop = FALSE]
  na_rows <- which(!complete.cases(raw))
  if (length(na_rows)) {
    message(sprintf("fLPPanel: dropping %d row(s) with NA.", length(na_rows)))
    raw <- raw[-na_rows, , drop = FALSE]
  }

  # ---- numeric checks on all vars used in the model matrix -------------
  for (nm in unique(c(lhs_vars, rhs_base_vars))) {
    if (!is.numeric(raw[[nm]]))
      stop(sprintf("fLPPanel: column '%s' must be numeric.", nm))
  }

  # ---- build the RHS design matrix from the formula --------------------
  # One column per term label (interactions computed as elementwise
  # products). No intercept — HDFE absorbs it.
  design_formula <- as.formula(paste("~", rhs_str, "- 1"), env = env)
  design <- stats::model.matrix(design_formula, data = raw)

  # model.matrix may rename term labels (e.g. reorder factors). Since all
  # our terms are numeric, colnames should match rhs_labels one-for-one,
  # in order.
  if (ncol(design) != length(rhs_labels)) {
    stop("fLPPanel: internal error \u2014 design matrix has ",
         ncol(design), " columns but ", length(rhs_labels), " term labels.")
  }
  colnames(design) <- rhs_labels

  shock_cols   <- match(shock, rhs_labels)
  control_cols <- setdiff(seq_along(rhs_labels), shock_cols)

  y_mat <- matrix(as.double(raw[[lhs_vars]]), ncol = 1L)
  W_mat <- if (length(control_cols)) design[, control_cols, drop = FALSE]
           else matrix(0, nrow(raw), 0L)
  storage.mode(W_mat) <- "double"

  # ---- shock / heterogeneity split -----------------------------------
  # Convention:
  #   * Bare shock term ("shock")  → X = data[["shock"]], s = ones
  #   * Single interaction "A:B"   → X = data[[A]], s = data[[B]]
  #     (engine forms s ⊗ X internally; mathematically identical to
  #      passing X = A * B and s = ones for point estimates and the
  #      time-clustered LAHR SE, and matches the paper's Imbens-Kolesar
  #      refinement — which uses the heterogeneity vector s as the
  #      per-period weighting variable)
  #   * Multiple shock terms       → X = design's shock cols, s = ones
  is_pairwise <- function(lbl) {
    p <- strsplit(lbl, ":", fixed = TRUE)[[1]]
    length(p) == 2L && all(p %in% names(raw)) &&
      all(vapply(p, \(nm) is.numeric(raw[[nm]]), logical(1L)))
  }

  if (length(shock) == 1L && is_pairwise(shock)) {
    parts  <- strsplit(shock, ":", fixed = TRUE)[[1]]
    X_mat  <- matrix(as.double(raw[[parts[1L]]]), ncol = 1L)
    s_mat  <- matrix(as.double(raw[[parts[2L]]]), ncol = 1L)
  } else {
    X_mat  <- design[, shock_cols, drop = FALSE]
    storage.mode(X_mat) <- "double"
    s_mat  <- matrix(1.0, nrow(raw), 1L)
  }

  fe_mat <- if (length(fe_vars)) {
    m <- matrix(0L, nrow(raw), length(fe_vars))
    for (j in seq_along(fe_vars))
      m[, j] <- as.integer(factor(raw[[fe_vars[j]]]))
    m
  } else {
    matrix(0L, nrow(raw), 0L)
  }

  # Time was already normalized above (before formula expansion) so gaps
  # are preserved. Rebuild the aligned integer time from the same data
  # subset after NA-drop and off-grid drop.
  i_index <- as.integer(factor(raw[[id_var]]))
  t_index <- .fLPPanel_normalize_time(raw[[time_var]])$t_norm

  # ---- call C++ engine --------------------------------------------------
  res <- fLPPanel_cpp(
    y            = y_mat,
    s            = s_mat,
    X            = X_mat,
    W            = W_mat,
    FE           = fe_mat,
    i_index      = i_index,
    t_index      = t_index,
    H            = H,
    p_max        = p_max,
    ssc_mode     = ssc_mode,
    cumulative   = isTRUE(cumulative),
    cluster_mode = cluster_spec$mode,
    n_threads    = n_threads,
    verbose      = FALSE
  )

  # ---- build fLP-shaped result -----------------------------------------
  h_seq <- 0L:H

  # One IRF column per shock term (label = term label from the formula).
  irf_names <- shock

  dimnames(res$estimate) <- list(as.character(h_seq), irf_names)
  dimnames(res$SE)       <- dimnames(res$estimate)
  dimnames(res$df)       <- dimnames(res$estimate)

  # Bands at each requested conf level using per-cell t quantiles. The C++
  # engine no longer computes CI90/95/99 or p-values (they were the only
  # R-API calls inside the OpenMP loop); they are rebuilt here.
  t_crit <- function(cl)
    array(stats::qt(0.5 * (1 + cl / 100), df = res$df), dim(res$df))
  fixed_ci <- function(cl) {
    b <- .fLP_bands(res$estimate, res$SE, cl, t_crit)
    a <- array(NA_real_, c(H + 1L, length(irf_names), 2L))
    a[, , 1L] <- b$lower
    a[, , 2L] <- b$upper
    a
  }
  CI90 <- fixed_ci(90)
  CI95 <- fixed_ci(95)
  CI99 <- fixed_ci(99)
  # Two-sided p-values from |t| under t(df).
  pval <- array(2 * stats::pt(abs(res$estimate / res$SE), df = res$df,
                              lower.tail = FALSE), dim(res$df))

  bands      <- .fLP_bands(res$estimate, res$SE, conf, t_crit)
  irfs_upper <- bands$upper
  irfs_lower <- bands$lower

  out <- list(
    irfs         = res$estimate,
    irfs_se      = res$SE,
    irfs_upper   = irfs_upper,
    irfs_lower   = irfs_lower,
    df           = res$df,
    pval         = pval,
    CI90         = CI90,
    CI95         = CI95,
    CI99         = CI99,
    status       = res$status,
    converged    = res$converged,

    lhs_vars     = lhs_vars,
    rhs_vars     = rhs_labels,
    fe_vars      = fe_vars,
    shock        = shock,
    id           = id_var,
    time         = time_var,
    panel_id     = panel_id,
    cluster      = cluster_spec$formula,
    cluster_vars = cluster_spec$vars,
    cluster_type = cluster_spec$label,
    horizons     = h_seq,
    conf         = conf,
    ssc          = ssc,
    small_sample = identical(ssc, "IK"),   # kept for backward compatibility
    cumulative   = isTRUE(cumulative),
    p_max        = p_max,
    n_threads    = as.integer(n_threads),
    nobs         = nrow(raw),
    formula      = main_formula,
    formula_orig = formula,
    call         = match.call()
  )

  class(out) <- c("fLPPanel", "fLP")
  out
}


# =======================================================================
# S3 methods
# =======================================================================

#' @export
print.fLPPanel <- function(x, digits = 4L, ...) {
  cat("\nPanel Local Projections (fLPPanel)\n")
  cat(strrep("-", 45L), "\n")
  cat("Original formula : ", deparse(x$formula_orig), "\n")
  cat("Panel index      :  id = ", x$id, ", time = ", x$time, "\n", sep = "")
  cat("Shock term(s)    : ", paste(x$shock, collapse = ", "), "\n")
  cat("Response         : ", x$lhs_vars, "\n")
  if (length(x$fe_vars))
    cat("Fixed effects    : ", paste(x$fe_vars, collapse = " + "), "\n")
  else
    cat("Fixed effects    :  (none)\n")
  cat("Horizons         :  0 to ", max(x$horizons), "\n", sep = "")
  cat("Cumulative       : ", isTRUE(x$cumulative), "\n")
  cat("Clustering       : ", paste(x$cluster_vars, collapse = " + "), "\n")
  cat(
    "SE               : ",
    if (identical(x$ssc, "IK")) {
      "Imbens-Kolesar (2016) time-clustered small-sample"
    } else if (identical(x$ssc, "CR1")) {
      paste0("CR1 (Stata / fixest) ", x$cluster_type, " cluster-robust")
    } else if (identical(x$cluster_type, "time")) {
      "asymptotic time-clustered (LAHR)"
    } else {
      paste0("asymptotic ", x$cluster_type, " cluster-robust")
    },
    "\n"
  )
  cat("Confidence       : ",
      paste0(format(x$conf, trim = TRUE), "%", collapse = ", "), "\n")
  cat("Observations     : ", x$nobs, "\n")
  cat("Lag order p_max  : ", x$p_max, "\n")
  cat("\nIRF (shock = '", x$shock, "'):\n", sep = "")
  print(round(x$irfs, digits))
  invisible(x)
}

#' @export
coef.fLPPanel <- function(object, ...) object$irfs

#' Tidy Panel Local Projection Estimates
#'
#' @param x An object returned by \code{\link{fLPPanel}}.
#' @param ... Reserved for additional arguments.
#' @return A data frame with horizon, shock, estimate, standard error,
#'   degrees of freedom, p-value, and confidence bounds. Multiple confidence
#'   levels use columns named \code{lower_<level>} and \code{upper_<level>}.
#' @export
tidy.fLPPanel <- function(x, ...) {
  h_seq <- 0L:max(x$horizons)
  base <- data.frame(
    horizon  = rep(h_seq, times = ncol(x$irfs)),
    shock    = rep(colnames(x$irfs), each = length(h_seq)),
    estimate = as.vector(x$irfs),
    se       = as.vector(x$irfs_se),
    df       = as.vector(x$df),
    pval     = as.vector(x$pval),
    stringsAsFactors = FALSE
  )
  .fLP_tidy_bands(base, x)
}
