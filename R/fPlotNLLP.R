# -----------------------------------------------------------------------
# fPlotNLLP.R — plot nonlinear / state-dependent LP responses
#
#   fPlotNLLP(fit)                       one fit, one panel per LHS variable
#   fPlotNLLP(fit_gdp, fit_p, fit_ffr)   several fits side by side
#
#   what = "components" : both shock coefficients, one colour each
#   what = "difference" : c1 - c2 with its HAC band (state and sign only)
#
# Author: Dr. Muhsin Ciftci
# -----------------------------------------------------------------------

#' @export
fPlotNLLP <- function(...,
                      variables       = NULL,
                      what            = c("components", "difference"),
                      scale           = 1,
                      return_data     = FALSE,
                      colors          = c("#910048", "#407EC9"),
                      ribbon_alpha    = 0.2,
                      zero_line_color = "#707372",
                      facet_scales    = "free_y",
                      facet_ncol      = NULL) {
  fits <- list(...)
  if (length(fits) == 0 || !all(vapply(fits, inherits, logical(1), "fNLLP")))
    stop("fPlotNLLP: pass one or more objects returned by fNLLP().", call. = FALSE)
  what <- match.arg(what)

  x <- fits[[1]]
  for (f in fits[-1]) {
    if (!identical(f$type, x$type) || !identical(f$components, x$components) ||
        !identical(f$conf, x$conf))
      stop("fPlotNLLP: all fits must have the same type, components and conf.",
           call. = FALSE)
  }
  if (what == "difference" && identical(x$type, "cubic"))
    stop("fPlotNLLP: type = \"cubic\" has no difference; plot the components ",
         "or use predict() for responses to shocks of given sizes.", call. = FALSE)
  if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale))
    stop("fPlotNLLP: 'scale' must be a finite numeric scalar.", call. = FALSE)

  all_vars <- unlist(lapply(fits, `[[`, "lhs_vars"))
  if (anyDuplicated(all_vars))
    stop("fPlotNLLP: the same LHS variable appears in more than one fit.",
         call. = FALSE)

  varnames <- all_vars
  if (!is.null(variables)) {
    missing_vars <- setdiff(variables, varnames)
    if (length(missing_vars) > 0)
      stop(sprintf("fPlotNLLP: variable(s) not found in the fit(s): %s",
                   paste(missing_vars, collapse = ", ")), call. = FALSE)
    varnames <- variables
  }
  if (is.null(facet_ncol)) {
    facet_ncol <- if (length(varnames) <= 4) length(varnames) else
      ceiling(sqrt(length(varnames)))
  }

  td <- do.call(rbind, lapply(fits, tidy.fNLLP,
                              difference = (what == "difference")))
  td <- td[td$lhs %in% varnames, , drop = FALSE]
  td <- if (what == "difference") td[td$component == "difference", , drop = FALSE] else
    td[td$component != "difference", , drop = FALSE]
  td$component <- droplevels(td$component)

  # Long format over confidence levels, widest first.
  keys  <- format(x$conf, trim = TRUE)
  multi <- length(keys) > 1
  bands <- lapply(keys, function(k) {
    lo <- if (multi) paste0("lower_", k) else "lower"
    up <- if (multi) paste0("upper_", k) else "upper"
    tibble::tibble(
      variable  = factor(td$lhs, levels = varnames),
      component = td$component,
      horizon   = td$horizon,
      conf      = as.numeric(k),
      point     = td$estimate * scale,
      # A negative scale flips the sign; pmin/pmax keep lower <= upper.
      lower     = pmin(td[[lo]] * scale, td[[up]] * scale),
      upper     = pmax(td[[lo]] * scale, td[[up]] * scale)
    )
  })
  plot_data <- do.call(rbind, bands)

  if (return_data) {
    if (!multi) plot_data$conf <- NULL
    return(plot_data)
  }

  comp_levels <- levels(plot_data$component)
  pal <- rep_len(colors, length(comp_levels))
  names(pal) <- comp_levels

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$horizon))
  for (i in seq_along(keys)) {
    band_k <- plot_data[plot_data$conf == as.numeric(keys[i]), , drop = FALSE]
    p <- p + ggplot2::geom_ribbon(
      data = band_k,
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper,
                   fill = .data$component),
      alpha = ribbon_alpha
    )
  }

  p +
    ggplot2::geom_line(
      data = plot_data[plot_data$conf == as.numeric(keys[1]), , drop = FALSE],
      ggplot2::aes(y = .data$point, color = .data$component),
      linewidth = 0.8
    ) +
    ggplot2::geom_hline(yintercept = 0, color = zero_line_color,
                        linetype = "dashed", linewidth = 0.6) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::facet_wrap(~ variable, scales = facet_scales, ncol = facet_ncol) +
    ggplot2::labs(x = "Horizon", y = NULL, color = NULL, fill = NULL) +
    ggplot2::theme(legend.position = if (what == "difference") "none" else "bottom")
}
