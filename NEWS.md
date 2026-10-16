# tidyMacro 1.1.0

## Breaking changes

* One plot function replaces the per-scheme plotters. `fPlotIRF()` draws the
  impulse responses of every identification scheme, and `fIRF()` converts any
  result into one tidy format (`model`, `shock`, `variable`, `horizon`,
  `band`, `point`, `lower`, `upper`). Several results can be overlaid in one
  figure, which was not possible before:

  ```r
  fPlotIRF(point_irf, bands = boot_chol, varnames = vn, shock = 1)
  fPlotIRF(Cholesky = fIRF(point_irf, bands = boot_chol),
           `External IV` = fIRF(boot_iv), varnames = vn, shock = 1)
  ```

  `fPlotIRFChol()`, `fPlotIRFBQ()`, `fPlotIRFLR()`, `fPlotIRFIV()`,
  `fPlotIRFHetero()`, `fPlotIRFMSW()`, `fPlotIRFSign()` and `fPlotLP()` have
  been removed. The plotted numbers are unchanged.

* `fBootstrapIVMBB()` now returns `point`, the impulse responses estimated on
  the original sample (the Wold coefficients times the proxy-SVAR impact
  vector), and this is what `fPlotIRF()` draws. It previously returned
  `meanirf` and `medianirf`, summaries of the bootstrap draws, and plotted
  the mean. The mean inherits the small-sample bias of the resampled VARs and
  can move visibly under extreme draws, so it is not the estimate reported in
  the source papers. `meanirf` and `medianirf` have been removed here and in
  `fBootstrapHetero()`, which already reported `point`.

* `fLPPanel(small_sample = TRUE)` is deprecated in favour of `ssc = "IK"` and
  will be removed in the next release.

## New features

* `fLPPanel(ssc = )` selects the small-sample correction of the
  cluster-robust variance:
  * `"none"` (default, unchanged): the asymptotic sandwich with normal
    critical values, the convention of the Almuzara-Sancibrian theory.
  * `"CR1"`: the correction used by Stata and **fixest**. Each one-way
    component of the meat is scaled by `G/(G-1) * (N-1)/(N-K)` with its own
    cluster count, and inference uses t with `min(G) - 1` degrees of
    freedom. `K` counts the regressors, the intercept, and `L - 1`
    parameters for every absorbed factor with `L` levels that is not nested
    in a clustering dimension. Available for time, unit and two-way
    clustering.
  * `"IK"`: the Imbens-Kolesar (2016) refinement, as before, for time
    clustering only.

* Every bootstrap function takes `verbose`, which defaults to `FALSE`. Thread
  and bias-correction progress is no longer printed unless asked for.

* `fBootstrapIVInvertible()` and `fBootstrapIVRecover()` take `n_threads`;
  they previously always used all cores minus one.

## Robustness

* A bootstrap draw that fails — a rank-deficient resample, a singular
  covariance, or a degenerate instrument — is now recorded and dropped after
  the parallel region instead of throwing from inside it, which could abort
  the R session. The bands are computed from the surviving draws and a
  warning reports how many were dropped; an error is raised only when every
  draw fails.

* `fVAR()` validates its arguments: `p`, `c`, a finite and long enough `y`,
  matching `exog` rows, and positive residual degrees of freedom. The
  degrees-of-freedom check in the bootstrap entry points now uses the actual
  design width, which it previously computed ignoring the intercept and any
  exogenous regressors.

* `nboot`, coverage levels, `policyvar` and the dimensions of `instr`,
  `noise` and `delta` are checked before any memory is allocated.

## Internal

* `clean_instr_inv()` and `clean_instr_rec()`, which were identical, are now
  one shared helper.
* `fSR_cpp()` builds its result list in two steps, so the package compiles
  with Rcpp versions that cap `List::create()` at 20 elements.
