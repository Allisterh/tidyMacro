#ifndef FBOOTSTRAPVAR_H
#define FBOOTSTRAPVAR_H

#include <RcppArmadillo.h>
#include <vector>
#include <R_ext/Print.h>
#include "fVAR.h"  // Include to access VARResult struct
#include "fWoldIRF.h"
#include "omp_threads.h"

// Resampling plan for the residual/wild VAR bootstrap. It is drawn on the R
// main thread before any OpenMP loop (RcppArmadillo's RNG is R's RNG, which
// is main-thread only), so workers only read it and results do not depend on
// the thread count.
//   residual: idx.col(b)   = T-p residual row indices for draw b
//   wild    : signs.col(b) = Rademacher sign for each residual row
struct BootstrapPlan {
  bool       wild;
  arma::umat idx;
  arma::mat  signs;
};

BootstrapPlan fBootstrapPlan_cpp(const std::string& bootscheme,
                                 int T_iter, int T_resid, int nboot);

// Simulate T x N data from the fitted VAR, driven by a (T-p) x N shock path.
// The first p rows are copied from y. When exog_ptr is not null, the
// exogenous contribution (rows after the lag block of beta) is added.
arma::mat fSimulateVAR_cpp(const arma::mat& y,
                           const VARResult& var_result,
                           const arma::mat& shocks,
                           const arma::mat* exog_ptr);

// Bootstrap sample b of a plan: resampled shocks fed to fSimulateVAR_cpp.
arma::mat fBootstrapVAR_cpp(const arma::mat& y,
                            const VARResult& var_result,
                            const BootstrapPlan& plan, int b,
                            const arma::mat* exog_ptr);

// Mean bootstrap coefficient matrix (N x n_coef) over nboot residual/wild
// draws. Pass 1 of the bias-corrected bootstraps. Draws are averaged in
// replication order, so the result does not depend on the thread count.
arma::mat fBootstrapMeanBeta_cpp(const arma::mat& y,
                                 const VARResult& var_result,
                                 int nboot, const std::string& bootscheme,
                                 const arma::mat* exog_ptr, int n_threads);

// Bias-corrected VAR for the two-pass ("Corrected") bootstraps: pass 1 mean
// of nboot1 draws, then fRemoveBias_cpp. Prints the pass-1 summary.
VARResult fBiasCorrectedVAR_cpp(const arma::mat& y,
                                const VARResult& var_result,
                                int nboot1, const std::string& bootscheme,
                                const arma::mat* exog_ptr, int n_threads,
                                bool verbose);

// Lower Cholesky factor of a bootstrap covariance. A draw that is only
// numerically indefinite is rescued with a ridge that grows until the
// factorisation succeeds; `ok` reports whether it ever did, so the caller
// can drop the draw. Never throws, so it is safe inside an OpenMP region.
inline arma::mat fLowerChol_cpp(const arma::mat& sigma, bool* ok = nullptr) {
  arma::mat S;
  if (arma::chol(S, sigma, "lower")) {
    if (ok) *ok = true;
    return S;
  }
  arma::mat reg = 0.5 * (sigma + sigma.t());
  const double scale = std::max(arma::trace(reg) / reg.n_rows, 1e-12);
  for (double eps = 1e-8; eps <= 1e-2; eps *= 100.0) {
    arma::mat tryreg = reg;
    tryreg.diag() += eps * scale;
    if (arma::chol(S, tryreg, "lower")) {
      if (ok) *ok = true;
      return S;
    }
  }
  if (ok) *ok = false;
  S.zeros(sigma.n_rows, sigma.n_cols);
  return S;
}

// Argument checks shared by the identified-VAR bootstraps (main thread).
void fCheckBootArgs_cpp(const arma::mat& y, const VARResult& var_result,
                        int nboot, int horizon, double conf, double conf2);

// Minimal checks for the bootstraps that do not re-estimate a VAR first.
void fCheckDrawArgs_cpp(int nboot, int horizon, double conf, double conf2);

// Keeps the columns of `flat` whose draw both succeeded and is finite.
// A draw can fail on resampled data when the instrument is weak or the VAR
// is near-singular; workers record that in `ok` without touching the R API,
// and this runs after the parallel region. Stops when no draw survives and
// warns when some are dropped, so the bands are never built on NaN.
arma::mat fKeepValidDraws_cpp(const arma::mat& flat,
                              const std::vector<unsigned char>& ok,
                              const char* label);

// Same rule for bootstraps that store one N x H slice per draw. A failed
// draw would otherwise leave a slice of zeros and pull the bands towards
// zero, which is worse than an obvious NaN.
arma::cube fKeepValidSlices_cpp(const arma::cube& draws,
                                const std::vector<unsigned char>& ok,
                                const char* label);

// Residual of the instrument after projecting out p of its own lags and p
// lags of y: the pre-whitening step shared by the invertibility and
// recoverability bootstraps.
arma::vec fCleanInstrument_cpp(const arma::vec& instr_raw, int p,
                               const arma::mat& y);

// ------------------------------------------------------------------------- //
//  Identified-VAR bootstrap driver (Cholesky, BQ, LR-max, Uhlig max-share)  //
// ------------------------------------------------------------------------- //
struct SVARBootResult {
  arma::mat  flat;       // stat_size x nboot identified responses per draw
  arma::mat  bands;      // stat_size x 5: upper, lower, upper2, lower2, median
  arma::cube boot_beta;  // N x n_coef x nboot bootstrapped coefficients
  arma::vec  conf;       // band levels in percent: conf, conf2
  std::vector<unsigned char> ok;  // per-draw success; see fKeepValidDraws_cpp
};

// Bands of an SVARBootResult (fills r.bands from the draws that succeeded).
void fSVARBootBands_cpp(SVARBootResult& r, double conf, double conf2,
                        int n_threads, const char* label);

// R list {<draws_name>, upper, lower, upper2, lower2, boot_beta}; `dims`
// gives the dimensions of one draw (e.g. N, N, H or N, H).
Rcpp::List fSVARBootList_cpp(const SVARBootResult& r, const char* draws_name,
                             const std::vector<int>& dims);

// For each draw: resample -> simulate -> re-estimate the VAR -> Wold IRFs ->
// identify(var_b, wold_b, out), where `out` is the draw's stat_size column.
// `identify` runs on OpenMP workers and must not touch the R API.
template <class Identify>
SVARBootResult fBootstrapSVAR_cpp(const arma::mat& y,
                                  const VARResult& var_result,
                                  int nboot, int horizon,
                                  const std::string& bootscheme,
                                  const arma::mat* exog_ptr, int n_threads,
                                  int stat_size, const char* label,
                                  Identify identify, bool verbose = false) {
  const int p      = var_result.p;
  const int c      = var_result.c;
  const int N      = static_cast<int>(y.n_cols);
  const int n_coef = static_cast<int>(var_result.beta.n_rows);
  const int T_iter = static_cast<int>(y.n_rows) - p;

  const BootstrapPlan plan = fBootstrapPlan_cpp(
      bootscheme, T_iter, var_result.residuals.n_rows, nboot);
  const int n_thr = fResolveThreads_cpp(n_threads);
  if (verbose)
    Rprintf("Using %d thread(s) for %s bootstrap...\n", n_thr, label);

  SVARBootResult r;
  r.flat.set_size(stat_size, nboot);
  r.boot_beta.set_size(N, n_coef, nboot);
  r.ok.assign(nboot, 0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_thr)
#endif
  for (int b = 0; b < nboot; ++b) {
    arma::mat yb = fBootstrapVAR_cpp(y, var_result, plan, b, exog_ptr);
    VARResult var_b = fVAR_cpp(yb, p, c, exog_ptr);
    r.boot_beta.slice(b) = var_b.beta.t();
    if (!var_b.ok) continue;                 // rank-deficient resample
    WoldIRFResult wold_b = fWoldIRF_cpp(var_b, horizon);
    arma::vec out(r.flat.colptr(b), stat_size, false, true);
    // `identify` runs on workers: it reports failure instead of throwing.
    if (!identify(var_b, wold_b.irfwold, out)) continue;
    r.ok[b] = 1;
  }
  return r;
}

// R wrapper function (for calling from R)
Rcpp::List fBootstrapVAR(const arma::mat& y,
                         const Rcpp::List& fVAR_result,
                         const std::string& bootscheme);

#endif
