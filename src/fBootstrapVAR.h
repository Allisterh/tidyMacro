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
                                const arma::mat* exog_ptr, int n_threads);

// Lower Cholesky factor of a bootstrap covariance, with a small ridge when
// the draw is numerically not positive definite (never throws).
inline arma::mat fLowerChol_cpp(const arma::mat& sigma) {
  arma::mat S;
  if (!arma::chol(S, sigma, "lower")) {
    arma::mat reg = 0.5 * (sigma + sigma.t());
    reg.diag() += 1e-8 * arma::trace(reg) / reg.n_rows;
    arma::chol(S, reg, "lower");
  }
  return S;
}

// Argument checks shared by the identified-VAR bootstraps (main thread).
void fCheckBootArgs_cpp(const arma::mat& y, const VARResult& var_result,
                        int nboot, int horizon, double conf, double conf2);

// ------------------------------------------------------------------------- //
//  Identified-VAR bootstrap driver (Cholesky, BQ, LR-max, Uhlig max-share)  //
// ------------------------------------------------------------------------- //
struct SVARBootResult {
  arma::mat  flat;       // stat_size x nboot identified responses per draw
  arma::mat  bands;      // stat_size x 5: upper, lower, upper2, lower2, median
  arma::cube boot_beta;  // N x n_coef x nboot bootstrapped coefficients
  arma::vec  conf;       // band levels in percent: conf, conf2
};

// Bands of an SVARBootResult (fills r.bands from r.flat).
void fSVARBootBands_cpp(SVARBootResult& r, double conf, double conf2,
                        int n_threads);

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
                                  Identify identify) {
  const int p      = var_result.p;
  const int c      = var_result.c;
  const int N      = static_cast<int>(y.n_cols);
  const int n_coef = static_cast<int>(var_result.beta.n_rows);
  const int T_iter = static_cast<int>(y.n_rows) - p;

  const BootstrapPlan plan = fBootstrapPlan_cpp(
      bootscheme, T_iter, var_result.residuals.n_rows, nboot);
  const int n_thr = fResolveThreads_cpp(n_threads);
  Rprintf("Using %d thread(s) for %s bootstrap...\n", n_thr, label);

  SVARBootResult r;
  r.flat.set_size(stat_size, nboot);
  r.boot_beta.set_size(N, n_coef, nboot);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_thr)
#endif
  for (int b = 0; b < nboot; ++b) {
    arma::mat yb = fBootstrapVAR_cpp(y, var_result, plan, b, exog_ptr);
    VARResult var_b = fVAR_cpp(yb, p, c, exog_ptr);
    r.boot_beta.slice(b) = var_b.beta.t();
    WoldIRFResult wold_b = fWoldIRF_cpp(var_b, horizon);
    arma::vec out(r.flat.colptr(b), stat_size, false, true);
    identify(var_b, wold_b.irfwold, out);
  }
  return r;
}

// R wrapper function (for calling from R)
Rcpp::List fBootstrapVAR(const arma::mat& y,
                         const Rcpp::List& fVAR_result,
                         const std::string& bootscheme);

#endif
