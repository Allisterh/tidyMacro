#ifndef FLPPANEL_H
#define FLPPANEL_H

#include <RcppArmadillo.h>
#include <vector>

// -----------------------------------------------------------------------
// Panel Local Projections — C++ engine
//
// Model at horizon h (per s-component j = 1..n_s):
//     y_{i,t+h} = alpha_{FE} + sum_j (s_{it,j} * X_it) * beta_h^j
//                 + W'_{it} * G^h
//                 + sum_{k=1..p} lag_k(y, s⊗X)' * D^h_k + u_{i,t+h}
//   (cumulative: replaces y_{i,t+h} by sum_{r=0..h} y_{i,t+r})
//
// Sample: y, s, X and W may contain NaN. Each horizon keeps the rows whose
// led (or cumulated) y, s⊗X, lags and W are all finite, so the sample is
// chosen per horizon, as with leads built before the regression in Stata.
//
// Inference:
//   * ssc_mode = 0 (none) → asymptotic cluster-robust sandwich SEs:
//       cluster_mode = 0: time, 1: unit, 2: unit + time (two-way).
//     Two-way clustering uses the Cameron-Gelbach-Miller inclusion-exclusion
//     meat M_i + M_t - M_it, with df = Inf (normal critical values).
//   * ssc_mode = 1 (CR1)  → the Stata / fixest finite-sample correction:
//     each one-way component of the meat is scaled by
//     G/(G-1) * (N-1)/(N-K) with its own cluster count, and inference uses
//     t with min(G) - 1 degrees of freedom. K counts the regressors, the
//     intercept, and L - 1 parameters for every absorbed factor with L
//     levels that is not nested in a clustering dimension.
//   * ssc_mode = 2 (IK)   → Imbens-Kolesar (2016, REStat) refinement,
//     per s-component, using the within-period projection X_t and the
//     hat-matrix P0 = I_T - X_t * (X_t'X_t)^{-1} * X_t'. This refinement
//     is defined only for cluster_mode = 0 (time clustering).
//
// The horizon loop is parallelized with OpenMP; y_h and the lag columns
// used across horizons are precomputed once, so each horizon runs
// independently on a shared read-only cache.
//
// Author: Dr. Muhsin Ciftci
// -----------------------------------------------------------------------

struct LPPanelResult {
  arma::mat  estimate;   // (H+1) x n_s
  arma::mat  SE;         // (H+1) x n_s
  arma::mat  df;         // (H+1) x n_s
  arma::ivec status;     // (H+1)   0 = OK, 1 = no complete rows,
                         //         2 = HDFE non-convergence
  arma::ivec converged;  // (H+1)   1 = HDFE converged, 0 = did not
  arma::ivec nobs;       // (H+1)   rows used at each horizon
};

LPPanelResult fLPPanel_internal(
    const arma::mat& y,
    const arma::mat& s,
    const arma::mat& X,
    const arma::mat& W,
    const arma::imat& FE,
    const arma::ivec& i_index,
    const arma::ivec& t_index,
    int  H,
    int  p_max,
    int  ssc_mode,
    bool cumulative,
    int  cluster_mode,
    int  n_threads,
    bool verbose
);

Rcpp::List fLPPanel_cpp(
    const arma::mat& y,
    const arma::mat& s,
    const arma::mat& X,
    const arma::mat& W,
    const arma::imat& FE,
    const arma::ivec& i_index,
    const arma::ivec& t_index,
    int  H,
    int  p_max,
    int  ssc_mode,
    bool cumulative,
    int  cluster_mode,
    int  n_threads,
    bool verbose
);

#endif // FLPPANEL_H
