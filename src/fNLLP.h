#ifndef FNLLP_H
#define FNLLP_H

#include <RcppArmadillo.h>
#include <vector>

// -----------------------------------------------------------------------
// Nonlinear and state-dependent Local Projections — C++ engine
//
// The caller passes the RHS exactly as written in the formula (X, no
// constant). The engine builds the regression design once from X and the
// specification, then runs one OLS per horizon. Two shock coefficients are
// tracked: c1 and c2.
//
//   specification 0 (sign):  [1, ..., max(s,0), min(s,0), ...]
//                            c1 = beta+, c2 = beta-
//   specification 1 (cubic): [1, ..., s, s^3, ...]
//                            c1 = beta_1, c2 = beta_3
//   specification 2 (state): [M, M*X_nc, (1-M), (1-M)*X_nc, X_c]
//                            c1 = coefficient on M*s, c2 = on (1-M)*s
//
// In the sign and cubic designs the shock column of X is replaced in place
// by its two transforms. In the state design the intercept is split into M
// and (1-M) (a common constant would be collinear with them); X_nc are the
// RHS columns interacted with the state, X_c the columns listed in
// common_cols, which enter once with a single coefficient (e.g. a trend).
//
// Sample:
//   balanced = false: horizon h uses rows t = t0..T-1-h (sample shrinks
//                     with h, as in fLP). nrow(Y) == nrow(X).
//   balanced = true:  every horizon uses rows t = t0..T-1 of X and reads the
//                     outcome from Y rows t+h. nrow(Y) >= nrow(X) + H.
//                     This is the fixed regressor sample of Tenreyro and
//                     Thwaites (2016, stlpm.m).
//
// Variance: Newey-West HAC sandwich per horizon and equation, Bartlett
// kernel, bandwidth min(max(nw_lags_base + h + nw_offset, 0), T_h - 1),
// identical to fLP. The 2 x 2 block of (c1, c2) is returned so the
// difference c1 - c2 and any response a*c1 + b*c2 have valid SEs.
//
// Author: Dr. Muhsin Ciftci
// -----------------------------------------------------------------------

struct NLLPResult {
  arma::cube irfs;       // (H+1) x n_y x 2 — c1 and c2
  arma::cube irfs_se;    // (H+1) x n_y x 2 — HAC SE of c1 and c2
  arma::mat  irfs_cov;   // (H+1) x n_y     — HAC Cov(c1, c2)
  arma::mat  diff;       // (H+1) x n_y     — c1 - c2 (NaN for cubic)
  arma::mat  diff_se;    // (H+1) x n_y     — SE of c1 - c2 (NaN for cubic)
  std::vector<int> nobs; // H+1             — observations per horizon

  bool rank_deficient = false;
  int  rank_fail_h    = -1;

  // Only populated when store_full = true:
  std::vector<arma::mat> betas;  // H+1 matrices (kd x n_y)
  std::vector<arma::mat> ses;    // H+1 matrices (kd x n_y)
};

// Regression design implied by X and the specification. c1 / c2 receive
// the 0-indexed design columns of the two shock coefficients.
arma::mat fNLLPDesign_cpp(
    const arma::mat&  X,
    int               shock_col,
    int               specification,
    const arma::vec&  state,
    const arma::uvec& common_cols,
    int&              c1,
    int&              c2
);

// Internal C++ function (callable from other translation units).
// Defaults live in the .cpp definitions (project convention).
NLLPResult fNLLP_internal(
    const arma::mat&  Y,
    const arma::mat&  X,
    int               H,
    int               shock_col,
    int               specification,
    const arma::vec&  state,
    const arma::uvec& common_cols,
    int               nw_lags_base,
    bool              store_full,
    bool              cumulative,
    bool              balanced,
    int               n_threads,
    int               nw_offset,
    bool              verbose,
    const arma::mat&  Y_pre
);

// R-callable wrapper
Rcpp::List fNLLP_cpp(
    const arma::mat&    Y,
    const arma::mat&    X,
    int                 H,
    int                 shock_col,
    int                 specification,
    const arma::vec&    state,
    Rcpp::IntegerVector common_cols,
    int                 nw_lags_base,
    bool                store_full,
    bool                cumulative,
    bool                balanced,
    int                 n_threads,
    int                 nw_offset,
    bool                verbose,
    Rcpp::Nullable<arma::mat> Y_pre
);

#endif // FNLLP_H
