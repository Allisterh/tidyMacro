#ifndef FVAR_H
#define FVAR_H

#include <RcppArmadillo.h>

// Struct to hold VAR estimation results
struct VARResult {
  arma::mat beta;        // Coefficient matrix (Np + c + M) x N
  arma::mat residuals;   // (T-p) x N matrix of OLS residuals
  arma::mat sigma;       // N x N residual covariance matrix (DF-corrected)
  int p;                 // Lag order
  int c;                 // Intercept indicator
  int n_exog;            // Number of exogenous variables
  bool ok = true;        // false when the design was rank deficient, so beta
                         // comes from a pseudo-inverse and is not unique.
                         // Set without touching the R API, so bootstrap
                         // workers can record a failed draw and the caller
                         // can reject it after the parallel region.
};

// Internal C++ function — no-exog overload. Use this inside parallel regions
// to avoid constructing Rcpp objects on worker threads.
VARResult fVAR_cpp(const arma::mat& y, int p, int c);

// Internal C++ function (for use in other C++ code) — Nullable-exog overload
VARResult fVAR_cpp(const arma::mat& y, int p, int c,
                   Rcpp::Nullable<arma::mat> exog);

// Internal C++ function — pre-converted exog overload (thread-safe, no Rcpp::as)
// Use this from bootstrap loops where exog has already been converted once
// outside the parallel region.
VARResult fVAR_cpp(const arma::mat& y, int p, int c,
                   const arma::mat& exog);

// Optional pre-converted exog: nullptr means no exogenous variables.
VARResult fVAR_cpp(const arma::mat& y, int p, int c,
                   const arma::mat* exog_ptr);

// VARResult from the list returned by fVAR() (n_exog defaults to 0).
VARResult fVARResultFromList_cpp(const Rcpp::List& var_result);

// Validates the arguments of an exported VAR entry point: p, c, a finite and
// long enough y, matching exogenous rows, and positive residual degrees of
// freedom. Calls Rcpp::stop, so it runs on the main thread only, once per
// call — never inside a bootstrap draw.
void fCheckVARArgs_cpp(const arma::mat& y, int p, int c,
                       const arma::mat* exog_ptr);

// Checks that `exog` is supplied exactly when the fitted VAR used exogenous
// variables. Call on the main thread, before any OpenMP region.
void fCheckExogArg_cpp(int n_exog, bool has_exog);

// For routines that do not support exogenous variables.
void fRequireNoExog_cpp(int n_exog, const char* fun);


// R wrapper function (for calling from R)
Rcpp::List fVAR(const arma::mat& y, int p, int c,
                Rcpp::Nullable<arma::mat> exog);

#endif // FVAR_H
