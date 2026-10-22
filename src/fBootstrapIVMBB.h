#ifndef FBOOTSTRAPIV_MBB_H
#define FBOOTSTRAPIV_MBB_H

#include "fVAR.h" // Include to access VARResult struct
#include <RcppArmadillo.h>

// Struct to hold bootstrap IV MBB results
struct BootstrapIVMBBResult {
  // Bands are percentiles recentred on `point` (q - median + point)
  arma::mat upper;   // N x (horizon+1) matrix of upper confidence bands (conf)
  arma::mat lower;   // N x (horizon+1) matrix of lower confidence bands (conf)
  arma::mat upper2;  // N x (horizon+1) matrix of upper confidence bands (conf2)
  arma::mat lower2;  // N x (horizon+1) matrix of lower confidence bands (conf2)
  arma::mat point;   // N x (horizon+1) impulse responses on the original
                     // sample: Wold coefficients times the proxy-SVAR impact
                     // vector, as in fHeteroIRF() and fMSW(). The bootstrap
                     // supplies the bands only, never the point estimate.
};

// Proxy-SVAR impact vector (Gertler-Karadi / Mertens-Ravn, unit effect on
// the policy variable): first stage of u(:, pv) on Z and a constant, then the
// remaining residuals on the fitted values and a constant. u and Z are row-aligned; pv is 0-based.
// Returns false when the second stage is singular (a weak or degenerate
// instrument on this sample) instead of throwing, so it is safe to call
// from an OpenMP worker; `s` is then left unusable and the draw is dropped.
bool fProxyImpact_cpp(const arma::mat& u, const arma::mat& Z, int pv,
                      arma::vec& s);

// Full-sample residuals with rows u_start.. replaced by an MBB draw.
arma::mat fMBBPadResiduals_cpp(const arma::mat& residuals, const arma::mat& eps_boot,
                               int u_start);

BootstrapIVMBBResult
fBootstrapIVMBB_cpp(const arma::mat &y, const VARResult &var_result,
                     const arma::mat &Z, int nboot, int blocksize,
                     const arma::ivec &adjustZ, const arma::ivec &adjustu,
                     int policyvar, int horizon, double conf, double conf2,
                     const arma::mat* exog_ptr, int n_threads, bool verbose);

#endif
