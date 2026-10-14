#ifndef FBOOTSTRAPIV_MBB_H
#define FBOOTSTRAPIV_MBB_H

#include "fVAR.h" // Include to access VARResult struct
#include <RcppArmadillo.h>

// Struct to hold bootstrap IV MBB results
struct BootstrapIVMBBResult {
  arma::mat upper;     // N x (horizon+1) matrix of upper confidence bands (conf)
  arma::mat lower;     // N x (horizon+1) matrix of lower confidence bands (conf)
  arma::mat upper2;    // N x (horizon+1) matrix of upper confidence bands (conf2)
  arma::mat lower2;    // N x (horizon+1) matrix of lower confidence bands (conf2)
  arma::mat meanirf;   // N x (horizon+1) matrix of mean impulse responses
  arma::mat medianirf; // N x (horizon+1) matrix of median impulse responses
};

// Proxy-SVAR impact vector (Gertler-Karadi / Mertens-Ravn, unit effect on
// the policy variable): first stage of u(:, pv) on Z, then the remaining
// residuals on the fitted values. u and Z are row-aligned; pv is 0-based.
arma::vec fProxyImpact_cpp(const arma::mat& u, const arma::mat& Z, int pv);

// Full-sample residuals with rows u_start.. replaced by an MBB draw.
arma::mat fMBBPadResiduals_cpp(const arma::mat& residuals, const arma::mat& eps_boot,
                               int u_start);

BootstrapIVMBBResult
fBootstrapIVMBB_cpp(const arma::mat &y, const VARResult &var_result,
                     const arma::mat &Z, int nboot, int blocksize,
                     const arma::ivec &adjustZ, const arma::ivec &adjustu,
                     int policyvar, int horizon, double conf, double conf2,
                     const arma::mat* exog_ptr, int n_threads);

#endif
