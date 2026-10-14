#ifndef FREMOVE_BIAS_H
#define FREMOVE_BIAS_H

#include <RcppArmadillo.h>

struct RemoveBiasResult {
  arma::mat Beta;       // N x n_coef bias-corrected coefficients (transposed)
  int corrections;      // 1 = full correction; > 1 = shrinkage steps used
};

// Kilian (1998) bias correction. beta is the n_coef x N VAR estimate and
// boot_mean the N x n_coef mean of first-pass bootstrap estimates. The bias
// is shrunk in steps of 0.01 until the corrected VAR is stable; an unstable
// original estimate is returned uncorrected.
RemoveBiasResult fRemoveBias_cpp(const arma::mat& beta, int c, int p,
                                 const arma::mat& boot_mean);

#endif // FREMOVE_BIAS_H
