#ifndef FBOOTSTRAPCHOL_H
#define FBOOTSTRAPCHOL_H

#include <RcppArmadillo.h>
#include "fBootstrapVAR.h"

// Recursive (Cholesky) bootstrap. Each draw stores the N x N x (horizon+1)
// structural IRF cube, vectorised.
SVARBootResult fBootstrapChol_cpp(const arma::mat& y, const VARResult& var_result,
                                  int nboot, int horizon, double conf, double conf2,
                                  const std::string& bootscheme,
                                  const arma::mat* exog_ptr, int n_threads, bool verbose);

#endif
