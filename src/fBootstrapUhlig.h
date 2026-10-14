#ifndef FBOOTSTRAPUHLIG_H
#define FBOOTSTRAPUHLIG_H

#include <RcppArmadillo.h>
#include "fBootstrapVAR.h"

// Uhlig max-share residual bootstrap. Each draw stores the N x (horizon+1)
// response to the identified shock, vectorised. idx and cumulate are 0-based.
SVARBootResult fBootstrapUhlig_cpp(const arma::mat& y, const VARResult& var_result,
                                   int nboot, int horizon, int idx,
                                   double conf, double conf2,
                                   const arma::uvec& cumulate,
                                   const arma::mat* exog_ptr, int n_threads);

#endif
