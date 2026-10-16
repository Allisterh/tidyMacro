#ifndef FBOOTSTRAPMAX_H
#define FBOOTSTRAPMAX_H

#include <RcppArmadillo.h>
#include "fBootstrapVAR.h"

// LR-max residual bootstrap. Each draw stores the N x (horizon+1) response
// to the identified shock, vectorised. var_idx and cumulate are 0-based.
SVARBootResult fBootstrapMax_cpp(const arma::mat& y, const VARResult& var_result,
                                 int nboot, int horizon, int var_idx,
                                 double conf, double conf2,
                                 const arma::uvec& cumulate,
                                 bool has_scaling, int scale_idx, double scale_size,
                                 const arma::mat* exog_ptr, int n_threads, bool verbose);

#endif
