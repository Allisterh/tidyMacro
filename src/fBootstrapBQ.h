#ifndef FBOOTSTRAPBQ_H
#define FBOOTSTRAPBQ_H

#include <RcppArmadillo.h>
#include "fBootstrapVAR.h"

// Blanchard-Quah bootstrap. Each draw stores the N x N x (horizon+1) BQ IRF
// cube (cumulated rows in `cumulate`, 0-based), vectorised.
SVARBootResult fBootstrapBQ_cpp(const arma::mat& y, const VARResult& var_result,
                                int nboot, int horizon, double conf, double conf2,
                                const std::string& bootscheme,
                                const arma::uvec& cumulate,
                                bool has_scaling, int scale_idx, double scale_size,
                                const arma::mat* exog_ptr, int n_threads, bool verbose);

#endif
