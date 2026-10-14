#include <RcppArmadillo.h>
#include "fCompanionMatrix.h"
#include "fRemoveBias.h"
// [[Rcpp::depends(RcppArmadillo)]]

static double max_modulus(const arma::mat& beta, int c, int p) {
    CompanionMatrixResult comp = fCompanionMatrix_cpp(beta, c, p);
    return arma::max(arma::abs(arma::eig_gen(comp.comp)));
}

RemoveBiasResult fRemoveBias_cpp(const arma::mat& beta, int c, int p,
                                 const arma::mat& boot_mean) {
    const arma::mat beta_t = beta.t();
    const arma::mat bias   = boot_mean - beta_t;

    RemoveBiasResult result;
    result.Beta        = beta_t;
    result.corrections = 1;

    // Original estimates already unstable: no correction
    if (max_modulus(beta, c, p) >= 1.0) return result;

    arma::mat Beta_t = beta_t - bias;
    double max_ev = max_modulus(Beta_t.t(), c, p);

    // Shrink the bias towards zero until the corrected VAR is stable
    double delta = 1.0;
    while (max_ev >= 1.0) {
        delta -= 0.01;
        result.corrections += 1;
        if (delta < 0.0 || result.corrections > 200) {
            Beta_t = beta_t;
            break;
        }
        Beta_t = beta_t - bias * delta;
        max_ev = max_modulus(Beta_t.t(), c, p);
    }

    result.Beta = Beta_t;
    return result;
}
