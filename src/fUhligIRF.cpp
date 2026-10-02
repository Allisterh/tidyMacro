// [[Rcpp::depends(RcppArmadillo)]]
#include "fUhligMaxShare.h"
#include <RcppArmadillo.h>

// Uhlig max-share structural IRF (point estimate companion to fBootstrapUhlig).
// Uses fUhligMaxShare_cpp to find the FEV-maximising h2, then applies the
// same sign convention as the bootstrap: first variable's last-horizon
// non-structural response is non-negative.

//' Uhlig Maximum-Share Impulse Responses
//'
//' Computes impulse responses using the direction returned by
//' \code{\link{fUhligMaxShare}}. The sign is chosen so the first variable's
//' last-horizon response is non-negative.
//'
//' @inheritParams fMaxIRF
//' @param idx Index (1-based) of the variable whose forecast error variance
//'   contribution is maximised.
//' @return An N x (horizon+1) matrix of structural impulse responses.
//' @seealso \code{\link{fBootstrapUhlig}}
//' @export
// [[Rcpp::export]]
arma::mat fUhligIRF(const arma::cube& wold, const arma::mat& S, int idx) {
    const int H = static_cast<int>(wold.n_slices);
    const int N = static_cast<int>(wold.n_rows);

    if (H < 1)
        Rcpp::stop("'wold' must contain at least one horizon slice.");
    if (N < 2 || wold.n_cols != wold.n_rows)
        Rcpp::stop("'wold' must be a square N x N x H array with N >= 2.");
    if (S.n_rows != wold.n_rows || S.n_cols != wold.n_cols)
        Rcpp::stop("'S' must be an N x N matrix conformable with 'wold'.");
    if (idx < 1 || idx > N)
        Rcpp::stop("'idx' is out of range.");

    arma::vec h2 = fUhligMaxShare_cpp(wold, S, idx - 1);  // 1-based -> 0-based

    // Precompute impact = S * h2 once (used for both sign check and IRF loop)
    arma::vec impact = S * h2;

    // Sign: first variable's last-horizon non-structural response non-negative
    if (arma::dot(wold.slice(H - 1).row(0), impact) < 0.0) impact = -impact;

    // Structural IRFs: N x H matrix
    arma::mat irf(N, H, arma::fill::none);
    for (int h = 0; h < H; ++h)
        irf.col(h) = wold.slice(h) * impact;

    return irf;
}
