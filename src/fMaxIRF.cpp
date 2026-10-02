// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Closed-form LR-max structural IRF (point estimate companion to fBootstrapMax).
// Maximises wold.slice(H-1).row(var_idx-1) * S * h2 subject to h2(0)=0, ||h2||=1.
// Sign: ensures the selected variable's last-horizon response is non-negative.

//' Long-Run Maximum-Response Impulse Responses
//'
//' Selects the shock direction that maximises the response of a selected
//' variable at the last supplied horizon, subject to a zero first shock
//' coordinate and unit norm. The sign is chosen so the selected variable's
//' last-horizon response is non-negative.
//'
//' @param wold N x N x (horizon+1) Wold impulse response array, including
//'   horizon zero in the first slice, as returned by \code{\link{fWoldIRF}}.
//' @param S N x N lower triangular Cholesky factor of the VAR residual
//'   covariance matrix. At least two variables are required.
//' @param var_idx Index (1-based) of the variable whose response is maximised.
//' @return An N x (horizon+1) matrix of structural impulse responses.
//' @seealso \code{\link{fBootstrapMax}}, \code{\link{fUhligIRF}}
//' @export
// [[Rcpp::export]]
arma::mat fMaxIRF(const arma::cube& wold, const arma::mat& S, int var_idx) {
    const int H = static_cast<int>(wold.n_slices);
    const int N = static_cast<int>(wold.n_rows);

    if (H < 1)
        Rcpp::stop("'wold' must contain at least one horizon slice.");
    if (N < 2 || wold.n_cols != wold.n_rows)
        Rcpp::stop("'wold' must be a square N x N x H array with N >= 2.");
    if (S.n_rows != wold.n_rows || S.n_cols != wold.n_cols)
        Rcpp::stop("'S' must be an N x N matrix conformable with 'wold'.");
    if (var_idx < 1 || var_idx > N)
        Rcpp::stop("'var_idx' is out of range.");

    // Optimal h2: normalised last column of (wold(H-1, var_idx-1, :) * S).cols(1..N-1)
    arma::rowvec M_full = wold.slice(H - 1).row(var_idx - 1) * S;
    arma::rowvec M_sub  = M_full.cols(1, N - 1);
    double M_norm = arma::norm(M_sub);

    arma::vec h2(N, arma::fill::zeros);
    if (M_norm > 1e-14)
        h2.rows(1, N - 1) = M_sub.t() / M_norm;

    // Precompute impact = S * h2 once (used for both sign check and IRF loop)
    arma::vec impact = S * h2;

    // Sign: selected variable's last-horizon response non-negative
    // Use a row-dot product with the precomputed impact instead of forming
    // the full (wold.slice(H-1) * S) matrix.
    if (arma::dot(wold.slice(H - 1).row(var_idx - 1), impact) < 0.0)
        impact = -impact;

    // Structural IRFs: N x H matrix
    arma::mat irf(N, H, arma::fill::none);
    for (int h = 0; h < H; ++h)
        irf.col(h) = wold.slice(h) * impact;

    return irf;
}
