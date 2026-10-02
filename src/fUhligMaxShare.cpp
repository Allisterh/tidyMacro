// [[Rcpp::depends(RcppArmadillo)]]
#include "fUhligMaxShare.h"
#include <RcppArmadillo.h>

// Replicates uhlig_maxshare.m (Barsky & Sims 2012 / Uhlig 2004 max-share criterion).
//
// Builds Omega = sum_{h=0}^{H-1} cumsum_h  where cumsum_h = sum_{k=0}^{h} D_k'D_k
// and D_k = wold.slice(k).row(idx) * S  (1 x N).
// This equals sum_{k=0}^{H-1} (H-k) * D_k'*D_k — a horizon-weighted FEV matrix.
// The leading eigenvector of the lower-right (N-1)x(N-1) block gives the free
// parameters; h2 = [0; eigenvec] imposes the zero-impact restriction.

arma::vec fUhligMaxShare_cpp(const arma::cube& wold, const arma::mat& S, int idx) {
    const int H = static_cast<int>(wold.n_slices);
    const int N = static_cast<int>(wold.n_rows);

    arma::mat omega(N, N, arma::fill::zeros);
    arma::mat temp(N, N, arma::fill::zeros);

    for (int h = 0; h < H; ++h) {
        arma::rowvec D = wold.slice(h).row(idx) * S;  // 1 x N
        temp += D.t() * D;                             // N x N, cumulative
        omega += temp;                                  // sum of cumulative sums
    }

    // Leading eigenvector of lower-right (N-1) x (N-1) submatrix
    arma::mat sub = omega.submat(1, 1, N - 1, N - 1);
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    arma::eig_sym(eigenvalues, eigenvectors, sub);

    // eig_sym sorts ascending — last column is the leading eigenvector
    arma::vec h2(N, arma::fill::zeros);
    h2.rows(1, N - 1) = eigenvectors.col(N - 2);

    return h2;
}

//' Uhlig Maximum-Share Shock Direction
//'
//' Finds the unit-length shock direction with a zero first coordinate that
//' maximises the selected variable's forecast error variance contribution,
//' summed over the supplied horizons. For H slices, horizon h receives
//' weight H-h, with h starting at zero. The sign of the direction is not
//' normalised; \code{\link{fUhligIRF}} applies the IRF sign convention.
//'
//' @inheritParams fMaxIRF
//' @inheritParams fUhligIRF
//' @return An N x 1 matrix containing the shock direction, with first
//'   element zero.
//' @export
// [[Rcpp::export]]
arma::vec fUhligMaxShare(const arma::cube& wold, const arma::mat& S, int idx) {
    const int N = static_cast<int>(wold.n_rows);
    if (wold.n_slices < 1)
        Rcpp::stop("'wold' must contain at least one horizon slice.");
    if (N < 2 || wold.n_cols != wold.n_rows)
        Rcpp::stop("'wold' must be a square N x N x H array with N >= 2.");
    if (S.n_rows != wold.n_rows || S.n_cols != wold.n_cols)
        Rcpp::stop("'S' must be an N x N matrix conformable with 'wold'.");
    if (idx < 1 || idx > N)
        Rcpp::stop("'idx' is out of range.");
    return fUhligMaxShare_cpp(wold, S, idx - 1);  // 1-based -> 0-based
}
