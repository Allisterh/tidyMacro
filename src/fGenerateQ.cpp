#include <RcppArmadillo.h>
#include "fGenerateQ.h"
#include <stdexcept>
// [[Rcpp::depends(RcppArmadillo)]]

void fGenerateQ_inplace(arma::mat& Q,
                        arma::mat& R,
                        arma::mat& G,
                        const arma::uword N) {
    if (N == 0u) {
        Rcpp::stop("N must be positive.");
    }

    G.randn(N, N);
    arma::qr_econ(Q, R, G);

    for (arma::uword j = 0; j < N; ++j) {
        if (R(j, j) < 0.0) {
            Q.col(j) *= -1.0;
        }
    }
}

void fGenerateQ_inplace(arma::mat& Q,
                        const arma::uword N,
                        tidymacro::RNG& rng) {
    if (N == 0u) {
        throw std::invalid_argument("N must be positive.");
    }

    // Gram-Schmidt on i.i.d. Gaussian columns yields the Q factor of the QR
    // decomposition with a positive diagonal in R, i.e. the Haar draw above.
    // Columns are filled in the same order as a column-major Gaussian matrix,
    // and a second orthogonalisation pass (CGS2) keeps Q orthonormal to
    // machine precision.
    Q.set_size(N, N);
    for (arma::uword j = 0; j < N; ++j) {
        double* q = Q.colptr(j);
        double nrm2 = 0.0;
        do {
            for (arma::uword i = 0; i < N; ++i) q[i] = rng.norm();
            for (int pass = 0; pass < 2; ++pass) {
                for (arma::uword l = 0; l < j; ++l) {
                    const double* ql = Q.colptr(l);
                    double d = 0.0;
                    for (arma::uword i = 0; i < N; ++i) d += ql[i] * q[i];
                    for (arma::uword i = 0; i < N; ++i) q[i] -= d * ql[i];
                }
            }
            nrm2 = 0.0;
            for (arma::uword i = 0; i < N; ++i) nrm2 += q[i] * q[i];
        } while (!(nrm2 > 1e-24));   // probability-zero degenerate draw
        const double inv = 1.0 / std::sqrt(nrm2);
        for (arma::uword i = 0; i < N; ++i) q[i] *= inv;
    }
}

arma::mat fGenerateQ_cpp(int N) {
    if (N <= 0) {
        Rcpp::stop("N must be positive.");
    }

    arma::mat Q;
    arma::mat R;
    arma::mat G;
    fGenerateQ_inplace(Q, R, G, static_cast<arma::uword>(N));
    return Q;
}

//' Draw a Random Orthonormal Matrix (RZW 2010)
//'
//' Generates an N x N orthonormal matrix Q (QQ' = Q'Q = I) by applying QR
//' decomposition to a matrix of i.i.d. standard normals and normalising the
//' diagonal of R to be positive, yielding a unique decomposition.
//'
//' @param N Integer. Dimension of the square matrix.
//'
//' @return An N x N orthonormal matrix.
//'
//' @details
//' Draw an N x N matrix M of independent standard normals, compute the
//' QR decomposition M = QR, then set
//' Q[,i] = -Q[,i] whenever R[i,i] < 0.  The resulting Q is uniformly
//' distributed on the Stiefel manifold (Haar measure), as required for
//' sign-restriction identification.
//'
//' @references
//' Rubio-Ramirez, J. F., Waggoner, D. F., & Zha, T. (2010). Structural vector
//' autoregressions: Theory of identification and algorithms for inference.
//' \emph{Review of Economic Studies}, 77(2), 665--696.
//'
//' @examples
//' \dontrun{
//' Q <- fGenerateQ(3)
//' round(t(Q) %*% Q, 10)  # should be identity
//' }
//'
//' @export
// [[Rcpp::export]]
arma::mat fGenerateQ(int N) {
    return fGenerateQ_cpp(N);
}
