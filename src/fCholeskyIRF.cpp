
#include "fCholeskyIRF.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Compute Cholesky Impulse Response Functions
//'
//' @param wold Wold representation cube (N x N x horizon)
//' @param S Cholesky factor matrix (N x N), lower triangular
//'
//' @return A cube containing:
//'   Cholesky structural impulse response functions (N x N x horizon)
//'
//' @details
//' This function computes structural impulse response functions using the
//' Cholesky decomposition identification scheme. For each horizon h, the
//' structural IRF is computed as:
//' \deqn{IRF_h = \Psi_h \cdot S}
//' where \eqn{\Psi_h} is the Wold representation at horizon h and S is the
//' Cholesky factor (lower triangular) of the covariance matrix of reduced-form
//' residuals.
//'
//' The Cholesky identification imposes a recursive structure on the
//' contemporaneous relationships between variables, with the ordering
//' of variables determining the causal structure.
//'
//' @examples
//' set.seed(123)
//' y     <- matrix(rnorm(400), ncol = 2)
//' VAR   <- fVAR(y, p = 2, c = 1)
//' wold  <- fWoldIRF(VAR, horizon = 10)
//' S     <- t(chol(VAR$sigma))        # lower-triangular factor
//' chol_irf <- fCholeskyIRF(wold, S)
//' chol_irf[, 1, 1]                   # impact responses to shock 1
//'
//' @export
// [[Rcpp::export]]
arma::cube fCholeskyIRF(const arma::cube &wold, const arma::mat &S) {

  int N = wold.n_rows;
  int horizon = wold.n_slices;

  // Pre-allocate output cube
  arma::cube chol(N, N, horizon);

  // Compute Cholesky IRF for each horizon
  for (int h = 0; h < horizon; h++) {
    chol.slice(h) = wold.slice(h) * S;
  }

  return chol;
}
