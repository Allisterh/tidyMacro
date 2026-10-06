#include <RcppArmadillo.h>
#include "fWoldIRF.h"
#include "fVAR.h"
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

// Internal C++ function (called from other C++ code)
// Block recursion Phi_h = sum_{l=1..min(h,p)} Phi_{h-l} * A_l, evaluated as
// one product per horizon: cube slices are contiguous, so slices h-L..h-1 form
// the N x (N L) matrix [Phi_{h-L} ... Phi_{h-1}], which multiplies the stacked
// coefficients [A_L; ...; A_1].  O(N^3 * p) per step vs O((N*p)^3) for
// companion matrix powers.
void fWoldIRF_into_cpp(const arma::mat& beta, int c, int p, int horizon,
                       arma::cube& irfwold) {
  const arma::uword N = beta.n_cols;

  // S = [A_p; ...; A_1], with A_l = beta block l transposed.
  arma::mat S(N * p, N, arma::fill::none);
  for (int l = 1; l <= p; ++l)
    S.rows((p - l) * N, (p - l + 1) * N - 1) =
      beta.rows(c + (l - 1) * N, c + l * N - 1).t();

  if (irfwold.n_rows != N || irfwold.n_cols != N ||
      irfwold.n_slices != static_cast<arma::uword>(horizon + 1))
    irfwold.set_size(N, N, horizon + 1);
  irfwold.slice(0).eye();

  for (int h = 1; h <= horizon; ++h) {
    const int L = std::min(h, p);
    const arma::mat lagged(irfwold.slice_memptr(h - L), N, N * L, false, true);
    arma::mat out(irfwold.slice_memptr(h), N, N, false, true);
    if (L == p) out = lagged * S;
    else        out = lagged * S.tail_rows(N * L);
  }
}

WoldIRFResult fWoldIRF_cpp(const VARResult& var_result, int horizon) {
  WoldIRFResult result;
  fWoldIRF_into_cpp(var_result.beta, var_result.c, var_result.p, horizon,
                    result.irfwold);
  return result;
}

//' Compute Wold Impulse Response Functions for VAR Model
//'
//' @param fVAR A list containing VAR estimation results with elements:
//'   \itemize{
//'     \item beta: Coefficient matrix
//'     \item c: Integer indicator for intercept (1 if intercept, 0 otherwise)
//'     \item p: Integer lag order
//'     \item n_exog: Number of exogenous variables (optional)
//'   }
//' @param horizon Integer, the IRF horizon (number of periods ahead)
//'
//' @return A 3D array (cube) of Wold impulse response functions with dimensions
//'   N x N x (horizon+1), where irfwold[,,h] contains the response at horizon h
//'
//' @details
//' This function computes the Wold (reduced-form) impulse response functions
//' for a VAR(p) model with optional exogenous variables by calculating powers 
//' of the companion matrix. The Wold IRF shows the effect of a one-unit shock 
//' to each variable on all variables in the system over time, without imposing 
//' any identifying restrictions.
//'
//' Exogenous variables do not contribute to the dynamic propagation of shocks
//' and are therefore excluded from the companion matrix representation.
//'
//' The IRFs are computed as:
//' \deqn{\Phi_h = A^h[1:N, 1:N]}
//' where A is the companion matrix and h is the horizon.
//'
//' @examples
//' \dontrun{
//' # VAR(2) model with 3 variables and 1 exogenous variable
//' y <- matrix(rnorm(200), ncol = 2)
//' exog <- matrix(rnorm(100), ncol = 1)
//' VAR <- fVAR(y, p = 2, c = 1, exog = exog)
//' irf <- fWoldIRF(VAR, horizon = 10)
//' # irf is a 2x2x11 array
//' 
//' # Plot impulse response
//' plot(irf[1, 2, ], type = "l", 
//'      main = "Response of Variable 2 to Shock in Variable 1")
//' }
//'
//' @export
// [[Rcpp::export]]
arma::cube fWoldIRF(const Rcpp::List& fVAR, int horizon) {
  // Extract elements from R list and construct VARResult struct
  VARResult var_result;
  var_result.beta = Rcpp::as<arma::mat>(fVAR["beta"]);
  var_result.residuals = Rcpp::as<arma::mat>(fVAR["residuals"]);
  var_result.sigma = Rcpp::as<arma::mat>(fVAR["sigma"]);
  var_result.p = Rcpp::as<int>(fVAR["p"]);
  var_result.c = Rcpp::as<int>(fVAR["c"]);
  var_result.n_exog = Rcpp::as<int>(fVAR["n_exog"]);
  
  // Call the C++ function with the struct
  WoldIRFResult result = fWoldIRF_cpp(var_result, horizon);
  
  // Return the IRF cube for R
  return result.irfwold;
}
