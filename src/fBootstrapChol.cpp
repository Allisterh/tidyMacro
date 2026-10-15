// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include <RcppArmadillo.h>
#include "fBootstrapChol.h"
#include "fBootstrapVAR.h"
#include "fCholeskyIRF.h"

SVARBootResult fBootstrapChol_cpp(const arma::mat& y, const VARResult& var_result,
                                  int nboot, int horizon, double conf, double conf2,
                                  const std::string& bootscheme,
                                  const arma::mat* exog_ptr, int n_threads) {
  fCheckBootArgs_cpp(y, var_result, nboot, horizon, conf, conf2);
  const int N = static_cast<int>(y.n_cols);

  SVARBootResult r = fBootstrapSVAR_cpp(
      y, var_result, nboot, horizon, bootscheme, exog_ptr, n_threads,
      N * N * (horizon + 1), "Cholesky",
      [](const VARResult& var_b, const arma::cube& wold_b, arma::vec& out) {
        out = arma::vectorise(fCholeskyIRF(wold_b, fLowerChol_cpp(var_b.sigma)));
      });
  fSVARBootBands_cpp(r, conf, conf2, n_threads);
  return r;
}

static Rcpp::List chol_boot_list(const SVARBootResult& r, int N, int horizon) {
  return fSVARBootList_cpp(r, "bootchol", {N, N, horizon + 1});
}

//' Bootstrap Cholesky Identified Impulse Response Functions
//'
//' @param y T x N matrix of original endogenous variables.
//' @param var_result List from \code{fVAR()}.
//' @param nboot Number of bootstrap replications.
//' @param horizon Maximum IRF horizon.
//' @param conf Outer confidence level in percent (e.g. 90).
//' @param conf2 Inner confidence level in percent (e.g. 68).
//' @param bootscheme \code{"residual"} or \code{"wild"}.
//' @param exog Optional T x M matrix of exogenous variables. Must be provided
//'   if the original VAR was estimated with exogenous variables; must be
//'   omitted otherwise.
//' @param n_threads Number of OpenMP threads. 0 = all cores minus one.
//'   Default 0.
//'
//' @return List with elements:
//'   \describe{
//'     \item{bootchol}{N x N x (horizon+1) x nboot array of bootstrapped Cholesky IRFs}
//'     \item{upper}{N x N x (horizon+1) outer upper confidence bands}
//'     \item{lower}{N x N x (horizon+1) outer lower confidence bands}
//'     \item{upper2}{N x N x (horizon+1) inner upper confidence bands}
//'     \item{lower2}{N x N x (horizon+1) inner lower confidence bands}
//'     \item{boot_beta}{N x (Np+c) x nboot array of bootstrapped coefficients}
//'     \item{conf}{Band levels in percent, \code{c(conf, conf2)}}
//'   }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fBootstrapChol(const arma::mat &y, const Rcpp::List &var_result,
                          int nboot, int horizon, double conf = 90.0,
                          double conf2 = 68.0,
                          const std::string &bootscheme = "residual",
                          Rcpp::Nullable<arma::mat> exog = R_NilValue,
                          int n_threads = 0) {
  VARResult vr = fVARResultFromList_cpp(var_result);
  fCheckExogArg_cpp(vr.n_exog, exog.isNotNull());
  arma::mat exog_mat;
  if (exog.isNotNull()) exog_mat = Rcpp::as<arma::mat>(exog);

  SVARBootResult r = fBootstrapChol_cpp(y, vr, nboot, horizon, conf, conf2,
                                        bootscheme,
                                        exog.isNotNull() ? &exog_mat : nullptr,
                                        n_threads);
  return chol_boot_list(r, y.n_cols, horizon);
}

//' Bootstrap Bias-Corrected Cholesky Identified Impulse Response Functions
//'
//' @param y T x N matrix of original endogenous variables.
//' @param var_result List from \code{fVAR()}.
//' @param nboot1 Number of first-pass replications for bias estimation.
//' @param nboot2 Number of second-pass replications for confidence bands.
//' @param horizon Maximum IRF horizon.
//' @param conf Outer confidence level in percent (e.g. 90).
//' @param conf2 Inner confidence level in percent (e.g. 68).
//' @param bootscheme \code{"residual"} or \code{"wild"}.
//' @param exog Optional T x M matrix of exogenous variables. Must be provided
//'   if the original VAR was estimated with exogenous variables; must be
//'   omitted otherwise.
//' @param n_threads Number of OpenMP threads. 0 = all cores minus one.
//'   Default 0.
//'
//' @return List with elements:
//'   \describe{
//'     \item{bootchol}{N x N x (horizon+1) x nboot2 array of bootstrapped Cholesky IRFs}
//'     \item{upper}{N x N x (horizon+1) outer upper confidence bands}
//'     \item{lower}{N x N x (horizon+1) outer lower confidence bands}
//'     \item{upper2}{N x N x (horizon+1) inner upper confidence bands}
//'     \item{lower2}{N x N x (horizon+1) inner lower confidence bands}
//'     \item{boot_beta}{N x (Np+c) x nboot2 array of bootstrapped coefficients}
//'     \item{conf}{Band levels in percent, \code{c(conf, conf2)}}
//'   }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fBootstrapCholCorrected(const arma::mat &y, const Rcpp::List &var_result,
                                    int nboot1, int nboot2, int horizon,
                                    double conf = 90.0, double conf2 = 68.0,
                                    const std::string &bootscheme = "residual",
                                    Rcpp::Nullable<arma::mat> exog = R_NilValue,
                                    int n_threads = 0) {
  VARResult vr = fVARResultFromList_cpp(var_result);
  fCheckExogArg_cpp(vr.n_exog, exog.isNotNull());
  arma::mat exog_mat;
  if (exog.isNotNull()) exog_mat = Rcpp::as<arma::mat>(exog);
  const arma::mat* exog_ptr = exog.isNotNull() ? &exog_mat : nullptr;

  VARResult corrected = fBiasCorrectedVAR_cpp(y, vr, nboot1, bootscheme,
                                              exog_ptr, n_threads);
  SVARBootResult r = fBootstrapChol_cpp(y, corrected, nboot2, horizon, conf,
                                        conf2, bootscheme, exog_ptr, n_threads);
  return chol_boot_list(r, y.n_cols, horizon);
}
