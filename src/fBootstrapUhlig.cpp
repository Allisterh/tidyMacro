// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include <RcppArmadillo.h>
#include "fBootstrapUhlig.h"
#include "fBootstrapVAR.h"
#include "fUhligMaxShare.h"

SVARBootResult fBootstrapUhlig_cpp(const arma::mat& y, const VARResult& var_result,
                                   int nboot, int horizon, int idx,
                                   double conf, double conf2,
                                   const arma::uvec& cumulate,
                                   const arma::mat* exog_ptr, int n_threads, bool verbose) {
    fCheckBootArgs_cpp(y, var_result, nboot, horizon, conf, conf2);
    const int N = static_cast<int>(y.n_cols);
    const int H = horizon + 1;

    SVARBootResult r = fBootstrapSVAR_cpp(
        y, var_result, nboot, horizon, "residual", exog_ptr, n_threads,
        N * H, "Uhlig max-share",
        [&](const VARResult& var_b, const arma::cube& wold_b, arma::vec& out) {
            bool chol_ok = true;
            const arma::mat S = fLowerChol_cpp(var_b.sigma, &chol_ok);
            if (!chol_ok) return false;
            arma::vec impact = S * fUhligMaxShare_cpp(wold_b, S, idx);

            // Sign: last-horizon response of variable 1 is non-negative
            if (arma::dot(wold_b.slice(H - 1).row(0), impact) < 0.0)
                impact = -impact;

            arma::mat irf(N, H, arma::fill::none);
            for (int h = 0; h < H; ++h) irf.col(h) = wold_b.slice(h) * impact;
            for (arma::uword ci = 0; ci < cumulate.n_elem; ++ci)
                irf.row(cumulate(ci)) = arma::cumsum(irf.row(cumulate(ci)));
            out = arma::vectorise(irf);
            return true;
        }, verbose);
    fSVARBootBands_cpp(r, conf, conf2, n_threads, "fBootstrapUhlig");
    return r;
}

static arma::uvec uhlig_args(int N, int idx, Rcpp::IntegerVector cumulate) {
    if (N < 2)
        Rcpp::stop("Uhlig max-share identification requires at least two variables.");
    if (idx < 1 || idx > N)
        Rcpp::stop("'idx' is out of range.");
    arma::uvec cumulate_cpp = (cumulate.size() == 0)
        ? arma::uvec() : arma::uvec(Rcpp::as<arma::uvec>(cumulate) - 1);
    if (!cumulate_cpp.is_empty() && cumulate_cpp.max() >= static_cast<arma::uword>(N))
        Rcpp::stop("'cumulate' contains an out-of-range variable index.");
    return cumulate_cpp;
}

//' Bootstrap Uhlig Maximum-Share Impulse Responses
//'
//' Computes residual-bootstrap confidence bands for the identification
//' used by \code{\link{fUhligIRF}}. The corrected variant first estimates
//' coefficient bias, shrinking the correction if needed for VAR stability,
//' and then bootstraps the bias-corrected VAR.
//'
//' @inheritParams fBootstrapMax
//' @param idx Index (1-based) of the variable whose forecast error variance
//'   contribution is maximised.
//' @return A list with \code{bootuhlig}, an N x (horizon+1) x nboot array;
//'   \code{upper}, \code{lower}, \code{upper2}, and \code{lower2}, each
//'   an N x (horizon+1) matrix; and \code{boot_beta}, an N x K x nboot
//'   coefficient array, where K is the number of regressors; and
//'   \code{conf}, the band levels \code{c(conf, conf2)}. The corrected
//'   variant stores \code{nboot2} draws.
//' @export
// [[Rcpp::export]]
Rcpp::List fBootstrapUhlig(const arma::mat& y, const Rcpp::List& var_result,
                            int nboot, int horizon, int idx,
                            double conf = 90.0, double conf2 = 68.0,
                            Rcpp::IntegerVector cumulate = Rcpp::IntegerVector::create(),
                            Rcpp::Nullable<arma::mat> exog = R_NilValue,
                            int n_threads = 0, bool verbose = false) {
  VARResult vr = fVARResultFromList_cpp(var_result);
  fCheckExogArg_cpp(vr.n_exog, exog.isNotNull());
  arma::mat exog_mat;
  if (exog.isNotNull()) exog_mat = Rcpp::as<arma::mat>(exog);
  const arma::mat* exog_ptr = exog.isNotNull() ? &exog_mat : nullptr;
  const int N = static_cast<int>(y.n_cols);
  arma::uvec cumulate_cpp = uhlig_args(N, idx, cumulate);

  SVARBootResult r = fBootstrapUhlig_cpp(y, vr, nboot, horizon, idx - 1, conf,
                                         conf2, cumulate_cpp, exog_ptr, n_threads, verbose);
  return fSVARBootList_cpp(r, "bootuhlig", {N, horizon + 1});
}

//' @rdname fBootstrapUhlig
//' @inheritParams fBootstrapMax
//' @export
// [[Rcpp::export]]
Rcpp::List fBootstrapUhligCorrected(const arma::mat& y, const Rcpp::List& var_result,
                                     int nboot1, int nboot2, int horizon, int idx,
                                     double conf = 90.0, double conf2 = 68.0,
                                     Rcpp::IntegerVector cumulate = Rcpp::IntegerVector::create(),
                                     Rcpp::Nullable<arma::mat> exog = R_NilValue,
                                     int n_threads = 0, bool verbose = false) {
  VARResult vr = fVARResultFromList_cpp(var_result);
  fCheckExogArg_cpp(vr.n_exog, exog.isNotNull());
  arma::mat exog_mat;
  if (exog.isNotNull()) exog_mat = Rcpp::as<arma::mat>(exog);
  const arma::mat* exog_ptr = exog.isNotNull() ? &exog_mat : nullptr;
  const int N = static_cast<int>(y.n_cols);
  arma::uvec cumulate_cpp = uhlig_args(N, idx, cumulate);
  fCheckBootArgs_cpp(y, vr, nboot2, horizon, conf, conf2);

  VARResult corrected = fBiasCorrectedVAR_cpp(y, vr, nboot1, "residual",
                                              exog_ptr, n_threads, verbose);
  SVARBootResult r = fBootstrapUhlig_cpp(y, corrected, nboot2, horizon, idx - 1,
                                         conf, conf2, cumulate_cpp, exog_ptr,
                                         n_threads, verbose);
  return fSVARBootList_cpp(r, "bootuhlig", {N, horizon + 1});
}
