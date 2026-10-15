// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include <RcppArmadillo.h>
#include "fBootstrapMax.h"
#include "fBootstrapVAR.h"
#include <cmath>

// Closed-form LR_max solution (replaces MATLAB fminsearch on LR_max.m).
// Maximises wold.slice(H-1).row(var_idx) * S * h2  (last-horizon response only,
// matching MATLAB irfs(var,:,end)) subject to h2(0) = 0, ||h2|| = 1.
// Optimal: h2 = [0; M_sub / ||M_sub||] where
// M_sub = (wold.slice(H-1).row(var_idx) * S).cols(1, N-1).
static arma::vec lr_max_solve(const arma::cube& wold, const arma::mat& S,
                               int var_idx) {
    const int H = static_cast<int>(wold.n_slices);
    const int N = static_cast<int>(wold.n_rows);

    arma::rowvec M_full = wold.slice(H - 1).row(var_idx) * S;
    arma::rowvec M_sub = M_full.cols(1, N - 1);
    double M_norm = arma::norm(M_sub);

    arma::vec h2(N, arma::fill::zeros);
    if (M_norm > 1e-14) {
        h2.rows(1, N - 1) = M_sub.t() / M_norm;
    }
    return h2;
}

SVARBootResult fBootstrapMax_cpp(const arma::mat& y, const VARResult& var_result,
                                 int nboot, int horizon, int var_idx,
                                 double conf, double conf2,
                                 const arma::uvec& cumulate,
                                 bool has_scaling, int scale_idx, double scale_size,
                                 const arma::mat* exog_ptr, int n_threads) {
    fCheckBootArgs_cpp(y, var_result, nboot, horizon, conf, conf2);
    const int N = static_cast<int>(y.n_cols);
    const int H = horizon + 1;

    SVARBootResult r = fBootstrapSVAR_cpp(
        y, var_result, nboot, horizon, "residual", exog_ptr, n_threads,
        N * H, "LR-Max",
        [&](const VARResult& var_b, const arma::cube& wold_b, arma::vec& out) {
            arma::mat S = fLowerChol_cpp(var_b.sigma);
            arma::vec impact = S * lr_max_solve(wold_b, S, var_idx);

            // Orient the shock so the selected variable's maximised response
            // is non-negative.
            if (arma::dot(wold_b.slice(H - 1).row(var_idx), impact) < 0.0)
                impact = -impact;

            arma::mat irf(N, H, arma::fill::none);
            for (int h = 0; h < H; ++h) irf.col(h) = wold_b.slice(h) * impact;

            if (has_scaling) {
                double scale_val = irf(scale_idx, 0) * scale_size;
                if (std::abs(scale_val) > 1e-14) irf /= scale_val;
            }
            for (arma::uword ci = 0; ci < cumulate.n_elem; ++ci)
                irf.row(cumulate(ci)) = arma::cumsum(irf.row(cumulate(ci)));
            out = arma::vectorise(irf);
        });
    fSVARBootBands_cpp(r, conf, conf2, n_threads);
    return r;
}

// R-side argument handling shared by fBootstrapMax and its corrected variant.
static void max_args(int N, int var_idx, Rcpp::IntegerVector cumulate,
                     Rcpp::Nullable<arma::vec> scaling, arma::uvec& cumulate_cpp,
                     bool& has_scaling, int& scale_idx, double& scale_size) {
    if (N < 2)
        Rcpp::stop("LR-Max identification requires at least two variables.");
    if (var_idx < 1 || var_idx > N)
        Rcpp::stop("'var_idx' is out of range.");
    cumulate_cpp = (cumulate.size() == 0) ? arma::uvec()
                                          : Rcpp::as<arma::uvec>(cumulate) - 1;
    if (!cumulate_cpp.is_empty() && cumulate_cpp.max() >= static_cast<arma::uword>(N))
        Rcpp::stop("'cumulate' contains an out-of-range variable index.");

    has_scaling = scaling.isNotNull();
    scale_idx = 0;
    scale_size = 0.0;
    if (has_scaling) {
        arma::vec s = Rcpp::as<arma::vec>(scaling);
        if (s.n_elem != 2)
            Rcpp::stop("'scaling' must contain a variable index and shock size.");
        scale_idx = static_cast<int>(s(0)) - 1;
        scale_size = s(1);
        if (scale_idx < 0 || scale_idx >= N)
            Rcpp::stop("The variable index in 'scaling' is out of range.");
        if (!std::isfinite(scale_size) || scale_size == 0.0)
            Rcpp::stop("The shock size in 'scaling' must be finite and nonzero.");
    }
}

//' Bootstrap Long-Run Maximum-Response Impulse Responses
//'
//' Computes residual-bootstrap confidence bands for the identification
//' used by \code{\link{fMaxIRF}}. The corrected variant first estimates
//' coefficient bias, shrinking the correction if needed for VAR stability,
//' and then bootstraps the bias-corrected VAR.
//'
//' @param y T x N matrix of endogenous variables.
//' @param var_result Fitted VAR returned by \code{\link{fVAR}}.
//' @param nboot Number of bootstrap replications.
//' @param horizon Maximum response horizon; horizon zero is included.
//' @param var_idx Index (1-based) of the variable whose last-horizon
//'   response is maximised.
//' @param conf Confidence level in percent (default 90).
//' @param conf2 Secondary confidence level in percent (default 68).
//' @param cumulate Integer vector of variable indices (1-based) whose
//'   responses should be cumulated. Defaults to no cumulation.
//' @param scaling Optional length-two numeric vector. Responses are divided
//'   by the impact response of variable \code{scaling[1]} multiplied by
//'   \code{scaling[2]}, when that product is nonzero. Default \code{NULL}.
//' @param exog T x M matrix of exogenous variables, required when the
//'   original VAR included them; otherwise \code{NULL}.
//' @param n_threads Number of OpenMP threads. Zero uses all available
//'   cores minus one, with a minimum of one. Without OpenMP, uses one.
//' @return A list with \code{bootmax}, an N x (horizon+1) x nboot array;
//'   \code{upper}, \code{lower}, \code{upper2}, and \code{lower2}, each
//'   an N x (horizon+1) matrix; and \code{boot_beta}, an N x K x nboot
//'   coefficient array, where K is the number of regressors; and
//'   \code{conf}, the band levels \code{c(conf, conf2)}. The corrected
//'   variant stores \code{nboot2} draws.
//' @export
// [[Rcpp::export]]
Rcpp::List fBootstrapMax(const arma::mat& y, const Rcpp::List& var_result,
                         int nboot, int horizon, int var_idx,
                         double conf = 90.0, double conf2 = 68.0,
                         Rcpp::IntegerVector cumulate = Rcpp::IntegerVector::create(),
                         Rcpp::Nullable<arma::vec> scaling = R_NilValue,
                         Rcpp::Nullable<arma::mat> exog    = R_NilValue,
                         int n_threads = 0) {
  VARResult vr = fVARResultFromList_cpp(var_result);
  fCheckExogArg_cpp(vr.n_exog, exog.isNotNull());
  arma::mat exog_mat;
  if (exog.isNotNull()) exog_mat = Rcpp::as<arma::mat>(exog);
  const arma::mat* exog_ptr = exog.isNotNull() ? &exog_mat : nullptr;
  const int N = static_cast<int>(y.n_cols);
  arma::uvec cumulate_cpp; bool has_scaling; int scale_idx; double scale_size;
  max_args(N, var_idx, cumulate, scaling, cumulate_cpp, has_scaling,
           scale_idx, scale_size);

  SVARBootResult r = fBootstrapMax_cpp(y, vr, nboot, horizon, var_idx - 1,
                                       conf, conf2, cumulate_cpp, has_scaling,
                                       scale_idx, scale_size, exog_ptr, n_threads);
  return fSVARBootList_cpp(r, "bootmax", {N, horizon + 1});
}

//' @rdname fBootstrapMax
//' @param nboot1 Number of first-pass replications used to estimate bias.
//' @param nboot2 Number of second-pass replications used for the bands.
//' @export
// [[Rcpp::export]]
Rcpp::List fBootstrapMaxCorrected(const arma::mat& y, const Rcpp::List& var_result,
                                   int nboot1, int nboot2, int horizon, int var_idx,
                                   double conf = 90.0, double conf2 = 68.0,
                                   Rcpp::IntegerVector cumulate = Rcpp::IntegerVector::create(),
                                   Rcpp::Nullable<arma::vec> scaling = R_NilValue,
                                   Rcpp::Nullable<arma::mat> exog    = R_NilValue,
                                   int n_threads = 0) {
  VARResult vr = fVARResultFromList_cpp(var_result);
  fCheckExogArg_cpp(vr.n_exog, exog.isNotNull());
  arma::mat exog_mat;
  if (exog.isNotNull()) exog_mat = Rcpp::as<arma::mat>(exog);
  const arma::mat* exog_ptr = exog.isNotNull() ? &exog_mat : nullptr;
  const int N = static_cast<int>(y.n_cols);
  arma::uvec cumulate_cpp; bool has_scaling; int scale_idx; double scale_size;
  max_args(N, var_idx, cumulate, scaling, cumulate_cpp, has_scaling,
           scale_idx, scale_size);
  fCheckBootArgs_cpp(y, vr, nboot2, horizon, conf, conf2);

  VARResult corrected = fBiasCorrectedVAR_cpp(y, vr, nboot1, "residual",
                                              exog_ptr, n_threads);
  SVARBootResult r = fBootstrapMax_cpp(y, corrected, nboot2, horizon, var_idx - 1,
                                       conf, conf2, cumulate_cpp, has_scaling,
                                       scale_idx, scale_size, exog_ptr, n_threads);
  return fSVARBootList_cpp(r, "bootmax", {N, horizon + 1});
}
