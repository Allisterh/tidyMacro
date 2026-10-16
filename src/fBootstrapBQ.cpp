// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include <RcppArmadillo.h>
#include "fBootstrapBQ.h"
#include "fBootstrapVAR.h"
#include "fBQIRF.h"

SVARBootResult fBootstrapBQ_cpp(const arma::mat& y, const VARResult& var_result,
                                int nboot, int horizon, double conf, double conf2,
                                const std::string& bootscheme,
                                const arma::uvec& cumulate,
                                bool has_scaling, int scale_idx, double scale_size,
                                const arma::mat* exog_ptr, int n_threads, bool verbose) {
  fCheckBootArgs_cpp(y, var_result, nboot, horizon, conf, conf2);
  const int N = static_cast<int>(y.n_cols);
  const int H = horizon + 1;

  SVARBootResult r = fBootstrapSVAR_cpp(
      y, var_result, nboot, horizon, bootscheme, exog_ptr, n_threads,
      N * N * H, "Blanchard-Quah",
      [&](const VARResult& var_b, const arma::cube& wold_b, arma::vec& out) {
        // BQ identification: C1 = sum of Wold slices, D1 = lower chol of the
        // long-run covariance, K = C1 \ D1
        arma::mat C1(N, N, arma::fill::zeros);
        for (int h = 0; h < H; ++h) C1 += wold_b.slice(h);
        bool chol_ok = true;
        const arma::mat D1 = fLowerChol_cpp(C1 * var_b.sigma * C1.t(), &chol_ok);
        if (!chol_ok) return false;
        arma::mat K;
        if (!arma::solve(K, C1, D1)) return false;

        arma::cube irf = fBQIRF_cpp(wold_b, K, has_scaling, scale_idx, scale_size);

        // Cumulate selected variables along the horizon
        for (arma::uword ci = 0; ci < cumulate.n_elem; ++ci)
          for (int j = 0; j < N; ++j) {
            arma::vec tmp = irf.tube(cumulate(ci), j);
            irf.tube(cumulate(ci), j) = arma::cumsum(tmp);
          }
        out = arma::vectorise(irf);
        return true;
      }, verbose);
  fSVARBootBands_cpp(r, conf, conf2, n_threads, "fBootstrapBQ");
  return r;
}

//' Bootstrap Blanchard-Quah Long-Run Identified Impulse Response Functions
//'
//' @param y T x N matrix of original endogenous variables.
//' @param var_result List from \code{fVAR()}.
//' @param nboot Number of bootstrap replications.
//' @param horizon Maximum IRF horizon.
//' @param conf Confidence level in percent (e.g. 68).
//' @param conf2 Secondary confidence level in percent (default 68).
//' @param bootscheme \code{"residual"} or \code{"wild"}.
//' @param cumulate Integer vector (1-based) of variable indices whose IRFs
//'   should be cumulated along the horizon. Typically used when the VAR is
//'   estimated in first differences and level responses are required.
//' @param scaling Optional numeric vector of length 2. First element is the
//'   variable index (1-based) for normalisation; second is the shock size.
//'   Default NULL (no normalisation).
//' @param exog Optional T x M matrix of exogenous variables. Must be provided
//'   if the original VAR was estimated with exogenous variables; must be
//'   omitted otherwise.
//' @param verbose Logical. If \code{TRUE}, report the number of threads
//'   used. Default \code{FALSE}, so the function is silent.
//' @param n_threads Number of OpenMP threads. 0 = all cores minus one.
//'   Default 0.
//'
//' @return List with elements:
//'   \describe{
//'     \item{bootbq}{N x N x (horizon+1) x nboot array of bootstrapped BQ IRFs}
//'     \item{upper}{N x N x (horizon+1) upper confidence bands}
//'     \item{lower}{N x N x (horizon+1) lower confidence bands}
//'     \item{boot_beta}{N x (Np+c) x nboot array of bootstrapped coefficients}
//'     \item{conf}{Band levels in percent, \code{c(conf, conf2)}}
//'   }
//'
//' @details
//' Implements the residual or wild bootstrap for BQ long-run identified VARs.
//' In each replication the full BQ identification is re-computed:
//' \eqn{C_1 = \sum_h \Psi_h}, \eqn{D_1 = \text{chol}(C_1 \Sigma C_1')^\top},
//' \eqn{K = C_1^{-1} D_1}. Selected variables are then cumulated if the VAR
//' is estimated in differences.
//'
//' @references
//' Blanchard, O. J., & Quah, D. (1989). The dynamic effects of aggregate
//' demand and supply disturbances. \emph{American Economic Review}, 79(4),
//' 655--673.
//'
//' Gali, J. (1999). Technology, employment, and the business cycle: Do
//' technology shocks explain aggregate fluctuations?
//' \emph{American Economic Review}, 89(1), 249--271.
//'
//' @examples
//' set.seed(123)
//' y          <- matrix(rnorm(400), ncol = 2)
//' var_result <- fVAR(y, p = 2, c = 1)
//' wold       <- fWoldIRF(var_result, horizon = 12)
//' C1         <- apply(wold, c(1, 2), sum)
//' D1         <- t(chol(C1 %*% var_result$sigma %*% t(C1)))
//' point_irf  <- fBQIRF(wold, solve(C1, D1))
//'
//' # Few draws here to keep the example quick; use 1000 or more in practice.
//' boot <- fBootstrapBQ(y, var_result, nboot = 50, horizon = 12,
//'                      conf = 68, bootscheme = "residual",
//'                      cumulate = c(1L, 2L), n_threads = 2)
//' fPlotIRF(point_irf, bands = boot, shock = 1)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fBootstrapBQ(const arma::mat &y, const Rcpp::List &var_result,
                        int nboot, int horizon, double conf = 90.0, double conf2 = 68.0,
                        const std::string &bootscheme = "residual",
                        Rcpp::IntegerVector cumulate = Rcpp::IntegerVector::create(),
                        Rcpp::Nullable<arma::vec> scaling = R_NilValue,
                        Rcpp::Nullable<arma::mat> exog = R_NilValue,
                        int n_threads = 0, bool verbose = false) {
  VARResult vr = fVARResultFromList_cpp(var_result);
  fCheckExogArg_cpp(vr.n_exog, exog.isNotNull());
  arma::mat exog_mat;
  if (exog.isNotNull()) exog_mat = Rcpp::as<arma::mat>(exog);

  const int N = static_cast<int>(y.n_cols);
  arma::uvec cumulate_cpp = (cumulate.size() == 0)
      ? arma::uvec()
      : Rcpp::as<arma::uvec>(cumulate) - 1;
  if (!cumulate_cpp.is_empty() && cumulate_cpp.max() >= static_cast<arma::uword>(N))
    Rcpp::stop("'cumulate' contains an out-of-range variable index.");

  // Parse scaling once on the main thread (no Rcpp::as inside OpenMP)
  const bool has_scaling = scaling.isNotNull();
  int    scale_idx  = 0;
  double scale_size = 0.0;
  if (has_scaling) {
    arma::vec s = Rcpp::as<arma::vec>(scaling);
    scale_idx  = static_cast<int>(s(0)) - 1;
    scale_size = s(1);
  }

  SVARBootResult r = fBootstrapBQ_cpp(y, vr, nboot, horizon, conf, conf2,
                                      bootscheme, cumulate_cpp, has_scaling,
                                      scale_idx, scale_size,
                                      exog.isNotNull() ? &exog_mat : nullptr,
                                      n_threads, verbose);
  return fSVARBootList_cpp(r, "bootbq", {N, N, horizon + 1});
}
