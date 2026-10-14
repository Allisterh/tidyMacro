// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include <RcppArmadillo.h>
#include "fBootstrapIVMBB.h"
#include "fBootstrapVAR.h"
#include "fGetBands.h"
#include "fMBBVAR.h"
#include "fOLS.h"

arma::vec fProxyImpact_cpp(const arma::mat& u, const arma::mat& Z, int pv) {
  const int N = static_cast<int>(u.n_cols);

  // First stage: policy residual on the instrument
  arma::vec u_p = u.col(pv);
  arma::mat u_q = u;
  u_q.shed_col(pv);
  OLSResult ols1 = fOLS_cpp(u_p, Z, 0, 0, 0, 0);
  const arma::mat& uhat = ols1.fitted_partial;

  // Second stage: non-policy residuals on the fitted policy residual
  arma::mat sq_sp = arma::solve(uhat.t() * uhat, uhat.t() * u_q);

  arma::vec s(N, arma::fill::none);
  int counter = 0;
  for (int i = 0; i < N; ++i)
    s(i) = (i == pv) ? 1.0 : sq_sp(counter++);
  return s;
}

arma::mat fMBBPadResiduals_cpp(const arma::mat& residuals, const arma::mat& eps_boot,
                               int u_start) {
  arma::mat res_boot = residuals;
  res_boot.rows(u_start, u_start + eps_boot.n_rows - 1) = eps_boot;
  return res_boot;
}

BootstrapIVMBBResult
fBootstrapIVMBB_cpp(const arma::mat &y, const VARResult &var_result,
                     const arma::mat &Z, int nboot, int blocksize,
                     const arma::ivec &adjustZ, const arma::ivec &adjustu,
                     int policyvar, int horizon, double conf, double conf2,
                     const arma::mat* exog_ptr, int n_threads) {

  const arma::mat &residuals = var_result.residuals;
  const int p  = var_result.p;
  const int c  = var_result.c;
  const int N  = y.n_cols;
  const int H  = horizon + 1;
  const int pv = policyvar - 1;

  // Balance IV and residuals to overlap (1-based R indices to 0-based)
  const int Z_start = adjustZ(0) - 1;
  const int Z_end   = adjustZ(1) - 1;
  const int u_start = adjustu(0) - 1;
  const int u_end   = adjustu(1) - 1;
  const arma::mat Zstar   = Z.rows(Z_start, Z_end);
  const arma::mat Epsstar = residuals.rows(u_start, u_end);

  // Blocks and centering are fixed; block indices are drawn on the main thread
  const MBBSamplerState mbb = fmbb_precompute_cpp(Epsstar, p, blocksize, &Zstar);
  const arma::imat index    = fmbb_draw_cpp(mbb, nboot);

  const int n_thr = fResolveThreads_cpp(n_threads);
  Rprintf("Using %d thread(s) for parallel bootstrap computation...\n", n_thr);

  // Draw b's N x H responses are column b of the flat (N*H) x nboot matrix
  arma::mat flat(N * H, nboot, arma::fill::none);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_thr)
#endif
  for (int b = 0; b < nboot; ++b) {
    MBBVARResult draw = fmbb_sample_cpp(mbb, index.col(b));

    // Resampled residuals inside the proxy window, original ones outside
    arma::mat res_boot = fMBBPadResiduals_cpp(residuals, draw.eps_boot, u_start);
    arma::mat yb = fSimulateVAR_cpp(y, var_result, res_boot, exog_ptr);
    VARResult var_b = fVAR_cpp(yb, p, c, exog_ptr);

    arma::vec s = fProxyImpact_cpp(var_b.residuals.rows(u_start, u_end),
                                   draw.M_boot, pv);
    WoldIRFResult wold_b = fWoldIRF_cpp(var_b, horizon);
    for (int h = 0; h < H; ++h)
      flat.col(b).rows(h * N, h * N + N - 1) = wold_b.irfwold.slice(h) * s;
  }

  arma::mat q = fRowPercentiles_cpp(flat, fBandLevels_cpp(conf, conf2), n_threads);

  BootstrapIVMBBResult result;
  result.upper     = arma::reshape(q.col(0), N, H);
  result.lower     = arma::reshape(q.col(1), N, H);
  result.upper2    = arma::reshape(q.col(2), N, H);
  result.lower2    = arma::reshape(q.col(3), N, H);
  result.medianirf = arma::reshape(q.col(4), N, H);
  result.meanirf   = arma::reshape(arma::mean(flat, 1), N, H);
  return result;
}

//' Bootstrap IV Impulse Response Functions using Moving Block Bootstrap
//'
//' Computes bootstrap confidence bands for instrumental variable identified
//' impulse response functions using the moving block bootstrap to preserve
//' temporal dependence in residuals and instruments.
//'
//' @param y T x N matrix of endogenous variables
//' @param var_result List output from fVAR containing beta, residuals, p, c
//' @param Z (T-p) x K matrix of instrumental variables
//' @param nboot Integer number of bootstrap replications
//' @param blocksize Integer block size for moving block bootstrap
//' @param adjustZ Integer vector of length 2: start and end row indices
//'   (1-based, inclusive) for instrument alignment.
//' @param adjustu Integer vector of length 2: start and end row indices
//'   (1-based, inclusive) for residual alignment.
//' @param policyvar Integer index (1-based) of the policy variable.
//' @param horizon Integer maximum impulse response horizon.
//' @param conf Confidence level in percent (default 90).
//' @param conf2 Secondary confidence level in percent (default 68).
//' @param exog Optional matrix of exogenous variables (T x M). Default NULL.
//' @param n_threads Integer number of threads for parallel computation.
//'   Default is 0 (uses all available cores minus one, with a minimum of one).
//'   Set to 1 for single-threaded execution. If OpenMP is not available,
//'   automatically falls back to single-threaded execution.
//'
//' @return A list containing:
//'   \itemize{
//'     \item upper: N x (horizon+1) matrix of upper confidence bands
//'     \item lower: N x (horizon+1) matrix of lower confidence bands
//'     \item upper2: N x (horizon+1) matrix of secondary upper bands
//'     \item lower2: N x (horizon+1) matrix of secondary lower bands
//'     \item meanirf: N x (horizon+1) matrix of mean impulse responses
//'     \item medianirf: N x (horizon+1) matrix of median impulse responses
//'   }
//'
//' @details
//' This function implements the moving block bootstrap for IV-identified SVARs.
//' The first stage regresses the policy variable residual on the instrument(s),
//' and the second stage recovers the structural impact matrix. The
//' normalization sets the policy variable shock to have unit impact on itself.
//'
//' The function uses OpenMP for parallel computation when available,
//' speeding up bootstrap iterations. Each bootstrap replication and percentile
//' computation is independent and can be parallelized. If OpenMP is not
//' available, the function automatically falls back to single-threaded execution.
//'
//' @examples
//' \dontrun{
//' var_result <- fVAR(y, p = 2, c = 1)
//' result <- fBootstrapIVMBB(y, var_result, Z,
//'                             nboot = 1000, blocksize = 10,
//'                             adjustZ = c(1, 100), adjustu = c(1, 100),
//'                             policyvar = 1, horizon = 20, conf = 68)
//'
//' # Use 4 threads for parallel computation
//' result <- fBootstrapIVMBB(y, var_result, Z,
//'                             nboot = 1000, blocksize = 10,
//'                             adjustZ = c(1, 100), adjustu = c(1, 100),
//'                             policyvar = 1, horizon = 20, conf = 68,
//'                             n_threads = 4)
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fBootstrapIVMBB(const arma::mat &y, const Rcpp::List &var_result,
                            const arma::mat &Z, int nboot, int blocksize,
                            const arma::ivec &adjustZ,
                            const arma::ivec &adjustu, int policyvar,
                            int horizon, double conf = 90.0, double conf2 = 68.0,
                            Rcpp::Nullable<arma::mat> exog = R_NilValue,
                            int n_threads = 0) {
  VARResult vr = fVARResultFromList_cpp(var_result);
  fCheckExogArg_cpp(vr.n_exog, exog.isNotNull());
  arma::mat exog_mat;
  if (exog.isNotNull()) exog_mat = Rcpp::as<arma::mat>(exog);

  BootstrapIVMBBResult result =
      fBootstrapIVMBB_cpp(y, vr, Z, nboot, blocksize, adjustZ, adjustu,
                          policyvar, horizon, conf, conf2,
                          exog.isNotNull() ? &exog_mat : nullptr, n_threads);

  return Rcpp::List::create(Rcpp::Named("upper")     = result.upper,
                            Rcpp::Named("lower")     = result.lower,
                            Rcpp::Named("upper2")    = result.upper2,
                            Rcpp::Named("lower2")    = result.lower2,
                            Rcpp::Named("meanirf")   = result.meanirf,
                            Rcpp::Named("medianirf") = result.medianirf);
}
