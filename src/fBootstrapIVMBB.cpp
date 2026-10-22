// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include <RcppArmadillo.h>
#include "fBootstrapIVMBB.h"
#include "fBootstrapVAR.h"
#include "fGetBands.h"
#include "fMBBVAR.h"
#include "fOLS.h"

bool fProxyImpact_cpp(const arma::mat& u, const arma::mat& Z, int pv,
                      arma::vec& s) {
  const int N = static_cast<int>(u.n_cols);

  // First stage: policy residual on the instrument and a constant. The
  // demeaned instrument gives the same slopes, and its fitted values are
  // demeaned, so the second stage needs no constant either.
  arma::vec u_p = u.col(pv);
  arma::mat u_q = u;
  u_q.shed_col(pv);
  const arma::mat Zc = Z.each_row() - arma::mean(Z, 0);
  OLSResult ols1 = fOLS_cpp(u_p, Zc, 0, 0, 0, 0);
  const arma::mat& uhat = ols1.fitted_partial;

  // Second stage: non-policy residuals on the fitted policy residual
  arma::mat sq_sp;
  if (!arma::solve(sq_sp, uhat.t() * uhat, uhat.t() * u_q,
                   arma::solve_opts::fast))
    return false;
  if (!sq_sp.is_finite()) return false;

  s.set_size(N);
  int counter = 0;
  for (int i = 0; i < N; ++i)
    s(i) = (i == pv) ? 1.0 : sq_sp(counter++);
  return true;
}

arma::mat fMBBPadResiduals_cpp(const arma::mat& residuals, const arma::mat& eps_boot,
                               int u_start) {
  arma::mat res_boot = residuals;
  res_boot.rows(u_start, u_start + eps_boot.n_rows - 1) = eps_boot;
  return res_boot;
}

// runProxyVAR.m redraws a resample whose identification window holds fewer
// than this many non-zero proxy values
static const int kMinProxyObs = 15;

// Block indices for nboot accepted draws (nBlock x nboot). A draw is redrawn
// when every instrument column has fewer than kMinProxyObs non-zero values
// in rows u_start..u_end. Uses R's RNG: main thread only.
static arma::imat fIVMBBDraw_cpp(const MBBSamplerState& mbb, int nboot,
                                 int u_start, int u_end) {
  arma::imat index(mbb.nBlock, nboot, arma::fill::none);
  const long long max_tries = std::max(1000LL, 100LL * nboot);
  long long tries = 0;
  int b = 0;
  while (b < nboot) {
    if (++tries > max_tries)
      Rcpp::stop("fBootstrapIVMBB: too few resamples keep %d non-zero "
                 "instrument values in the identification window.",
                 kMinProxyObs);
    arma::ivec idx = fmbb_draw_cpp(mbb, 1).col(0);
    const arma::mat Mb = fmbb_sample_cpp(mbb, idx).M_boot.rows(u_start, u_end);
    bool enough = false;
    for (arma::uword col = 0; col < Mb.n_cols; ++col)
      if (arma::accu(Mb.col(col) != 0.0) >= kMinProxyObs) enough = true;
    if (enough) index.col(b++) = idx;
  }
  return index;
}

BootstrapIVMBBResult
fBootstrapIVMBB_cpp(const arma::mat &y, const VARResult &var_result,
                     const arma::mat &Z, int nboot, int blocksize,
                     const arma::ivec &adjustZ, const arma::ivec &adjustu,
                     int policyvar, int horizon, double conf, double conf2,
                     const arma::mat* exog_ptr, int n_threads, bool verbose) {

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
  if (Z_end - Z_start != u_end - u_start)
    Rcpp::stop("fBootstrapIVMBB: 'adjustZ' and 'adjustu' must select the same "
               "number of rows.");
  const arma::mat Zstar   = Z.rows(Z_start, Z_end);
  const arma::mat Epsstar = residuals.rows(u_start, u_end);

  // Point estimate on the original sample, before any resampling.
  arma::vec s_point;
  if (!fProxyImpact_cpp(Epsstar, Zstar, pv, s_point))
    Rcpp::stop("fBootstrapIVMBB: the proxy-SVAR second stage is singular on "
               "the original sample; the instrument carries no usable "
               "variation for this policy variable.");
  const WoldIRFResult wold_pt = fWoldIRF_cpp(var_result, horizon);
  arma::mat point(N, H, arma::fill::none);
  for (int h = 0; h < H; ++h) point.col(h) = wold_pt.irfwold.slice(h) * s_point;

  // Kaenzig (2021): blocks run over all residuals, with the instrument set to
  // zero outside the identification window. Blocks and centering are fixed;
  // block indices are drawn on the main thread.
  arma::mat Zlong(residuals.n_rows, Zstar.n_cols, arma::fill::zeros);
  Zlong.rows(u_start, u_end) = Zstar;
  const MBBSamplerState mbb = fmbb_precompute_cpp(residuals, p, blocksize, &Zlong);
  const arma::imat index    = fIVMBBDraw_cpp(mbb, nboot, u_start, u_end);

  const int n_thr = fResolveThreads_cpp(n_threads);
  if (verbose)
    Rprintf("Using %d thread(s) for parallel bootstrap computation...\n", n_thr);

  // Draw b's N x H responses are column b of the flat (N*H) x nboot matrix
  arma::mat flat(N * H, nboot, arma::fill::none);
  std::vector<unsigned char> ok(nboot, 0);   // per-draw success, no R API

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_thr)
#endif
  for (int b = 0; b < nboot; ++b) {
    MBBVARResult draw = fmbb_sample_cpp(mbb, index.col(b));

    arma::mat yb = fSimulateVAR_cpp(y, var_result, draw.eps_boot, exog_ptr);
    VARResult var_b = fVAR_cpp(yb, p, c, exog_ptr);
    if (!var_b.ok) continue;                 // rank-deficient resample

    arma::vec s;
    if (!fProxyImpact_cpp(var_b.residuals.rows(u_start, u_end),
                          draw.M_boot.rows(u_start, u_end), pv, s)) continue;

    WoldIRFResult wold_b = fWoldIRF_cpp(var_b, horizon);
    for (int h = 0; h < H; ++h)
      flat.col(b).rows(h * N, h * N + N - 1) = wold_b.irfwold.slice(h) * s;
    ok[b] = 1;
  }

  const arma::mat kept = fKeepValidDraws_cpp(flat, ok, "fBootstrapIVMBB");
  arma::mat q = fRowPercentiles_cpp(kept, fBandLevels_cpp(conf, conf2), n_threads);

  // Percentiles recentred on the point estimate (Kaenzig 2021):
  // band = q(pct) - median(boot) + point
  const arma::mat shift = point - arma::reshape(q.col(4), N, H);

  BootstrapIVMBBResult result;
  result.upper  = arma::reshape(q.col(0), N, H) + shift;
  result.lower  = arma::reshape(q.col(1), N, H) + shift;
  result.upper2 = arma::reshape(q.col(2), N, H) + shift;
  result.lower2 = arma::reshape(q.col(3), N, H) + shift;
  result.point  = point;
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
//' @param Z Matrix of instrumental variables with K columns; rows
//'   \code{adjustZ} are used. Zeros are treated as months without an
//'   instrument observation.
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
//' @param verbose Logical. If \code{TRUE}, report the number of threads used.
//'   Default \code{FALSE}, so the function is silent.
//'
//' @return A list containing:
//'   \itemize{
//'     \item upper: N x (horizon+1) matrix of upper confidence bands
//'     \item lower: N x (horizon+1) matrix of lower confidence bands
//'     \item upper2: N x (horizon+1) matrix of secondary upper bands
//'     \item lower2: N x (horizon+1) matrix of secondary lower bands
//'
//'     All four bands are bootstrap percentiles recentred on the point
//'     estimate, percentile minus bootstrap median plus \code{point}, so
//'     \code{point} always lies inside them.
//'     \item point: N x (horizon+1) matrix of impulse responses on the
//'       original sample, computed as the Wold coefficients times the
//'       proxy-SVAR impact vector, exactly as in \code{fHeteroIRF()} and
//'       \code{fMSW()}. The bootstrap supplies the bands only, never the
//'       point estimate.
//'     \item conf: band levels in percent, \code{c(conf, conf2)}
//'   }
//'
//' @details
//' This function implements the moving block bootstrap for IV-identified
//' SVARs, following the replication code of Kaenzig (2021). The first stage
//' regresses the policy variable residual on the instrument(s) and a
//' constant, and the second stage regresses the other residuals on the
//' fitted values and a constant. The normalization sets the policy variable
//' shock to have unit impact on itself.
//'
//' Blocks are drawn from all VAR residuals jointly with the instrument, which
//' is set to zero outside the rows selected by \code{adjustu}; the
//' identification step uses those rows only. Non-zero instrument values are
//' centred so that they keep their sample mean. A resample with fewer than 15
//' non-zero instrument values in the identification window is drawn again.
//'
//' The function uses OpenMP for parallel computation when available,
//' speeding up bootstrap iterations. Each bootstrap replication and percentile
//' computation is independent and can be parallelized. If OpenMP is not
//' available, the function automatically falls back to single-threaded execution.
//'
//' @examples
//' set.seed(123)
//' y <- matrix(rnorm(400), ncol = 2)
//' var_result <- fVAR(y, p = 2, c = 1)
//'
//' # An instrument correlated with the first reduced-form residual
//' u <- var_result$residuals
//' Z <- matrix(u[, 1] + rnorm(nrow(u), sd = 0.5), ncol = 1)
//' rows <- c(1, nrow(u))
//'
//' # Few draws here to keep the example quick; use 1000 or more in practice.
//' result <- fBootstrapIVMBB(y, var_result, Z, nboot = 50, blocksize = 10,
//'                           adjustZ = rows, adjustu = rows,
//'                           policyvar = 1, horizon = 12, conf = 68,
//'                           n_threads = 2)
//' fPlotIRF(result, varnames = c("y1", "y2"))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fBootstrapIVMBB(const arma::mat &y, const Rcpp::List &var_result,
                            const arma::mat &Z, int nboot, int blocksize,
                            const arma::ivec &adjustZ,
                            const arma::ivec &adjustu, int policyvar,
                            int horizon, double conf = 90.0, double conf2 = 68.0,
                            Rcpp::Nullable<arma::mat> exog = R_NilValue,
                            int n_threads = 0, bool verbose = false) {
  VARResult vr = fVARResultFromList_cpp(var_result);
  fCheckExogArg_cpp(vr.n_exog, exog.isNotNull());
  fCheckBootArgs_cpp(y, vr, nboot, horizon, conf, conf2);
  if (blocksize < 0)
    Rcpp::stop("'blocksize' must be non-negative (0 selects it automatically).");
  if (policyvar < 1 || policyvar > static_cast<int>(y.n_cols))
    Rcpp::stop("'policyvar' must be between 1 and %d.",
               static_cast<int>(y.n_cols));
  arma::mat exog_mat;
  if (exog.isNotNull()) exog_mat = Rcpp::as<arma::mat>(exog);

  BootstrapIVMBBResult result =
      fBootstrapIVMBB_cpp(y, vr, Z, nboot, blocksize, adjustZ, adjustu,
                          policyvar, horizon, conf, conf2,
                          exog.isNotNull() ? &exog_mat : nullptr, n_threads,
                          verbose);

  return Rcpp::List::create(Rcpp::Named("upper")  = result.upper,
                            Rcpp::Named("lower")  = result.lower,
                            Rcpp::Named("upper2") = result.upper2,
                            Rcpp::Named("lower2") = result.lower2,
                            Rcpp::Named("point")  = result.point,
                            Rcpp::Named("conf")   = Rcpp::NumericVector::create(conf, conf2));
}
