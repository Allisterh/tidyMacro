// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// Prevent Armadillo from spawning its own OpenMP threads inside our parallel region
#define ARMA_DONT_USE_OPENMP

#include "fHDIV.h"
#include "fGetShock.h"
#include "fCompanionMatrix.h"
#include "fOLS.h"
#include "fVAR.h"
#include "fBootstrapVAR.h"
#include "fBootstrapIVMBB.h"
#include "fGetBands.h"
#include "fMBBVAR.h"
#include <RcppArmadillo.h>

// =========================================================================
// Internal helper: companion-form HD recursion
//
// Given unit-normalised structural impact vector s and VAR parameters,
// recovers the structural shock via Stock-Watson (2018) and propagates it
// through the companion form to obtain the T_eff x N HD matrix.
//
// HD_big[:, t] = F * HD_big[:, t-1]  +  (s * V[t-1]) padded with zeros
// HDshock      = HD_big[0:N-1, 1:T_eff].t()   →  T_eff x N
// =========================================================================
static arma::mat hd_recursion_iv(const arma::mat& residuals,
                                  const arma::mat& sigma,
                                  const arma::vec& s,
                                  const arma::mat& beta,
                                  int c, int p) {

    const int N       = static_cast<int>(residuals.n_cols);
    const int T_eff   = static_cast<int>(residuals.n_rows);
    const int nvarXeq = N * p;

    // Structural shock series (unit normalisation, shockSize = 1)
    arma::vec V = fGetShock_cpp(residuals, sigma, s, 1.0, "unit");

    // Companion matrix F  (nvarXeq x nvarXeq)
    CompanionMatrixResult comp = fCompanionMatrix_cpp(beta, c, p);
    const arma::mat& F = comp.comp;

    // Companion-form recursion
    arma::mat HD_big(nvarXeq, T_eff + 1, arma::fill::zeros);

    for (int i = 1; i <= T_eff; ++i) {
        HD_big.col(i)         = F * HD_big.col(i - 1);
        HD_big.col(i).head(N) += s * V(i - 1);
    }

    // Return T_eff x N  (drop the t=0 initialisation column, keep top N rows)
    return HD_big.head_rows(N).cols(1, T_eff).t();
}


// =========================================================================
// fHDIV_cpp  —  point-estimate companion-form HD for IV-identified shock
// =========================================================================
HDIVResult fHDIV_cpp(const arma::mat& residuals,
                     const arma::mat& sigma,
                     const arma::vec& s,
                     const arma::mat& beta,
                     int c, int p) {
    HDIVResult result;
    result.HDshock = hd_recursion_iv(residuals, sigma, s, beta, c, p);
    return result;
}


// =========================================================================
// fBootstrapHDIV_cpp  —  MBB bootstrap uncertainty bands for the IV HD
//
// For each bootstrap draw:
//   1. Resample residuals and instrument via MBB
//   2. Generate and re-estimate bootstrap VAR
//   3. Re-identify the IV shock (first-stage OLS, second-stage projection)
//   4. Run companion-form HD recursion
//
// Quantiles are recentered around the point estimate following:
//   Kaenzig (2021) / Gertler-Karadi MATLAB convention:
//   upper = q(1-a/2) - median(boot) + HDshock_point
// =========================================================================
BootstrapHDIVResult fBootstrapHDIV_cpp(const arma::mat& y,
                                        const VARResult& var_result,
                                        const arma::mat& Z,
                                        const arma::vec& s,
                                        int nboot, int blocksize,
                                        const arma::ivec& adjustZ,
                                        const arma::ivec& adjustu,
                                        int policyvar,
                                        double conf,
                                        int n_threads,
                                        bool verbose) {

    const int N      = static_cast<int>(y.n_cols);
    const int T_eff  = static_cast<int>(var_result.residuals.n_rows);
    const int p      = var_result.p;
    const int c      = var_result.c;
    const int pv_idx = policyvar - 1;   // 0-based

    // Convert 1-based R indices to 0-based C++ indices
    const int Z_start = adjustZ(0) - 1;
    const int Z_end   = adjustZ(1) - 1;
    const int u_start = adjustu(0) - 1;
    const int u_end   = adjustu(1) - 1;

    const arma::mat Zstar   = Z.rows(Z_start, Z_end);
    const arma::mat Epsstar = var_result.residuals.rows(u_start, u_end);

    // Point-estimate HD (used for recentering)
    arma::mat HDshock_pt = hd_recursion_iv(var_result.residuals,
                                            var_result.sigma,
                                            s, var_result.beta, c, p);

    // Blocks and centering are fixed; block indices are drawn on the main thread
    const MBBSamplerState mbb = fmbb_precompute_cpp(Epsstar, p, blocksize, &Zstar);
    const arma::imat index    = fmbb_draw_cpp(mbb, nboot);

    const int n_thr = fResolveThreads_cpp(n_threads);
    if (verbose)
        Rprintf("Using %d thread(s) for bootstrap HD computation...\n", n_thr);

    // Draw b's T_eff x N HD is column b of the flat (T_eff*N) x nboot matrix
    arma::mat flat(T_eff * N, nboot, arma::fill::none);
    std::vector<unsigned char> ok(nboot, 0);   // per-draw success, no R API

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_thr)
#endif
    for (int b = 0; b < nboot; ++b) {
        MBBVARResult draw = fmbb_sample_cpp(mbb, index.col(b));

        // Resampled residuals inside the proxy window, original ones outside
        arma::mat res_boot = fMBBPadResiduals_cpp(var_result.residuals,
                                                  draw.eps_boot, u_start);
        arma::mat yb = fSimulateVAR_cpp(y, var_result, res_boot, nullptr);
        VARResult var_b = fVAR_cpp(yb, p, c);
        if (!var_b.ok) continue;               // rank-deficient resample

        arma::vec s_b;
        if (!fProxyImpact_cpp(var_b.residuals.rows(u_start, u_end),
                              draw.M_boot, pv_idx, s_b)) continue;

        flat.col(b) = arma::vectorise(hd_recursion_iv(var_b.residuals, var_b.sigma,
                                                      s_b, var_b.beta, c, p));
        ok[b] = 1;
    }

    const arma::mat kept = fKeepValidDraws_cpp(flat, ok, "fBootstrapHDIV");

    // Percentiles recentred on the point estimate (Kaenzig / Gertler-Karadi):
    // band = q(pct) - median(boot) + HDshock_point
    arma::mat q = fRowPercentiles_cpp(
        kept, {50.0 + conf * 0.5, 50.0 - conf * 0.5, 50.0}, n_threads);
    arma::mat shift = arma::reshape(q.col(2), T_eff, N);

    BootstrapHDIVResult result;
    result.HDshock = HDshock_pt;
    result.upper   = arma::reshape(q.col(0), T_eff, N) - shift + HDshock_pt;
    result.lower   = arma::reshape(q.col(1), T_eff, N) - shift + HDshock_pt;
    return result;
}


// =========================================================================
// R wrappers
// =========================================================================

//' Companion-Form Historical Decomposition for IV-Identified Shock
//'
//' Computes the companion-form historical decomposition attributing each
//' observation of all VAR variables to the IV-identified structural shock,
//' following Stock and Watson (2018).
//'
//' @param residuals (T-p) x N matrix of VAR residuals from \code{fVAR}.
//' @param sigma N x N residual covariance matrix (e.g. \code{var_result$sigma}).
//' @param s N x 1 unit-normalised structural impact vector (first element = 1).
//' @param beta (Np+c) x N VAR coefficient matrix from \code{fVAR}.
//' @param c Integer intercept indicator (1 = include, 0 = exclude).
//' @param p Integer VAR lag order.
//'
//' @return A list with element:
//'   \item{HDshock}{(T-p) x N matrix. Column \code{j} is the contribution
//'     of the IV shock to variable \code{j} at each point in time.}
//'
//' @examples
//' set.seed(123)
//' y          <- matrix(rnorm(400), ncol = 2)
//' var_result <- fVAR(y, p = 2, c = 1)
//'
//' S  <- t(chol(var_result$sigma))
//' s  <- S[, 1] / S[1, 1]
//' hd <- fHDIV(var_result$residuals, var_result$sigma, s,
//'             var_result$beta, c = 1, p = 2)
//' dim(hd$HDshock)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fHDIV(const arma::mat& residuals,
                 const arma::mat& sigma,
                 const arma::vec& s,
                 const arma::mat& beta,
                 int c, int p) {

    HDIVResult result = fHDIV_cpp(residuals, sigma, s, beta, c, p);
    return Rcpp::List::create(Rcpp::Named("HDshock") = result.HDshock);
}


//' Bootstrap Confidence Bands for the IV Historical Decomposition
//'
//' Computes MBB bootstrap uncertainty bands for the companion-form historical
//' decomposition of an IV-identified structural shock. In each replication
//' the VAR is re-estimated, the IV identification is re-run, and the
//' companion-form HD recursion is applied to the bootstrap draw. Bands are
//' recentered around the point estimate following Kaenzig (2021).
//'
//' @param y T x N matrix of original endogenous variables.
//' @param var_result List from \code{fVAR}.
//' @param Z Instrument matrix aligned to the proxy sample.
//' @param s N x 1 unit-normalised point-estimate structural impact vector.
//' @param nboot Number of bootstrap replications.
//' @param blocksize MBB block size.
//' @param adjustZ Integer vector \code{c(start, end)} selecting the proxy-sample
//'   rows of \code{Z} (1-based).
//' @param adjustu Integer vector \code{c(start, end)} selecting the proxy-sample
//'   rows of the residuals (1-based).
//' @param policyvar Integer (1-based) index of the IV policy variable. Default 1.
//' @param conf Confidence level in percent (e.g. 90 for 90\% CI). Default 90.
//' @param n_threads OpenMP threads. 0 = all cores minus one. Default 0.
//' @param verbose Logical. If \code{TRUE}, report the number of threads used.
//'   Default \code{FALSE}, so the function is silent.
//'
//' @return A list with elements:
//'   \item{HDshock}{(T-p) x N point-estimate HD matrix.}
//'   \item{upper}{(T-p) x N upper confidence bands.}
//'   \item{lower}{(T-p) x N lower confidence bands.}
//'
//' @examples
//' set.seed(123)
//' y          <- matrix(rnorm(400), ncol = 2)
//' var_result <- fVAR(y, p = 2, c = 1)
//'
//' # An instrument correlated with the first reduced-form residual
//' u    <- var_result$residuals
//' Z    <- matrix(u[, 1] + rnorm(nrow(u), sd = 0.5), ncol = 1)
//' rows <- c(1, nrow(u))
//' s    <- c(1, sum(Z * u[, 2]) / sum(Z * u[, 1]))
//'
//' # Few draws here to keep the example quick; use 1000 or more in practice.
//' boot <- fBootstrapHDIV(y, var_result, Z, s = s, nboot = 50, blocksize = 10,
//'                        adjustZ = rows, adjustu = rows, policyvar = 1,
//'                        n_threads = 2)
//' names(boot)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fBootstrapHDIV(const arma::mat& y,
                           const Rcpp::List& var_result,
                           const arma::mat& Z,
                           const arma::vec& s,
                           int nboot, int blocksize,
                           const arma::ivec& adjustZ,
                           const arma::ivec& adjustu,
                           int policyvar = 1,
                           double conf = 90.0,
                           int n_threads = 0,
                           bool verbose = false) {

    VARResult var_struct = fVARResultFromList_cpp(var_result);
    fRequireNoExog_cpp(var_struct.n_exog, "fBootstrapHDIV");
    fCheckDrawArgs_cpp(nboot, 0, conf, conf);
    if (policyvar < 1 || policyvar > static_cast<int>(y.n_cols))
        Rcpp::stop("'policyvar' must be between 1 and %d.",
                   static_cast<int>(y.n_cols));

    BootstrapHDIVResult result = fBootstrapHDIV_cpp(
        y, var_struct, Z, s, nboot, blocksize, adjustZ, adjustu,
        policyvar, conf, n_threads, verbose);

    return Rcpp::List::create(
        Rcpp::Named("HDshock") = result.HDshock,
        Rcpp::Named("upper")   = result.upper,
        Rcpp::Named("lower")   = result.lower);
}
