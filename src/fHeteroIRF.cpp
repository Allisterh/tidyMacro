// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include "fHeteroIRF.h"
#include "fVAR.h"
#include "fBootstrapVAR.h"
#include "fBootstrapIVMBB.h"
#include "fGetBands.h"
#include "fMBBVAR.h"
#include "fWoldIRF.h"
#include <RcppArmadillo.h>

// =========================================================================
// Internal helper: proxy-IV structural impact vector
// Matches MATLAB: b1 = [uhat] \ yiIV  (treatment group only)
//   XrIV = ZrIV = proxy(indsR1) - mean(proxy(indsR1))   [centred proxy]
//   yiIV = U(indsR1) - mean(U(indsR1))                  [centred residuals]
//   b1   = (XrIV' XrIV)^{-1} XrIV' yiIV                [OLS, no intercept]
//   then normalise so b1[nvar-1] = scale
// Returns an error message (nullptr on success) instead of throwing, so the
// bootstrap can call it on OpenMP workers.
// =========================================================================
static const char* hetero_iv_b1(const arma::mat& eta,    // T1 x N
                                const arma::mat& Z,      // T1 x 1
                                const arma::ivec& indsR1,
                                int nvar, double scale, arma::vec& b1) {
    arma::uvec idx = arma::find(indsR1 == 1);   // OPEC months
    if (idx.is_empty())
        return "fHeteroIRF: indsR1 contains no treatment observations.";

    arma::mat U_opec = eta.rows(idx);           // T_opec x N
    const arma::vec z_full = Z.col(0);
    arma::vec z_opec = z_full.elem(idx);        // T_opec x 1

    // Centre (MATLAB: proxy - mean(proxy), U - mean(U))
    z_opec -= arma::mean(z_opec);
    U_opec.each_row() -= arma::mean(U_opec, 0);

    // b1 = (z'z)^{-1} z'U  ->  N-vector
    const double zz = arma::dot(z_opec, z_opec);
    if (zz < 1e-14)
        return "fHeteroIRF: proxy has zero variance in treatment months.";
    b1 = (U_opec.t() * z_opec) / zz;

    // Normalise: b1[nvar-1] = scale
    const double nf = b1(nvar - 1);
    if (std::abs(nf) < 1e-14)
        return "fHeteroIRF: normalization element is numerically zero.";
    b1 *= (scale / nf);
    return nullptr;
}

// IRF.col(h) = Wold_h * b1
static arma::mat hetero_irf(const arma::mat& beta, int c, int p, int hor,
                            const arma::vec& b1) {
    arma::cube wold;
    fWoldIRF_into_cpp(beta, c, p, hor, wold);
    arma::mat IRF(b1.n_elem, hor + 1, arma::fill::none);
    for (int h = 0; h <= hor; ++h) IRF.col(h) = wold.slice(h) * b1;
    return IRF;
}

// =========================================================================
// fHeteroIRF_cpp: point estimate
// =========================================================================
HeteroIRFResult fHeteroIRF_cpp(const arma::mat& beta, int c,
                                const arma::mat& eta,
                                const arma::mat& Z,
                                const arma::ivec& indsR1,
                                int p, int hor, int nvar, double scale) {
    HeteroIRFResult res;
    const char* err = hetero_iv_b1(eta, Z, indsR1, nvar, scale, res.b1);
    if (err != nullptr) Rcpp::stop(err);
    res.IRF = hetero_irf(beta, c, p, hor, res.b1);
    return res;
}

// =========================================================================
// fBootstrapHetero_cpp: MBB bootstrap
// Resamples (residuals, Z) jointly in the proxy window — same as
// fBootstrapIVMBB — but identification uses OPEC months only.
// =========================================================================
BootHeteroResult fBootstrapHetero_cpp(const arma::mat& y,
                                      const VARResult& var_result,
                                      const arma::mat& Z,
                                      const arma::ivec& indsR1,
                                      const arma::ivec& adjustu,
                                      int nboot, int blocksize,
                                      int hor, int nvar, double scale,
                                      double conf, double conf2,
                                      int n_threads, bool verbose) {
    const arma::mat& residuals = var_result.residuals;
    const int p = var_result.p;
    const int c = var_result.c;
    const int N = (int)y.n_cols;
    const int H = hor + 1;

    const int u_start = adjustu(0) - 1;
    const int u_end   = adjustu(1) - 1;

    // Proxy-window residuals and proxy (Z is already aligned to the window)
    const arma::mat Epsstar = residuals.rows(u_start, u_end);   // T1 x N

    // Point estimate (centre of the Hall bands)
    const arma::mat irf_pt = fHeteroIRF_cpp(var_result.beta, c, Epsstar, Z,
                                            indsR1, p, hor, nvar, scale).IRF;

    // Blocks and centering are fixed; block indices are drawn on the main thread
    const MBBSamplerState mbb = fmbb_precompute_cpp(Epsstar, p, blocksize, &Z);
    const arma::imat index    = fmbb_draw_cpp(mbb, nboot);

    const int n_thr = fResolveThreads_cpp(n_threads);
    if (verbose)
        Rprintf("Using %d thread(s) for heteroskedasticity bootstrap...\n", n_thr);

    // Draw b's N x H responses are column b of the flat (N*H) x nboot matrix
    arma::mat flat(N * H, nboot, arma::fill::zeros);
    // A draw is kept if it has >= 15 treatment months and identifies
    std::vector<int> valid(nboot, 0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_thr)
#endif
    for (int b = 0; b < nboot; ++b) {
        MBBVARResult draw = fmbb_sample_cpp(mbb, index.col(b));
        const arma::mat& Z_boot = draw.M_boot;

        // Treatment months of the bootstrapped proxy (non-zero = OPEC month)
        arma::ivec indsR1_b = arma::conv_to<arma::ivec>::from(arma::abs(Z_boot.col(0)) > 1e-14);
        if (arma::accu(indsR1_b) < 15) continue;

        // Regenerate VAR data and re-estimate
        arma::mat res_boot = fMBBPadResiduals_cpp(residuals, draw.eps_boot, u_start);
        arma::mat yb = fSimulateVAR_cpp(y, var_result, res_boot, nullptr);
        VARResult var_b = fVAR_cpp(yb, p, c);

        // Identify using proxy IV in bootstrap OPEC months
        arma::vec b1;
        if (hetero_iv_b1(var_b.residuals.rows(u_start, u_end), Z_boot, indsR1_b,
                         nvar, scale, b1) != nullptr) continue;
        flat.col(b) = arma::vectorise(hetero_irf(var_b.beta, c, p, hor, b1));
        valid[b] = 1;
    }

    arma::uvec keep = arma::find(arma::conv_to<arma::ivec>::from(valid) == 1);
    if (keep.n_elem < 10)
        Rcpp::stop("fBootstrapHetero: fewer than 10 valid bootstrap iterations. "
                   "Check proxy and adjustu settings.");
    const arma::mat kept = flat.cols(keep);

    arma::mat q = fRowPercentiles_cpp(kept, fBandLevels_cpp(conf, conf2), n_threads);
    const arma::mat bmed  = arma::reshape(q.col(4), N, H);
    const arma::mat shift = irf_pt - bmed;

    // Hall (1992) recentered percentile bands
    BootHeteroResult res;
    res.upper     = arma::reshape(q.col(0), N, H) + shift;
    res.lower     = arma::reshape(q.col(1), N, H) + shift;
    res.upper2    = arma::reshape(q.col(2), N, H) + shift;
    res.lower2    = arma::reshape(q.col(3), N, H) + shift;
    res.point     = irf_pt;
    return res;
}

// =========================================================================
// R wrappers
// =========================================================================

//' Heteroskedasticity-Based VAR Identification (Proxy IV, Treatment Months)
//'
//' Identifies a structural shock using the external instrument restricted to
//' high-variance (treatment) months, following the proxy-SVAR approach of
//' Kaenzig (2021). The structural impact vector is estimated via OLS of the
//' centred residuals on the centred proxy, restricted to OPEC announcement
//' months (\code{indsR1 == 1}).
//'
//' @param var_result List returned by \code{fVAR()}.
//' @param Z          Proxy matrix (T1 x 1), aligned to the identification window
//'   defined by \code{adjustu}.
//' @param adjustu    Integer vector \code{c(start, end)} (1-based) giving the
//'   rows of \code{var_result$residuals} corresponding to the proxy window.
//' @param indsR1     Integer 0/1 vector of length \code{adjustu[2]-adjustu[1]+1};
//'   1 = treatment month (e.g. OPEC announcement), 0 otherwise.
//' @param hor        IRF horizon.
//' @param nvar       1-based normalization variable index.
//' @param scale      Shock size for normalization (default 10).
//'
//' @return A list with \code{IRF} (N x hor+1) and \code{b1} (N-vector).
//'
//' @examples
//' set.seed(123)
//' y          <- matrix(rnorm(400), ncol = 2)
//' var_result <- fVAR(y, p = 2, c = 1)
//'
//' # A proxy and a 0/1 indicator for the high-variance (treatment) periods
//' u      <- var_result$residuals
//' Z      <- matrix(u[, 1] + rnorm(nrow(u), sd = 0.5), ncol = 1)
//' rows   <- c(1, nrow(u))
//' indsR1 <- as.integer(abs(Z[, 1]) > stats::median(abs(Z[, 1])))
//'
//' irf <- fHeteroIRF(var_result, Z, rows, indsR1, hor = 12, nvar = 1,
//'                   scale = 1)
//' dim(irf$IRF)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fHeteroIRF(const Rcpp::List& var_result,
                      const arma::mat& Z,
                      const arma::ivec& adjustu,
                      const arma::ivec& indsR1,
                      int hor,
                      int nvar,
                      double scale = 10.0) {

    VARResult vr  = fVARResultFromList_cpp(var_result);
    arma::mat eta = vr.residuals.rows(adjustu(0) - 1, adjustu(1) - 1);

    HeteroIRFResult res = fHeteroIRF_cpp(vr.beta, vr.c, eta, Z, indsR1, vr.p,
                                         hor, nvar, scale);

    return Rcpp::List::create(
        Rcpp::Named("IRF") = res.IRF,
        Rcpp::Named("b1")  = res.b1
    );
}

//' MBB Bootstrap for Heteroskedasticity-Based VAR Identification
//'
//' Moving-block bootstrap confidence bands for structural IRFs identified via
//' proxy IV restricted to OPEC announcement months (Kaenzig 2021, treatment
//' group only). Residuals and proxy are resampled jointly.
//'
//' @param y          T x N data matrix.
//' @param var_result List returned by \code{fVAR()}.
//' @param Z          Proxy matrix (T1 x 1), aligned to \code{adjustu}.
//' @param indsR1     Integer 0/1 vector of length T1 (1 = treatment month).
//' @param adjustu    Integer vector \code{c(start, end)} (1-based).
//' @param nboot      Bootstrap replications.
//' @param blocksize  MBB block size (0 = automatic \code{5.03*T^0.25}).
//' @param hor        IRF horizon.
//' @param nvar       1-based normalization variable.
//' @param scale      Shock size.
//' @param conf        Primary confidence level (e.g. 90).
//' @param conf2       Secondary confidence level (e.g. 68).
//' @param n_threads  OpenMP threads (0 = all cores minus one).
//' @param verbose    Logical. If \code{TRUE}, report the number of threads
//'   used. Default \code{FALSE}, so the function is silent.
//'
//' @return A list with \code{upper}, \code{lower}, \code{upper2},
//'   \code{lower2}, \code{point} (each N x hor+1) and \code{conf}, the band
//'   levels \code{c(conf, conf2)}. \code{point} is the plug-in estimate on
//'   the original sample; the bootstrap supplies the bands only.
//'
//' @examples
//' set.seed(123)
//' y          <- matrix(rnorm(400), ncol = 2)
//' var_result <- fVAR(y, p = 2, c = 1)
//'
//' u      <- var_result$residuals
//' Z      <- matrix(u[, 1] + rnorm(nrow(u), sd = 0.5), ncol = 1)
//' rows   <- c(1, nrow(u))
//' indsR1 <- as.integer(abs(Z[, 1]) > stats::median(abs(Z[, 1])))
//'
//' # Few draws here to keep the example quick; use 1000 or more in practice.
//' boot <- fBootstrapHetero(y, var_result, Z, indsR1, rows, nboot = 50,
//'                          blocksize = 10, hor = 12, nvar = 1, scale = 1,
//'                          n_threads = 2)
//' dim(boot$upper)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fBootstrapHetero(const arma::mat& y,
                             const Rcpp::List& var_result,
                             const arma::mat& Z,
                             const arma::ivec& indsR1,
                             const arma::ivec& adjustu,
                             int nboot      = 1000,
                             int blocksize  = 0,
                             int hor        = 48,
                             int nvar       = 1,
                             double scale   = 10.0,
                             double conf     = 90.0,
                             double conf2    = 68.0,
                             int n_threads  = 0,
                             bool verbose   = false) {

    VARResult vr = fVARResultFromList_cpp(var_result);
    fRequireNoExog_cpp(vr.n_exog, "fBootstrapHetero");
    fCheckDrawArgs_cpp(nboot, hor, conf, conf2);
    if (nvar < 1 || nvar > static_cast<int>(y.n_cols))
        Rcpp::stop("'nvar' must be between 1 and %d.",
                   static_cast<int>(y.n_cols));

    BootHeteroResult res = fBootstrapHetero_cpp(y, vr, Z, indsR1, adjustu,
                                                 nboot, blocksize, hor,
                                                 nvar, scale, conf, conf2,
                                                 n_threads, verbose);

    return Rcpp::List::create(
        Rcpp::Named("upper")  = res.upper,
        Rcpp::Named("lower")  = res.lower,
        Rcpp::Named("upper2") = res.upper2,
        Rcpp::Named("lower2") = res.lower2,
        Rcpp::Named("point")  = res.point,
        Rcpp::Named("conf")   = Rcpp::NumericVector::create(conf, conf2)
    );
}
