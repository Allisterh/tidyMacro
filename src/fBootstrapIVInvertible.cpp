#include "fBootstrapIVInvertible.h"
#include "fGetBands.h"
#include "fOLS.h"
#include "fVAR.h"
#include "fGenerateVARData.h"
#include "fLagMakerMatrix.h"
#include "fWoldIRF.h"
#include "fBootstrapVAR.h"
#include "omp_threads.h"
#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

arma::cube fBootstrapIVInvertible_cpp(
    const arma::mat& y,
    const arma::vec& instr,
    const VARResult& var_result,
    int nboot, int p, int c, int hor,
    const arma::ivec& cumu, int n_threads,
    std::vector<unsigned char>& ok)
{
    const arma::mat& beta = var_result.beta;
    const arma::mat& eps  = var_result.residuals;

    const int N = static_cast<int>(y.n_cols);
    const int T = static_cast<int>(y.n_rows);
    const int m = static_cast<int>(eps.n_rows);
    const int H = hor + 1;

    arma::cube bootirf(N, H, nboot, arma::fill::zeros);
    ok.assign(nboot, 0);   // per-draw success, recorded without the R API

    // Pre-generate all Rademacher signs sequentially (Armadillo RNG is not thread-safe)
    arma::mat all_rad(m, nboot);
    for (int b = 0; b < nboot; ++b) {
        arma::vec u = arma::randu<arma::vec>(m);
        for (int i = 0; i < m; ++i) all_rad(i, b) = (u(i) > 0.5) ? 1.0 : -1.0;
    }

    const int n_thr = fResolveThreads_cpp(n_threads);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_thr)
#endif
    for (int b = 0; b < nboot; ++b) {

        // 1. Wild bootstrap: apply pre-generated Rademacher signs to VAR residuals
        const arma::vec rad = all_rad.col(b);

        arma::mat eps_wild(m, N, arma::fill::none);
        for (int i = 0; i < m; ++i) eps_wild.row(i) = eps.row(i) * rad(i);
        arma::mat y_new = fGenerateVARData(y, p, c, beta, eps_wild);

        // 2. Apply same Rademacher signs to instrument (indices p..T-1)
        arma::vec z_new = instr;
        for (int i = p; i < T; ++i) z_new(i) = instr(i) * rad(i - p);

        // 3. Clean instrument
        arma::vec z_b = fCleanInstrument_cpp(z_new, p, y_new);

        // 4. Estimate VAR and Wold IRF. A resample can be rank deficient or
        //    produce a singular covariance; drop the draw instead of throwing,
        //    which would abort the R session from inside the parallel region.
        VARResult vr_b = fVAR_cpp(y_new, p, c);
        if (!vr_b.ok) continue;
        arma::mat Sigma_inv;
        if (!arma::inv_sympd(Sigma_inv, vr_b.sigma)) continue;
        WoldIRFResult wr_b = fWoldIRF_cpp(vr_b, hor);
        const arma::mat& eps_b   = vr_b.residuals;
        const arma::cube& wold_b = wr_b.irfwold;

        // 5. Re-estimate structural impact vector via IV (lightweight OLS: beta only)
        arma::mat z_b_mat(z_b.n_elem, 1);
        z_b_mat.col(0) = z_b;
        arma::mat s0 = fOLS_cpp(arma::mat(eps_b), z_b_mat, 0, 0, 0, 1).beta; // 1 x N
        const double s0_1 = s0(0, 0);
        if (!std::isfinite(s0_1) || std::fabs(s0_1) < 1e-12) continue;
        arma::vec s_b = s0.t() / s0_1;
        const double scaler_b = std::sqrt(arma::dot(s_b, Sigma_inv * s_b));
        if (!std::isfinite(scaler_b) || scaler_b < 1e-12) continue;

        // 6. Compute IRF, cumulate, sign-normalise
        arma::mat irf_b(N, H, arma::fill::none);
        for (int h = 0; h < H; ++h)
            irf_b.col(h) = wold_b.slice(h) * s_b / scaler_b;

        for (int i = 0; i < static_cast<int>(cumu.n_elem); ++i) {
            int v = cumu(i) - 1;
            for (int h = 1; h < H; ++h) irf_b(v, h) += irf_b(v, h - 1);
        }

        double sgn = (irf_b(0, 0) >= 0.0) ? 1.0 : -1.0;
        bootirf.slice(b) = irf_b * sgn;
        ok[b] = 1;
    }

    return bootirf;
}

//' Bootstrap IRFs for the Invertible Case (Wild Bootstrap)
//'
//' Implements the wild bootstrap with Rademacher signs for standard SVAR-IV.
//' The same Rademacher vector is applied to both the VAR residuals and the
//' raw instrument; at each draw the structural impact vector is re-estimated
//' by IV regression.
//'
//' @param y T x N matrix of endogenous variables.
//' @param instr Length-T raw instrument vector (before cleaning).
//' @param var_result List from \code{fVAR}.
//' @param nboot Number of bootstrap replications.
//' @param p VAR lag order.
//' @param c Integer (0/1): include intercept.
//' @param hor Impulse-response horizon.
//' @param cumu 1-indexed integer vector of variable positions to cumulate.
//' @param conf Primary confidence level (e.g. 90). Default 90.
//' @param conf2 Secondary confidence level (e.g. 68). Default 68.
//' @param n_threads OpenMP threads. 0 (default) uses all cores minus one.
//'
//' @return List with N x (hor+1) matrices: upper, lower (at \code{conf}),
//'   upper2, lower2 (at \code{conf2}), median.
//'
//' @examples
//' set.seed(123)
//' y          <- matrix(rnorm(400), ncol = 2)
//' var_result <- fVAR(y, p = 1, c = 1)
//'
//' # A length-T instrument (the first observation is lost to the lag)
//' u     <- var_result$residuals
//' instr <- c(0, u[, 1] + rnorm(nrow(u), sd = 0.5))
//'
//' # Few draws here to keep the example quick; use 1000 or more in practice.
//' boot <- fBootstrapIVInvertible(y, instr, var_result, nboot = 50, p = 1,
//'                                c = 1, hor = 12, cumu = integer(),
//'                                n_threads = 2)
//' dim(boot$median)
//'
//' @seealso \code{\link{fBootstrapIVRecover}}, \code{\link{fGetBands}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fBootstrapIVInvertible(
    const arma::mat& y,
    const arma::vec& instr,
    const Rcpp::List& var_result,
    int nboot, int p, int c, int hor,
    const arma::ivec& cumu,
    double conf = 90.0, double conf2 = 68.0,
    int n_threads = 0)
{
    VARResult vr = fVARResultFromList_cpp(var_result);
    fRequireNoExog_cpp(vr.n_exog, "fBootstrapIVInvertible");
    fCheckDrawArgs_cpp(nboot, hor, conf, conf2);
    if (static_cast<int>(instr.n_elem) != static_cast<int>(y.n_rows))
        Rcpp::stop("'instr' must have one value per row of 'y' (%d vs %d).",
                   static_cast<int>(instr.n_elem), static_cast<int>(y.n_rows));

    std::vector<unsigned char> ok;
    arma::cube bootirf = fBootstrapIVInvertible_cpp(
        y, instr, vr, nboot, p, c, hor, cumu, n_threads, ok);
    bootirf = fKeepValidSlices_cpp(bootirf, ok, "fBootstrapIVInvertible");

    // Sort once, extract both confidence bands in a single pass
    FGetBands2Result bands = fGetBands2_cpp(bootirf, conf, conf2);

    return Rcpp::List::create(
        Rcpp::Named("upper")   = bands.upper,
        Rcpp::Named("lower")   = bands.lower,
        Rcpp::Named("upper2")  = bands.upper2,
        Rcpp::Named("lower2")  = bands.lower2,
        Rcpp::Named("median")  = bands.median
    );
}
