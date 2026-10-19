// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include <RcppArmadillo.h>
#include "fBootstrapVAR.h"
#include "fVAR.h"
#include "omp_threads.h"
#include "fGetBands.h"
#include "fRemoveBias.h"
#include "fLagMakerMatrix.h"
#include "fOLS.h"
#include <cstdio>
#include <string>

BootstrapPlan fBootstrapPlan_cpp(const std::string& bootscheme,
                                 int T_iter, int T_resid, int nboot) {
    BootstrapPlan plan;
    if (bootscheme == "residual") {
        // iid resampling of residual rows
        plan.wild = false;
        plan.idx.set_size(T_iter, nboot);
        for (int b = 0; b < nboot; ++b)
            plan.idx.col(b) = arma::randi<arma::uvec>(
                T_iter, arma::distr_param(0, T_resid - 1));
    } else if (bootscheme == "wild") {
        // Rademacher signs: +1 or -1 with equal probability
        plan.wild = true;
        plan.signs.set_size(T_resid, nboot);
        for (int b = 0; b < nboot; ++b) {
            arma::vec u = arma::randu<arma::vec>(T_resid);
            plan.signs.col(b) = arma::sign(u - 0.5);
        }
    } else {
        Rcpp::stop("'bootscheme' must be \"residual\" or \"wild\".");
    }
    return plan;
}

arma::mat fSimulateVAR_cpp(const arma::mat& y,
                           const VARResult& var_result,
                           const arma::mat& shocks,
                           const arma::mat* exog_ptr) {
    const arma::mat& beta = var_result.beta;
    const int p = var_result.p;
    const int c = var_result.c;
    const int T = static_cast<int>(y.n_rows);
    const int N = static_cast<int>(y.n_cols);
    const int T_iter = T - p;

    // beta layout (rows): [intercept (if c=1)] | [A_1 ; ... ; A_p] | [exog]
    const int lag_end_row = c + N * p - 1;
    const arma::rowvec intercept = (c == 1) ? arma::rowvec(beta.row(0))
                                            : arma::zeros<arma::rowvec>(N);
    const arma::mat Pi = beta.rows(c, lag_end_row);
    const bool add_exog = var_result.n_exog > 0 && exog_ptr != nullptr;
    arma::mat B_exog;
    if (add_exog) B_exog = beta.rows(lag_end_row + 1, beta.n_rows - 1);

    arma::mat ynext(T, N, arma::fill::none);
    ynext.rows(0, p - 1) = y.rows(0, p - 1);   // seed with original data

    // lag vector: [y_{t-1}' , y_{t-2}' , ... , y_{t-p}']  (column, length N*p)
    arma::vec ylag(N * p, arma::fill::none);
    for (int lag = 0; lag < p; ++lag)
        ylag.rows(lag * N, lag * N + N - 1) = y.row(p - 1 - lag).t();

    for (int i = 0; i < T_iter; ++i) {
        arma::rowvec yt = intercept + (Pi.t() * ylag).t() + shocks.row(i);
        if (add_exog)
            yt += (B_exog.t() * exog_ptr->row(p + i).t()).t();
        ynext.row(p + i) = yt;

        // Shift the lag vector down and insert the new observation on top
        if (p > 1) ylag.rows(N, N * p - 1) = ylag.rows(0, N * (p - 1) - 1);
        ylag.rows(0, N - 1) = yt.t();
    }
    return ynext;
}

arma::mat fBootstrapVAR_cpp(const arma::mat& y,
                            const VARResult& var_result,
                            const BootstrapPlan& plan, int b,
                            const arma::mat* exog_ptr) {
    const arma::mat& residuals = var_result.residuals;
    const int T_iter  = static_cast<int>(y.n_rows) - var_result.p;
    const int T_resid = static_cast<int>(residuals.n_rows);

    arma::mat shocks(T_iter, residuals.n_cols, arma::fill::none);
    if (plan.wild) {
        for (int i = 0; i < T_iter; ++i) {
            const int ri = i % T_resid;
            shocks.row(i) = residuals.row(ri) * plan.signs(ri, b);
        }
    } else {
        for (int i = 0; i < T_iter; ++i)
            shocks.row(i) = residuals.row(plan.idx(i, b));
    }
    return fSimulateVAR_cpp(y, var_result, shocks, exog_ptr);
}

arma::mat fBootstrapMeanBeta_cpp(const arma::mat& y,
                                 const VARResult& var_result,
                                 int nboot, const std::string& bootscheme,
                                 const arma::mat* exog_ptr, int n_threads) {
    const int p      = var_result.p;
    const int c      = var_result.c;
    const int N      = static_cast<int>(y.n_cols);
    const int n_coef = static_cast<int>(var_result.beta.n_rows);
    const int T_iter = static_cast<int>(y.n_rows) - p;
    const int T_resid = static_cast<int>(var_result.residuals.n_rows);

    const BootstrapPlan plan = fBootstrapPlan_cpp(bootscheme, T_iter, T_resid, nboot);
    const int n_thr = fResolveThreads_cpp(n_threads);

    // Store each draw and average in replication order afterwards; a
    // thread-local reduction would make the sum depend on the thread count.
    arma::cube draws(N, n_coef, nboot, arma::fill::none);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_thr)
#endif
    for (int b = 0; b < nboot; ++b) {
        arma::mat yb = fBootstrapVAR_cpp(y, var_result, plan, b, exog_ptr);
        draws.slice(b) = fVAR_cpp(yb, p, c, exog_ptr).beta.t();
    }

    arma::mat boot_mean(N, n_coef, arma::fill::zeros);
    for (int b = 0; b < nboot; ++b) boot_mean += draws.slice(b);
    return boot_mean / static_cast<double>(nboot);
}

VARResult fBiasCorrectedVAR_cpp(const arma::mat& y,
                                const VARResult& var_result,
                                int nboot1, const std::string& bootscheme,
                                const arma::mat* exog_ptr, int n_threads,
                                bool verbose) {
    if (nboot1 <= 0) Rcpp::stop("'nboot1' must be positive.");
    if (verbose)
        Rprintf("[Pass 1] Bias estimation: %d reps, %d thread(s)\n",
                nboot1, fResolveThreads_cpp(n_threads));
    arma::mat boot_mean = fBootstrapMeanBeta_cpp(y, var_result, nboot1,
                                                 bootscheme, exog_ptr, n_threads);
    RemoveBiasResult rb = fRemoveBias_cpp(var_result.beta, var_result.c,
                                          var_result.p, boot_mean);
    if (verbose) {
        if (rb.corrections > 1)
            Rprintf("[Bias correction] %d shrinkage iteration(s)\n", rb.corrections);
        else
            Rprintf("[Bias correction] Full correction applied\n");
    }

    VARResult corrected = var_result;
    corrected.beta = rb.Beta.t();
    return corrected;
}

void fCheckBootArgs_cpp(const arma::mat& y, const VARResult& var_result,
                        int nboot, int horizon, double conf, double conf2) {
    if (nboot <= 0)
        Rcpp::stop("'nboot' must be positive.");
    if (horizon < 0)
        Rcpp::stop("'horizon' must be non-negative.");
    if (conf <= 0.0 || conf > 100.0 || conf2 <= 0.0 || conf2 > 100.0)
        Rcpp::stop("'conf' and 'conf2' must be in (0, 100].");
    if (!y.is_finite())
        Rcpp::stop("'y' contains NA, NaN or Inf values.");

    const int T = static_cast<int>(y.n_rows);
    const int N = static_cast<int>(y.n_cols);
    const int p = var_result.p;
    if (T <= p)
        Rcpp::stop("'y' has %d rows, which is not more than p = %d.", T, p);

    // Degrees of freedom from the actual design width, which includes the
    // intercept and any exogenous regressors.
    const int k  = var_result.c + N * p + var_result.n_exog;
    const int df = (T - p) - k;
    if (df <= 0)
        Rcpp::stop("%d usable observations for %d regressors: no residual "
                   "degrees of freedom left.", T - p, k);
    if (static_cast<int>(var_result.beta.n_rows) != k)
        Rcpp::stop("'var_result' does not match 'y': beta has %d rows but the "
                   "design needs %d (c + N*p + n_exog).",
                   static_cast<int>(var_result.beta.n_rows), k);
}

void fCheckDrawArgs_cpp(int nboot, int horizon, double conf, double conf2) {
    if (nboot <= 0)
        Rcpp::stop("'nboot' must be positive, got %d.", nboot);
    if (horizon < 0)
        Rcpp::stop("'horizon' must be non-negative, got %d.", horizon);
    if (conf <= 0.0 || conf > 100.0 || conf2 <= 0.0 || conf2 > 100.0)
        Rcpp::stop("'conf' and 'conf2' must be in (0, 100].");
}

arma::mat fKeepValidDraws_cpp(const arma::mat& flat,
                              const std::vector<unsigned char>& ok,
                              const char* label) {
    const arma::uword nboot = ok.size();
    arma::uvec keep(nboot);
    arma::uword n_keep = 0;
    for (arma::uword b = 0; b < nboot; ++b)
        if (ok[b] && flat.col(b).is_finite()) keep(n_keep++) = b;

    if (n_keep == 0)
        Rcpp::stop("%s: every bootstrap draw failed. The instrument is "
                   "probably too weak, or the VAR is near-singular on "
                   "resampled data.", label);
    if (n_keep == nboot) return flat;

    char msg[256];
    std::snprintf(msg, sizeof(msg),
                  "%s: %d of %d bootstrap draws failed and were dropped "
                  "before computing the bands.", label,
                  static_cast<int>(nboot - n_keep), static_cast<int>(nboot));
    Rcpp::warning(std::string(msg));
    return flat.cols(keep.head(n_keep));
}

arma::cube fKeepValidSlices_cpp(const arma::cube& draws,
                                const std::vector<unsigned char>& ok,
                                const char* label) {
    const arma::uword nboot = ok.size();
    arma::uvec keep(nboot);
    arma::uword n_keep = 0;
    for (arma::uword b = 0; b < nboot; ++b)
        if (ok[b] && draws.slice(b).is_finite()) keep(n_keep++) = b;

    if (n_keep == 0)
        Rcpp::stop("%s: every bootstrap draw failed. The instrument is "
                   "probably too weak, or the VAR is near-singular on "
                   "resampled data.", label);
    if (n_keep == nboot) return draws;

    char msg[256];
    std::snprintf(msg, sizeof(msg),
                  "%s: %d of %d bootstrap draws failed and were dropped "
                  "before computing the bands.", label,
                  static_cast<int>(nboot - n_keep), static_cast<int>(nboot));
    Rcpp::warning(std::string(msg));

    arma::cube kept(draws.n_rows, draws.n_cols, n_keep, arma::fill::none);
    for (arma::uword j = 0; j < n_keep; ++j) kept.slice(j) = draws.slice(keep(j));
    return kept;
}

arma::vec fCleanInstrument_cpp(const arma::vec& instr_raw, int p,
                               const arma::mat& y) {
    const int T = static_cast<int>(y.n_rows);
    arma::mat instr_mat(instr_raw.n_elem, 1);
    instr_mat.col(0) = instr_raw;
    arma::vec z0 = instr_raw.rows(p, T - 1);
    arma::mat z0_mat(z0.n_elem, 1);
    z0_mat.col(0) = z0;
    arma::mat lag_i   = fLagMakerMatrix(instr_mat, p);
    arma::mat lag_X   = fLagMakerMatrix(y, p);
    arma::mat X_clean = arma::join_horiz(lag_i, lag_X);
    return fOLS_cpp(z0_mat, X_clean, 1, 0, 0, 1).err.col(0);
}

void fSVARBootBands_cpp(SVARBootResult& r, double conf, double conf2,
                        int n_threads, const char* label) {
    const arma::mat kept = fKeepValidDraws_cpp(r.flat, r.ok, label);
    r.bands = fRowPercentiles_cpp(kept, fBandLevels_cpp(conf, conf2), n_threads);
    r.conf  = {conf, conf2};
}

Rcpp::List fSVARBootList_cpp(const SVARBootResult& r, const char* draws_name,
                             const std::vector<int>& dims) {
    Rcpp::IntegerVector band_dim(dims.begin(), dims.end());
    Rcpp::IntegerVector draw_dim(dims.begin(), dims.end());
    draw_dim.push_back(static_cast<int>(r.flat.n_cols));

    Rcpp::NumericVector draws(r.flat.begin(), r.flat.end());
    draws.attr("dim") = draw_dim;

    auto band = [&](int j) {
        Rcpp::NumericVector v(r.bands.colptr(j), r.bands.colptr(j) + r.bands.n_rows);
        v.attr("dim") = band_dim;
        return v;
    };

    return Rcpp::List::create(
        Rcpp::Named(draws_name)  = draws,
        Rcpp::Named("upper")     = band(0),
        Rcpp::Named("lower")     = band(1),
        Rcpp::Named("upper2")    = band(2),
        Rcpp::Named("lower2")    = band(3),
        Rcpp::Named("boot_beta") = r.boot_beta,
        Rcpp::Named("conf")      = Rcpp::NumericVector(r.conf.begin(), r.conf.end()));
}

//' Bootstrap VAR Model
//'
//' Generates a bootstrap pseudo-sample from a fitted VAR model using either
//' residual (iid) resampling or a wild (Rademacher) bootstrap scheme.
//'
//' @param y T x N matrix of original endogenous variables.
//' @param fVAR_result List from \code{fVAR} containing beta, residuals, p,
//'   and c.
//' @param bootscheme String indicating the bootstrap method: "residual" for
//'   iid resampling of residuals, or "wild" for a Rademacher wild bootstrap.
//'
//' @return A list with two elements: ynext — T x N bootstrapped data matrix
//'   (first p rows copied from \code{y}); rademacher — T_eff-length vector of
//'   Rademacher signs used in the wild scheme (all zeros for residual scheme).
//'
//' @details
//' The recursion is \deqn{y_t = c + A_1 y_{t-1} + \cdots + A_p y_{t-p} + e_t^*}
//' where \eqn{e_t^*} is either a randomly resampled residual (residual scheme)
//' or the t-th residual multiplied by a Rademacher sign (wild scheme).
//'
//' @examples
//' set.seed(123)
//' y          <- matrix(rnorm(400), ncol = 2)
//' var_result <- fVAR(y, p = 2, c = 1)
//'
//' # One artificial sample drawn from the estimated VAR
//' y_star <- fBootstrapVAR(y, var_result, bootscheme = "wild")
//' dim(y_star$ynext)
//'
//' @seealso \code{\link{fBootstrapChol}}, \code{\link{fVAR}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fBootstrapVAR(const arma::mat&    y,
                         const Rcpp::List&   fVAR_result,
                         const std::string&  bootscheme = "residual") {

    VARResult var_result = fVARResultFromList_cpp(fVAR_result);
    const int T_iter = static_cast<int>(y.n_rows) - var_result.p;
    BootstrapPlan plan = fBootstrapPlan_cpp(bootscheme, T_iter,
                                            var_result.residuals.n_rows, 1);

    arma::mat ynext = fBootstrapVAR_cpp(y, var_result, plan, 0, nullptr);
    arma::vec rademacher = plan.wild ? arma::vec(plan.signs.col(0))
                                     : arma::vec(T_iter, arma::fill::zeros);

    return Rcpp::List::create(
        Rcpp::Named("ynext")      = ynext,
        Rcpp::Named("rademacher") = rademacher
    );
}
