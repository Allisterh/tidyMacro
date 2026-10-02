// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include "fBootstrapUhlig.h"
#include "fUhligMaxShare.h"
#include "fVAR.h"
#include "fBootstrapVAR.h"
#include "fWoldIRF.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <R_ext/Print.h>
#include <algorithm>
#include <cmath>


BootstrapUhligResult
fBootstrapUhlig_cpp(const arma::mat& y, const VARResult& var_result,
                    int nboot, int horizon, int idx, double conf, double conf2,
                    const arma::uvec& cumulate,
                    Rcpp::Nullable<arma::mat> exog,
                    int n_threads) {

    const int p        = var_result.p;
    const int c        = var_result.c;
    const int n_exog   = var_result.n_exog;
    const int T        = static_cast<int>(y.n_rows);
    const int N        = static_cast<int>(y.n_cols);
    const int H        = horizon + 1;
    const int n_coef   = static_cast<int>(var_result.beta.n_rows);
    const int slice_sz = N * H;

    if (nboot <= 0)
        Rcpp::stop("'nboot' must be positive.");
    if (horizon < 0)
        Rcpp::stop("'horizon' must be non-negative.");
    if (N < 2)
        Rcpp::stop("Uhlig max-share identification requires at least two variables.");
    if (idx < 0 || idx >= N)
        Rcpp::stop("'idx' is out of range.");
    if (conf <= 0.0 || conf > 100.0 || conf2 <= 0.0 || conf2 > 100.0)
        Rcpp::stop("'conf' and 'conf2' must be in (0, 100].");
    if (!cumulate.is_empty() && cumulate.max() >= static_cast<arma::uword>(N))
        Rcpp::stop("'cumulate' contains an out-of-range variable index.");

    if (n_exog > 0 && exog.isNull())
        Rcpp::stop("Original VAR used exogenous variables. You must provide the 'exog' parameter.");
    if (n_exog == 0 && exog.isNotNull())
        Rcpp::stop("Original VAR did not use exogenous variables. Do not provide the 'exog' parameter.");

    const double df = static_cast<double>(T - 1 - p - N * p);
    if (df <= 0.0)
        Rcpp::stop("Degrees of freedom T-1-p-N*p = %.0f <= 0.", df);

    arma::mat bootuhlig_flat(slice_sz, nboot, arma::fill::zeros);
    arma::cube boot_beta(N, n_coef, nboot, arma::fill::zeros);

    const bool has_exog = exog.isNotNull();
    arma::mat exog_mat;
    if (has_exog) exog_mat = Rcpp::as<arma::mat>(exog);
    const arma::mat* exog_ptr = has_exog ? &exog_mat : nullptr;

    // Pre-generate residual indices on R's main thread. The OpenMP loop then
    // contains no RNG/R API calls and is reproducible across thread counts.
    const int T_iter = T - p;
    const int T_resid = static_cast<int>(var_result.residuals.n_rows);
    arma::umat resample_indices(T_iter, nboot, arma::fill::none);
    for (int b = 0; b < nboot; ++b)
        resample_indices.col(b) = arma::randi<arma::uvec>(
            T_iter, arma::distr_param(0, T_resid - 1));

    int actual_threads = 1;
#ifdef _OPENMP
    actual_threads = (n_threads <= 0)
                         ? std::max(1, omp_get_max_threads() - 1)
                         : n_threads;
    omp_set_num_threads(actual_threads);
    Rprintf("Using %d thread(s) for Uhlig max-share bootstrap...\n", actual_threads);
#else
    Rprintf("OpenMP not available. Running single-threaded Uhlig bootstrap.\n");
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(actual_threads)
#endif
    for (int b = 0; b < nboot; ++b) {

        arma::uvec indices = resample_indices.col(b);
        BootstrapVARResult boot_data =
            fBootstrapVAR_cpp(y, var_result, indices, exog_ptr);

        VARResult var_loop = has_exog
                                 ? fVAR_cpp(boot_data.ynext, p, c, exog_mat)
                                 : fVAR_cpp(boot_data.ynext, p, c);
        boot_beta.slice(b) = var_loop.beta.t();

        WoldIRFResult wold_res = fWoldIRF_cpp(var_loop, horizon);
        const arma::cube& wold_loop = wold_res.irfwold;

        arma::mat S_loop;
        if (!arma::chol(S_loop, var_loop.sigma, "lower")) {
            arma::mat reg = var_loop.sigma;
            reg.diag() += 1e-8 * arma::trace(reg) / N;
            arma::chol(S_loop, reg, "lower");
        }

        arma::vec h2 = fUhligMaxShare_cpp(wold_loop, S_loop, idx);

        // Precompute impact = S_loop * h2 once per draw
        arma::vec impact = S_loop * h2;

        // Sign: ensure last-horizon non-structural response of var 0 is non-negative
        if (arma::dot(wold_loop.slice(H - 1).row(0), impact) < 0.0)
            impact = -impact;

        arma::mat struct_irf(N, H, arma::fill::none);
        for (int hh = 0; hh < H; ++hh) {
            struct_irf.col(hh) = wold_loop.slice(hh) * impact;
        }

        for (arma::uword ci = 0; ci < cumulate.n_elem; ++ci) {
            arma::uword ri = cumulate(ci);
            struct_irf.row(ri) = arma::cumsum(struct_irf.row(ri));
        }

        bootuhlig_flat.col(b) = arma::vectorise(struct_irf);
    }

    const double up_pct   = 50.0 + conf  * 0.5;
    const double low_pct  = 50.0 - conf  * 0.5;
    const double up_pct2  = 50.0 + conf2 * 0.5;
    const double low_pct2 = 50.0 - conf2 * 0.5;

    arma::mat upper (N, H, arma::fill::zeros);
    arma::mat lower (N, H, arma::fill::zeros);
    arma::mat upper2(N, H, arma::fill::zeros);
    arma::mat lower2(N, H, arma::fill::zeros);

    auto nth_pct = [](std::vector<double>& v, double pct) -> double {
        const int n = static_cast<int>(v.size());
        // Match MATLAB prctile's default exact/midpoint interpolation.
        const double raw = (pct / 100.0) * n - 0.5;
        if (raw <= 0.0) return *std::min_element(v.begin(), v.end());
        if (raw >= n - 1.0) return *std::max_element(v.begin(), v.end());
        const int lo = static_cast<int>(std::floor(raw));
        const double frac = raw - lo;
        std::nth_element(v.begin(), v.begin() + lo, v.end());
        const double lo_val = v[lo];
        if (frac < 1e-12 || lo + 1 >= n) return lo_val;
        return lo_val * (1.0 - frac) +
               *std::min_element(v.begin() + lo + 1, v.end()) * frac;
    };

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < slice_sz; ++i) {
        arma::rowvec row_rv = bootuhlig_flat.row(i);
        std::vector<double> sv(row_rv.begin(), row_rv.end());
        int row_i = i % N;
        int col_i = i / N;
        upper (row_i, col_i) = nth_pct(sv, up_pct);
        lower (row_i, col_i) = nth_pct(sv, low_pct);
        upper2(row_i, col_i) = nth_pct(sv, up_pct2);
        lower2(row_i, col_i) = nth_pct(sv, low_pct2);
    }

    BootstrapUhligResult result;
    result.bootuhlig_flat = bootuhlig_flat;
    result.upper          = upper;
    result.lower          = lower;
    result.upper2         = upper2;
    result.lower2         = lower2;
    result.boot_beta      = boot_beta;
    result.N              = N;
    result.H              = H;
    return result;
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
//'   coefficient array, where K is the number of regressors. The corrected
//'   variant stores \code{nboot2} draws.
//' @export
// [[Rcpp::export]]
Rcpp::List fBootstrapUhlig(const arma::mat& y, const Rcpp::List& var_result,
                            int nboot, int horizon, int idx,
                            double conf = 90.0, double conf2 = 68.0,
                            Rcpp::IntegerVector cumulate = Rcpp::IntegerVector::create(),
                            Rcpp::Nullable<arma::mat> exog = R_NilValue,
                            int n_threads = 0) {
    VARResult vr;
    vr.beta       = Rcpp::as<arma::mat>(var_result["beta"]);
    vr.residuals  = Rcpp::as<arma::mat>(var_result["residuals"]);
    vr.sigma = Rcpp::as<arma::mat>(var_result["sigma"]);
    vr.p          = Rcpp::as<int>(var_result["p"]);
    vr.c          = Rcpp::as<int>(var_result["c"]);
    vr.n_exog     = var_result.containsElementNamed("n_exog")
                        ? Rcpp::as<int>(var_result["n_exog"]) : 0;

    arma::uvec cumulate_cpp = (cumulate.size() == 0)
        ? arma::uvec()
        : Rcpp::as<arma::uvec>(cumulate) - 1;
    BootstrapUhligResult res = fBootstrapUhlig_cpp(y, vr, nboot, horizon,
                                                   idx - 1, conf, conf2,
                                                   cumulate_cpp, exog, n_threads);

    Rcpp::NumericVector bootuhlig_out(res.bootuhlig_flat.begin(), res.bootuhlig_flat.end());
    bootuhlig_out.attr("dim") = Rcpp::IntegerVector::create(res.N, res.H, nboot);

    return Rcpp::List::create(
        Rcpp::Named("bootuhlig") = bootuhlig_out,
        Rcpp::Named("upper")     = res.upper,
        Rcpp::Named("lower")     = res.lower,
        Rcpp::Named("upper2")    = res.upper2,
        Rcpp::Named("lower2")    = res.lower2,
        Rcpp::Named("boot_beta") = res.boot_beta
    );
}
