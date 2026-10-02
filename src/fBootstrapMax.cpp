// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include "fBootstrapMax.h"
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

// Closed-form LR_max solution (replaces MATLAB fminsearch on LR_max.m).
// Maximises wold.slice(H-1).row(var_idx) * S * h2  (last-horizon response only,
// matching MATLAB irfs(var,:,end)) subject to h2(0) = 0, ||h2|| = 1.
// Optimal: h2 = [0; M_sub / ||M_sub||] where
// M_sub = (wold.slice(H-1).row(var_idx) * S).cols(1, N-1).
static arma::vec lr_max_solve(const arma::cube& wold, const arma::mat& S,
                               int var_idx) {
    const int H = static_cast<int>(wold.n_slices);
    const int N = static_cast<int>(wold.n_rows);

    // Maximise the last-horizon response of var_idx (MATLAB LR_max.m uses irfs(var,:,end))
    arma::rowvec M_full = wold.slice(H - 1).row(var_idx) * S;

    arma::rowvec M_sub = M_full.cols(1, N - 1);
    double M_norm = arma::norm(M_sub);

    arma::vec h2(N, arma::fill::zeros);
    if (M_norm > 1e-14) {
        h2.rows(1, N - 1) = M_sub.t() / M_norm;
    }
    return h2;
}


BootstrapMaxResult
fBootstrapMax_cpp(const arma::mat& y, const VARResult& var_result,
                  int nboot, int horizon, int var_idx, double conf, double conf2,
                  const arma::uvec& cumulate,
                  Rcpp::Nullable<arma::vec> scaling,
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
        Rcpp::stop("LR-Max identification requires at least two variables.");
    if (var_idx < 0 || var_idx >= N)
        Rcpp::stop("'var_idx' is out of range.");
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

    arma::mat bootmax_flat(slice_sz, nboot, arma::fill::zeros);
    arma::cube boot_beta(N, n_coef, nboot, arma::fill::zeros);

    const bool has_scaling = scaling.isNotNull();
    arma::vec scaling_vec;
    if (has_scaling) {
        scaling_vec = Rcpp::as<arma::vec>(scaling);
        if (scaling_vec.n_elem != 2)
            Rcpp::stop("'scaling' must contain a variable index and shock size.");
        const int s_idx = static_cast<int>(scaling_vec(0)) - 1;
        if (s_idx < 0 || s_idx >= N)
            Rcpp::stop("The variable index in 'scaling' is out of range.");
        if (!std::isfinite(scaling_vec(1)) || scaling_vec(1) == 0.0)
            Rcpp::stop("The shock size in 'scaling' must be finite and nonzero.");
    }

    const bool has_exog = exog.isNotNull();
    arma::mat exog_mat;
    if (has_exog) exog_mat = Rcpp::as<arma::mat>(exog);
    const arma::mat* exog_ptr = has_exog ? &exog_mat : nullptr;

    // RcppArmadillo random draws use R's main-thread-only RNG. Generate the
    // complete resampling plan before the OpenMP loop so workers perform no R
    // API calls and results are invariant to thread count for a fixed seed.
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
    Rprintf("Using %d thread(s) for LR-Max bootstrap...\n", actual_threads);
#else
    Rprintf("OpenMP not available. Running single-threaded LR-Max bootstrap.\n");
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

        arma::vec h2 = lr_max_solve(wold_loop, S_loop, var_idx);

        // Precompute impact = S_loop * h2 once per draw — reused for sign
        // check and every horizon below.
        arma::vec impact = S_loop * h2;

        // Orient the shock so the selected variable's maximised response is
        // non-negative. For the news-shock application var_idx is TFP (0).
        if (arma::dot(wold_loop.slice(H - 1).row(var_idx), impact) < 0.0)
            impact = -impact;

        arma::mat struct_irf(N, H, arma::fill::none);
        for (int hh = 0; hh < H; ++hh) {
            struct_irf.col(hh) = wold_loop.slice(hh) * impact;
        }

        if (has_scaling) {
            int s_idx = static_cast<int>(scaling_vec(0)) - 1;
            double scale_val = struct_irf(s_idx, 0) * scaling_vec(1);
            if (std::abs(scale_val) > 1e-14) struct_irf /= scale_val;
        }

        for (arma::uword ci = 0; ci < cumulate.n_elem; ++ci) {
            arma::uword ri = cumulate(ci);
            struct_irf.row(ri) = arma::cumsum(struct_irf.row(ri));
        }

        bootmax_flat.col(b) = arma::vectorise(struct_irf);
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
        // MATLAB prctile's default exact/midpoint rule assigns sorted
        // observation i the probability (i - 0.5) / n.
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
        arma::rowvec row_rv = bootmax_flat.row(i);
        std::vector<double> sv(row_rv.begin(), row_rv.end());
        int row_i = i % N;
        int col_i = i / N;
        upper (row_i, col_i) = nth_pct(sv, up_pct);
        lower (row_i, col_i) = nth_pct(sv, low_pct);
        upper2(row_i, col_i) = nth_pct(sv, up_pct2);
        lower2(row_i, col_i) = nth_pct(sv, low_pct2);
    }

    BootstrapMaxResult result;
    result.bootmax_flat = bootmax_flat;
    result.upper        = upper;
    result.lower        = lower;
    result.upper2       = upper2;
    result.lower2       = lower2;
    result.boot_beta    = boot_beta;
    result.N            = N;
    result.H            = H;
    return result;
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
//'   coefficient array, where K is the number of regressors. The corrected
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
    BootstrapMaxResult res = fBootstrapMax_cpp(y, vr, nboot, horizon,
                                               var_idx - 1, conf, conf2,
                                               cumulate_cpp, scaling,
                                               exog, n_threads);

    Rcpp::NumericVector bootmax_out(res.bootmax_flat.begin(), res.bootmax_flat.end());
    bootmax_out.attr("dim") = Rcpp::IntegerVector::create(res.N, res.H, nboot);

    return Rcpp::List::create(
        Rcpp::Named("bootmax")   = bootmax_out,
        Rcpp::Named("upper")     = res.upper,
        Rcpp::Named("lower")     = res.lower,
        Rcpp::Named("upper2")    = res.upper2,
        Rcpp::Named("lower2")    = res.lower2,
        Rcpp::Named("boot_beta") = res.boot_beta
    );
}
