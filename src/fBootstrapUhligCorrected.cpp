// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include "fBootstrapUhligCorrected.h"
#include "fBootstrapUhlig.h"
#include "fVAR.h"
#include "fBootstrapVAR.h"
#include "fCompanionMatrix.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <R_ext/Print.h>
#include <algorithm>

static void bias_correct_uhlig(const arma::mat& beta, int c, int p,
                                const arma::mat& boot_mean,
                                arma::mat& Beta_t, int& corrections) {
    const arma::mat beta_t = beta.t();
    arma::mat bias = boot_mean - beta_t;

    CompanionMatrixResult comp0 = fCompanionMatrix_cpp(beta, c, p);
    arma::cx_vec ev0 = arma::eig_gen(comp0.comp);
    double max_ev = arma::max(arma::abs(ev0));

    corrections = 1;

    if (max_ev >= 1.0) {
        Beta_t = beta_t;
        return;
    }

    Beta_t = beta_t - bias;

    {
        CompanionMatrixResult comp1 = fCompanionMatrix_cpp(Beta_t.t(), c, p);
        arma::cx_vec ev1 = arma::eig_gen(comp1.comp);
        max_ev = arma::max(arma::abs(ev1));
    }

    double delta = 1.0;
    while (max_ev >= 1.0) {
        delta -= 0.01;
        corrections += 1;
        if (delta < 0.0 || corrections > 200) {
            Beta_t = beta_t;
            break;
        }
        Beta_t = beta_t - bias * delta;
        CompanionMatrixResult comp_loop = fCompanionMatrix_cpp(Beta_t.t(), c, p);
        arma::cx_vec ev_loop = arma::eig_gen(comp_loop.comp);
        max_ev = arma::max(arma::abs(ev_loop));
    }
}

BootstrapUhligResult
fBootstrapUhligCorrected_cpp(const arma::mat& y, const VARResult& var_result,
                              int nboot1, int nboot2, int horizon, int idx,
                              double conf, double conf2, const arma::uvec& cumulate,
                              Rcpp::Nullable<arma::mat> exog,
                              int n_threads) {

    const int p      = var_result.p;
    const int c      = var_result.c;
    const int n_exog = var_result.n_exog;
    const int N      = static_cast<int>(y.n_cols);
    const int n_coef = static_cast<int>(var_result.beta.n_rows);
    const int T_iter = static_cast<int>(y.n_rows) - p;
    const int T_resid = static_cast<int>(var_result.residuals.n_rows);

    if (nboot1 <= 0 || nboot2 <= 0)
        Rcpp::stop("'nboot1' and 'nboot2' must be positive.");

    if (n_exog > 0 && exog.isNull())
        Rcpp::stop("Original VAR used exogenous variables. You must provide the 'exog' parameter.");
    if (n_exog == 0 && exog.isNotNull())
        Rcpp::stop("Original VAR did not use exogenous variables. Do not provide the 'exog' parameter.");

    int actual_threads = 1;
#ifdef _OPENMP
    actual_threads = (n_threads <= 0)
                         ? std::max(1, omp_get_max_threads() - 1)
                         : n_threads;
    omp_set_num_threads(actual_threads);
    Rprintf("[Pass 1] Bias estimation: %d reps, %d thread(s)\n",
                nboot1, actual_threads);
#else
    Rprintf("[Pass 1] Bias estimation: %d reps, single-threaded\n", nboot1);
#endif

    const bool has_exog = exog.isNotNull();
    arma::mat exog_mat;
    if (has_exog) exog_mat = Rcpp::as<arma::mat>(exog);
    const arma::mat* exog_ptr = has_exog ? &exog_mat : nullptr;

    // Generate all residual resampling indices on R's main thread. Workers
    // then consume fixed columns without touching R's RNG, and a serial mean
    // below preserves the same floating-point order for every thread count.
    arma::umat resample_indices(T_iter, nboot1, arma::fill::none);
    for (int b = 0; b < nboot1; ++b)
        resample_indices.col(b) = arma::randi<arma::uvec>(
            T_iter, arma::distr_param(0, T_resid - 1));

    arma::cube pass1_beta(N, n_coef, nboot1, arma::fill::none);
    arma::mat boot_mean(N, n_coef, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(actual_threads)
#endif
    for (int b = 0; b < nboot1; ++b) {
        arma::uvec indices = resample_indices.col(b);
        BootstrapVARResult boot_data =
            fBootstrapVAR_cpp(y, var_result, indices, exog_ptr);
        VARResult var_loop = has_exog
                                 ? fVAR_cpp(boot_data.ynext, p, c, exog_mat)
                                 : fVAR_cpp(boot_data.ynext, p, c);
        pass1_beta.slice(b) = var_loop.beta.t();
    }

    for (int b = 0; b < nboot1; ++b)
        boot_mean += pass1_beta.slice(b);
    boot_mean /= static_cast<double>(nboot1);

    arma::mat Beta_t;
    int corrections = 1;
    bias_correct_uhlig(var_result.beta, c, p, boot_mean, Beta_t, corrections);

    if (corrections > 1)
        Rprintf("[Bias correction] %d shrinkage iteration(s)\n", corrections);
    else
        Rprintf("[Bias correction] Full correction applied\n");

#ifdef _OPENMP
    Rprintf("[Pass 2] Bias-corrected bootstrap: %d reps, %d thread(s)\n",
                nboot2, actual_threads);
#else
    Rprintf("[Pass 2] Bias-corrected bootstrap: %d reps, single-threaded\n", nboot2);
#endif

    VARResult corrected_var = var_result;
    corrected_var.beta      = Beta_t.t();

    return fBootstrapUhlig_cpp(y, corrected_var, nboot2, horizon, idx,
                               conf, conf2, cumulate, exog, actual_threads);
}

//' @rdname fBootstrapUhlig
//' @inheritParams fBootstrapMax
//' @export
// [[Rcpp::export]]
Rcpp::List fBootstrapUhligCorrected(const arma::mat& y, const Rcpp::List& var_result,
                                     int nboot1, int nboot2, int horizon, int idx,
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
    BootstrapUhligResult res = fBootstrapUhligCorrected_cpp(y, vr, nboot1, nboot2,
                                                            horizon, idx - 1,
                                                            conf, conf2, cumulate_cpp,
                                                            exog, n_threads);

    Rcpp::NumericVector bootuhlig_out(res.bootuhlig_flat.begin(), res.bootuhlig_flat.end());
    bootuhlig_out.attr("dim") = Rcpp::IntegerVector::create(res.N, res.H, nboot2);

    return Rcpp::List::create(
        Rcpp::Named("bootuhlig") = bootuhlig_out,
        Rcpp::Named("upper")     = res.upper,
        Rcpp::Named("lower")     = res.lower,
        Rcpp::Named("upper2")    = res.upper2,
        Rcpp::Named("lower2")    = res.lower2,
        Rcpp::Named("boot_beta") = res.boot_beta
    );
}
