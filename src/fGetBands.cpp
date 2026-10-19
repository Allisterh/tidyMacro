#include "fGetBands.h"
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(openmp)]]
#include "omp_threads.h"

double fPercentileSorted_cpp(const double* sorted, int n, double pct) {
    const double raw  = (pct / 100.0) * (n - 1);
    const int    lo   = static_cast<int>(raw);
    const double frac = raw - lo;
    if (frac < 1e-12 || lo + 1 >= n) return sorted[lo];
    return sorted[lo] * (1.0 - frac) + sorted[lo + 1] * frac;
}

arma::mat fRowPercentiles_cpp(const arma::mat& flat,
                              const std::vector<double>& pcts,
                              int n_threads) {
    const int n_rows = static_cast<int>(flat.n_rows);
    const int B      = static_cast<int>(flat.n_cols);
    const int n_pct  = static_cast<int>(pcts.size());
    arma::mat out(n_rows, n_pct, arma::fill::none);
    const int n_thr = fResolveThreads_cpp(n_threads);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_thr)
#endif
    for (int i = 0; i < n_rows; ++i) {
        std::vector<double> sv(B);
        for (int b = 0; b < B; ++b) sv[b] = flat(i, b);
        std::sort(sv.begin(), sv.end());
        for (int j = 0; j < n_pct; ++j)
            out(i, j) = fPercentileSorted_cpp(sv.data(), B, pcts[j]);
    }
    return out;
}

std::vector<double> fBandLevels_cpp(double conf, double conf2) {
    return {50.0 + conf * 0.5, 50.0 - conf * 0.5,
            50.0 + conf2 * 0.5, 50.0 - conf2 * 0.5, 50.0};
}

// bootirf is N x H x B; its memory is the flat (N*H) x B matrix.
static arma::mat cube_flat(const arma::cube& bootirf) {
    return arma::mat(const_cast<double*>(bootirf.memptr()),
                     bootirf.n_rows * bootirf.n_cols, bootirf.n_slices,
                     false, true);
}

static arma::mat as_NH(const arma::mat& q, int col, int N, int H) {
    return arma::reshape(q.col(col), N, H);
}

FGetBandsResult fGetBands_cpp(const arma::cube& bootirf, double conf) {
    const int N = static_cast<int>(bootirf.n_rows);
    const int H = static_cast<int>(bootirf.n_cols);
    arma::mat q = fRowPercentiles_cpp(cube_flat(bootirf),
                                      {50.0 + conf / 2.0, 50.0 - conf / 2.0, 50.0},
                                      1);
    FGetBandsResult result;
    result.upper  = as_NH(q, 0, N, H);
    result.lower  = as_NH(q, 1, N, H);
    result.median = as_NH(q, 2, N, H);
    return result;
}

FGetBands2Result fGetBands2_cpp(const arma::cube& bootirf, double conf, double conf2) {
    const int N = static_cast<int>(bootirf.n_rows);
    const int H = static_cast<int>(bootirf.n_cols);
    arma::mat q = fRowPercentiles_cpp(cube_flat(bootirf),
                                      fBandLevels_cpp(conf, conf2), 1);
    FGetBands2Result res;
    res.upper  = as_NH(q, 0, N, H);
    res.lower  = as_NH(q, 1, N, H);
    res.upper2 = as_NH(q, 2, N, H);
    res.lower2 = as_NH(q, 3, N, H);
    res.median = as_NH(q, 4, N, H);
    return res;
}

//' Compute Confidence Bands from Bootstrap IRF Cube
//'
//' Computes upper, lower, and median bands from a 3-D bootstrap IRF array
//' along the third (bootstrap) dimension. Percentiles follow R's
//' \code{quantile(type = 7)}.
//'
//' @param bootirf N x (hor+1) x nboot cube of bootstrapped IRFs.
//' @param conf Confidence level (e.g. 68 for 68\% band). Upper and lower
//'   quantiles are \eqn{(50 + conf/2)}\% and \eqn{(50 - conf/2)}\%.
//'   Default 90.
//'
//' @return List with three N x (hor+1) matrices: upper, lower, median.
//'
//' @examples
//' set.seed(123)
//' # 200 bootstrap draws of a 2-variable IRF over horizons 0 to 12
//' bootirf <- array(rnorm(2 * 13 * 200), c(2, 13, 200))
//' bands   <- fGetBands(bootirf, conf = 68)
//' bands$median[, 1:3]
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fGetBands(const arma::cube& bootirf, double conf = 90.0) {
    FGetBandsResult r = fGetBands_cpp(bootirf, conf);
    return Rcpp::List::create(
        Rcpp::Named("upper")  = r.upper,
        Rcpp::Named("lower")  = r.lower,
        Rcpp::Named("median") = r.median
    );
}
