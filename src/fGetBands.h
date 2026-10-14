#ifndef FGETBANDS_H
#define FGETBANDS_H

#include <RcppArmadillo.h>
#include <vector>

struct FGetBandsResult {
    arma::mat upper;
    arma::mat lower;
    arma::mat median;
};

struct FGetBands2Result {
    arma::mat upper;
    arma::mat lower;
    arma::mat upper2;
    arma::mat lower2;
    arma::mat median;
};

// Percentile (R quantile type 7) of a sorted array of length n.
double fPercentileSorted_cpp(const double* sorted, int n, double pct);

// Row-wise percentiles of a flat bootstrap matrix (rows = statistics,
// columns = draws). Column j of the result holds percentile pcts[j].
arma::mat fRowPercentiles_cpp(const arma::mat& flat,
                              const std::vector<double>& pcts,
                              int n_threads);

// {50 + conf/2, 50 - conf/2, 50 + conf2/2, 50 - conf2/2, 50}: upper, lower,
// upper2, lower2 and median columns for fRowPercentiles_cpp.
std::vector<double> fBandLevels_cpp(double conf, double conf2);

FGetBandsResult  fGetBands_cpp(const arma::cube& bootirf, double conf);
FGetBands2Result fGetBands2_cpp(const arma::cube& bootirf, double conf, double conf2);

#endif
