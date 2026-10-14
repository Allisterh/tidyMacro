#ifndef SE_NEWEY_WEST_H
#define SE_NEWEY_WEST_H

// Newey-West (Bartlett kernel) long-run moments shared by fLP, fNLLP, fLPIV,
// fOLS and fMSW. All sums are unnormalised; callers divide by T where their
// estimator does. Weights: w_l = (L + 1 - l) / (L + 1), l = 0..L.

#include <RcppArmadillo.h>

double se_bartlett_weight_cpp(int l, int L);

// Long-run variance of one series a of length n:
//   sum_t a_t^2 + 2 sum_{l=1..L} w_l sum_{t>=l} a_t a_{t-l}
double se_newey_west_lrv_cpp(const double* a, int n, int L);

// Long-run variances v11, v22 and covariance v12 of two series a, b:
//   v12 = sum_t a_t b_t + sum_{l=1..L} w_l sum_{t>=l} (a_t b_{t-l} + b_t a_{t-l})
void se_newey_west_lrv2_cpp(const double* a, const double* b, int n, int L,
                            double& v11, double& v22, double& v12);

// HAC "meat" of a score matrix g (T x k, row t = g_t'):
//   Gamma_0 + sum_{l=1..L} w_l (Gamma_l + Gamma_l'),  Gamma_l = sum_{t>=l} g_t g_{t-l}'
// prewhite = true applies Andrews & Monahan (1992) VAR(1) prewhitening:
// fit g_t = A g_{t-1} + r_t, take the meat of r with the same L and recolor
// with D = (I - A)^{-1}. A is shrunk (singular values capped at 0.97) when
// its spectral radius exceeds 0.97. Never calls the R API (OpenMP-safe).
arma::mat se_newey_west_cpp(const arma::mat& g, int L, bool prewhite);

#endif // SE_NEWEY_WEST_H
