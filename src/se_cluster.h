#ifndef SE_CLUSTER_H
#define SE_CLUSTER_H

// Cluster-robust long-run moments shared by fLPPanel and fLPDID. All sums
// are unnormalised; finite-sample factors (CR1, reghdfe) stay with the
// caller. Never calls the R API (OpenMP-safe).

#include <RcppArmadillo.h>
#include <vector>

// Score sums within clusters (G x k): row g = sum of the rows r of `scores`
// with cluster(r) equal to the g-th label in ascending order.
arma::mat se_cluster_sums_cpp(const arma::mat& scores, const arma::ivec& cluster);

// Same with compact ids already assigned: gid[r] in 0..G-1, row g = cluster g.
arma::mat se_cluster_sums_cpp(const arma::mat& scores,
                              const std::vector<arma::uword>& gid,
                              arma::uword G);

// One-way cluster meat S'S, S = se_cluster_sums_cpp(scores, cluster).
arma::mat se_cluster_cpp(const arma::mat& scores, const arma::ivec& cluster);

// Cameron-Gelbach-Miller two-way meat M_a + M_b - M_(a,b), where M_(a,b)
// clusters on the intersections of a and b.
arma::mat se_cluster_twoway_cpp(const arma::mat& scores,
                                const arma::ivec& cluster_a,
                                const arma::ivec& cluster_b);

// Meat of the intersection of two cluster dimensions: the term the
// Cameron-Gelbach-Miller two-way estimator subtracts.
arma::mat se_cluster_pair_cpp(const arma::mat& scores,
                              const arma::ivec& cluster_a,
                              const arma::ivec& cluster_b);

// Number of distinct clusters in a dimension.
arma::uword se_n_clusters_cpp(const arma::ivec& cluster);

// CR1 finite-sample factor, as applied by Stata and fixest:
//   G / (G - 1) * (N - 1) / (N - K)
// G clusters, N observations, K estimated parameters. Returns 1 when the
// factor is not defined (one cluster, or no residual degrees of freedom).
double se_cr1_factor_cpp(arma::uword G, arma::uword N, arma::uword K);

// True when every level of `fe` sits inside a single cluster, i.e. the fixed
// effect is nested in the clustering dimension. Nested fixed effects do not
// count towards K in the CR1 factor (fixest's fixef.K = "nested" rule).
bool se_fe_nested_cpp(const arma::ivec& fe, const arma::ivec& cluster);

// Number of distinct levels in a fixed-effect column.
arma::uword se_n_levels_cpp(const arma::ivec& fe);

#endif // SE_CLUSTER_H
