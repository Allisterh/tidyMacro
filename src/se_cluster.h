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

#endif // SE_CLUSTER_H
