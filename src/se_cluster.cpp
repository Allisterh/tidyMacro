#include <RcppArmadillo.h>
#include "se_cluster.h"
#include <cstdint>
#include <unordered_map>
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat se_cluster_sums_cpp(const arma::mat& scores,
                              const std::vector<arma::uword>& gid,
                              arma::uword G) {
  arma::mat sums(G, scores.n_cols, arma::fill::zeros);
  for (arma::uword r = 0; r < scores.n_rows; ++r) {
    sums.row(gid[r]) += scores.row(r);
  }
  return sums;
}

arma::mat se_cluster_sums_cpp(const arma::mat& scores, const arma::ivec& cluster) {
  const arma::ivec labels = arma::sort(arma::unique(cluster));
  std::unordered_map<int, arma::uword> row_of;
  row_of.reserve(labels.n_elem * 2);
  for (arma::uword g = 0; g < labels.n_elem; ++g) {
    row_of.emplace(labels(g), g);
  }

  std::vector<arma::uword> gid(scores.n_rows);
  for (arma::uword r = 0; r < scores.n_rows; ++r) gid[r] = row_of[cluster(r)];
  return se_cluster_sums_cpp(scores, gid, labels.n_elem);
}

arma::mat se_cluster_cpp(const arma::mat& scores, const arma::ivec& cluster) {
  const arma::mat S = se_cluster_sums_cpp(scores, cluster);
  return S.t() * S;
}

// Score sums within intersections of two cluster dimensions, rows in order
// of first appearance. Duplicate (a, b) cells are summed.
static arma::mat se_cluster_pair_sums(
    const arma::mat& scores,
    const arma::ivec& cluster_a,
    const arma::ivec& cluster_b
) {
  std::unordered_map<std::uint64_t, arma::uword> row_of;
  row_of.reserve(scores.n_rows * 2);
  arma::mat sums(scores.n_rows, scores.n_cols, arma::fill::zeros);
  arma::uword n_groups = 0;

  for (arma::uword r = 0; r < scores.n_rows; ++r) {
    const std::uint64_t a = static_cast<std::uint32_t>(cluster_a(r));
    const std::uint64_t b = static_cast<std::uint32_t>(cluster_b(r));
    const std::uint64_t key = (a << 32) | b;
    auto it = row_of.find(key);
    if (it == row_of.end()) {
      row_of.emplace(key, n_groups);
      sums.row(n_groups) = scores.row(r);
      ++n_groups;
    } else {
      sums.row(it->second) += scores.row(r);
    }
  }

  sums.resize(n_groups, scores.n_cols);
  return sums;
}

arma::mat se_cluster_twoway_cpp(const arma::mat& scores,
                                const arma::ivec& cluster_a,
                                const arma::ivec& cluster_b) {
  const arma::mat Sa  = se_cluster_sums_cpp(scores, cluster_a);
  const arma::mat Sb  = se_cluster_sums_cpp(scores, cluster_b);
  const arma::mat Sab = se_cluster_pair_sums(scores, cluster_a, cluster_b);
  const arma::mat Ma  = Sa.t() * Sa;
  return Ma + Sb.t() * Sb - Sab.t() * Sab;
}
