#ifndef LP_SAMPLE_H
#define LP_SAMPLE_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>

// -----------------------------------------------------------------------
// Horizon samples of the time-series LP engines (fLP.cpp, fLPIV.cpp and
// the unbalanced path of fNLLP.cpp).
//
// Row r of Y and row r of the regressors refer to the same date; rows of Y
// past the last regressor row hold the outcomes' leads. Y may contain NaN
// (outcome not observed). Y_pre is the outcome one period before row 0
// (1 x n_y, may hold NaN, or empty when no earlier date exists).
//
// For horizon h and outcome j the sample is every regressor row t with
// y_{t+h} observed and, for a cumulative LP, y_{t-1} observed: Stata's
// F.h y followed by casewise deletion, so the outcomes never shorten the
// regressor sample and one outcome's NaN never touches another's sample.
//
// Each outcome must be observed on one block of consecutive rows over the
// rows the engine reads (leading and trailing NaN are fine). Every sample
// is then a block [lo, hi] of consecutive regressor rows, so the Newey-West
// autocovariances keep their spacing. With [elo, ehi] the observed block of
// outcome j (elo = -1 when the cumulative base Y_pre is observed and
// adjacent):
//   level      : t in [max(0, elo - h), min(T - 1, ehi - h)]
//   cumulative : t in [max(0, elo + 1), min(T - 1, ehi - h)]
// -----------------------------------------------------------------------

struct LPSamples {
  arma::imat lo;     // (H+1) x n_y — first regressor row of the sample
  arma::imat hi;     // (H+1) x n_y — last regressor row (hi < lo: empty)
  arma::imat nobs;   // (H+1) x n_y — observations, hi - lo + 1 floored at 0
};

// T = number of regressor rows. `fn` prefixes the error messages. Called
// before any parallel region, so Rcpp::stop is safe here.
inline LPSamples lp_horizon_samples_cpp(const arma::mat& Y,
                                        const arma::mat& Y_pre,
                                        int T, int H, bool cumulative,
                                        const char* fn) {
  const int ny = static_cast<int>(Y.n_cols);
  const int TY = static_cast<int>(Y.n_rows);
  const int R  = std::min(TY, T + H);            // rows ever read as y_{t+h}
  const bool has_pre = cumulative && Y_pre.n_rows == 1;

  LPSamples s;
  s.lo.set_size(H + 1, ny);
  s.hi.set_size(H + 1, ny);
  s.nobs.set_size(H + 1, ny);

  for (int j = 0; j < ny; ++j) {
    // Extended column: index -1 is the cumulative base, 0..R-1 are rows of Y.
    auto observed = [&](int r) {
      return r < 0 ? (has_pre && !std::isnan(Y_pre(0, j)))
                   : !std::isnan(Y(r, j));
    };
    const int r0 = cumulative ? -1 : 0;
    int elo = R, ehi = -2;                       // empty unless found below
    for (int r = r0; r < R; ++r) {
      if (!observed(r)) continue;
      if (elo == R) elo = r;
      ehi = r;
    }
    for (int r = elo + 1; r < ehi; ++r) {
      if (!observed(r))
        Rcpp::stop("%s: outcome column %d is missing at row %d, inside the "
                   "rows where it is observed. Each horizon sample must be "
                   "a block of consecutive dates; fill the gap or subset "
                   "the data to a contiguous span.", fn, j + 1, r + 1);
    }

    for (int h = 0; h <= H; ++h) {
      const int lo = cumulative ? std::max(0, elo + 1) : std::max(0, elo - h);
      const int hi = std::min(T - 1, ehi - h);
      s.lo(h, j)   = lo;
      s.hi(h, j)   = hi;
      s.nobs(h, j) = std::max(0, hi - lo + 1);
    }
  }
  return s;
}

// Outcomes of one horizon grouped by identical samples, so outcomes that
// share [lo, hi] share one factorisation of the design.
struct LPGroup {
  int lo, hi;
  arma::uvec cols;   // outcome columns in this group
};

inline std::vector<LPGroup> lp_groups_cpp(const LPSamples& s, int h) {
  std::vector<LPGroup> g;
  std::vector<std::vector<arma::uword>> cols;
  for (arma::uword j = 0; j < s.lo.n_cols; ++j) {
    const int lo = s.lo(h, j), hi = s.hi(h, j);
    size_t k = 0;
    while (k < g.size() && (g[k].lo != lo || g[k].hi != hi)) ++k;
    if (k == g.size()) {
      g.push_back({lo, hi, arma::uvec()});
      cols.emplace_back();
    }
    cols[k].push_back(j);
  }
  for (size_t k = 0; k < g.size(); ++k)
    g[k].cols = arma::conv_to<arma::uvec>::from(cols[k]);
  return g;
}

// Rows r0..r1 of M restricted to a group's columns. A group holding every
// column (the usual case) is a contiguous slice, so no gather is needed.
inline arma::mat lp_block_cpp(const arma::mat& M, int r0, int r1,
                              const arma::uvec& cols) {
  if (cols.n_elem == M.n_cols) return M.rows(r0, r1);
  return M.submat(arma::regspace<arma::uvec>(r0, r1), cols);
}

#endif // LP_SAMPLE_H
