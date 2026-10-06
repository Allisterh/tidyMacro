// -----------------------------------------------------------------------
// fNLLP.cpp — Nonlinear and state-dependent Local Projections (C++ engine)
//
//   sign  : y_{t+h} = a_h + b+_h max(s_t,0) + b-_h min(s_t,0) + g_h' x_t + u
//   cubic : y_{t+h} = a_h + b1_h s_t + b3_h s_t^3 + g_h' x_t + u
//   state : y_{t+h} = M_t (aH_h + bH_h s_t + gH_h' x_t)
//                   + (1-M_t)(aL_h + bL_h s_t + gL_h' x_t) + d_h' c_t + u
//
//   x_t are the RHS terms of the formula, c_t the terms declared common
//   (state design only). Cumulative LHS: y_{t+h} - y_{t-1}.
//
//   Newey-West HAC — full sandwich, same estimator and bandwidth as fLP:
//     V_h = (X'X)^{-1} G (X'X)^{-1}
//     G   = Gamma_0 + sum_{a=1..nwL} w_a (Gamma_a + Gamma_a'),
//     w_a = (nwL + 1 - a) / (nwL + 1)
//
//   Horizon loop is parallelized with OpenMP (each horizon is independent).
//
// Author: Dr. Muhsin Ciftci
// -----------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// CRITICAL: Define this BEFORE including RcppArmadillo to prevent conflicts
#define ARMA_DONT_USE_OPENMP

#include "fNLLP.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <R_ext/Print.h>

using namespace arma;


// =======================================================================
// Design matrix
// =======================================================================

arma::mat fNLLPDesign_cpp(
    const arma::mat&  X,
    int               shock_col,
    int               specification,
    const arma::vec&  state,
    const arma::uvec& common_cols,
    int&              c1,
    int&              c2
) {
  const int T = static_cast<int>(X.n_rows);
  const int k = static_cast<int>(X.n_cols);
  const arma::vec s = X.col(shock_col);

  // ---- sign / cubic: [1, X] with the shock column replaced by two -------
  if (specification == 0 || specification == 1) {
    arma::mat D(T, k + 2);
    D.col(0).ones();
    int j = 1;
    for (int c = 0; c < k; ++c) {
      if (c == shock_col) {
        if (specification == 0) {
          D.col(j)     = arma::clamp(s, 0.0, arma::datum::inf);
          D.col(j + 1) = arma::clamp(s, -arma::datum::inf, 0.0);
        } else {
          D.col(j)     = s;
          D.col(j + 1) = arma::pow(s, 3);
        }
        c1 = j;
        c2 = j + 1;
        j += 2;
      } else {
        D.col(j++) = X.col(c);
      }
    }
    return D;
  }

  // ---- state: [M, M*X_nc, (1-M), (1-M)*X_nc, X_c] ------------------------
  std::vector<char> is_common(k, 0);
  for (arma::uword i = 0; i < common_cols.n_elem; ++i) is_common[common_cols(i)] = 1;

  std::vector<int> nc;
  nc.reserve(k);
  for (int c = 0; c < k; ++c) if (!is_common[c]) nc.push_back(c);

  const int knc = static_cast<int>(nc.size());
  const int kc  = k - knc;
  const int kb  = 1 + knc;                       // columns per regime block

  const arma::vec lo = 1.0 - state;
  arma::mat D(T, 2 * kb + kc);
  D.col(0)  = state;
  D.col(kb) = lo;
  for (int i = 0; i < knc; ++i) {
    D.col(1 + i)      = state % X.col(nc[i]);
    D.col(kb + 1 + i) = lo % X.col(nc[i]);
    if (nc[i] == shock_col) {
      c1 = 1 + i;
      c2 = kb + 1 + i;
    }
  }
  int j = 2 * kb;
  for (int c = 0; c < k; ++c) if (is_common[c]) D.col(j++) = X.col(c);

  return D;
}


// =======================================================================
// Helpers
// =======================================================================

namespace {

// (X'X)^{-1} from X'X by Cholesky of the column-equilibrated matrix
// A = S X'X S, S = diag(1/sqrt(diag(X'X))). Equilibration makes the test
// scale-free (a trend next to a shock does not trigger it). Returns false
// when A is not numerically positive definite or rcond(A) < 1e-8, i.e.
// when normal equations could lose more than ~8 digits; the caller then
// uses the QR path, which also detects rank deficiency.
bool inv_xpx_chol(const arma::mat& XtX, arma::mat& out) {
  const arma::vec d = XtX.diag();
  if (!d.is_finite() || d.min() <= 0.0) return false;
  const arma::vec s  = 1.0 / arma::sqrt(d);
  const arma::mat SS = s * s.t();
  const arma::mat A  = XtX % SS;
  arma::mat L;
  if (!arma::chol(L, A, "lower")) return false;
  if (arma::rcond(A) < 1e-8) return false;
  const arma::mat Linv = arma::inv(arma::trimatl(L));
  out = (Linv.t() * Linv) % SS;
  return true;
}

// Least-squares solver for one design X, reused for every Y regressed on X.
//   Cholesky : beta = (X'X)^{-1} X'Y, (X'X)^{-1} from inv_xpx_chol.
//   QR       : when the Cholesky test fails. On the column-normalised design
//              Xs = X C^{-1}, C = diag(||x_j||), Xs = Q R:
//                beta      = C^{-1} R^{-1} Q'Y
//                (X'X)^{-1} = C^{-1} R^{-1} R^{-T} C^{-1}.
//              Normalising makes the rank test independent of the units of
//              the regressors; rank is judged from the singular values of the
//              kd x kd factor R. X'X is never formed on this path.
struct LSSolver {
  bool      qr = false;
  bool      ok = false;   // false: X is numerically rank deficient
  arma::mat inv;          // (X'X)^{-1}
  arma::mat Q, Rinv;      // QR path only
  arma::vec c;            // QR path only: column norms of X

  arma::mat beta(const arma::mat& X, const arma::mat& Y) const {
    if (!qr) return inv * (X.t() * Y);
    arma::mat b = Rinv * (Q.t() * Y);
    b.each_col() /= c;
    return b;
  }
};

LSSolver make_solver(const arma::mat& X, const arma::mat& XtX) {
  LSSolver S;
  if (inv_xpx_chol(XtX, S.inv)) {
    S.ok = true;
    return S;
  }
  S.qr = true;
  S.c  = arma::sqrt(arma::sum(arma::square(X), 0)).t();
  if (!S.c.is_finite() || S.c.min() <= 0.0) return S;     // zero column
  arma::mat Xs = X;
  Xs.each_row() /= S.c.t();
  arma::mat R;
  if (!arma::qr_econ(S.Q, R, Xs)) return S;
  const arma::vec sv = arma::svd(R);
  if (sv.min() <= sv.max() * std::max(X.n_rows, X.n_cols) * arma::datum::eps)
    return S;
  S.Rinv = arma::inv(arma::trimatu(R));
  S.inv  = (S.Rinv * S.Rinv.t()) / (S.c * S.c.t());
  S.ok   = true;
  return S;
}

// Bartlett-weighted long-run moments of two series a, b of length n:
//   v11 = sum_t a_t^2 + 2 sum_{l=1..L} w_l sum_{t>=l} a_t a_{t-l}
//   v22 = same for b
//   v12 = sum_t a_t b_t + sum_{l=1..L} w_l sum_{t>=l} (a_t b_{t-l} + b_t a_{t-l})
// w_l = (L + 1 - l) / (L + 1). With a_t = u_t P(t, j) these are the HAC
// variances and covariance of coefficients j (see the horizon loop).
// One pass per lag over raw pointers: no temporaries, vectorisable.
inline void bartlett2(const double* a, const double* b, int n, int L,
                      double& v11, double& v22, double& v12) {
  double s11 = 0.0, s22 = 0.0, s12 = 0.0;
  for (int t = 0; t < n; ++t) {
    s11 += a[t] * a[t];
    s22 += b[t] * b[t];
    s12 += a[t] * b[t];
  }
  for (int l = 1; l <= L; ++l) {
    const double w = static_cast<double>(L + 1 - l) / static_cast<double>(L + 1);
    double g11 = 0.0, g22 = 0.0, g12 = 0.0;
    for (int t = l; t < n; ++t) {
      g11 += a[t] * a[t - l];
      g22 += b[t] * b[t - l];
      g12 += a[t] * b[t - l] + b[t] * a[t - l];
    }
    s11 += 2.0 * w * g11;
    s22 += 2.0 * w * g22;
    s12 += w * g12;
  }
  v11 = s11; v22 = s22; v12 = s12;
}

// Bartlett-weighted long-run variance of one series (v11 of bartlett2).
inline double bartlett1(const double* a, int n, int L) {
  double s = 0.0;
  for (int t = 0; t < n; ++t) s += a[t] * a[t];
  for (int l = 1; l <= L; ++l) {
    const double w = static_cast<double>(L + 1 - l) / static_cast<double>(L + 1);
    double g = 0.0;
    for (int t = l; t < n; ++t) g += a[t] * a[t - l];
    s += 2.0 * w * g;
  }
  return s;
}

} // namespace


// =======================================================================
// Internal C++ implementation
// =======================================================================

NLLPResult fNLLP_internal(
    const arma::mat&  Y,
    const arma::mat&  X,
    int               H,
    int               shock_col,
    int               specification,
    const arma::vec&  state,
    const arma::uvec& common_cols,
    int               nw_lags_base,
    bool              store_full,
    bool              cumulative,
    bool              balanced,
    int               n_threads,
    int               nw_offset,
    bool              verbose,
    const arma::mat&  Y_pre
) {
  // ---------- dimensions -----------------------------------------------
  const int T  = static_cast<int>(X.n_rows);
  const int TY = static_cast<int>(Y.n_rows);
  const int ny = static_cast<int>(Y.n_cols);

  if (specification < 0 || specification > 2)
    Rcpp::stop("fNLLP: specification must be 0 (sign), 1 (cubic) or 2 (state).");
  if (shock_col < 0 || shock_col >= static_cast<int>(X.n_cols))
    Rcpp::stop("fNLLP: shock_col must index a column of X.");
  if (balanced) {
    if (TY < T + H)
      Rcpp::stop("fNLLP: balanced = TRUE needs nrow(Y) >= nrow(X) + H.");
  } else if (TY != T) {
    Rcpp::stop("fNLLP: Y and X must have the same number of rows.");
  }
  if (specification == 2) {
    if (static_cast<int>(state.n_elem) != T)
      Rcpp::stop("fNLLP: state must have nrow(X) elements.");
    for (arma::uword i = 0; i < common_cols.n_elem; ++i) {
      if (static_cast<int>(common_cols(i)) >= static_cast<int>(X.n_cols))
        Rcpp::stop("fNLLP: common_cols must index columns of X.");
      if (static_cast<int>(common_cols(i)) == shock_col)
        Rcpp::stop("fNLLP: the shock cannot be a common column.");
    }
  }

  // ---------- design (built once, sliced per horizon) -------------------
  int c1 = -1, c2 = -1;
  const arma::mat D  = fNLLPDesign_cpp(X, shock_col, specification, state,
                                       common_cols, c1, c2);
  const int       kd = static_cast<int>(D.n_cols);

  // ---------- cumulative long-difference base ---------------------------
  if (cumulative && Y_pre.n_rows > 0 &&
      (Y_pre.n_rows != 1 || static_cast<int>(Y_pre.n_cols) != ny)) {
    Rcpp::stop("fNLLP: Y_pre must be a 1 x ncol(Y) matrix, or empty.");
  }
  const bool has_pre = cumulative && (Y_pre.n_rows == 1);
  const int  t0      = (cumulative && !has_pre) ? 1 : 0;

  const int T_min = balanced ? T - t0 : T - H - t0;
  if (T_min <= kd)
    Rcpp::stop("fNLLP: not enough observations for the largest horizon.");

  arma::mat Ylag;
  if (cumulative) {
    Ylag.set_size(T, ny);
    if (has_pre) Ylag.row(0) = Y_pre.row(0);
    else         Ylag.row(0).zeros();
    if (T > 1) Ylag.rows(1, T - 1) = Y.rows(0, T - 2);
  }

  // ---------- X'X and solvers, before the parallel loop ----------------
  // balanced  : one design for all horizons, so the design rows, the solver
  //             (Cholesky or QR) and the projection columns are built once
  //             here and shared by every horizon.
  // shrinking : X'X over rows t0..T-1-h. Start from the smallest sample
  //             (h = H) and add one row outer product per horizon, so no
  //             horizon re-forms X'X from scratch. Each horizon builds its
  //             solver inside the parallel loop.
  const int n_inv = balanced ? 1 : H + 1;
  std::vector<arma::mat> XtX(n_inv);
  arma::mat Xb;                                  // balanced design rows
  {
    const int t_last_H = balanced ? T - 1 : T - 1 - H;
    if (balanced) {
      Xb = D.rows(t0, t_last_H);
      XtX[0] = Xb.t() * Xb;
    } else {
      const arma::mat Dh = D.rows(t0, t_last_H);
      XtX[n_inv - 1] = Dh.t() * Dh;
      for (int i = n_inv - 2; i >= 0; --i) {
        const arma::rowvec r = D.row(T - 1 - i);
        XtX[i] = XtX[i + 1] + r.t() * r;
      }
    }
  }
  const arma::uvec tgt = {static_cast<arma::uword>(c1), static_cast<arma::uword>(c2)};
  LSSolver  Sb;                                  // balanced solver
  arma::mat Pb;                                  // balanced projection columns
  if (balanced) {
    Sb = make_solver(Xb, XtX[0]);
    if (Sb.ok) Pb = store_full ? arma::mat(Xb * Sb.inv) : arma::mat(Xb * Sb.inv.cols(tgt));
  }

  // ---------- output storage -------------------------------------------
  NLLPResult out;
  out.irfs.zeros(H + 1, ny, 2);
  out.irfs_se.zeros(H + 1, ny, 2);
  out.irfs_cov.zeros(H + 1, ny);
  out.diff.zeros(H + 1, ny);
  out.diff_se.zeros(H + 1, ny);
  out.nobs.assign(H + 1, 0);
  if (store_full) {
    out.betas.resize(H + 1);
    out.ses.resize(H + 1);
  }

  std::vector<char> rank_ok(H + 1, 1);

  // ---------- threading setup ------------------------------------------
  int actual_threads = 1;
#ifdef _OPENMP
  actual_threads = (n_threads <= 0)
                       ? std::max(1, omp_get_max_threads())
                       : n_threads;
  actual_threads = std::min(actual_threads, H + 1);
  if (verbose) {
    Rprintf("fNLLP: using %d thread(s) for parallel horizon loop...\n",
            actual_threads);
  }
#else
  (void) n_threads;
  if (verbose) {
    Rprintf("fNLLP: OpenMP not available. Running single-threaded.\n");
  }
#endif

  // ---------- horizon loop (parallel over h) ---------------------------
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(actual_threads)
#endif
  for (int h = 0; h <= H; h++) {

    // ---- align data: regressors at t, outcome at t + h ---------------
    const int t_last = balanced ? T - 1 : T - 1 - h;
    arma::mat Yh = Y.rows(t0 + h, t_last + h);
    if (cumulative) Yh -= Ylag.rows(t0, t_last);
    const int Th = static_cast<int>(Yh.n_rows);

    // ---- design, solver and projection: shared when balanced ----------
    arma::mat Xloc, Ploc;
    LSSolver  Sloc;
    const arma::mat* Xp = &Xb;
    const LSSolver*  S  = &Sb;
    const arma::mat* Pp = &Pb;
    if (!balanced) {
      Xloc = D.rows(t0, t_last);
      Sloc = make_solver(Xloc, XtX[h]);
      Xp = &Xloc;
      S  = &Sloc;
    }
    out.nobs[h] = Th;
    if (!S->ok) {                  // not identified: the wrapper stops anyway
      rank_ok[h] = 0;
      continue;
    }
    if (!balanced) {
      Ploc = store_full ? arma::mat(Xloc * Sloc.inv) : arma::mat(Xloc * Sloc.inv.cols(tgt));
      Pp = &Ploc;
    }
    const arma::mat& Xreg = *Xp;
    const arma::mat& P    = *Pp;

    // ---- OLS for all equations ----------------------------------------
    const arma::mat Beta = S->beta(Xreg, Yh);              // kd x ny
    const arma::mat U    = Yh - Xreg * Beta;               // Th x ny

    const int nwL = std::min(std::max(nw_lags_base + h + nw_offset, 0), Th - 1);

    // ---- Newey-West HAC in projected form ------------------------------
    // With P = X (X'X)^{-1}, V(j,j') = sum_a w_a sum_t a_j(t) a_j'(t-a)
    // where a_j(t) = u_t P(t,j). Only the columns that are reported are
    // formed: c1 and c2 on the fast path, every column for store_full.
    // P = X (X'X)^{-1}: all kd columns for store_full, else only c1, c2.
    arma::mat Se_h;
    if (store_full) Se_h.set_size(kd, ny);
    const arma::uword j1 = store_full ? static_cast<arma::uword>(c1) : 0;
    const arma::uword j2 = store_full ? static_cast<arma::uword>(c2) : 1;
    arma::vec a1(Th), a2(Th), aj(Th);

    for (int eq = 0; eq < ny; eq++) {
      const arma::vec u = U.col(eq);
      a1 = u % P.col(j1);
      a2 = u % P.col(j2);
      double v11, v22, v12;
      bartlett2(a1.memptr(), a2.memptr(), Th, nwL, v11, v22, v12);

      if (store_full) {
        for (int j = 0; j < kd; j++) {
          double vj;
          if (j == c1)      vj = v11;
          else if (j == c2) vj = v22;
          else {
            aj = u % P.col(j);
            vj = bartlett1(aj.memptr(), Th, nwL);
          }
          Se_h(j, eq) = std::sqrt(std::max(0.0, vj));
        }
      }

      out.irfs(h, eq, 0)    = Beta(c1, eq);
      out.irfs(h, eq, 1)    = Beta(c2, eq);
      out.irfs_se(h, eq, 0) = std::sqrt(std::max(0.0, v11));
      out.irfs_se(h, eq, 1) = std::sqrt(std::max(0.0, v22));
      out.irfs_cov(h, eq)   = v12;
      if (specification == 1) {
        out.diff(h, eq)    = arma::datum::nan;
        out.diff_se(h, eq) = arma::datum::nan;
      } else {
        out.diff(h, eq)    = Beta(c1, eq) - Beta(c2, eq);
        out.diff_se(h, eq) = std::sqrt(std::max(0.0, v11 + v22 - 2.0 * v12));
      }
    }

    if (store_full) {
      out.betas[h] = Beta;
      out.ses[h]   = Se_h;
    }

  } // end horizon loop

  for (int h = 0; h <= H; h++) {
    if (!rank_ok[h]) {
      out.rank_deficient = true;
      out.rank_fail_h    = h;
      break;
    }
  }

  return out;
}


// =======================================================================
// R-callable wrapper — exported to R via Rcpp
// =======================================================================

// [[Rcpp::export]]
Rcpp::List fNLLP_cpp(
    const arma::mat&    Y,
    const arma::mat&    X,
    int                 H,
    int                 shock_col,
    int                 specification,
    const arma::vec&    state,
    Rcpp::IntegerVector common_cols  = Rcpp::IntegerVector::create(),
    int                 nw_lags_base = 0,
    bool                store_full   = false,
    bool                cumulative   = false,
    bool                balanced     = false,
    int                 n_threads    = 0,
    int                 nw_offset    = 1,
    bool                verbose      = false,
    Rcpp::Nullable<arma::mat> Y_pre  = R_NilValue
) {
  if (Y.n_cols == 0) Rcpp::stop("fNLLP_cpp: Y must have at least one column.");
  if (X.n_cols == 0) Rcpp::stop("fNLLP_cpp: X must have at least one column.");
  if (H < 0)         Rcpp::stop("fNLLP_cpp: H must be non-negative.");
  if (nw_lags_base < 0)
    Rcpp::stop("fNLLP_cpp: nw_lags_base must be non-negative.");
  if (!Y.is_finite() || !X.is_finite())
    Rcpp::stop("fNLLP_cpp: Y and X must be finite (no NA/NaN/Inf).");
  if (specification == 2 && !state.is_finite())
    Rcpp::stop("fNLLP_cpp: state must be finite (no NA/NaN/Inf).");

  arma::uvec common(common_cols.size());
  for (int i = 0; i < common_cols.size(); ++i) {
    if (common_cols[i] == NA_INTEGER || common_cols[i] < 0)
      Rcpp::stop("fNLLP_cpp: common_cols must be non-negative (0-indexed).");
    common(i) = static_cast<arma::uword>(common_cols[i]);
  }

  arma::mat Y_pre_mat;
  if (Y_pre.isNotNull()) {
    Y_pre_mat = Rcpp::as<arma::mat>(Y_pre.get());
    if (!Y_pre_mat.is_finite())
      Rcpp::stop("fNLLP_cpp: Y_pre must be finite (no NA/NaN/Inf).");
  }

  NLLPResult res = fNLLP_internal(Y, X, H, shock_col, specification, state,
                                  common, nw_lags_base, store_full, cumulative,
                                  balanced, n_threads, nw_offset, verbose,
                                  Y_pre_mat);

  if (res.rank_deficient) {
    Rcpp::stop(
      "fNLLP_cpp: the regression design is rank deficient at horizon %d. "
      "The shock coefficients are not identified. Drop collinear columns "
      "(in the state design the intercept is already split into the state "
      "and one minus the state, so do not add a constant-like column).",
      res.rank_fail_h);
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("irfs")     = res.irfs,
    Rcpp::Named("irfs_se")  = res.irfs_se,
    Rcpp::Named("irfs_cov") = res.irfs_cov,
    Rcpp::Named("diff")     = res.diff,
    Rcpp::Named("diff_se")  = res.diff_se,
    Rcpp::Named("nobs")     = Rcpp::wrap(res.nobs)
  );

  if (store_full) {
    Rcpp::List betas_list(res.betas.size());
    Rcpp::List ses_list(res.ses.size());
    for (size_t h = 0; h < res.betas.size(); ++h) {
      betas_list[h] = res.betas[h];
      ses_list[h]   = res.ses[h];
    }
    out["betas"] = betas_list;
    out["ses"]   = ses_list;
  }

  return out;
}
