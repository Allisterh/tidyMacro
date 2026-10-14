#include <RcppArmadillo.h>
#include "se_newey_west.h"
// [[Rcpp::depends(RcppArmadillo)]]

double se_bartlett_weight_cpp(int l, int L) {
  return static_cast<double>(L + 1 - l) / static_cast<double>(L + 1);
}

double se_newey_west_lrv_cpp(const double* a, int n, int L) {
  double s = 0.0;
  for (int t = 0; t < n; ++t) s += a[t] * a[t];
  for (int l = 1; l <= L; ++l) {
    double g = 0.0;
    for (int t = l; t < n; ++t) g += a[t] * a[t - l];
    s += 2.0 * se_bartlett_weight_cpp(l, L) * g;
  }
  return s;
}

// One pass per lag over raw pointers: no temporaries, vectorisable.
void se_newey_west_lrv2_cpp(const double* a, const double* b, int n, int L,
                            double& v11, double& v22, double& v12) {
  double s11 = 0.0, s22 = 0.0, s12 = 0.0;
  for (int t = 0; t < n; ++t) {
    s11 += a[t] * a[t];
    s22 += b[t] * b[t];
    s12 += a[t] * b[t];
  }
  for (int l = 1; l <= L; ++l) {
    const double w = se_bartlett_weight_cpp(l, L);
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

static arma::mat nw_meat(const arma::mat& g, int L) {
  const arma::uword T = g.n_rows;
  arma::mat S = g.t() * g;
  for (int l = 1; l <= L && static_cast<arma::uword>(l) < T; ++l) {
    const arma::mat Gl = g.rows(l, T - 1).t() * g.rows(0, T - 1 - l);
    S += se_bartlett_weight_cpp(l, L) * (Gl + Gl.t());
  }
  return S;
}

arma::mat se_newey_west_cpp(const arma::mat& g, int L, bool prewhite) {
  const arma::uword T = g.n_rows;
  const arma::uword k = g.n_cols;
  if (!prewhite || T < k + 2) return nw_meat(g, L);

  // VAR(1) without intercept: g_t' = g_{t-1}' A' + r_t'
  const arma::mat G0 = g.rows(1, T - 1);
  const arma::mat G1 = g.rows(0, T - 2);
  arma::mat G1_pinv;
  if (!arma::pinv(G1_pinv, G1)) return nw_meat(g, L);
  arma::mat A = (G1_pinv * G0).t();

  // Keep the recoloring finite: cap the singular values at 0.97 when the
  // spectral radius exceeds 0.97 (Andrews and Monahan, 1992)
  const double cap = 0.97;
  arma::cx_vec ev;
  if (!arma::eig_gen(ev, A)) return nw_meat(g, L);
  if (arma::max(arma::abs(ev)) > cap) {
    arma::mat U, V;
    arma::vec s;
    if (!arma::svd(U, s, V, A)) return nw_meat(g, L);
    A = U * arma::diagmat(arma::clamp(s, 0.0, cap)) * V.t();
  }

  arma::mat D;
  if (!arma::inv(D, arma::eye(k, k) - A)) return nw_meat(g, L);
  const arma::mat R = G0 - G1 * A.t();
  return D * nw_meat(R, L) * D.t();
}
