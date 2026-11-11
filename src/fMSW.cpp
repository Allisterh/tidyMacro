// [[Rcpp::depends(RcppArmadillo)]]

#include "fMSW.h"
#include "fVAR.h"
#include "fWoldIRF.h"
#include "se_newey_west.h"
#include "fLagMakerMatrix.h"
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <algorithm>

// =========================================================================
// Internal helpers
// =========================================================================

// Score derivatives of the IRF numerator, for every horizon and variable.
//
// MSW need the 1 x n^2 p row  kron(Gamma', e_j') * G_h,  where
//   G_h = sum_{k=0}^{h-1} kron( (Alut^k J')', C_{h-1-k} )
// is the derivative of vec(C_h) wrt vec(AL) (Alut = companion matrix,
// J = [I_n | 0]). Since kron(Gamma', e_j') picks row j of each n x n block,
//   kron(Gamma', e_j') G_h = row j of Y_h,
//   Y_h = sum_{k=0}^{h-1} kron(gamma_k, C_{h-1-k})        (n x n^2 p),
//   gamma_k = (Alut^k J' Gamma)'                            (1 x np),
// so the n^2 x n^2 p matrices G_h are never formed. Column-major,
//   Y_h = reshape( sum_c vec(C_c) gamma_{h-1-c}, n, n^2 p ),
// one small matrix product per horizon.
// Returns Y: n*(hor+1) x n^2 p, rows h*n .. h*n+n-1 hold Y_h (Y_0 = 0).
static arma::mat irf_score_rows(const arma::mat& AL,
                                const arma::vec& Gamma,
                                const std::vector<arma::mat>& C_vec,  // C_0..C_hor
                                int p, int hor, int n) {
    const int np  = n * p;
    const int n2p = n * np;
    arma::mat Y(n * (hor + 1), n2p, arma::fill::zeros);
    if (hor < 1) return Y;

    // Companion matrix Alut (np x np)
    arma::mat Alut(np, np, arma::fill::zeros);
    Alut.rows(0, n - 1) = AL;
    if (p > 1) {
        Alut.submat(n, 0, np - 1, n * (p - 1) - 1) =
            arma::eye(n * (p - 1), n * (p - 1));
    }

    // Gam row k = gamma_k = (Alut^k J' Gamma)', k = 0..hor-1  (hor x np)
    arma::mat Gam(hor, np);
    arma::vec s(np, arma::fill::zeros);
    s.head(n) = Gamma;
    for (int k = 0; k < hor; ++k) {
        Gam.row(k) = s.t();
        if (k < hor - 1) s = Alut * s;
    }

    // Cmat column c = vec(C_c), c = 0..hor-1  (n^2 x hor)
    arma::mat Cmat(n * n, hor);
    for (int c = 0; c < hor; ++c) Cmat.col(c) = arma::vectorise(C_vec[c]);

    for (int h = 1; h <= hor; ++h) {
        // M = sum_{c=0}^{h-1} vec(C_c) gamma_{h-1-c}  (n^2 x np)
        const arma::mat M = Cmat.cols(0, h - 1) *
                            arma::flipud(Gam.rows(0, h - 1));
        Y.rows(h * n, h * n + n - 1) = arma::reshape(M, n, n2p);
    }
    return Y;
}

// =========================================================================
// Scores of [vec(Ahat); vec(Gammahat)]  (CovAhat_Sigmahat_Gamma equivalent)
//
// MSW define WHat, the asymptotic variance of [vec(Ahat)', vec(Gammahat)']',
// as the [vec(A), vec(Gamma)] rows and columns of Shat * HAC(g_t) * Shat',
// with g_t = vec(eta_t [x_t; eta_t; z_t]'). The HAC estimator is a quadratic
// form in the demeaned scores, so WHat = HAC(h_t) / T1, where h_t stacks the
// [vec(A), vec(Gamma)] rows of Shat g_t. h_t has a closed form:
//   vec(A) rows    : vec(eta_t (P x_t)'),            P = rows m..m+np-1 of Q1^{-1}
//   vec(Gamma) rows: vec(eta_t (z_t - Q2 Q1^{-1} x_t)')
// i.e. h_t = vec(eta_t w_t'), w_t = [P x_t; z_t - Q2 Q1^{-1} x_t].
// Returns the demeaned T1 x (n^2 p + nk) matrix H with row t = h_t'.
// fMSW_cpp never forms WHat: every quadratic form a' WHat b it needs is the
// NW long-run covariance of the series H a and H b, divided by T1.
// =========================================================================
static arma::mat msw_scores(
        int p,
        const arma::mat& X,      // T1 x (m + n*p) design matrix
        const arma::mat& Z,      // T1 x k instrument
        const arma::mat& eta) {  // n x T1 residuals (transposed)

    const int n  = (int)eta.n_rows;
    const int T1 = (int)eta.n_cols;
    const int k  = (int)Z.n_cols;
    const int np = n * p;
    const int m  = (int)X.n_cols - np;   // deterministic regressors (intercept count)
    const int nw = np + k;               // length of w_t

    // Q1 = X'X / T1, Q2 = Z'X / T1
    const arma::mat Q1inv = arma::inv_sympd((X.t() * X) / T1);
    const arma::mat Q2    = (Z.t() * X) / T1;

    // W: T1 x (np + k), row t = w_t'
    arma::mat W(T1, nw);
    W.cols(0, np - 1)  = X * Q1inv.rows(m, m + np - 1).t();
    W.cols(np, nw - 1) = Z - X * (Q2 * Q1inv).t();

    // H: T1 x n*(np + k), column q*n + i = eta_i,t * w_q,t  (row t = vec(eta_t w_t')')
    arma::mat H(T1, n * nw);
    for (int q = 0; q < nw; ++q) {
        const double* wq = W.colptr(q);
        for (int i = 0; i < n; ++i) {
            double* hc = H.colptr(q * n + i);
            for (int t = 0; t < T1; ++t) hc[t] = eta(i, t) * wq[t];
        }
    }

    H.each_row() -= arma::mean(H, 0);
    return H;
}

// =========================================================================
// fMSW_cpp  — main Montiel-Olea, Stock, Watson (2021) inference
// =========================================================================
MSWResult fMSW_cpp(const arma::mat& AL,
                   const arma::mat& Sigma_sub,
                   const arma::mat& eta,
                   const arma::mat& X,
                   const arma::mat& Z,
                   int p, int hor,
                   int nvar,
                   double scale,
                   double confidence,
                   int NWlags) {

    int n  = (int)AL.n_rows;
    int T1 = (int)eta.n_cols;

    // Critical value (two-sided chi^2 with 1 df)
    // norminv(1-(1-conf)/2)^2 = qnorm(0.5+conf/2)^2
    // Use the approximation / exact value via R's qnorm via Rcpp
    double alpha_tail = (1.0 - confidence) / 2.0;
    // Beasley-Springer-Moro approximation for the normal quantile
    auto qnorm_approx = [](double p_val) -> double {
        // Rational approximation for tail quantiles
        const double a[] = {-3.969683028665376e+01,  2.209460984245205e+02,
                            -2.759285104469687e+02,  1.383577518672690e+02,
                            -3.066479806614716e+01,  2.506628277459239e+00};
        const double b[] = {-5.447609879822406e+01,  1.615858368580409e+02,
                            -1.556989798598866e+02,  6.680131188771972e+01,
                            -1.328068155288572e+01};
        const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
                            -2.400758277161838e+00, -2.549732539343734e+00,
                             4.374664141464968e+00,  2.938163982698783e+00};
        const double d[] = { 7.784695709041462e-03,  3.224671290700398e-01,
                              2.445134137142996e+00,  3.754408661907416e+00};
        double q, r;
        if (p_val < 0.02425) {
            q = std::sqrt(-2.0 * std::log(p_val));
            return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                   ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
        } else if (p_val <= 0.97575) {
            q = p_val - 0.5;
            r = q * q;
            return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
                   (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
        } else {
            q = std::sqrt(-2.0 * std::log(1.0 - p_val));
            return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                    ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
        }
    };
    double z_crit  = qnorm_approx(1.0 - alpha_tail);
    double critval = z_crit * z_crit;   // chi^2(1) critical value

    // Gamma = eta * Z / T1   (n x k)
    arma::mat Gamma = (eta * Z) / T1;

    // Demeaned scores of [vec(Ahat); vec(Gammahat)]: T1 x (n^2 p + nk)
    const arma::mat H = msw_scores(p, X, Z, eta);

    const int n2p = n * n * p;
    const int N   = (int)H.n_cols;         // n^2 p + n k
    const int nv0 = nvar - 1;              // 0-based normalization variable
    const int e   = n2p + nv0;             // its position in [vec(A); Gamma]
    const double* he = H.colptr(e);

    // WHat(e, e): asymptotic variance of Gammahat_nvar
    const double W_ee = se_newey_west_lrv_cpp(he, T1, NWlags) / T1;

    // Waldstat for instrument relevance (HAC robust); reported as F as well
    double Gamma_nvar = Gamma(nv0, 0);
    double Waldstat = ((double)T1 * Gamma_nvar * Gamma_nvar) / W_ee;
    double Fstat = Waldstat;

    // MA representation: C_vec[h] = C_h for h=0..hor (Wold multipliers;
    // AL' is the lag block of beta without intercept)
    arma::cube wold;
    fWoldIRF_into_cpp(AL.t(), 0, p, hor, wold);
    std::vector<arma::mat> C_vec(hor + 1);
    for (int h = 0; h <= hor; ++h) C_vec[h] = wold.slice(h);

    // Stack every (j, h) pair, row r = h*n + j. With
    //   sf_r = scale * [kron(Gamma', e_j') G_h,  e_j' C_h]   (1 x N)
    // the delta-method gradient is d_r = sf_r - lambda_jh * e', where e is
    // the unit vector of Gammahat_nvar. D holds the rows d_r.
    const int R = n * (hor + 1);
    arma::mat D(R, N, arma::fill::zeros);
    D.cols(0, n2p - 1) = scale * irf_score_rows(AL, Gamma.col(0), C_vec, p, hor, n);
    arma::vec CG(R);                        // e_j' C_h Gamma
    for (int h = 0; h <= hor; ++h) {
        D.submat(h * n, n2p, h * n + n - 1, n2p + n - 1) = scale * C_vec[h];
        CG.subvec(h * n, h * n + n - 1) = C_vec[h] * Gamma.col(0);
    }
    const arma::vec lambda = scale * CG / Gamma_nvar;   // plug-in IRFs
    D.col(e) -= lambda;

    // Projected scores, one series per (j, h): column r of HD = H d_r', and
    // H sf_r' = HD.col(r) + lambda_jh * H.col(e). Then
    //   d_r WHat d_r'  = lrv(H d_r') / T1
    //   sf_r WHat sf_r' = lrv(H sf_r') / T1,  sf_r WHat e = lrcov(H sf_r', H e) / T1
    const arma::mat HD = H * D.t();         // T1 x R
    arma::vec hsf(T1);

    // Allocate output matrices: n x (hor+1)
    arma::mat IRF(n, hor + 1, arma::fill::zeros);
    arma::mat IRFstderr(n, hor + 1, arma::fill::zeros);
    arma::mat Dmlbound(n, hor + 1, arma::fill::zeros);
    arma::mat Dmubound(n, hor + 1, arma::fill::zeros);
    arma::mat MSWlbound(n, hor + 1, arma::fill::zeros);
    arma::mat MSWubound(n, hor + 1, arma::fill::zeros);

    // ahat is the same for all (j, h)
    double ahat = (double)T1 * Gamma_nvar * Gamma_nvar - critval * W_ee;

    const double sqrt_T1 = std::sqrt((double)T1);
    const double abs_Gnv = std::abs(Gamma_nvar);

    for (int h = 0; h <= hor; ++h) {
        for (int j = 0; j < n; ++j) {
            const int r = h * n + j;
            const double lambda_jh = lambda(r);
            IRF(j, h) = lambda_jh;

            // Delta-method variance d' WHat d
            const double* hd = HD.colptr(r);
            const double DmVar = se_newey_west_lrv_cpp(hd, T1, NWlags) / T1;
            IRFstderr(j, h)   = std::sqrt(DmVar) / (sqrt_T1 * abs_Gnv);
            double half_width = std::sqrt(critval / T1) * std::sqrt(DmVar) / abs_Gnv;
            Dmlbound(j, h)    = lambda_jh - half_width;
            Dmubound(j, h)    = lambda_jh + half_width;

            // MSW (Anderson-Rubin) bounds.  With f = sf_r / scale,
            //   b2 + b3      = 2 critval scale   (f WHat)_e = 2 critval sf_r WHat e
            //   c2 + c3 + c4 =  critval scale^2  f WHat f'  = critval sf_r WHat sf_r'
            for (int t = 0; t < T1; ++t) hsf[t] = hd[t] + lambda_jh * he[t];
            double v_sf, v_e, v_sfe;
            se_newey_west_lrv2_cpp(hsf.memptr(), he, T1, NWlags, v_sf, v_e, v_sfe);
            const double ej_Ch_Gamma = CG(r);
            double bhat = -2.0 * T1 * scale * ej_Ch_Gamma * Gamma_nvar
                          + 2.0 * critval * v_sfe / T1;
            double chat = std::pow(sqrt_T1 * scale * ej_Ch_Gamma, 2.0)
                          - critval * v_sf / T1;

            double Delta = bhat * bhat - 4.0 * ahat * chat;

            double lb, ub;
            if (ahat > 0 && Delta > 0) {
                lb = (-bhat - std::sqrt(Delta)) / (2.0 * ahat);
                ub = (-bhat + std::sqrt(Delta)) / (2.0 * ahat);
            } else if (ahat < 0 && Delta > 0) {
                lb = (-bhat + std::sqrt(Delta)) / (2.0 * ahat);
                ub = (-bhat - std::sqrt(Delta)) / (2.0 * ahat);
            } else if (ahat > 0 && Delta <= 0) {
                lb = arma::datum::nan;
                ub = arma::datum::nan;
            } else {
                lb = -arma::datum::inf;
                ub =  arma::datum::inf;
            }
            MSWlbound(j, h) = lb;
            MSWubound(j, h) = ub;
        }
    }

    // Force normalization variable at h=0 to exactly scale
    MSWlbound(nv0, 0) = scale;
    MSWubound(nv0, 0) = scale;

    MSWResult result;
    result.IRF       = IRF;
    result.IRFstderr = IRFstderr;
    result.Dmlbound  = Dmlbound;
    result.Dmubound  = Dmubound;
    result.MSWlbound = MSWlbound;
    result.MSWubound = MSWubound;
    result.Waldstat  = Waldstat;
    result.Fstat     = Fstat;
    return result;
}

// =========================================================================
// R wrapper
// =========================================================================

//' Weak-IV Robust Impulse Response Inference (Montiel-Olea, Stock, Watson 2021)
//'
//' Computes both plug-in (delta method) and Anderson-Rubin weak-IV robust
//' confidence sets for IV-identified structural IRFs following
//' Montiel-Olea, Stock and Watson (2021).
//'
//' @param var_result List from \code{fVAR}.
//' @param Z One-column instrument matrix (proxy sample, no missing values).
//' @param finaldata T x N matrix of original endogenous variables (full sample).
//' @param adjustu Integer vector \code{c(start, end)} (1-based) selecting the
//'   proxy-sample rows of the VAR residuals.
//' @param hor Forecast horizon. Default 48.
//' @param nvar 1-based index of the normalization variable. Default 1.
//' @param scale Scale of the structural shock. Default 1.
//' @param confidence Confidence level (0-1). Default 0.9.
//' @param NWlags Number of Newey-West lags for HAC correction: 0 = White
//'   heteroskedasticity-consistent only (no autocorrelation), k > 0 = Newey-West
//'   with k lags. Must be >= 0. Default 0.
//'
//' @return A list with elements:
//'   \item{IRF}{n x (hor+1) plug-in IRF point estimates.}
//'   \item{IRFstderr}{n x (hor+1) delta-method standard errors.}
//'   \item{Dmlbound}{n x (hor+1) delta-method lower confidence bounds.}
//'   \item{Dmubound}{n x (hor+1) delta-method upper confidence bounds.}
//'   \item{MSWlbound}{n x (hor+1) Anderson-Rubin lower bounds.}
//'   \item{MSWubound}{n x (hor+1) Anderson-Rubin upper bounds.}
//'   \item{Waldstat}{Wald statistic for instrument relevance (HAC-robust F-stat).}
//'   \item{Fstat}{Same as Waldstat (reported for compatibility).}
//'
//' @examples
//' set.seed(123)
//' y          <- matrix(rnorm(400), ncol = 2)
//' var_result <- fVAR(y, p = 2, c = 1)
//'
//' u    <- var_result$residuals
//' Z    <- matrix(u[, 1] + rnorm(nrow(u), sd = 0.5), ncol = 1)
//' rows <- c(1, nrow(u))
//'
//' msw <- fMSW(var_result, Z, y, rows, hor = 12, nvar = 1, NWlags = 2)
//' msw$Fstat
//' dim(msw$MSWlbound)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fMSW(const Rcpp::List& var_result,
                const arma::mat& Z,
                const arma::mat& finaldata,
                const arma::ivec& adjustu,
                int hor        = 48,
                int nvar       = 1,
                double scale   = 1.0,
                double confidence = 0.9,
                int NWlags     = 0) {

    if (NWlags < 0)
        Rcpp::stop("NWlags must be >= 0. Use 0 for White HC-only, or a positive integer for Newey-West with that many lags.");
    if (Z.n_cols != 1)
        Rcpp::stop("fMSW supports a single instrument: Z must have exactly one column.");

    // Extract VAR components
    VARResult vr = fVARResultFromList_cpp(var_result);
    fRequireNoExog_cpp(vr.n_exog, "fMSW");
    const arma::mat& beta      = vr.beta;
    const arma::mat& residuals = vr.residuals;
    const int p = vr.p;
    const int c = vr.c;

    int T = (int)finaldata.n_rows;

    // Proxy-sample residuals (0-based)
    int u_start = adjustu(0) - 1;
    int u_end   = adjustu(1) - 1;
    arma::mat u = residuals.rows(u_start, u_end);  // T1 x N
    int T1 = (int)u.n_rows;

    // AL: n x n*p (strip constant row from beta, then transpose)
    // beta is (c + N*p) x N; first c rows are deterministic
    arma::mat AL = beta.rows(c, (int)beta.n_rows - 1).t();  // N x N*p

    // Sigma_sub = cov(u) with denominator T1-1 (standard sample covariance)
    arma::mat u_dm = u.each_row() - arma::mean(u, 0);
    arma::mat Sigma_sub = (u_dm.t() * u_dm) / (T1 - 1);

    // eta = u' (N x T1)
    arma::mat eta = u.t();

    // Build full design matrix for VAR sample, then restrict to proxy sample
    // X_full: (T-p) x (c + N*p)
    arma::mat lags_full = fLagMakerMatrix(finaldata, p);   // (T-p) x N*p
    arma::mat X_full;
    if (c == 1) {
        X_full = arma::join_rows(arma::ones<arma::mat>(T - p, 1), lags_full);
    } else {
        X_full = lags_full;
    }
    // Restrict to proxy sample (same rows as residuals[u_start:u_end])
    arma::mat X = X_full.rows(u_start, u_end);  // T1 x (c + N*p)

    MSWResult res = fMSW_cpp(AL, Sigma_sub, eta, X, Z,
                             p, hor, nvar, scale, confidence, NWlags);

    return Rcpp::List::create(
        Rcpp::Named("IRF")       = res.IRF,
        Rcpp::Named("IRFstderr") = res.IRFstderr,
        Rcpp::Named("Dmlbound")  = res.Dmlbound,
        Rcpp::Named("Dmubound")  = res.Dmubound,
        Rcpp::Named("MSWlbound") = res.MSWlbound,
        Rcpp::Named("MSWubound") = res.MSWubound,
        Rcpp::Named("Waldstat")  = res.Waldstat,
        Rcpp::Named("Fstat")     = res.Fstat
    );
}
