#include "fVAR.h"
#include "fLagMakerMatrix.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// Shared implementation: takes either an empty mat (no exog) or a T x M exog mat
// Kept static so the public _cpp overloads stay as the only external entry points.
static VARResult fVAR_cpp_impl(const arma::mat &y, int p, int c,
                               const arma::mat *exog_ptr) {
  int T = y.n_rows;

  arma::mat yfinal = y.rows(p, T - 1);
  arma::mat lags = fLagMakerMatrix(y, p);

  arma::mat X;
  if (c == 1) {
    X = arma::join_rows(arma::ones<arma::mat>(T - p, 1), lags);
  } else {
    X = lags;
  }

  int n_exog = 0;
  if (exog_ptr != nullptr) {
    const arma::mat &exog_mat = *exog_ptr;
    arma::mat exog_eff = exog_mat.rows(p, T - 1);
    X = arma::join_rows(X, exog_eff);
    n_exog = exog_mat.n_cols;
  }

  arma::mat XtX = X.t() * X;
  arma::mat Xty = X.t() * yfinal;

  // A rank-deficient design must not throw: this runs inside OpenMP workers.
  // Fall back to a pseudo-inverse and flag the draw instead.
  arma::mat beta;
  bool ok = arma::solve(beta, XtX, Xty, arma::solve_opts::fast);
  if (!ok) {
    arma::mat XtX_pinv;
    if (arma::pinv(XtX_pinv, XtX)) beta = XtX_pinv * Xty;
    else                           beta.zeros(X.n_cols, y.n_cols);
  }

  arma::mat residuals = yfinal - X * beta;
  int n_obs = residuals.n_rows;   // T_full - p (residual sample)
  int N     = y.n_cols;
  const int df = n_obs - c - N * p - n_exog;   // residual degrees of freedom
  arma::mat sigma = (df > 0)
      ? arma::mat(residuals.t() * residuals / df)
      : arma::mat(N, N, arma::fill::value(arma::datum::nan));
  if (df <= 0) ok = false;

  VARResult result;
  result.beta = beta;
  result.residuals = residuals;
  result.sigma = sigma;
  result.p = p;
  result.c = c;
  result.n_exog = n_exog;
  result.ok = ok;
  return result;
}

void fCheckVARArgs_cpp(const arma::mat &y, int p, int c,
                       const arma::mat *exog_ptr) {
  if (y.n_rows == 0 || y.n_cols == 0)
    Rcpp::stop("fVAR: 'y' must be a non-empty T x N matrix.");
  if (p < 1)
    Rcpp::stop("fVAR: 'p' must be at least 1, got %d.", p);
  if (c != 0 && c != 1)
    Rcpp::stop("fVAR: 'c' must be 0 (no intercept) or 1 (intercept), got %d.", c);
  if (!y.is_finite())
    Rcpp::stop("fVAR: 'y' contains NA, NaN or Inf values.");

  const int T = static_cast<int>(y.n_rows);
  const int N = static_cast<int>(y.n_cols);
  if (T <= p)
    Rcpp::stop("fVAR: 'y' has %d rows, which is not more than p = %d.", T, p);

  int n_exog = 0;
  if (exog_ptr != nullptr) {
    if (static_cast<int>(exog_ptr->n_rows) != T)
      Rcpp::stop("fVAR: 'exog' has %d rows but 'y' has %d; they must match.",
                 static_cast<int>(exog_ptr->n_rows), T);
    if (!exog_ptr->is_finite())
      Rcpp::stop("fVAR: 'exog' contains NA, NaN or Inf values.");
    n_exog = static_cast<int>(exog_ptr->n_cols);
  }

  const int n_obs = T - p;                 // rows of the regression
  const int k     = c + N * p + n_exog;    // columns of the design matrix
  if (n_obs <= k)
    Rcpp::stop("fVAR: %d usable observations for %d regressors "
               "(c + N*p + n_exog); no residual degrees of freedom left. "
               "Reduce 'p' or use a longer sample.", n_obs, k);
}

// No-exog overload for C++ callers, including OpenMP worker threads.
VARResult fVAR_cpp(const arma::mat &y, int p, int c) {
  return fVAR_cpp_impl(y, p, c, nullptr);
}

// Internal C++ function (called from other C++ code) — Nullable overload
VARResult fVAR_cpp(const arma::mat &y, int p, int c,
                   Rcpp::Nullable<arma::mat> exog) {
  if (exog.isNotNull()) {
    arma::mat exog_mat = Rcpp::as<arma::mat>(exog);
    return fVAR_cpp_impl(y, p, c, &exog_mat);
  }
  return fVAR_cpp_impl(y, p, c, nullptr);
}

// Pre-converted exog overload: avoids Rcpp::as in hot loops
VARResult fVAR_cpp(const arma::mat &y, int p, int c,
                   const arma::mat &exog) {
  return fVAR_cpp_impl(y, p, c, &exog);
}

VARResult fVAR_cpp(const arma::mat &y, int p, int c,
                   const arma::mat *exog_ptr) {
  return fVAR_cpp_impl(y, p, c, exog_ptr);
}

VARResult fVARResultFromList_cpp(const Rcpp::List &var_result) {
  VARResult vr;
  vr.beta      = Rcpp::as<arma::mat>(var_result["beta"]);
  vr.residuals = Rcpp::as<arma::mat>(var_result["residuals"]);
  vr.sigma     = var_result.containsElementNamed("sigma")
                     ? Rcpp::as<arma::mat>(var_result["sigma"])
                     : arma::mat();
  vr.p         = Rcpp::as<int>(var_result["p"]);
  vr.c         = Rcpp::as<int>(var_result["c"]);
  vr.n_exog    = var_result.containsElementNamed("n_exog")
                     ? Rcpp::as<int>(var_result["n_exog"])
                     : 0;
  return vr;
}

void fCheckExogArg_cpp(int n_exog, bool has_exog) {
  if (n_exog > 0 && !has_exog)
    Rcpp::stop("Original VAR used exogenous variables. You must provide the 'exog' parameter.");
  if (n_exog == 0 && has_exog)
    Rcpp::stop("Original VAR did not use exogenous variables. Do not provide the 'exog' parameter.");
}

void fRequireNoExog_cpp(int n_exog, const char *fun) {
  if (n_exog > 0)
    Rcpp::stop("%s does not support VARs with exogenous variables.", fun);
}


//' Vector Autoregression (VAR) Model Estimation
//'
//' Estimates a Vector Autoregression model with optional exogenous variables
//' using ordinary least squares (OLS). Optimised for maximum speed and memory
//' efficiency via RcppArmadillo.
//'
//' @param y A numeric matrix of time series data with T rows (time periods)
//'   and N columns (endogenous variables).
//' @param p An integer specifying the number of lags to include in the VAR.
//' @param c An integer indicator for including a constant term (1 = include
//'   intercept, 0 = no intercept).
//' @param exog An optional T x M matrix of exogenous variables. If provided,
//'   these variables enter contemporaneously (not lagged). Default is NULL.
//'
//' @return A list with elements: beta — coefficient matrix of dimensions
//'   (Np + c + M) x N, where the first row is the intercept (if c = 1),
//'   followed by N*p lag-coefficient rows, then M exogenous-coefficient rows;
//'   residuals — (T-p) x N matrix of OLS residuals;
//'   sigma — N x N residual covariance matrix with degrees-of-freedom
//'   correction: denominator is (T - p) - c - N*p - n_exog, matching
//'   Bloom (2009) and Kaenzig (2021);
//'   p — lag order (echoed from input);
//'   c — intercept indicator (echoed from input);
//'   n_exog — number of exogenous variables (0 if none provided).
//'
//' @details
//' The model is \deqn{Y_t = c + A_1 Y_{t-1} + \cdots + A_p Y_{t-p} + B X_t + e_t}
//' where \eqn{X_t} are contemporaneous exogenous variables. Without exogenous
//' variables the model reduces to a standard VAR(p). Estimation uses
//' \code{arma::solve} for numerical stability.
//'
//' @examples
//' set.seed(123)
//' y      <- matrix(rnorm(400), ncol = 2)
//' result <- fVAR(y, p = 2, c = 1)
//' result$beta
//'
//' # With one exogenous regressor
//' exog   <- matrix(rnorm(200), ncol = 1)
//' result <- fVAR(y, p = 2, c = 1, exog = exog)
//' dim(result$residuals)
//'
//' @seealso \code{\link{fWoldIRF}}, \code{\link{fBootstrapChol}},
//'   \code{\link{fAICBIC}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fVAR(const arma::mat &y, int p, int c,
                Rcpp::Nullable<arma::mat> exog = R_NilValue) {
  arma::mat exog_mat;
  if (exog.isNotNull()) exog_mat = Rcpp::as<arma::mat>(exog);
  const arma::mat *exog_ptr = exog.isNotNull() ? &exog_mat : nullptr;
  fCheckVARArgs_cpp(y, p, c, exog_ptr);

  VARResult result = fVAR_cpp(y, p, c, exog_ptr);
  if (!result.ok)
    Rcpp::stop("fVAR: the regressor matrix is rank deficient (collinear "
               "variables, or a constant series with c = 1).");

  return Rcpp::List::create(Rcpp::Named("beta") = result.beta,
                            Rcpp::Named("residuals") = result.residuals,
                            Rcpp::Named("sigma") = result.sigma,
                            Rcpp::Named("p") = result.p,
                            Rcpp::Named("c") = result.c,
                            Rcpp::Named("n_exog") = result.n_exog);
}
