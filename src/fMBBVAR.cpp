#include <RcppArmadillo.h>
#include "fMBBVAR.h"
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

// Tile a blockSize-row centering matrix down to final_size rows
static arma::mat make_centering_full(const arma::mat& centering,
                                     int nBlock, int actualBlockSize, int final_size) {
  arma::mat cf(final_size, centering.n_cols, arma::fill::none);
  for (int i = 0; i < nBlock; ++i) {
    int start_row = i * actualBlockSize;
    int end_row   = std::min(start_row + actualBlockSize - 1, final_size - 1);
    cf.rows(start_row, end_row) = centering.rows(0, end_row - start_row);
  }
  return cf;
}

MBBSamplerState fmbb_precompute_cpp(const arma::mat& eps, int lags, int BlockSize,
                                    const arma::mat* M) {
  const int Tp = eps.n_rows;
  const int N  = eps.n_cols;
  const int T  = Tp + lags;
  const bool has_M = (M != nullptr);
  const int k  = has_M ? static_cast<int>(M->n_cols) : 0;

  int actualBlockSize = BlockSize;
  if (actualBlockSize <= 0)
    actualBlockSize = static_cast<int>(std::round(5.03 * std::pow(Tp, 0.25)));

  const int nBlock        = static_cast<int>(std::ceil(static_cast<double>(T - lags) / actualBlockSize));
  const int nBlocks_avail = T - lags - actualBlockSize + 1;
  const int final_size    = T - lags;
  const int center_end    = T - lags - actualBlockSize;

  MBBSamplerState s;
  s.Blocks.set_size(actualBlockSize, N, nBlocks_avail);
  for (int j = 0; j < nBlocks_avail; ++j)
    s.Blocks.slice(j) = eps.rows(j, j + actualBlockSize - 1);

  arma::mat centering(actualBlockSize, N, arma::fill::none);
  for (int j = 0; j < actualBlockSize; ++j)
    centering.row(j) = arma::mean(eps.rows(j, j + center_end), 0);
  s.centering_full = make_centering_full(centering, nBlock, actualBlockSize, final_size);

  if (has_M) {
    s.MBlocks.set_size(actualBlockSize, k, nBlocks_avail);
    for (int j = 0; j < nBlocks_avail; ++j)
      s.MBlocks.slice(j) = M->rows(j, j + actualBlockSize - 1);

    // Instruments are centred on their non-zero values only, keeping the
    // sample mean of the non-zero values (Kaenzig 2021, runProxyVAR.m)
    arma::rowvec Mmean(k, arma::fill::zeros);
    for (int col = 0; col < k; ++col) {
      arma::vec col_vec = M->col(col);
      arma::uvec nz = arma::find(col_vec != 0.0);
      if (nz.n_elem > 0) Mmean(col) = arma::mean(col_vec.elem(nz));
    }
    arma::mat Mcentering(actualBlockSize, k, arma::fill::zeros);
    for (int j = 0; j < actualBlockSize; ++j) {
      arma::mat subM = M->rows(j, j + center_end);
      for (int col = 0; col < k; ++col) {
        arma::vec col_vec = subM.col(col);
        arma::uvec nz = arma::find(col_vec != 0.0);
        if (nz.n_elem > 0)
          Mcentering(j, col) = arma::mean(col_vec.elem(nz)) - Mmean(col);
      }
    }
    s.Mcentering_full = make_centering_full(Mcentering, nBlock, actualBlockSize, final_size);
  }

  s.nBlock          = nBlock;
  s.nBlocks_avail   = nBlocks_avail;
  s.final_size      = final_size;
  s.actualBlockSize = actualBlockSize;
  s.N               = N;
  s.k               = k;
  s.has_M           = has_M;
  return s;
}

arma::imat fmbb_draw_cpp(const MBBSamplerState& s, int nboot) {
  arma::imat index(s.nBlock, nboot, arma::fill::none);
  for (int b = 0; b < nboot; ++b)
    index.col(b) = arma::randi<arma::ivec>(s.nBlock,
                       arma::distr_param(0, s.nBlocks_avail - 1));
  return index;
}

MBBVARResult fmbb_sample_cpp(const MBBSamplerState& s, const arma::ivec& index) {
  arma::mat U_boot(s.final_size, s.N, arma::fill::none);
  arma::mat M_boot;
  if (s.has_M) M_boot.set_size(s.final_size, s.k);

  for (int j = 0; j < s.nBlock; ++j) {
    int start_row    = j * s.actualBlockSize;
    int end_row      = std::min(start_row + s.actualBlockSize - 1, s.final_size - 1);
    int rows_to_copy = end_row - start_row + 1;
    U_boot.rows(start_row, end_row) =
        s.Blocks.slice(index(j)).rows(0, rows_to_copy - 1);
    if (s.has_M)
      M_boot.rows(start_row, end_row) =
          s.MBlocks.slice(index(j)).rows(0, rows_to_copy - 1);
  }

  MBBVARResult result;
  result.eps_boot = U_boot - s.centering_full;

  if (s.has_M) {
    for (int col = 0; col < s.k; ++col)
      for (int i = 0; i < s.final_size; ++i)
        if (M_boot(i, col) != 0.0)
          M_boot(i, col) -= s.Mcentering_full(i, col);
  }
  result.M_boot = M_boot;
  return result;
}

//' Moving Block Bootstrap for VAR Residuals and Instruments
//' 
//' Computes a moving-block bootstrap sample for the reduced form errors
//' of a VAR(p) model. Optionally bootstraps a matrix of instruments M.
//' 
//' @param eps (T-p) x N matrix of VAR residuals
//' @param lags Integer lag order p of the VAR model
//' @param BlockSize Integer block size for moving block bootstrap. 
//'   If 0 or negative, automatically computed as round(5.03 * T^0.25)
//' @param M Optional (T-p) x K matrix of instruments to bootstrap (default = NULL)
//' 
//' @return A list containing:
//'   \itemize{
//'     \item eps_boot: (T-p) x N matrix of bootstrapped residuals (centered)
//'     \item M_boot: (T-p) x K matrix of bootstrapped instruments (centered), 
//'                   or empty matrix if M not provided
//'   }
//' 
//' @details
//' The moving block bootstrap preserves the temporal dependence structure
//' in the residuals. The default block size follows the recommendation from
//' "Proxy SVARs: Asymptotic Theory, Bootstrap Inference, and the Effects 
//' of Income Tax Changes in the United States".
//' 
//' Both residuals and instruments are centered to preserve their sample properties.
//' For instruments, only non-zero values are centered (to handle sparse
//' instruments), and the centering keeps the sample mean of the non-zero
//' values, as in Kaenzig (2021).
//' 
//' @examples
//' set.seed(123)
//' eps    <- matrix(rnorm(200), 100, 2)
//' result <- fMBBVAR(eps, lags = 2, BlockSize = 0)
//' dim(result$eps_boot)
//'
//' # Bootstrap the residuals and an instrument jointly
//' M      <- matrix(rnorm(100), 100, 1)
//' result <- fMBBVAR(eps, lags = 2, BlockSize = 0, M)
//' dim(result$M_boot)
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List fMBBVAR(const arma::mat& eps, 
                    int lags, 
                    int BlockSize, 
                    Rcpp::Nullable<arma::mat> M = R_NilValue) {
  arma::mat M_mat;
  if (M.isNotNull()) M_mat = Rcpp::as<arma::mat>(M);
  MBBSamplerState state = fmbb_precompute_cpp(eps, lags, BlockSize,
                                              M.isNotNull() ? &M_mat : nullptr);
  MBBVARResult result = fmbb_sample_cpp(state, fmbb_draw_cpp(state, 1).col(0));

  return Rcpp::List::create(
    Rcpp::Named("eps_boot") = result.eps_boot,
    Rcpp::Named("M_boot") = result.M_boot
  );
}
