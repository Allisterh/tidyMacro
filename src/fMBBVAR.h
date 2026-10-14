#ifndef FMBB_VAR_H
#define FMBB_VAR_H

#include <RcppArmadillo.h>

// One moving-block bootstrap draw
struct MBBVARResult {
  arma::mat eps_boot;  // (T-p) x N bootstrapped residuals (centered)
  arma::mat M_boot;    // (T-p) x K bootstrapped instruments (centered); empty without M
};

// Blocks and centering terms, built once before the bootstrap loop.
struct MBBSamplerState {
  arma::cube Blocks;           // (blockSize x N x nBlocks_avail)
  arma::cube MBlocks;          // (blockSize x K x nBlocks_avail), empty if !has_M
  arma::mat  centering_full;   // (final_size x N)
  arma::mat  Mcentering_full;  // (final_size x K), empty if !has_M
  int nBlock;
  int nBlocks_avail;
  int final_size;
  int actualBlockSize;
  int N;
  int k;
  bool has_M;
};

// Build the sampler state from the residuals and optional instruments
// (M = nullptr for none). BlockSize <= 0 uses round(5.03 * T^(1/4)).
MBBSamplerState fmbb_precompute_cpp(const arma::mat& eps, int lags, int BlockSize,
                                    const arma::mat* M);

// Block indices for nboot draws (nBlock x nboot). Uses R's RNG: call on the
// main thread, before any OpenMP loop.
arma::imat fmbb_draw_cpp(const MBBSamplerState& state, int nboot);

// One draw from pre-drawn block indices (thread-safe).
MBBVARResult fmbb_sample_cpp(const MBBSamplerState& state, const arma::ivec& index);

// R wrapper function (for calling from R)
Rcpp::List fMBBVAR(const arma::mat& eps,
                    int lags,
                    int BlockSize,
                    Rcpp::Nullable<arma::mat> M);

#endif // FMBB_VAR_H
