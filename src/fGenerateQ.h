#ifndef FGENERATEQ_H
#define FGENERATEQ_H

#include <RcppArmadillo.h>
#include "rng_tidymacro.h"

// Fill Q with an N x N random orthonormal matrix using QR decomposition.
// R and G are caller-owned scratch buffers reused by hot loops.
void fGenerateQ_inplace(arma::mat& Q,
                        arma::mat& R,
                        arma::mat& G,
                        arma::uword N);

// Same distribution, driven by a caller-owned RNG.  Required inside OpenMP
// regions, where arma::randn (which reads R's global RNG state) must not be
// used.  Gram-Schmidt on Gaussian columns needs no scratch and avoids the
// LAPACK call overhead that dominates at the small N of a rotation search.
void fGenerateQ_inplace(arma::mat& Q,
                        arma::uword N,
                        tidymacro::RNG& rng);

// Generate an N x N random orthonormal matrix using QR decomposition.
// Sign normalisation: diagonal of R forced positive for a unique decomposition
// (Rubio-Ramirez, Waggoner & Zha 2010).
arma::mat fGenerateQ_cpp(int N);

// R wrapper
arma::mat fGenerateQ(int N);

#endif // FGENERATEQ_H
