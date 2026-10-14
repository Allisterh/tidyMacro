#ifndef OMP_THREADS_H
#define OMP_THREADS_H

#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>

// Package-wide OpenMP rule: n_threads <= 0 uses all cores minus one (at
// least one). Callers pass the result to a per-loop num_threads() clause and
// never call omp_set_num_threads(), so the global OpenMP setting (and hence
// omp_get_max_threads() in later calls) is left untouched.
inline int fResolveThreads_cpp(int n_threads) {
#ifdef _OPENMP
  return n_threads > 0 ? n_threads : std::max(1, omp_get_max_threads() - 1);
#else
  (void) n_threads;
  return 1;
#endif
}

#endif // OMP_THREADS_H
