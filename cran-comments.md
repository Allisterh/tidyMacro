## Resubmission

This replaces the pending submission of version 1.0.0, which I have asked to
withdraw. Version 1.1.0 fixes the issues found since then: console output is
now suppressible, bootstrap draws that fail no longer abort the R session,
`fVAR()` validates its arguments, and the DESCRIPTION describes the
implemented methods with references.

## Test environments

* local macOS, R 4.5.x
* Ubuntu 22.04, R 4.3.3 (`R CMD check --as-cran` with
  `_R_CHECK_LIMIT_CORES_=TRUE` and `OMP_THREAD_LIMIT=2`)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking installed package size ... NOTE
    installed size is 46.9Mb
    sub-directories of 1Mb or more:
      libs  46.1Mb

  The package is compiled against RcppArmadillo with OpenMP, and the shared
  object carries debug symbols on the check machine. There are no data or
  documentation files of any size; the libs directory is the compiled code
  alone.

## Parallelism

All functions that use OpenMP take `n_threads`, and the examples and tests
use at most two threads. Random draws are generated on the main thread before
any parallel region, so results do not depend on the thread count; this is
covered by tests.
