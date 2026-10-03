# This file is part of the standard testthat setup for CRAN.
# CRAN allows at most 2 cores in tests/examples.
Sys.setenv(OMP_THREAD_LIMIT = "2")

library(testthat)
library(tidyMacro)

test_check("tidyMacro")
