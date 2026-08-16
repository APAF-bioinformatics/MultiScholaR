library(testthat)
library(MultiScholaR)

# Repository-wide characterization tests run in isolated processes via the
# contract harness. R CMD check runs the installed-package contract smoke tests.
test_check("MultiScholaR", filter = "^quality-runtime-contracts$")
