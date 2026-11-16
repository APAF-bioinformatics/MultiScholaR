# KIT-130: show that the network probe used by check-mofa.R observes what it claims to. Run once
# with the network up ("online") and once with it cut ("offline"); each run asserts its state.
#
# Usage: <bundle>/.../Rscript --vanilla net-probe.R <bundle directory> <online|offline>

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), "common.R"))
expected <- commandArgs(trailingOnly = TRUE)[2L]
reachable <- canReachNetwork()
cat(sprintf("network reachable from R: %s (expected: %s)\n", reachable, expected))
stopifnot(expected %in% c("online", "offline"), reachable == (expected == "online"))
checkPassed(sprintf("network probe sees '%s'", expected))
