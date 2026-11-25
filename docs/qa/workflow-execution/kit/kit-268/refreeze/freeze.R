# KIT-268 O-3 (a) (VAL-314): re-freeze the GUI-entry oracle of the lane G0 listed (lipid_canonical) on the runner, at the fix tree.
# MULTISCHOLAR_GUI_ENTRY_FREEZE=1 makes guiEntryExpectLane rewrite tests/testdata/gui-entry/lipid_canonical.json (header
# included); the file is copied to kit268-freeze/ for the review of its diff. Run from the package root.
Sys.setenv(MULTISCHOLAR_GUI_ENTRY_FREEZE = "1")
pkg <- normalizePath(".")
out <- file.path(pkg, "kit268-freeze")
dir.create(out, showWarnings = FALSE)
suppressMessages(pkgload::load_all(pkg, quiet = TRUE, export_all = TRUE, helpers = TRUE))
testthat::test_file(file.path(pkg, "tests/testthat/test-gui-entry-lipidomics.R"), reporter = testthat::SummaryReporter$new(),
  stop_on_failure = FALSE)
file.copy(file.path(pkg, "tests/testdata/gui-entry/lipid_canonical.json"), out, overwrite = TRUE)
cat("copied:", file.exists(file.path(out, "lipid_canonical.json")), "\n")
