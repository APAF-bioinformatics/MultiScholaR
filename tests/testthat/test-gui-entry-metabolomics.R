# KIT-101 W3: the entry-point behaviour suite, metabolomics. `mod_metabolomics_server` under testServer from import to
# the summary on each committed metabolomics fixture (helper-gui-entry.R), compared step by step with its frozen oracle
# (tests/testdata/gui-entry/<lane>.json). No skip: a lane that cannot run fails.

test_that("KIT-101 entry lane metab_lc: LC-MS positive and negative assays (custom format)", {
  guiEntryExpectLane("metab_lc")
})

test_that("KIT-101 entry lane metab_gc: one GC-MS assay (custom format)", {
  guiEntryExpectLane("metab_gc")
})

test_that("KIT-101 entry lane metab_combined: LC-MS and GC-MS assays (custom format)", {
  guiEntryExpectLane("metab_combined")
})
