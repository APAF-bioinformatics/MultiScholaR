# KIT-101 W3: the entry-point behaviour suite, lipidomics. `mod_lipidomics_server` under testServer from import to the
# summary on the committed lipid fixture (helper-gui-entry.R), compared step by step with its frozen oracle
# (tests/testdata/gui-entry/<lane>.json). No skip: a lane that cannot run fails.

test_that("KIT-101 entry lane lipid_canonical: two LipidSearch assays", {
  guiEntryExpectLane("lipid_canonical")
})
