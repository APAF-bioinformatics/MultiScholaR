# KIT-101 W3: the entry-point behaviour suite, transcriptomics. `mod_transcriptomics_server` under testServer from
# import to the summary on the committed trx_screen fixture (helper-gui-entry.R), compared step by step with its frozen
# oracle (tests/testdata/gui-entry/<lane>.json). No skip: a lane that cannot run fails.

test_that("KIT-101 entry lane trx_screen: the KIT-040 count corpus, 12 samples, 300 genes", {
  guiEntryExpectLane("trx_screen")
})
