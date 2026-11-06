# KIT-101 W3: the entry-point behaviour suite, proteomics. `mod_proteomics_server` under testServer from import to the
# summary on each committed proteomics fixture (helper-gui-entry.R), compared step by step with its frozen oracle
# (tests/testdata/gui-entry/<lane>.json). No skip: a lane that cannot run fails.

test_that("KIT-101 entry lane prot_dia: DIA-NN, the DIA route, with enrichment through the doubles", {
  guiEntryExpectLane("prot_dia")
})

test_that("KIT-101 entry lane prot_dia_limpa: DIA-NN, the DIA_LIMPA route", {
  guiEntryExpectLane("prot_dia_limpa")
})

test_that("KIT-101 entry lane prot_tmt: Proteome Discoverer TMT, the protein route", {
  guiEntryExpectLane("prot_tmt")
})

test_that("KIT-101 entry lane prot_lfq: MaxQuant LFQ, the protein route", {
  guiEntryExpectLane("prot_lfq")
})

test_that("KIT-101 entry lane prot_lfq_fragpipe: FragPipe LFQ, the protein route", {
  guiEntryExpectLane("prot_lfq_fragpipe")
})
