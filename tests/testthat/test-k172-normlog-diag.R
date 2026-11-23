# KIT-172 diagnosis (VAL-68; run branch run-k172-win only, never merged): the four lipidomics and metabolomics lanes'
# `norm_log` text on the bundle, to read why its digest differs from the re-frozen oracle on Windows (bundle run
# 45151293709). The lane capture holds a log longer than 1000 characters as its digest only, so the capture function
# of the helpers is wrapped to keep any text holding "Applying log2 transformation" whole; the text at step `norm` is
# returned in failure messages of 500 characters (SUITEFAIL lines), with the digest the capture would have written and
# the oracle's. Diagnosis only: no package code, oracle or ticket test is changed.
library(testthat)

test_that("KIT-172 diagnosis: norm_log text per lane at step norm", {
  helpers <- environment(guiEntryRunLane)
  original <- get("guiEntryOutputValue", envir = helpers)
  assign("guiEntryOutputValue", function(value, roots) {
    if (is.character(value) && any(grepl("Applying log2 transformation", value, fixed = TRUE))) {
      return(paste("<fulltext>", guiEntryMask(paste(value, collapse = "\n"), roots)))
    }
    original(value, roots)
  }, envir = helpers)
  withr::defer(assign("guiEntryOutputValue", original, envir = helpers))
  # Mode `testthat`: as the suite runs it (C collation). Mode `system`: two lanes again under the runner's own locale
  # (LC_ALL ""), as the bundle's lane step runs the lanes outside testthat.
  # Rerun (VAL-68, the one mechanics rerun): the first run returned the text at step `norm`, where it equals the
  # oracle; the differing steps are `session_export` on, where the log passes 1000 characters and is held as a digest.
  runs <- list(list(mode = "testthat", lanes = c("lipid_canonical", "metab_gc")),
    list(mode = "system", lanes = "lipid_canonical"))
  for (run in runs) for (lane in run$lanes) {
    if (identical(run$mode, "system")) {
      saved <- Sys.getlocale()
      suppressWarnings(Sys.setlocale("LC_ALL", ""))
      record <- guiEntryRunLane(lane)
      locale <- Sys.getlocale()
      suppressWarnings(Sys.setlocale("LC_ALL", saved))
    } else {
      record <- guiEntryRunLane(lane)
      locale <- Sys.getlocale()
    }
    lane_label <- paste0(lane, "[", run$mode, "]")
    omic <- if (startsWith(lane, "lipid")) "lipidomics" else "metabolomics"
    key <- paste0(omic, "_workflow-norm-norm_log")
    oracle <- jsonlite::read_json(guiEntryOraclePath(lane), simplifyVector = FALSE)$steps$session_export$outputs[[key]]
    text <- sub("^<fulltext> ", "", record$steps$session_export$outputs[[key]])
    computed <- paste("<text>", digest::digest(text, algo = "sha256", serialize = FALSE))
    bytes <- charToRaw(text)
    fail(sprintf("K172DIAG %s step session_export chars %d bytes %d encoding %s nonascii %s computed %s oracle %s equal %s crlf %s locale %s",
      lane_label, nchar(text), length(bytes), paste(unique(Encoding(text)), collapse = ","),
      paste(as.character(bytes[bytes > as.raw(127)]), collapse = " "), computed, oracle, identical(computed, oracle),
      grepl("\r", text, fixed = TRUE), substr(locale, 1L, 120L)))
    flat <- gsub("\r", "<CR>", gsub("\n", " | ", text, fixed = TRUE), fixed = TRUE)
    starts <- seq(1L, nchar(flat), by = 500L)
    for (i in seq_along(starts)) fail(sprintf("K172TEXT %s %02d %s", lane_label, i, substr(flat, starts[[i]], starts[[i]] + 499L)))
  }
})
