# KIT-275 r2.1 (audit Q1 (b), pilot read): the cross-implementation identity check. The pilot built its V1a values with its own `pilotPlant()` (pilot/pilot-k275.R);
# the counted run builds them inside `generate()` (ticket 5.3's "add" branch). For r = 9601 (a pilot data set) at each design, the counted `makeVersion(design, r, "V1a",
# "none", v0)` must equal the pilot's values EXACTLY (`identical()` on y, on the values passed to the pipeline and on planted_delta): the pilot's control IS the
# counted control. A mismatch stops (exit 1). Run from the package root (it sources sim-k275.R as a library, `options(k275.library = TRUE)`). It runs no journey.
# Usage: Rscript crosscheck-k275.R          Output: one line per design and the verdict.
options(warnPartialMatchDollar = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])))
options(k275.library = TRUE)
source(file.path(here, "sim-k275.R"), local = FALSE)
# take the pilot's `PILOT_S` and `pilotPlant` from its file, unchanged (the pilot script itself is not sourced: it would run a job)
exprs <- parse(file = file.path(here, "pilot", "pilot-k275.R"))
for (e in exprs) {
    if (is.call(e) && identical(e[[1L]], as.name("<-")) && as.character(e[[2L]]) %in% c("PILOT_S", "pilotPlant")) eval(e, envir = globalenv())
}
stopifnot(exists("pilotPlant"), exists("PILOT_S"), identical(unname(PILOT_S[["V1a"]]), 4))
r <- 9601L
bad <- 0L
for (design in c("D3", "D4", "D6")) {
    p <- CELLS[[design]]
    a_run <- sprintf("R%03d", p[["groups"]] * p[["per"]])
    v0 <- makeVersion(design, r, "V0", "none")
    v1 <- makeVersion(design, r, "V1a", "none", v0)
    plant <- pilotPlant(design, r, v0[["sim"]], a_run)
    pa <- plant[["versions"]][["V1a"]]
    pilot_values <- maskedLog2(pa, "none", v0[["sim"]][["y"]][v0[["sim"]][["feat_ids"]], , drop = FALSE])
    same_y <- identical(v1[["sim"]][["y"]], pa[["y"]])
    same_values <- identical(v1[["values"]], pilot_values)
    same_delta <- identical(v1[["sim"]][["planted_delta"]], pa[["planted_delta"]])
    cat(sprintf("crosscheck %s r=%d: y identical %s; values identical %s; planted_delta identical %s\n", design, r, same_y, same_values, same_delta))
    if (!(same_y && same_values && same_delta)) bad <- bad + 1L
}
cat(sprintf("crosscheck: %s\n", if (bad) sprintf("MISMATCH in %d design(s): the pilot's control is not the counted control; STOP", bad) else "the counted V1a equals the pilot's, exactly, at every design"))
quit(save = "no", status = if (bad) 1L else 0L)
