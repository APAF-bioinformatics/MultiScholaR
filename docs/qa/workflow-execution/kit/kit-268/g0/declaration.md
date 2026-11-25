# KIT-268 G0: declaration (written before the run)

Work item (AGENTS.md): KIT-268 r1.2 (owner decision VAL-312), gate G0 only, as the ticket's section 10 defines it (RC07).
Outcome: an observation, with no expectation, of which comparisons fall back to Benjamini-Hochberg per lane and assay,
which stored `da` steps of the four GUI-entry oracles would gain a notification, and whether the keys of
`qvalue_warnings` map to `contrasts` entries (the mapping observation). Exit: the tables are returned and recorded.
No product code and no test expectation is written in this work; owner questions O-1, O-2, O-3 are open.

- Purpose: observe, on BASE only, through `trace(completeOmicDaRun)`, the tab's `results$qvalue_warnings` (assay, then comparison
  keys) and `workflowData$contrasts_tbl[c("contrasts","friendly_names")]` at every call, with the test name, lane and step.
- Base commit (bound): 4d588ba3ed5c6fee76a5a1450c028b409ba419ba (origin/v0.9-staging, 2026-11-26). The request's
  `package_commit` is the head of branch kit-268 (base plus only the files under `docs/qa/workflow-execution/kit/kit-268/g0/`;
  nothing under `R/`, `inst/`, `tests/`); `base_commit` is "none" (an observation, no comparison).
- Files observed, each in its own process, in this order: g0-control.R (the positive control), then
  `tests/testthat/test-gui-entry-lipidomics.R`, `test-gui-entry-metabolomics.R` (lanes lipid_canonical, metab_lc, metab_gc,
  metab_combined), `test-omic-summary-module.R`, `test-omic-da-module.R` (the W1 captures; the fourth file the ticket names).
- Positive control first: the `tied_pvalues` case of `helper-module-ci-lipid-da.R` (kernel `runLipidsDA` with the helper's stub,
  non-empty `qvalue_warnings`) through `completeOmicDaRun`. If it prints no fallback comparison the instrument is blind: that
  script exits 2, every later script exits 3 at once without running its file, and nothing else is read. A blind instrument is
  an instrument defect (VAL-54 rule, at most three re-sends), not a "no fallback" result.
- Recorded (per process, TSV, verbatim): `<mode>-calls.tsv` (one row per call: test, lane, step, number of assays and keys,
  failed contrasts, the contrasts table's `contrasts` and `friendly_names`), `<mode>-keys.tsv` (one row per assay and key, with
  the mapping by the ticket's rule: exact match on `contrasts`, else the text before the first "=", else the key),
  `<mode>-tests.tsv` (the testthat counts, context only), `console.log`.
- Instrument-liveness check (post-hoc, in outcome): each GUI-entry lane's stored `da` step holds the completion notification,
  so each lane must show at least one traced call at step `da`; a lane without one means the instrument missed that route.
- Stop rules: (1) control empty: stop, do not read the rest; (2) a script timeout or crash is recorded as such and not re-run
  unless it is an instrument defect under VAL-54 (four conditions, at most three); (3) no change to any threshold or file list
  to obtain a result; (4) no stored file is touched (freeze variable is not set).
- Ceiling 45 job-minutes (expected 30). Script timeouts (minutes): control 4, lipid GUI 5, metab GUI 8, da-module 15,
  summary 7 (sum 39). Runner: ubuntu-22.04, reference BLAS, extra packages sva, limpa, mixOmics (as checkpoint-20261126, whose
  per-file times are the basis).
- Decision it feeds: the scope and predictions of R1/R2, the O-3 table, and stop rule 3 of the ticket (mapping).
