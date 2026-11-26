# KIT-276 step 0 (ticket section 13), read at origin/v0.9-staging 38443af8 (branch kit-276)

Reading only; no R, no run. `data/benchmark/` neither opened nor searched. All three dependencies are on staging
(KIT-216 = VAL-377, KIT-273 = VAL-374/375, KIT-266 = VAL-379/380).

## Result: STOPPED before any test or code, one expectation of the ticket cannot hold (S-1)

### S-1 (blocking): T18's "`saved$replay$mismatches` is empty" cannot hold on a correct implementation

- T18 replays row 6 (K, K', K) with `kit129SaveReopenReplay`-style code, i.e. `kitReplayJourney`
  (`tools/profiling/workflow_kit.R:1966-1977`). Its `mismatches` compare, for EVERY state name of the ORIGIN store,
  the origin fingerprint with the replay's (`kitFingerprintStore` takes `names(head$references)`, and
  `kitReadHead` sets `references = snapshot$states`, the whole states map, `:1380-1391`; `kitReproduceMismatches`
  reports `<state>/missing` for an origin state the replay does not have, `:827-833`).
- The ticket's own design keeps the superseded state's record in the store ("the state records stay in
  `snapshot$states`", 4.1; "its saved result stays in the project", VAL-346). After K, K', K the origin's states
  map holds `ruv_corrected`; the history, the record and the replay do not list it, so the replay's map lacks it:
  `mismatches` = `"ruv_corrected/missing"` on the fixed code, not empty. (`saved$replay$states`, the replay map's
  names, does equal the history `cf, rc`: that half holds.) T17's replay has the same property (the origin keeps
  `ruv_covariates`) but T17 asserts no `mismatches`.
- Audit S03 replaced "the replayed rc entity equals the live one" by "`mismatches` empty" and cited
  `test-kit129-run-path.R:1531, :1582, :1607, :1648`; in all four the origin's map equals its history, so they never
  met a superseded record.
- I did not write the test: an expectation believed wrong is not weakened or reworded by the implementer.
  Needs a ruling. Options (all keep the intent, "the replay's `rc` is computed from `cf`"):
  (A) assert exactly `as.character(saved$replay$mismatches) == "ruv_corrected/missing"` plus the replayed
  `ruv_covariates` fingerprint equal to the origin's (no `ruv_covariates/*` entry); (B) compare only the states of the
  record's history: `kitReproduceMismatches(populations(origin)[history], populations(replay))` called in the test
  (the harness file is not edited); empty on the fixed code, and on the base the history `cf, rc, rk` makes the
  replay refuse (class captured). I recommend (B).

### Further findings (each changes how a row is built; none blocks by itself)

- S-2 `completeCorrelationWorkflowFn` stubs with FIXED formals exist: `test-prot-05b-norm-module-contracts.R:2302` and
  `:2369` (`function(finalS4ForDe, workflowData, output, normData, correlationThreshold, skipped, successNotification,
  completionMessage, messagePrefix)`). Passing a new argument from the three call sites would break them ("unused
  argument"). Plan: the flag is passed only when TRUE (a RUV head before the correlation, which those tests never
  have). No stored test changes.
- S-3 `resetProtNormOutputs` (`R/mod_prot_norm_session_helpers.R:475-483`) is pinned by
  `test-prot-05b-norm-module-contracts.R:2877-2891` (`expect_identical(ls(output), "correlation_filter_summary")`, a
  KIT-122 F10 rule: the RUV summaries are not reassigned there). Ticket 4.3 puts the clearing "at Reset
  (`resetProtNormOutputs`)". Plan: clear `ruv_covariates_summary` in `runProtNormResetWorkflow` (limpa route only, new
  trailing argument `renderTextFn`), after `resetOutputsFn`; `resetProtNormOutputs` untouched. Same user-visible
  result; a departure from the ticket's place.
- S-4 The leftover-file removal (`source_dir/ruv_optimization_results.RDS`) cannot live in
  `completeProtNormCorrelationWorkflow` (no paths there); it is done at the three call sites (they hold
  `experimentPaths`) and at the Skip, through one helper beside `removeDiaLimpaRuvMarkers`.
- S-5 T1 "leaves `ruvcov_k` 'not recorded'" through `workflowRecordReportParameters` needs the template's parsed
  study object; T2 (the template's own chunks on the record's RUV outcomes) states the same thing at the report; T1
  asserts the record's steps instead.
- S-6 Managed rows: the dispatcher's hand-over of `held$rebuild` into `finish`, and the managed Skip's dispatch, are
  driven only by a recording runner (the full managed chain needs the module test server of KIT-266's 27 to 44
  minute file). The managed Skip row counts actions (one `limpa_composite` run, no in-application build).
- S-7 T6's enrichment half: no enrichment is produced by the row (it needs KIT-216's product machinery); the row
  asserts `diaEnrichmentProduct` NULL and the DA half in full.

### Line check (ticket citation -> head)

Holds, moved: `commitWorkflowReferenceState` `R/utils_workflow_dia_state.R:853` (KIT-216's marker line at `:893`,
history write `:891`); `truncateDIAStepProvenance` `:846` (keeps `enrichment_results$state`);
`selectWorkflowStateReference` `:1699`; adapter `revertToState` `:1920-1940`; Skip branch
`R/mod_prot_norm_workflow_helpers.R:80-91` (unchanged); the correlation Apply in process `:788-808`, managed
`dispatchProtNormCorrelation` `:818-870` (`held$rebuild` `:851`), Skip `:872-900`; `protNormLimpaRebuildAtCorrelation`
`R/utils_workflow_dia_norm.R:1004`, `removeDiaLimpaRuvMarkers` `:1014`; `buildDiaLimpaCombinedQcFigure` `:949`,
`dispatchProtNormLimpaComposite` `R/mod_prot_norm_workflow_helpers.R:~325`; `completeProtNormCorrelationWorkflow`
`R/mod_prot_norm_correlation_helpers.R:382-413`; `resolveDIARuvRunDescription` `R/utils_workflow_dia_products.R:721`
(now also returns `norm_method`); archive `collectDIANativeSessionMetadata` `R/utils_workflow_dia_sessions.R:1972`
(`ruv_mode` `:2001-2008`, `ruv_k` `:2009`), the refresh `:1800-1830` (`ruv_mode` rewritten only when the description
is non-NULL, `:1816-1820`); replay RUV branches `R/utils_workflow_analysis_record.R:1058-1066`;
`commitDIADaArguments` `R/utils_workflow_dia_da.R:672` (no-commit test `:680`); `diaEnrichmentProduct`
`R/utils_workflow_dia_state.R:965`; `workflowDaResultStatus` `R/utils_workflow_dia_da.R:125`. The declared edit
`tests/testthat/test-workflow-integration-input.R:703` is at the same line. R2 grep at the branch point: 96 files
(ticket: 91; the five added are `test-kit216-prot-session-after-enrichment.R`, `test-kit266-limpa-combined-qc-figure.R`,
`test-kit271-limpa-table-name.R`, `test-kit273-limpa-ruv-delivery.R`, `test-kit274-trx-no-ruv-export.R` or already
counted; recount with the new file at the request), so the budget table needs restating at the request (priced at
the newest timing).
Note 1.5 (cancor): `removeDiaLimpaRuvMarkers` is called before both RUV runs (`:117` in process, `:395` managed);
the file is written by the covariates publication (`R/utils_workflow_dia_products.R:1843-1845`).
