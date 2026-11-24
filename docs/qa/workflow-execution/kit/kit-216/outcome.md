# KIT-216 outcome (preparation stage; the ticket is at r1.3)

This records the preparation stage only. The tests-only commit's run (R1) is not sent until KIT-215 has merged
(ticket 15a). No fix code is written yet.

## G0: the proteomics Summary calls in the R2 scope (run 45491444529, commit f4a94eb0 = staging f93dfd85 + script)

**Instrument.** `g0-inventory.R` wraps the five Summary functions in its own process. The wrapped functions are
`saveWorkflowSummaryParameters`, `renderWorkflowSummaryReport`, `writeWorkflowSummarySessionState`,
`protSummaryIntegrationWrite` and `resolveProtSummaryFinalS4State`. Each wrapper records the call site, the owner and
its head, store and export, and the design's prediction. It then runs 14 test files, each in its own process with a
15-minute cap. The tables are `g0-calls.tsv` (124 calls) and `g0-files.tsv`. The test counts under the probe are not
regression evidence.

### Observed

- **Calls recorded.** Calls were recorded from 7 files: test-workflow-integration-input, test-workflow-native-session,
  test-workflow-analysis-record, test-kit206-prot-norm-export-guard, test-omic-summary-module, test-omics-summary-core
  and test-kit172-downstream. The last of these made no proteomics call. No proteomics call was recorded from
  test-kit207-summary-head-guard, test-workflow-omics-application-reads or test-module-ci-prot-summary.
- **Files that ran nothing.** test-workflow-memory-mechanism-member skipped all 13 tests on the runner.
  test-prot-12-summary-module-contracts and test-prot-12a ran 0 tests.
- **File cut short.** test-kit129-run-path was killed at the probe's 15-minute cap (exit 124), so no call was recorded
  from it.
- **Calls that move under I10 or SD-1** (reference-state owner, head outside the four DA states). All are already
  declared:
  - `test-workflow-integration-input.R:481` (call `:288`): head `imputed`, no normalisation in the store. Refused with
    `summary_none` (SD-1). This is the 6.4 row added by R1-02 (a).
  - `test-kit206-prot-norm-export-guard.R:629` (E11b, call `:645`): head `protein_replicate_filtered`, a normalisation
    in the store, an export. Refused with T3 (I10). This is the 6.4 row "E11b".
- **Calls that move under SD-2.** The `prot_native` case of test-omic-summary-module (`:4636`, call `:624`) and its
  reader test-omics-summary-core (`:295`, call `:135`): a native reference owner at `correlation_filtered` with no
  export. These are 6.2 and R1-01.
- **The fake owners.** The class-only `DIAStateAdapter` calls (`test-workflow-analysis-record.R:139`, `:255`, `:367`)
  have `usesDIAReferenceStates` FALSE and pass the guard, as R1-02 (b) requires.
- **No move.** The native-session Save (`test-workflow-native-session.R:394`) passes the guard, because its owner is not
  a reference-state owner. Its refresh returns at the same head.
- **No other proteomics call moves.** Every other proteomics call is at one of the four DA states, or its owner is not a
  reference-state owner.

### Inferred (by reading, not executed)

- **A test that needs item 11's result-head clause.** `test-workflow-integration-input.R:437` saves at an
  `enrichment_completed` head built by direct `saveState` (`:473`, a store shape listed in 3.1 item 12). G0's naive
  prediction for that call was "refuse T3". Under 3.1 item 11, a result head is judged by the table it was made from
  (item 6's rule), here `correlation_filtered`, so the call passes. This test therefore depends on item 11's
  result-head clause being implemented. That is not a contradiction of the ticket.
- **The final object moves only where nothing reads it.** The final object (P06) changes at `:437` (`normalized`),
  `:586` (`ruv_corrected`) and `:631` (`correlation_filtered` after `ruv_covariates`). These tests stub
  `createWorkflowArgsFn` and assert only the integration input, which `protIntegrationInputSource` chooses, so no
  assertion moves.
- **test-kit129-run-path, not observed.** Its call (`:1664`) is a direct `protSummaryIntegrationWrite` on a non-native
  owner. No guard is involved and SD-2 does not apply, so no move is expected.
- **Gui-entry lanes.** They were not probed. Their stored heads at the Summary steps were read from the JSON: `prot_dia`
  is at `enrichment_completed`, which moves by design (6.1). The other four lanes are at `ruv_covariates` or
  `correlation_filtered`, which pass.

**Result against the ticket.** No contradiction. 6.4 needs no new row from G0. The `:437` dependency above is
recorded for the fix.

## Not checked

- Calls in test-kit129-run-path (cap), test-workflow-memory-mechanism-member (skipped on the runner) and the
  test-prot-12 files (no tests ran).
- Calls that reach the Summary only through the managed worker (not in this process).
- Every gate G1 to G10 and every E row: the draft test file is committed (5e91f1cc) but not run.

## Core-hours

G0 took 76.0 job-minutes, which is 5.07 core-hours. It is not in the ticket's budget table; it is charged against the
total ceiling of 561.33.
