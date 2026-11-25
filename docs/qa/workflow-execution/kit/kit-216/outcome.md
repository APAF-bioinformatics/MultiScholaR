# KIT-216 outcome (ticket r1.4: one correction pass after R2, per audits/AUDIT-KIT-216-r2-20261126.md)

Runs: `runs.tsv`. The comparison run after the correction (R2c) and its result are added below when in. (Commit
fa006ff9's message announced this file's update; it carried `runs.tsv` only. This file was updated in the next commit.)

## Commits on `kit-216`

| Commit | What |
|---|---|
| f4a94eb0 | G0 script (preparation) |
| 92607bf8 | merge of staging 16112917 |
| cc986d2a | tests-only commit (R1) |
| 4320e668 | tests-only commit 2: base-made archives (`tests/testdata/kit216/`), E9, instrument repairs to E15 and G16 |
| 854dcfc4 | the fix (R2's head) |
| 7c9a0850 | r1.4 tests first: G21, E9's restore check |
| 45b3026c | r1.4 correction pass: F04, F01, test edits F05 to F07, T5 and SD-1 text, F03's NEWS sentence |

Base-made archives (R1, run 45588146575, on cc986d2a through the application's commit, then
`refreshDIANativeSummaryArchive`; 1,812,480 bytes each): `base-enrichment-1-run.msdia.tar` sha256
`03ea8cc5d1574ceef87fb49df8d613eb26902139382a15c52e72f58747c02716`; `base-enrichment-2-runs.msdia.tar` sha256
`e9bcf254b85aa4169935e9aa18e8d31ac50f01ad047bdd902f57279c7b3fddcf`.

## Gates and failing-first (observed)

- **G0, G0b** (45491444529, 45587846430): the proteomics Summary calls move only where 6.2 and 6.4 declare it
  (`test-workflow-integration-input.R:557` SD-1; KIT-206 E11b I10; `prot_native` SD-2); every other passes. Section
  "Preparation: G0" below; tables `g0-calls.tsv`, `g0b-calls.tsv`.
- **R1** (45588146575, cc986d2a): every failing-first row failed at its stated base value (`r1-evidence.txt`,
  `r1-failures.txt`): E1 "native export needs component metadata; explicitly convert the legacy session first"; E2
  "Study parameters saved with warnings", planted input kept; E3 and E4 Open refuses (`saved_entity_corrupt`
  "component checkpoint is incomplete"; the premise holds); E5b head `enrichment_completed` and the step listed; E7
  `cyclicloess` (also the reopened selector); E7b and E7c `quantile` (the selector); E8 no input and the step's input
  itself; E9 (G9) render called and RDS written; E10 parameters and record written, `failed`; E10c written; E11
  `covariates`; E14 `protein_replicate_filtered`; G18 refused by the refresh; G20 refused and input removed; E9b Open
  refuses. G6's feasibility passed (E11 and E7c ran).
- **Instrument defects (VAL-54) in R1**, repaired in 4320e668 with no expected value changed, rerun on R2's base side:
  E15's setup (the correlation filter at 0.5 removes no protein on this recipe, so the DA tab kept the bus manifest; the
  bus is set to NULL as the tab does) and G16 (the record value was not decoded). R2 base: E9 and E15 failed as stated;
  G16 passed, as a guard row should.

## R2 (45591899662; 854dcfc4 against 4320e668; 80 files; check)

- Observed as predicted: the new file `fixed` (145/0, 1 skip); check identical; the declared moves of 6.1 and 6.2 at
  their leaves (`prot_dia`, W1, W113, W2 readers, `prot_native` and its reader); `test-workflow-integration-input.R`
  count-changed; G11's six record printouts identical between base and head (no record field moves; stop rule 5).
- Not predicted (stop rules 2 and 3), the reason for the one correction pass: G7 (F04), KIT-024 (F01), S10 (F05), four
  history rewrites (F06), A107-08's mock (F07), three W1/W113 readers (F08).

## Declared design moves (6.4): old and new values

| Test (file:line) | Old | New |
|---|---|---|
| `test-omics-enrich-core.R` F-16 test 4 | mock `saveDIAOwnerState`; `current_state` "enrichment_completed" | mock `commitWorkflowEnrichmentProduct`; `current_state` "next_state"; not in the history |
| `test-kit205-enrichment-background.R:306`, `:335`, `:343` | "enrichment_completed"; `replay$states` ending with it | "correlation_filtered"; `c("initial","correlation_filtered")` |
| `test-go-propagation.R:308` | three states | `c("initial","correlation_filtered")` |
| `test-workflow-analysis-record.R:1018`, `:1072` | "enrichment_completed"; three states | "correlation_filtered"; two states |
| `test-omic-enrich-module.R:3787` (W3 `committed`) | `"enrichment_completed" %in% history` | `!is.null(diaEnrichmentProduct(manager))` |
| `test-omic-enrich-module.R:3804` (W3 ×2), `:4255` (P-07), `:4398` (A107-04 no terms) | `c("initial","correlation_filtered","enrichment_completed")` | `"<unchanged>"` |
| `test-omic-enrich-module.R:4415` (A107-08) | mock `saveDIAOwnerState` | mock `commitWorkflowEnrichmentProduct` |
| `test-kit206-prot-norm-export-guard.R:664` (E11b) | "The session could not be saved: " + KIT-206's text | ... + `summary_prot` (T3) |
| `test-workflow-integration-input.R:557-568` | "not written (no final proteomics matrix yet)" | refused, class `omic_norm_correlation_input_stale`; earlier input removed |
| `test-enrichment-significance.R` S10 | (no mock) | `omicEnrichManifestCurrent` mocked TRUE |

Stored references (6.1, 6.2) are re-frozen only after the controls (pending).

## Unreachable branches kept (3.1 item 12)

KIT-206's result-class branch (`R/mod_prot_norm_session_helpers.R`), `DIA_QC_REFUSED_HEADS`' `enrichment_completed`,
KIT-212's result-head clause (`workflowDaTableState`), item 11's result-head reading, and item 4's history clause (r1.4:
only when the head is the enrichment state itself): reachable only by test-built old-shape stores or an archive before
its conversion.

## Hand-off: F03 (coordinator's ruling (b), 2026-11-26)

A base project with enrichment, then 'Run Normalization & RUV' at the enrichment head and the correlation filter, then
exported, could not be opened before and now opens. Its analysis record still lists the earlier enrichment as the last
step with the current table as input, and the re-run normalisation's input names the enrichment record (the D2 class:
a commit made at a result head). No number changes. Stated in NEWS (DRAFT sentence, also in `wording.md`); no code;
G21 observes the store shape and its record. For a separate ticket: how the record treats a step committed at a result
head.

## Acceptance

| # | State |
|---|---|
| A00 | met (VAL-179; re-check 1; R2 audit) |
| A0 | met (R1; G0, G0b; archives committed with their sha256) |
| A1 | R2 head 145/0; after the correction: R2c (pending) |
| A2 | R2 failed it (stop rules 2, 3); R2c (pending) |
| A3 | not run: the merge check rides a train (coordinator) |
| A4 | controls and re-freeze pending, after R2c |
| A5 | NEWS, ledger row committed; no version, `NAMESPACE` or `man/` change |

## Not checked

- No R locally; every result is a runner's. G21 and E9's restore check were not run on 854dcfc4 (added in r1.4).
- The managed route's adoption and item 13 under a managed session: by reading only.
- LFQ and TMT native owners (G12): DIA only. G17 (managed Save at E10's head): no proteomics managed-runner fixture;
  by reading, KIT-207's guard runs in the controller. G5/E12 (peptide checkpoint, then Summary with no normalisation
  export): not built; by reading, the refresh requires `output_policy = "native_first_v1"`, which the peptide
  checkpoint's metadata lacks.
- The running application, the browser, Windows, real data.

## For the change report

- **What was wrong.** A proteomics project's session file written after enrichment could not be reopened; Export after
  enrichment refused with an untrue "legacy session" message; any session from a project that had ever run enrichment
  could not be opened; without a session export in the current application session the Summary's Save wrote no
  integration input and no session file and said only "saved with warnings"; a second enrichment run removed the
  integration input and the record named the enrichment as its own input; the session file named the selector's
  normalisation method, not the run's, and on DIA with limpa an earlier RUV mode; the proteomics Summary did not refuse
  a table that was not a normalisation result, and its final object ignored the head; the enrichment Run could run on a
  DA no longer current.
- **Who.** Every proteomics user who ran enrichment and saved or reopened; every reopened project at the Summary; users
  who changed the normalisation selector after a run.
- **Code** (base 16112917; fix 854dcfc4 + r1.4 45b3026c): `R/utils_workflow_dia_state.R` (`commitWorkflowEnrichmentProduct`,
  `diaEnrichmentProduct`, marker cleared on Revert and on a new commit); `R/utils_workflow_enrichment.R` (commit,
  restore); `R/utils_workflow_analysis_record.R` (the product step; replay's check); `R/utils_workflow_dia_sessions.R`
  (Open's loop; `convertDIANativeEnrichmentHead`; T4; method and RUV statements; missing-export class);
  `R/utils_workflow_dia_products.R` (`resolveDIANormalisationMethod`); `R/utils_workflow_dia_summary.R` (SD-2, SD-4; final
  object); `R/utils_workflow_summary.R` (proteomics head rule); `R/mod_omic_enrich.R` (item 13); `R/mod_omic_norm.R`
  (texts); `R/mod_prot_norm_observer_helpers.R`, `R/mod_prot_da_handlers_load_session.R` (selector preset; T1).
- **Results.** No computed value changes; the head after enrichment is the table; the refusals as listed; the archive's
  two statements; the final object at `normalized`.
- **Evidence.** G0 45491444529, G0b 45587846430, R1 45588146575, R2 45591899662, R2c (pending). Synthetic fixtures;
  mechanism only.

## Preparation: G0

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

