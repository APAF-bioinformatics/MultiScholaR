# KIT-034: the series' own quarantine entries

Status 2026-10-14: **complete** (no counted run; no install: no package code changed). Authority: `tickets/kit/00-KIT-PLAN.md` (decision row 2026-10-14, "plan it, ticket it, do it"; run ledger KIT-034: none); protocol v1.4. Control `kit_018`; KIT-017's candidate: `kit_024` (`core/browser-smoke.70tjXp`); this ticket changes no package code.

## Entries closed

| Entry | What was wrong | Repair |
|---|---|---|
| Q-01 | The intensity and dedup shared fixtures' fake states had no `protein_id_column` slot; KIT-015's in-process apply step counts proteins through it, so the step errored inside the module's handler and the tests saw empty captures | The fake classes declare the slot (`Protein.Ids`); the copy helper carries it |
| Q-05 | `perf102Code`/`perf102Evidence` lived at the top of `test-perf-003.R`, so the PERF-102 file could not run alone; the sourced copy of the execution file also predates the KIT-007 registry and the KIT-024 dispatch refusal | Helpers moved to `helper-perf102.R`; the sourced copy gets a minimal `workflowActionSpec` and a pass-through `assertWorkflowOperationDispatchable`, as it already got a stubbed validator |
| Q-09 | Expected the pre-KIT-008 `gui`-mode default variant | Expects the registry's `default_variant` |
| Q-15 | Expected `aa_seq_tbl_final` in the global environment | Expects it on the workflow bus and absent from the global environment (KIT-013) |
| Q-16 | The coordinator keeps the workflow type in a reactive value set by an observer; the eager test mocks had neither | `reactiveVal` and `observe` added to the eager mocks |
| Q-20 | The GUI diagnostic trace lists named the capacity-session functions KIT-011 removed | Names dropped; the owner observation follows the import adoption |
| Q-04 | Quarantined as "needs the installed package" | `skipIfNoInstalledRuntime()` (`helper-installed-runtime.R`); the `perf-003` dispatch test additionally had a job fixture without the KIT-008 `monitor_borrowed` flag, fixed; run against the installed `kit_024` library: perf-003 173, native-session 541, saved-recovery 108, all passing ([installed-runtime-run.md](installed-runtime-run.md)) |
| Q-14 | Quarantined as "needs a browser" | `skipIfNoBrowserRuntime()`; run against the installed library on this host: ruv-pearson 326 and the passive wrapper test passing with `CHROMOTE_CHROME=/usr/bin/brave` (the harness's browser; chromote's discovery does not find Brave) and the library named ([installed-runtime-run.md](installed-runtime-run.md)) |

## Entries kept (pre-series), one reclassified

Q-02, Q-03, Q-07, Q-08, Q-10, Q-11, Q-13, Q-17, Q-18, Q-19 stay for the stage-2 pre-step. Q-12 was listed as a warning; it is a fixture (the test parses a fake path and expects fixture sequences, but the reader opens the file itself), so it stays quarantined with the corrected cause.

## Q-11, raised for the user (no change)

`fitProtDaBaseModel` builds the block as `group_replicates` and uses `duplicateCorrelation` whenever that block has duplicates. On a design where `replicates` is a biological replicate index, the block collapses to the group and limma (3.68.4) warns that the block is already encoded in the design and sets the intrablock correlation to zero; the block then does nothing. Recommended rule: technical replicates are declared explicitly (the `technical_replicate_id` slot), the DA design carries the biological factors only, the declared technical replicates are the block (`duplicateCorrelation`) or are collapsed before modelling, and a block spanned by the design means there are no technical replicates (standard `lmFit`, recorded) rather than a silent zero. RUV-III's replicate structure and the DA block must come from the same declaration. Changing the model moves the DA oracles, so it is a user decision for after the closure, with a re-recorded oracle.

## Tests

| Suite | Result |
|---|---|
| test-prot-02d-qc-protein-intensity-shared.R / test-prot-02l-qc-protein-dedup-shared.R | 52 / 32 pass (Q-01) |
| test-perf-102-audit-repairs.R / test-perf-003.R | 74 / 137 (1 installed-runtime skip) pass (Q-05, Q-04) |
| test-perf106-repair-v2.R | 75 pass (Q-09) |
| test-workflow-server-shared.R | 70 pass (Q-15) |
| test-prot-qc-coordinator-direct-shared.R | 21 pass (Q-16) |
| test-workflow-gui-diagnostic.R | 72 pass, 4 unrelated skips (Q-20) |
| Installed-runtime run ([installed-runtime-run.md](installed-runtime-run.md)) | perf-003 173, native-session 541, saved-recovery 108, ruv-pearson 326, passive wrapper test passing |
| Full suite | 343 files, 29,633 pass, 0 fail, 0 error, 139 skips (152 before this ticket) |

## Claim limits

- The installed-runtime tests are proven on this host's installed candidate library once; the skip condition keeps them out of `load_all` runs by design.
- No package code changed except the two harness trace lists under `tools/profiling`; no install for this ticket.
