# KIT series remediation record (2026-10-09)

Register: `KIT-SERIES-AUDIT-20261009.md` (unchanged). Authority: user, "fix all blockers" (2026-10-09). Scope of this pass: the ten blockers only. Should-fix and note items stay open.

| Finding | State | Closure evidence |
|---|---|---|
| A1 peptide intensity filter unit | closed | `global_statistics` pre-pass added to the declaration contract (KIT-004 item 1); KIT-009 declares the intensity filter's dataset-wide quantile as the pre-pass; plan/README unit tables moved it to feature group after the pre-pass; KIT-014 wording updated. |
| A2 protein missing-value filter unit | closed | Plan/README unit tables and KIT-015 move the `removeRowsWithMissingValuesPercent` family to whole-matrix; KIT-010 proves the row partitioner on the protein replicate filter and declares q-value (also closes A3's contradiction for q-value: row, KIT-010; KIT-014 depends on KIT-010). |
| B1 CORE-pinned harness files | closed | KIT-001 and KIT-005 file lists contain only new harness files (`workflow_kit.R`, `workflow_kit_driver.R`, `workflow_kit_coordinator.R`, `run_workflow_kit.sh`); regression check now verifies the four pinned files' sha256 against `core/bindings.json` and `core-018-floor --validate-only`. |
| B2 populations/protocol circularity | closed | KIT-001 proves snapshot-and-fingerprint on the named screen corpus and freezes protocol v1; KIT-005's counted control journey produces every post-rollup population and oracle, its second member reproduces them; plan ledger and protocol §3/§5 updated. |
| B3 off-ledger oracle runs | closed | No in-process full-size oracle runs remain; oracles come from the counted control journey (KIT-005 item 2; protocol §3). |
| B4 runner identity scheme | closed | KIT-004 item 2: the runner calls `diaComponentScientificIdentity` verbatim with existing operation names and factors; normalisation and RUV stay one declaration; rollback note corrected; `R/mod_prot_norm_ruv_helpers.R` and DESCRIPTION added. |
| B18 no worker in unmanaged launches | closed | KIT-008 item 1 states the dispatch rule (worker when `session$userData$workflow_execution` is present, in-process otherwise) and tests both. |
| B19 KIT-008 before the adapter | closed | Index and dependencies reordered: KIT-011 (deps KIT-004, KIT-007) precedes KIT-008 (dep KIT-011); KIT-012 depends on KIT-008; controls re-pointed (KIT-011 vs KIT-007 build, KIT-012 vs KIT-008 build). |
| B20 member cannot exercise the GUI | closed | KIT-008's member is a `gui`-mode managed launch driven through the rollup tab (shinytest2 against `managed_url`), 10 min; legacy-vs-P14 page identity checked before launch. |
| B29 entry fixture moved to a helper | closed | KIT-013 item 4 keeps a generic production entry hook injected by the managed launch; file list adds `utils_workflow_execution.R`, the KIT driver, `test-perf107-entry.R`; the non-existent test file names are replaced. |

Validation after repair: all 26 KIT tickets pass the schema and OpenAI validators; the selector returns KIT-001; stage-1 counted total recounted to 34.
