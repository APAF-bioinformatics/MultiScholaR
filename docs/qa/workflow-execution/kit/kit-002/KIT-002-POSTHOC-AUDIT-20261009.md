# KIT-002 post-hoc audit (2026-10-09)

Mode: `post_hoc`, audit only. No audited artifact was changed. Scope: the completed KIT-002 at `4865be8e` (implementation) and `1589d0f0` (closeout) against the ticket's acceptance criteria, protocol v1.1 and the plan ledger.

## Mechanical and provenance

| Check | Result |
|---|---|
| `audit_inventory.py` on `tickets/kit` KIT-002 | 0 findings, exit 0 |
| Commit scope `4865be8e` | exactly the declared files: `DESCRIPTION` (Collate), `R/utils_workflow_class_declarations.R`, `R/utils_workflow_dia_state.R`, the new test; 188 insertions, 15 deletions |
| `NAMESPACE` | unchanged (484 `export(` lines); all new functions `@noRd` |
| Install record | `install-kit002-start.txt` (commit `4865be8e`), `core/browser-smoke.401rH2`, runtime patch 0 bytes |
| Installed check | `gui-bulk-rollup:writer-r-preflight` member `core/browser-smoke.qGXvfo`: 174 passed, 0 failed, 2 pre-existing skips; the exact relational rollup cases (129 assertions) and both P14 batch cases pass; peak 0.78 GiB |
| R CMD check | candidate (`4865be8e`) and baseline (`ec68dfb3`) built from clean worktrees, proof library on `R_LIBS`, suggests not forced: logs identical apart from timestamps (4 warnings, 3 notes); kept as `kit-002/r-cmd-check*.log` |
| Regressions | class declarations 48, checkpoint state 238, saved QC rollup 76 (+1 pre-existing skip), bulk relational 45 (+3 pre-existing skips) |

## Implementation and semantics

- Subclasses: no class in `R/` extends `PeptideQuantitativeData` or `ProteinQuantitativeData` (`PeptideQuantitativeDataDiann` is a constructor), so `class(object)[[1L]]` lookups are equivalent to the previous `methods::is()` branches.
- Undeclared S4 states cannot reach the rewritten loops: `saveDIAStateRecord` fails first at `rememberDIAWorkers` (`no slot of name "args"`, probed on the current code) and the baseline code failed at the hard-coded `protein_quant_table` slot access, so there is no behaviour change; `checkpointDIAExecutionResult` now refuses with a typed class instead of a slot error (improvement).
- The additive `declaration` field is a list of two atomics at the record's top level; it is not under `slots` (identity digest tested unchanged), is not visited by `collectDIANativeReferences` call sites (they take components, support, products, dependencies, snapshots), and round-trips through `qs2`/`rds` (tested). `newDIAPeptideCheckpoint` and `convertDIANativeLegacyCheckpoints` (through `checkpointDIAExecutionResult`) also write it.
- The KIT-001 fingerprint summary digests `slots`, `dimensions`, `class`, `config` and `description`, not `declaration`, so KIT-001's attempt-4 fingerprints stay comparable with records written by this build.

## Findings

| ID | Severity | Observed | Impact | Smallest repair | Blocked |
|---|---|---|---|---|---|
| PB-01 | medium (attributable to KIT-001) | The fourth R CMD check warning, present in both logs, is a code/documentation mismatch for `peptideIntensityFilteringHelper`: the argument `require_complete_design_runs` added by `9075b361` is missing from the roxygen `@param` line and `man/peptideIntensityFilteringHelper.Rd` was not regenerated. CORE-018 recorded three warnings. | The series introduced a new check warning in KIT-001, which had no R CMD check requirement; KIT-002's "no new warning" is true against its own baseline but not against CORE-018. | Add the argument to the `@param` line and regenerate the Rd; re-check (expect 3 warnings). | none |
| PB-02 | low | The ticket's regression check asks for a tarball "built on `/home`"; the checks ran from the session scratch directory under `/tmp` (the earlier CORE-018 check used `~/.cache`). The check installed and passed, so the location did not matter here. | Deviation from the ticket wording. | Record; use `~/.cache` next time or amend the wording. | none |
| PB-03 | low | Two environment facts were needed to run the check at all (proof library on `R_LIBS` for the configure-time exact-sum artifacts; suggests not forced because GSEABase and MOFA2 are not installed). They are recorded in the outcome but not in the ticket or protocol. | Repeatability of the regression check depends on undocumented environment. | Add the recipe to the KIT README's checks section (or the next ticket's regression wording). | none |
| PB-04 | info | `materializeDIAStates(share_immutable = TRUE)` now calls `workflowClassDeclaration()` on every S4 state instead of intersecting a fixed list; unreachable difference (undeclared S4 never saves), but no test pins the undeclared-class behaviour of this loop or of `saveDIAStateRecord`. | Test coverage gap for a refusal path. | One test expecting the typed refusal from `checkpointDIAExecutionResult` for an undeclared S4. | none |
| PB-05 | info | `readComponentMetadata(classes = NULL)` accepts any declared class; today that is the same two classes, so `readDIAStateMetadata` and the generic are equivalent until KIT-003 adds entries. | None now; KIT-003 must decide which callers may see the new classes. | — | none |

## States

| State | Value |
|---|---|
| implementation_verified | yes |
| execution_complete | yes (installed check and R CMD check comparison) |
| decision_supported | yes; no empirical decision was claimed (no counted member) |
| promotion_authorized | not applicable |
| claim_scope | DIA behaviour unchanged for the two declared classes; identity digests unchanged; one additive record field; the candidate build is not a control |

## Next action

No blocker for KIT-003. PB-01 is a KIT-001 documentation defect worth fixing before the next R CMD check comparison (one roxygen line and a regenerated Rd); PB-03 is a one-paragraph README addition. Both are bounded repairs for a remediation record on request.
