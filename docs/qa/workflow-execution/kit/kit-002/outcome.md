# KIT-002: class declaration registry and declaration-driven checkpoints

Status 2026-10-09: **complete** at `4865be8e` (implementation) plus this record. No counted member; one recorded install check. Authority: `tickets/kit/00-KIT-PLAN.md`; protocol v1.1.

## What changed

- `R/utils_workflow_class_declarations.R` (new): `workflowClassDeclarations()`, `workflowClassDeclaration(class)` (refuses unknown classes with `workflow_class_undeclared` and validates every named slot against the S4 definition), `validateWorkflowClassDeclaration`, `workflowClassRequiredComponents`, `workflowDeclarationField`. Entries: `PeptideQuantitativeData` (assay `peptide_data`, long, tables `peptide_data` + `design_matrix`, matrix `peptide_matrix`, no describer) and `ProteinQuantitativeData` (assay `protein_quant_table`, wide, tables + `protein_id_table`, describer `describeDIAProteinState`). Collate entry added.
- `R/utils_workflow_dia_state.R`: new `readComponentMetadata(manager, reference, classes)` reads any declared class and takes the required components from the declaration; `readDIAStateMetadata` wraps it with the two DIA classes and the unchanged `saved_entity_corrupt` refusal (45 callers unchanged). `checkpointDIAExecutionResult` and `newDIAPeptideCheckpoint` take the assay and design slots and the describer from the declaration and write the additive top-level `declaration` field (`class`, `version`), never inside `slots`. `saveDIAStateRecord` and `materializeDIAStates` take their field lists from the declaration. `inheritDIAExecutionSource`, `resolveDIAOwnedValue` and `diaComponentScientificIdentity` are untouched (audit B8, KIT-004 dependency).
- `tests/testthat/test-workflow-class-declarations.R` (new, 48 assertions): declarations and refusals; both classes round-trip exactly (full and projected reads, zero-row prototypes); a record without the field reads identically and the scientific identity ignores the field; a foreign class refuses at the metadata boundary with `saved_entity_corrupt`.

## Evidence

| Check | Result |
|---|---|
| `test-workflow-class-declarations.R` | 48 pass |
| `test-workflow-checkpoint-dia-state.R` | 238 pass |
| `test-workflow-saved-qc-rollup.R` | 76 pass, 1 pre-existing skip |
| `test-workflow-bulk-relational.R` | 45 pass, 3 pre-existing skips |
| `NAMESPACE` | unchanged, 484 `export(` lines |
| Installed check: `install-corrected` of `4865be8e` (`core/browser-smoke.401rH2`, no runtime patch), then guarded `gui-bulk-rollup:writer-r-preflight` (`core/browser-smoke.qGXvfo`) | 174 checks passed, 0 failed, 2 pre-existing skips; the exact relational rollup cases (129 assertions) and the P14 batch cases pass; peak 0.78 GiB |
| `R CMD check --no-manual --no-tests --no-build-vignettes` on tarballs built from clean worktrees of the candidate (`4865be8e`) and the baseline (`ec68dfb3`), proof library on `R_LIBS`, suggests not forced ([r-cmd-check.log](r-cmd-check.log), [r-cmd-check-baseline.log](r-cmd-check-baseline.log)) | identical apart from timestamps: 4 warnings, 3 notes on both; no new warning or note |

Note on the check: a first pass from the working tree failed at package dependencies (GSEABase and MOFA2 are not installed here) and, with suggests not forced, at configure (the exact-sum artifacts are located through installed dependency packages); both are environment facts, resolved by the proof library on `R_LIBS`. A working-tree build also carries ignored local artifacts the baseline lacks, so both sides were built from worktrees.

## Scientific contract

None changed: table components, products and identity digests are identical; the checkpoint record differs only by the top-level `declaration` field, which the identity digest does not read (tested).

## Hand-off

The registry is the entry point for KIT-003 (list-of-tables slots and result classes) and KIT-004 (operation runner). The candidate build is not a control; controls change at KIT-006.
