# KIT-049: DIA replicate-mean imputation optional; the DIA-NN route named at the entry point

Status 2026-10-21: **completed** (ticket revision 1.1). Implemented at `c07494e8`; install `kit_049` = `core/browser-smoke.4oKgSd`; no counted member.

## One mechanism: a declared skip on the imputation step

- `diaPeptideImputationOperation` declares `skip` (logical, default `FALSE`) beside `proportion_missing_values` and `exclusion_column`; `[peptideMissingValueImputation] skip = false` in `inst/config/config.ini`; label "Skip replicate-mean imputation" (KIT-039 map); the worker validator (`validateWorkflowPeptideImputationParameters`) accepts a logical scalar and refuses anything else.
- The compute path (`computeDIAReferencePeptideGroups`) never hands `skip` to the kernel. With `skip = TRUE` each batch runs `peptideMissingValueImputationSkipped` (the shape of the kernel's own skipped branch: `Peptide.Imputed` equals the raw quantity, `is_imputed` FALSE on every row, no support table, summary `skipped`/`user_skip`, the section arguments as the kernel records them plus `skip = TRUE`) and the result is committed as the `imputed` state through the usual finishing step. `skip` enters the recorded arguments, the configuration and the saved-operation identity only when TRUE, so the default path's arguments, configuration, identity and fingerprints are unchanged by construction.
- The hand-off (`buildDIAReferencePeptideHandoff`) carries `skip` and `skipped` into the state's config with the description "Missing value imputation skipped by user (missing values kept for IQ)"; KIT-038's record builder lifts `skipped` into the step's `outcomes`.
- Module: a checkbox "Skip imputation (keep missing values for IQ)" bound to the parameter on the managed-worker path (`parameters$skip <- TRUE`) and the in-process path (`runPeptideImputationStep(skipImputation = TRUE)`, which uses the same skipped output); the status text comes from the audit summary (`skipped`, `user_skip`).
- Entry point: the DIA-NN radio reads "DIA-NN workflow (choose before import)" with the one choice "DIA, iq MaxLFQ (optional replicate-mean imputation)" (`iq`); `limpa` is no longer selectable (its removal from the code and the second route's label are KIT-046's).

## Evidence

- Screen check by record replay ([skip-check-screen.R](skip-check-screen.R); the KIT-028 replay machinery on reflinked copies of the kit-journey-028 project, cut to import and design): [skip-check-tree.json](skip-check-tree.json) on the working tree and [skip-check-installed.json](skip-check-installed.json) on the installed build.
  - **Default path unchanged**: the record replayed as recorded on this build reproduces the origin's **16 states with 0 mismatches** (every component, support entry and summary; DA included).
  - **Skip path**: with `skip = TRUE` added to the imputation step as `user_set`, the states before `imputed` are identical to the origin's, `imputed` differs, every later state differs (as it must), and the replay ran IQ through DA (`DA-contrast-output` replayed). On the replay store the skipped state has 224,359 rows (the origin's imputed state 279,240: the skipped output generates no rows), `Peptide.Imputed` equals the raw quantity on every row with the same 16,306 missing values, no row is flagged imputed, the section arguments carry `skip = TRUE` with summary `skipped`/`user_skip`, the state's config carries `skip` and `skipped`, and the record written from the replay store shows `skip = TRUE` as `user_set` with `outcomes$skipped = TRUE` and validates. The replayed record with the skip is kept beside ([analysis_record-skip-tree.json](analysis_record-skip-tree.json)).
- Tests: `test-workflow-operation-registry.R` 323 (the declared parameter, its defaults mapping and label), `test-prot-02k-qc-peptide-impute-shared.R` 63 (the skipped output and the validator), `test-workflow-execution.R` 64, `test-workflow-dia-products.R` 176, `test-prot-02h-peptide-audit-shared.R` 93, `test-module-ci-prot-peptide-qc.R` 43, `test-prot-qc-peptide-methods-direct-shared.R` 119; dependents grep (`diaPeptideImputationOperation`, `computeWorkflowPeptideImputationAction`, `dia_quantification_method`, `peptideMissingValueImputation`): 13 files, all green (native session 606, partitions 306, declarations 119, harness 382, checkpoint state 236, saved exports 276, fixtures 101, perf107 49 and 29, config contracts 52, memory-mechanism member skipped by design).
- `R CMD check --no-manual --no-tests` at `c07494e8`: 3 WARNINGs, 3 NOTEs, the baseline ([log](check-c07494e8/r-cmd-check.log)). Install recorded; bindings regenerated (`kit_049` bound; the control stays `kit_036`).

## Backward compatibility and rollback

- Stores and records written before this ticket read unchanged: `skip` absent reads as `FALSE` (the default replay above is that proof); the config key is optional; the project `config.ini` written for new projects gains `skip=false` (the head's `defaults_sha256` moves for new projects only; `kitFingerprintState` excludes the head).
- `git revert c07494e8`: a store written with `skip = TRUE` reopens on the previous build as an ordinary component checkpoint; its provenance names a parameter the previous registry does not declare, which the record builder reports as `unknown`.

## Hand-off

KIT-046 takes the relabelled entry point and the retired `limpa` value (the radio's server-side `limpa` branches, `runProteinLimpaRollupApplyStep`, `buildProtTestModeLimpaProteinObject` and the rollup tab's limpa radio) and adds the `DIA_LIMPA` route beside it.
