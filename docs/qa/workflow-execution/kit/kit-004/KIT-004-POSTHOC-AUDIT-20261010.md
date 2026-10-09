# KIT-004 post-hoc audit (2026-10-10)

Mode: `post_hoc`, audit only (no repairs applied). Scope: the KIT-004 implementation commit `ea2b678d` (runner `R/utils_workflow_operations.R`; declarations in `R/utils_workflow_dia_products.R` and `R/utils_workflow_dia_da.R`; module-edge moves in `R/mod_prot_norm_workflow_helpers.R`, `R/mod_prot_norm_ruv_helpers.R`, `R/mod_prot_da_handlers_run_analysis.R`; tests; the in-process comparison and identity-continuity evidence under this directory) against protocol v1.1 (sha256 `8998db8f…`) and the plan ledger (`tickets/kit/00-KIT-PLAN.md`, KIT-004: no counted member).

Mechanical inventory: `audit_inventory.py` exit 0, 0 findings (ticket schema valid, index consistent, protocol hash matches, no missing dependencies). Working tree clean at the audit.

## Evidence layers

| Layer | Status | Evidence |
|---|---|---|
| mechanical | pass | inventory; `DESCRIPTION` Collate entry; `NAMESPACE` and `man/` unchanged; grep of the three runner files clean |
| semantic | pass with notes | the identity scheme is called verbatim (`diaComponentScientificIdentity`, names, parameter order, DA factors); kernels unchanged and called through the legacy aliases; store-handle checkpoints go through `saveDIAOwnedValue` → declaration-driven `checkpointDIAExecutionResult` |
| empirical | pass with limits | two baseline passes agree except the composite QC PDF; baseline vs post identical on every state, 31 products, DA manifest, 28 holder fields and 58 of 59 files; the KIT-001 member's own products hit under the runner with its software identity substituted; 66 runner assertions; affected regressions green except the known process-loss artefact and two module-CI errors reproduced at baseline |
| promotion | not claimed | no counted member; `kit_001` stays the control |

## Findings

| ID | Severity | Layer | Observed | Claim affected / impact | Smallest repair | Blocked |
|---|---|---|---|---|---|---|
| PC-01 | medium | tests | `test-workflow-operations.R`, the callr fresh-session test, loads the package in the child with `pkgload::load_all(test_path("..", ".."))`. In the installed suites (`load_package = "installed"`, as `test-workflow-execution.R` runs them) and in an `R CMD check` tests directory there is no `DESCRIPTION` at that root, so the test errors instead of running or skipping. | Regression check (the test passes only from a source tree); not a product or science claim | In the child, use the installed namespace when the root has no `DESCRIPTION` (`.libPaths` + `library(MultiScholaR)`), or `skip` with a reason | none |
| PC-02 | low | provenance | `fingerprints-post.json` records `head_commit e2bd5680`: the post pass ran on the uncommitted KIT-004 tree before `ea2b678d`. The code that ran differs from the commit only by the `rootDIAProduct` loop-variable rename and test edits (the kernel alias fix was in). | Exactness evidence provenance (stated in part in the outcome) | Rerun `in-process-comparison.R run post` on `ea2b678d` and replace the file, or add the exact provenance sentence to the outcome | none |
| PC-03 | low | design | On the reference path the correlation completion (`completeProtNormCorrelationWorkflow`) now runs in the observers after the wrapper's `withDIAStateScope` has closed; the legacy branch runs its completion inside a scope. Its reads do not need a scope (the in-process run and the 425 native-session assertions pass) and its conditional commit is a no-op after the wrapper's commit. | Working-cache semantics only; no evidence of a behaviour change | Wrap the reference branch of the two observers in `withDIAStateScope(getDIAActiveManager(), …)` like the legacy branch | none |
| PC-04 | low | behaviour | The module edge writes the RUV audit trail into `config_list` on cache hits as well as misses, and for manual mode takes the percentage from the settings rather than the computed result (both documented in the outcome). No test and no fingerprint covers the `config_list` write (the comparison does not fingerprint `config_list`). | Side effect outside the exactness claim | One assertion on `config_list$ruvParameters` after the module-edge call in a normalisation module test, or accept as recorded | none |
| PC-05 | info | duplication | `resolveProtRuvParameters` (kernel-only) repeats the percentage-range and `ruv_k` defaults of the adapter `resolveProtNormRuvParameters`, which the ticket keeps unchanged. | Drift risk between the two paths until the legacy path is removed | Have the adapter delegate when KIT-011/012 retire the legacy path; record in the hand-off | none |
| PC-06 | info | process | The ticket budget names one recorded comparison. Three passes completed (pre-a, pre-b, post) and two aborted (a script crash before any fingerprint; a first post pass that exposed the kernel-alias argument names, followed by one bounded correction). pre-b was added to establish baseline determinism. | Within the one-correction rule; the extra baseline pass was a design decision taken while running | Declare "two baseline passes" as part of the comparison design next time, before running | none |
| PC-07 | info | coverage | The `workflow` holder fingerprint covers five fields; `ruv_normalised_for_da_analysis_obj`, `design_matrix` and `config_list` are not fingerprinted; the composite QC PDF is excluded because it differs between two baseline runs. | Claim scope ("28 holder fields, 58 of 59 files") is stated correctly | None required; extend the holder list if a later ticket relies on those fields | none |
| PC-08 | info | compatibility | The wrapper names are kept but their signatures changed (`notify`, `output`, `session` removed; `contrasts`, `annotations` added); all four internal call sites updated; `tools/profiling/workflow_checkpoint_gui_diagnostic.R` references the names only for tracing. | Internal (`@noRd`) functions; public API unchanged | None | none |

## Closeout states

- `implementation_verified`: true (code maps to every ticket outcome; grep clean; tests present).
- `execution_complete`: true (no counted run required; the recorded comparison and continuity check ran to completion).
- `decision_supported`: true for the exactness claim as scoped (states, products, DA manifest, 28 holder fields, 58 of 59 files; PDF excluded with baseline evidence).
- `promotion_authorized`: false (not applicable; no control change).
- `claim_scope`: products, checkpoints, publication files and holder fields hash-identical to the pre-ticket build on the screen fixture, in-process, with the stated exclusions; identity continuity proven with the member's software identity substituted (an installed build is required for a native hit).

Next action: no blocker for KIT-005. PC-01 should be fixed before the next installed suite run; PC-02 is a four-minute rerun.
