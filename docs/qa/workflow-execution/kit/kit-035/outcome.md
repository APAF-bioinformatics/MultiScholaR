# KIT-035: technical replicates declared once

Status 2026-10-16: **complete.** Install `kit_035` recorded (control stays `kit_036`); no counted member; the `journey_da` oracle unchanged.

## What changed (one mechanism)

- `resolveProtTechnicalReplicateBlock(design, group_id, technical_replicate_id)` (`R/func_prot_norm_optimization_helpers.R`) is the one rule: the block is the object's declared technical replicate column (`technical_replicate_id`, which the design step sets to `tech_rep_group` and the class defaults to `replicates`), prefixed by the group only where a value recurs across groups, so the block string equals the `group_replicates` block the DA model built before. A design carries technical replicates when the block has duplicates and is not spanned by the group; otherwise the reason is recorded: `spanned_by_design`, `no_duplicates`, `no_declaration` or `missing_values`.
- `fitProtDaBaseModel` and the DIA execution path (`fitDIAExecutionBaseModel`) take the block from that rule; `runTestsContrasts` carries `group_id` and `technical_replicate_id`, and the five S4 callers (proteomics DA methods, the wrapper, the compat path, metabolomics) pass the object's slots. A spanned block is standard `lmFit` with the reason recorded, never `duplicateCorrelation` with a zeroed correlation; a true technical replicate block runs `duplicateCorrelation` exactly as before.
- `resolveProtRuvReplicateBlocks` reads the same rule, so RUV-III and the DA model report the same block on the same object.
- Recorded: the fit artifact, the saved `da_fit` factors and the DA manifest's `model` carry `group_id`, `technical_replicate_id`, `block_reason`, `uses_duplicate_correlation`, `consensus_correlation` (provenance for KIT-025 and KIT-027; the factors and the manifest are not part of any identity or oracle fingerprint).
- Kernels untouched: `limma::lmFit` and `limma::duplicateCorrelation` are called as before.

## Evidence

| Check | Before | After |
|---|---|---|
| Screen project (KIT-008 member, 8 samples; [da-identity-check.R](da-identity-check.R) → [before](da-identity-before.json), [after](da-identity-after.json)) | limma warning "Block factor already encoded in the design matrix: setting intrablock correlation to zero"; fit factors `uses_duplicate_correlation TRUE`, consensus 0 | no warning; `block_reason spanned_by_design`, `uses_duplicate_correlation FALSE`; **contrast table digest identical** (`e346d5b893f0`, 4,790 rows) and the contrast identity identical (`ede3655a36fc`) |
| 300-sample design (the KIT-036 counted candidate's store; [full-design-check.txt](full-design-check.txt)) | the stored `da_fit` used `duplicateCorrelation`, consensus 0.383, over 102 blocks | the declared block equals the former `group_replicates` string for every sample and the design keeps true technical replicates (`technical_replicates`, 102 blocks), so the fit, its identity and the `journey_da` oracle are unchanged; **no oracle re-recording** |

Tests: `test-prot-07d-da-model-stats-shared.R` 37 pass (new: the declared block, the four reasons, the spanned fixture now standard), `test-prot-ruv-replicate-structure.R` 21 pass (new: both resolvers agree; `group` on a spanned design), `test-workflow-dia-products.R` 149, `test-workflow-saved-da.R` 109 (the exact base-fit resolver carries the provenance factors), norm module contracts 562, DA module 245, metabolomics DA 44.

`R CMD check` at the ticket's commit flagged a code/documentation mismatch for `runTestsContrasts` (its two new arguments); documented at `bd7c591e` ([r-cmd-check-bd7c591e.log](r-cmd-check-bd7c591e.log)): 3 warnings, 3 notes, nothing new. Full suite at `ae87566c` (first run): 343 files, 29,658 pass, 1 fail, 4 errors; besides the PERF-101 runner artefact, four files failed and were classified against the ticket's own change (all four green at the previous commit, `6c9c4901`): three mocked `runTestsContrasts` stand-ins (`test-prot-da-compat-shared.R`, `test-metab-02q-da-analysis-helper-characterization.R`, `test-metab-s4-da-methods-direct-shared.R`) had fixed signatures and refused the two new arguments (now `...`), and `test-tech-reps-limma.R` expected `duplicateCorrelation` on a design with one biological sample per group, which the design spans (the case this ticket makes standard); its fixture now has two biological samples per group, each run twice, so the expectation holds on an estimable design. Full suite at the final commit `5a0d32b3` ([logs](check-logs/)): 343 files, 29,689 pass, 0 fail, 140 skip; the one error is the PERF-101 runner artefact. `R CMD check` 3 warnings, 3 notes.

## Claim limits

The declaration is what the design step wrote; on the synthetic corpora `tech_rep_group` is `group_replicates`, so the 300-sample corpus is modelled with technical replicates as before, which is the declaration's statement, not this ticket's. Whether a study's `replicates` are technical is the user's design decision, now visible in the record.
