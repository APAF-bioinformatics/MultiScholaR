# KIT-028: replay a record, the completeness proof

Status 2026-10-18: **complete.** Install `kit_028` (`run_workflow_kit.sh install`): installed `source_commit` `e28c3b28…` verified, the replay function present; bindings regenerated (`aef82ea2`; control `kit_036` unchanged). Control `kit_036` (unchanged; no counted member).

## Origin (rev 1.2, user decision 2026-10-18)

No existing store recorded step provenance, so one single journey member ran on `kit_027` (`QY8ksX`): stage `kit-journey-028`, selector `kit:journey-028:candidate:1`, screen corpus, no launch record, import through DA, ceiling 900 s bound per ticket, one attempt.
- **Result:** complete in 649 s, whole-slice anonymous peak 2.11 GiB; decision `accepted` (exactness and completion, self-consistency on the screen corpus). It measures nothing.
- **Its record:** every step names its registered operation and variant, no declared parameter is `missing`, correlation shows `preceded_by` `pearson-pairs`, and `analyses$differential_abundance` is present. Only the standing external-resources line is incomplete. This is the first real-journey confirmation of the KIT-027 audit repairs.

## Replay scope: the fixture's states

17 journey steps, 16 distinct states (`session_export` writes no state; the DA step commits its arguments onto `correlation_filtered`):

| State | Operation | Role |
|---|---|---|
| `initial` | `import` | precondition |
| `raw_data_s4` | `design_builder` | precondition |
| `qvalue_filtered` | `peptide_qvalue_filter` | replayed |
| `protein_peptide_filtered` | `peptide_protein_evidence_filter` | replayed |
| `precursor_rollup` | `precursor_rollup` | replayed |
| `intensity_filtered` | `peptide_intensity_filter` | replayed |
| `sample_filtered` | `minimum_peptides_per_sample` | replayed |
| `replicate_filtered` | `peptide_replicate_filter` | replayed |
| `imputed` | `peptide_imputation` | replayed |
| `protein_s4_created` | `protein_iq_rollup` | replayed |
| `protein_accession_cleaned` | `protein_accession_cleanup` | replayed |
| `protein_intensity_filtered` | `protein_intensity_filter` | replayed |
| `duplicates_removed` | `protein_duplicate_removal` | replayed |
| `protein_replicate_filtered` | `protein_replicate_filter` | replayed |
| `ruv_corrected` | `normalization-ruv` | replayed |
| `correlation_filtered` | `pearson-selection` | replayed |

The DA step is replayed after the states through `runDIAReferenceDaAnalysis`, with the record's `analyses$differential_abundance` settings and contrasts (publication deferred). Its argument commit moves `correlation_filtered`, which is then compared in its final form.

## What changed (one mechanism: replay)

- **Replay entry point.** `replayWorkflowAnalysisRecord(record, project_root, inputs, paths, scratch)` (`R/utils_workflow_analysis_record.R`):
  - Refuses a record whose input files differ from `inputs[].sha256` before opening the store (`workflow_replay_inputs`).
  - Verifies the preconditions: the store holds exactly the record's `import` and `design_builder` states, with their entity ids (`workflow_replay_precondition`).
  - Runs every later step in `state_history` order with its recorded operation, variant and parameters:
    - worker declarations: registry handler, then the controller's completion on one workflow owner; each peptide step's filtering progress is advanced as the GUI display does (`updateProteinFiltering`), since the next step's progress product records it;
    - normalisation and correlation: `runSavedOperation` and the commit their reference wrappers make;
    - DA: as above.
  - Stops on `legacy`, `import` or `design_builder` out of position, or an unregistered operation (`workflow_replay_unsupported`, naming the step).
- **Harness.**
  - `kitPrepareReplayProject` cuts a copied project's store head to the precondition states. The reopened index (rebuilt from the head's roots) then finds no saved result of a later step, and the function refuses if any scientific product is still indexed. So no replayed step can reuse the origin's saved results.
  - `kitReplayJourney` replays and compares every state with the origin through `kitReproduceMismatches` (`head_sha256` excluded).
  - `kitMemberSeconds` reads a ticket's own journey ceiling.
- **Record fix found while replaying.** A product's `product_type` read a field products do not carry (always NA); it now reads the identity's selector.

## Evidence

| Check | Result |
|---|---|
| Replay check on the origin ([replay-check.R](replay-check.R) → [replay-check-screen.json](replay-check-screen.json); origin record [analysis_record-origin-screen.json](analysis_record-origin-screen.json)), installed `kit_027` library, in process | **Exact:** 16 states, 0 mismatches against the origin (319 s). **Edited** (minimum peptides per protein 1 → 16): identical before the filter; at the filter 281,343 rows, equal to the direct kernel call (`filterMinNumPeptidesPerProtein` on the origin's Q-value state), against the origin's 320,872; every later state differs (260 s). **Inputs:** a different FASTA refused before any step. **Legacy:** the Q-value step marked `legacy` stops with `workflow_replay_unsupported` naming `qvalue_filtered`. Total 594 s (ceiling 900 s). |
| `test-workflow-analysis-record.R` (new: input refusal before the store opens, a missing input, a legacy step named, a moved import refused as a precondition) | 69 pass |
| `test-workflow-kit-harness.R`, `test-workflow-operation-registry.R` | 340, 219 pass |
| Test rule (user, 2026-10-18: targeted files, dependents, check, install; full suite at checkpoints) | Dependents of the changed functions (`buildWorkflowAnalysisRecord`, `writeWorkflowAnalysisRecord`, `replayWorkflowAnalysisRecord`, `kitMemberSeconds`, `kitPrepareReplayProject`, `kitReplayJourney`): only the two files above. `R CMD check` at `e28c3b28`: 3 warnings, 3 notes (baseline). Checkpoint full suite at `e28c3b28` ([logs](check-logs/)): 344 files, 29,910 pass, 0 fail, 140 skip; the one error is the known PERF-101 runner artefact. This is the checkpoint for KIT-029 onward. |

**Attempts of the replay check** (in process, not counted; recorded, not retried for a result):
1. The exact replay passed. The edit then failed: the ticket's Q-value edit (0.01 → 0.001) is refused by the package ("DIA-NN identification q-value thresholds are fixed at 0.01"), a science-contract rule, so the ticket's edit cannot exist (rev 1.3).
2. Exact passed again. The edit at 3 peptides per protein removed nothing: the fixture has at least 14 peptides per protein (the direct kernel call agreed, 320,872 rows).
3. The recorded run: the edit at 16, with its effect checked by the kernel call before replay.

A bug found by the new tests (an input role absent from `inputs` raised a subscript error instead of the refusal) was fixed before the recorded run.

## Claim limits

- Replay proves the record complete for the steps it replays, on this fixture: every state is reproduced from the record and the hash-checked inputs.
- The import and the GUI-built design are verified preconditions, not replayed.
- Replay deliberately reproduces the GUI's progress bookkeeping (the filtering progress the display advances) so that the progress products match. Without it, those products would differ while the science would not.
- Enrichment is not part of this journey.
- A record from a store committed before KIT-025 is all `legacy` and is refused.
