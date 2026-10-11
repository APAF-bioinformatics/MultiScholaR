# KIT-006: complete-sample partitioner

Status 2026-10-11: **accepted and complete.** Stage `kit-006` decision `accepted` (both pairs exact against their oracles, each candidate within 1.02 × its control); `kit_006` (`core/browser-smoke.qBR5tm`, source `62fb9aec`) is the control for KIT-007 onward. Pre-execution audit: [KIT-006-PREEXEC-AUDIT-20261011.md](KIT-006-PREEXEC-AUDIT-20261011.md) (fixed before launch). Authority: `tickets/kit/00-KIT-PLAN.md` (run ledger KIT-006: four counted operation members, 10 min each); protocol v1.3. Control: `kit_001` (`core/browser-smoke.X9aj9P`, source `8e77d92d`), the series control since the user's decisions of 2026-10-09/10 (the ticket text's `core_018` is superseded by the plan).

## What was built

- `R/utils_workflow_partitions.R` (new): `partitionBySamples(groups, batches, stage, between, consume)` visits complete groups in the declared order, never splitting a group, with the step's own kernels staging and consuming each batch; `chooseCompleteGroupBatches(counts, estimate_for, worker_cap, policy)` picks the fewest batches whose estimate fits the admitted cap and refuses above it; `completeGroupBatchBound(values, batches)` is the per-batch upper bound (the largest `ceiling(n / b)` groups) both estimates use. No scientific table, column or operation is named in the file (test).
- **P14 on the partitioner.** `persistDIARelationalRollup`'s batch loop is the partitioner's `stage` (ordered aggregate to one parquet stage per batch), `between` (engine close) and `consume` (stage read into the writer); `estimateDIARollupWorkingMemory` keeps its CORE-010 calibration (`diaRollupRowBytes = 64`) and hands the choice to the generic admission. Phases (`ordered_result_stage`, `query_engine_close`, `staged_result_read_summary_write`, `final_writer_finish`) and every `execution` field are unchanged; `batch_sizes` comes from the partitioner. `diaPrecursorRollupOperation()` declares the step (unit `complete_sample`, three global statistics, merge rules; kernel = the reference entry point whose aggregate is the ordered-sum SQL).
- **Sample filter on the partitioner.** `prepareDIASampleFilterCounts` (physical and selection modes) replaces its fixed three batches with the adaptive admission (floor three, `sample-adaptive-complete-batches-and-expanded-transfer-v4`): per-sample identity state bounded by the largest samples a batch can hold; the identity counts run through `partitionBySamples`. `diaSampleFilterOperation()` declares the step (`minimum_peptides_per_sample`, unit `complete_sample`, the global design/inclusion statistic, merge rules; kernel = the reference entry point in physical mode). The store handle of the operation runner exposes `manager` so both declarations run through `runProjected` (test). The worker's production default (`legacy`) is untouched (KIT-008).
- Harness: a ticket binds one or two (operation, population) pairs (`kit-006`: P14 on `full_p14`, the sample filter on `full_primary`; members `control-1`, `candidate-1`, `control-2`, `candidate-2`); `kitOperationDecision` (every member eligible and verified; per pair the candidate's kernel peak and anonymous peak ≤ 1.02 × the control's; a reduction ≤ 0.95 reported); `verifyKitOperationMember` verifies the sample filter against CORE's native oracle record (`oracles.filter_primary`: page hashes, rows, samples, design, summary and support) with `kitFingerprintFilterResult`, the oracle's own recipe. Harness tests 160.

Divergence from the ticket's sketch (pre-execution audit PE-04): the partitioner is `partitionBySamples(groups, batches, stage, between, consume)` plus a separate admission chooser, not `partitionBySamples(scope, declaration, stage, consume, cap, phase_labels)`; the phase labels stay in the callers (P14 wraps its own), and the two declarations run through `runProjected` with the store handle's `manager` rather than through `runSavedOperation`, whose identity and cache path belongs to the whole-object steps. One mechanism less; no behaviour differs.

## Unit gate

`kit-001/step-map.json`: `rollup` and `sample` carry unit `complete_sample`, grouping `sample_id`, no open flag (bound by hash in the bindings).

## In-process exactness before the full members

[in-process-screen-check.R](in-process-screen-check.R) → [in-process-screen-check.json](in-process-screen-check.json) (rerun on the committed tree after audit PE-03): the KIT-001 screen journey's `intensity_filtered` state (8 samples, 226,991 rows; materialised as a snapshot store, the evidence store untouched), cutoff 3000: native mode and physical mode agree on every page hash, row and sample count, design, summary and support (7 of 8 samples retained); physical took 3 batches (2/3/3) at the normal cap and 4 (2/2/2/2) at a cap set between the 4- and 3-batch estimates, identical both times. Unit tests: the small fixture takes 3 batches at the normal cap and 4 under a reduced cap, exact against the native method either way; P14's 3/6/32-batch cases unchanged (`test-workflow-bulk-relational.R`).

## Tests (source tree, the qualified fastduckdb library on the path)

`test-workflow-partitions.R` 61; `test-perf-106.R` 193; `test-perf106-repair-v2.R` 75; `test-workflow-bulk-relational.R` 45 (3 explicit-preflight skips as before); `test-workflow-saved-qc-rollup.R` 76; `test-workflow-operations.R` 66; `test-workflow-kit-harness.R` 160.

## Results (stage `kit-006`; tracked summary [kit-006-stage.json](kit-006-stage.json))

Install: `kit_006` = `core/browser-smoke.qBR5tm` from `62fb9aec`, empty runtime patch; the installed `gui-bulk-rollup:writer-r-preflight` on it (`core/browser-smoke.095yCi`) passed 174 checks (0 failed, 2 pre-existing skips) with the five exact P14 cases, closing KIT-005 audit PH-10.

| Member | Build | Operation / population | Status | Elapsed | Slice peak | Anonymous peak | Max gap | Batches | Verified |
|---|---|---|---|---|---|---|---|---|---|
| control-1 | kit_001 | P14 on `full_p14` (24,000,000 rows, 300 samples) | complete | 75 s | 6.607 GiB | 5.254 GiB | 80 ms | 3 × 100 | exact (1,200 pages) |
| candidate-1 | kit_006 | same | complete | 73 s | 6.547 GiB | 5.232 GiB | 70 ms | 3 × 100 | exact (1,200 pages) |
| control-2 | kit_001 | sample filter, physical, on `full_primary` (6,225,000 rows, cutoff 3000, 75 included runs) | complete | 30 s | 1.584 GiB | 1.332 GiB | 60 ms | 3 × 100 (fixed) | exact: 615 pages, 6,150,000 rows, 225 samples, design, summary, support (after correction 1) |
| candidate-2 | kit_006 | same | complete | 32 s | 1.578 GiB | 1.324 GiB | 70 ms | 3 × 100 (adaptive) | exact, same record |

Pair ratios (candidate / control): P14 slice peak 0.991, anonymous 0.996, elapsed 0.98; sample filter slice peak 0.996, anonymous 0.994, elapsed 1.09. Both `not_worse`; no reduction claimed (none expected: the mechanism moved, the batches did not). Predeclared batch rules held: P14's count equals the control's (3); the sample filter's adaptive count is 3 at the admitted cap (14.3 GiB; estimate 4.07 GiB against the control's 3.37 GiB fixed-policy estimate). Zero swap, every gap under 250 ms, every slice stopped. Host: `low-power` profile, CPU idle at start.

### Attempts and corrections (no member measurement rerun for a correction)

- **Attempt 1** (`8d8efccc`): the P14 pair completed and verified. Both sample-filter members were refused in 3 s by the pinned action validator ("inclusion entries must be character"): the bindings are read without simplification and the 75-run inclusion list reached the action as a list of strings. No work had run on either build. Harness fix: case parameters are coerced to atomic vectors. The user amended the ledger for a rerun of those two members only.
- **Rerun 1** (`run_workflow_kit.sh rerun kit-006 … control-2,candidate-2`, new runner mode): inside the finished stage, the attempt 1 directories kept as `control-2.attempt1` / `candidate-2.attempt1`, the corrected harness snapshotted as `source.rerun1`, the P14 receipts untouched, the decision recomputed over all four. Both members completed; candidate-2 verified.
- **Correction 1** (verifier only, [control-2-verifier-correction.R](control-2-verifier-correction.R)): the new sample-filter verifier required `admission$batches`, which only the adaptive admission writes; the control's fixed-three policy records `count_batch_sizes`. control-2's failed verification is kept as `verification.harness-failed.rds`; the verifier reran from the corrected harness (fallback to the count of `count_batch_sizes`) and accepted; decision `accepted`.

### Claim limits

- Memory decided on the whole-slice peak and the sampled anonymous peak of one member per build per operation (the protocol's pair rule); the two operations run under 7 GiB and 2 GiB respectively, far from the 16 GiB cap, so this ticket proves the mechanism and exactness rather than a reduction.
- The P14 receipts of attempt 1 predate the `science_execution` field in the operation receipt; their batch counts come from the verification records.
- Times are reported, not judged.

## Declared runs (as authorised; recorded above)

1. **Install** (recorded, not counted): the KIT-006 tree at its commit into `<proof>/core` with `run_workflow_checkpoint_proof.sh --service install-corrected <kit_001 library> <proof>/core`, bound as `kit_006`; bindings regenerated with `tickets$006`; `run_workflow_kit.sh validate kit-006` must pass. The adapters installed preflight deferred from KIT-005 (audit PH-10) runs on this install through the recorded preflight service.
2. **Stage `kit-006`, four counted operation members, 600 s ceiling each**, in order `control-1` (kit_001, P14 on `full_p14`, 24,000,000 rows), `candidate-1` (kit_006, same), `control-2` (kit_001, sample filter physical on `full_primary`, 6,225,000 rows, cutoff 3000 with the bound inclusion list), `candidate-2` (kit_006, same). Each in its own slice (16 GiB cap, 15.5 GiB early stop, 12 GiB start headroom, 2 GiB host floor, 250 ms gap rule); the verifier runs outside the slice after each member.
3. **Decision (predeclared):** both pairs exact against their oracles (P14 pages; the native filter record) and each candidate's kernel peak and anonymous peak ≤ 1.02 × its control's → `accepted`, the kit_006 build becomes the control. Any rule failing → `rejected` or `invalid` (instrumentation), recorded; one correction rerun of a single member only for a demonstrated harness failure. Predeclared batch rules: P14's batch count at the normal cap equals the control's; the sample filter's is ≥ 3 (the control's fixed three).

Attempts: one stage launch; no retry of a valid unfavourable result.
