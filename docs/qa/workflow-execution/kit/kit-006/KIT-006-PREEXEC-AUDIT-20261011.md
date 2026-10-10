# KIT-006 pre-execution audit (2026-10-11)

Mode: `audit_only` (no repairs; the counted members have not run). Scope: the KIT-006 implementation commit `57ad8613` (partitioner `R/utils_workflow_partitions.R`; the P14 loop and admission in `R/utils_workflow_dia_relational_rollup.R`; the sample-filter admission, counts and the two declarations in `R/utils_workflow_dia_qc_rollup.R`; the store handle's `manager` in `R/utils_workflow_operations.R`; harness pairs, operation decision and filter verifier in `tools/profiling/workflow_kit.R` and the coordinator; generator and bindings; tests; the in-process screen check and the declared runs in `outcome.md`) against protocol v1.3 and the plan ledger (KIT-006: four counted members, 10 min each).

Mechanical inventory: `audit_inventory.py` exit 0, 0 findings. Working tree clean at the audit. `NAMESPACE` and `man/` unchanged (every new function is `@noRd`).

## What was verified

| Check | Result |
|---|---|
| P14 loop equivalence (diff of `persistDIARelationalRollup`) | the per-batch SQL, the stage file naming, `record_stage`, the row count, the engine close and the stage reads are byte-identical statements moved into the partitioner's `stage`, `between` and `consume`; the read scope is an environment, so `scope$connection <- NULL` inside `between` has the effect the original assignment had; `batch_sizes` equals the former `lengths(selections)`; an empty input still yields one empty batch (`FALSE` filter) as before |
| P14 admission equivalence | `completeGroupBatchBound` is the former `batch()` (sum of the largest `ceiling(n/b)` values); the candidate counts, the first-fit choice, the refusal class and message, and every returned field (names, order, the `estimates` list with NULL gaps) are unchanged; `test-workflow-bulk-relational.R`'s 3/6/32-batch cases and the cap-boundary cases pass |
| Sample filter | `present` is the design order filtered to observed samples as before; the per-batch identity SQL is unchanged; `counts$identities` filled by `<<-` from the stage closure; the floor of three is the candidate minimum; the refusal class and message unchanged; `test-perf-106.R` 193 and `test-perf106-repair-v2.R` 75 pass, including the test that asserts three count batches on the small fixture |
| Declarations | both validate against the class registry; `runProjected` runs each kernel through the reference entry point on the small fixture with results identical to the native method / the R rollup kernel (`test-workflow-partitions.R`) |
| Harness | pinned adapters test passes `request$parameters` to the action and checks the filter against them (lines 274, 419, 476), so the bound cutoff and inclusion list reach the member; `variant` `p14`/`physical` by operation; the oracle key resolves to `filter_primary`; `kitOperationDecision` reads `sampled_anon_bytes` and `execution$admission$batches`, both present in operation receipts; harness tests 160 |
| In-process screen check | native and physical identical on all seven fingerprint fields at 3 and at 4 batches on the 8-sample pre-filter state; the evidence store untouched (snapshot store) |
| Declared runs | install, four members, ceiling, order, slice policy, predeclared decision and batch rules, attempt budget (one launch) are stated before execution |

## Findings

| ID | Severity | Layer | Observed | Impact | Smallest repair | Blocked |
|---|---|---|---|---|---|---|
| PE-01 | low | ticket metadata | `files_to_modify` omits `R/utils_workflow_operations.R` (the store handle gained `manager`), the two harness files, the generator, `bindings.json` and `test-workflow-kit-harness.R`; the ticket's authority line still names the control as `core_018` while the plan, protocol and bindings make it `kit_001`. | Allowlist discipline and a stale control name; no execution impact. | Ticket revision 1.1 with the file list and the control. | none |
| PE-02 | low | tracked authority | `kit-001/step-map.json` (hashed into the bindings) still describes the sample step's reference path as "`prepareDIASampleFilterCounts` (:1121, fixed 3 batches)". | The unit gate text is stale after this ticket. | One text edit and a bindings regeneration (the step-map hash is not a protocol hash). | none |
| PE-03 | info | provenance | `in-process-screen-check.json` records `head_commit 976222c8`: the check ran on the uncommitted KIT-006 tree. The R/ files were final at that point (the only later edit was a test regex), so the code that ran is the code at `57ad8613`. | Provenance sentence only. | State it in the outcome, or rerun the 10-second check on `57ad8613`. | none |
| PE-04 | info | design divergence | The partitioner's shape differs from the ticket's sketch: phases stay in the callers (P14 wraps its own labels), the admission is a separate chooser (`chooseCompleteGroupBatches`), and the two declarations run through `runProjected` with the store handle's `manager` rather than through `runSavedOperation` (whose identity/cache path is the whole-object runner's). The outcome describes what was built but does not name the divergence. | Readers of the ticket expect `partitionBySamples(scope, declaration, stage, consume, cap, phase_labels)`. | One sentence in the outcome's implementation summary. | none |
| PE-05 | info | admission semantics | The sample filter's estimate at three batches is now the bound over the largest `ceiling(N/3)` samples instead of the state of the actual three batches, so it is never smaller than before (on the screen state: 1.32 GB at a 12 GiB cap). At `full_primary` (6.2 million rows, about 430 bytes of identity state per row) the three-batch estimate is about 3 GB, far under the cap, so the predeclared "three at the normal cap" rule holds. | None expected; recorded so a surprise batch count can be traced. | None. | none |
| PE-06 | low | verifier coverage | `kitFingerprintFilterResult` reproduces the CORE recipe (`verifyPerf106Result`) and was exercised by the screen check, but no test compares it against a known CORE record. Its first real check is `control-2` (kit_001, physical) against the native oracle; CORE's repair-v2 confirmation showed physical and native identical on that population, so a mismatch there would be a harness defect, which the ledger allows one correction for. | A verifier recipe error would surface on the control member, not the candidate. | Accept as the control's job, or add a fixture test later. | none |
| PE-07 | info | verifier cost | The sample-filter verifier reads the 6.2 million-row output twice (table fingerprint, then the oracle recipe), outside the measured slice. | Time only. | None. | none |

## Closeout states

- `implementation_verified`: true (loop and admission equivalence by diff and tests; declarations exercised; harness paths traced to the pinned adapters test).
- `execution_complete`: false (install and four members declared, not run).
- `decision_supported`: not yet (the decision rules are predeclared and tested on synthetic receipts).
- `promotion_authorized`: false.
- `claim_scope` so far: physical-mode sample filter exact against the native method on the 8-sample screen state at 3 and 4 batches; P14 unchanged on its unit-test populations; no full-size or memory claim.

Next action: none blocking. PE-01 and PE-02 are text fixes that can be done with the closeout; the declared runs may proceed on the user's go.
