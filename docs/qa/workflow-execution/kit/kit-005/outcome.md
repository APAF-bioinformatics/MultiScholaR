# KIT-005: harness, operation members, per-step oracles, baseline journey

Status 2026-10-11: stage decisions recorded. Authority: `tickets/kit/00-KIT-PLAN.md` (run ledger KIT-005 with the user's amendments of 2026-10-10); protocol v1.3 (`../protocol.md`, sha256 in `../bindings.json`). Control: `kit_001` = install `core/browser-smoke.X9aj9P` (source `8e77d92d` on branch `kit-control-ceiling`: the KIT-001 control source `9075b361` plus three user decisions of 2026-10-10, below).

## Results

| Stage | Member | Decision | Elapsed | Whole-slice peak | Anonymous peak | Resident end | Evidence |
|---|---|---|---|---|---|---|---|
| `kit-lifecycle-operation-complete` | precursor rollup on `screen_p14`, installed action dispatch | `lifecycle_complete_observed` (every P14 phase complete, pages equal to the P14 oracle) | 38 s | 1.42 GiB | — | — | [operation-complete-verifier-correction.R](operation-complete-verifier-correction.R) |
| `kit-lifecycle-operation-fail` | same, worker fails inside the ordered result stage | `lifecycle_failure_observed` (head unchanged on the input component) | 28 s | 1.36 GiB | — | — | [operation-fail-verifier-correction.R](operation-fail-verifier-correction.R) |
| `kit-journey` control-1 | 300-sample journey, `kit-journey-v1` record | complete, verified, 18 populations in-run | 3,847 s | 7.89 GiB | 5.67 GiB | 3.03 GiB | [journey-attempt6.json](journey-attempt6.json) |
| `kit-journey` control-2 | same | complete, verified, **reproduced every fingerprint and file of control-1** | 3,941 s | 7.61 GiB | 5.51 GiB | 2.98 GiB | same |
| `kit-journey` default-1 | production default (no record, in-memory `WorkflowState`) | **memory guard in imputation after 8 of 18 steps** (valid recorded result) | 643 s | 15.51 GiB | 14.32 GiB | not reached | same |

Stage decision `baseline_recorded` under protocol v1.3 ([journey-decision-v13.R](journey-decision-v13.R)); noise between the two control members: anonymous peak ratio 0.971, resident end 0.983, elapsed 1.02. Control-1 had one 280 ms sampling gap in 66,823 samples (RUV); control-2 none above 250 ms.

Per step (seconds; anonymous peak of the whole slice in GiB) — control-1 / control-2 / default-1:

| Step | control-1 | control-2 | default-1 |
|---|---|---|---|
| import | 31 / 2.87 | 50 / 2.65 | 29 / 2.73 |
| design | 129 / 3.29 | 139 / 2.91 | 113 / 3.52 |
| q-value | 288 / 3.70 | 259 / 3.05 | 47 / 6.39 |
| protein evidence | 169 / 2.58 | 180 / 2.54 | 8 / 6.63 |
| rollup | 98 / 2.70 | 150 / 2.56 | 70 / 9.75 |
| intensity | 113 / 2.66 | 194 / 2.60 | 59 / 9.25 |
| sample | 53 / 2.36 | 88 / 2.16 | 12 / 9.25 |
| replicate | 112 / 2.61 | 169 / 2.54 | 64 / 9.25 |
| imputation | 592 / 3.02 | 589 / 2.87 | 220 / 14.32 (guard) |
| IQ | 231 / 3.10 | 232 / 2.96 | — |
| accession, protein intensity, duplicates, protein replicate | 33, 34, 38, 34 / 2.7 | 30, 32, 35, 31 / 2.6 | — |
| normalisation + RUV | 1,508 / 3.36 | 1,395 / 3.13 | — |
| correlation | 40 / 2.79 | 36 / 2.67 | — |
| session export | 40 / 2.37 | 45 / 2.24 | — |
| DA | 160 / 5.67 | 147 / 5.51 | — |

The in-memory default's per-step peaks (q-value 6.4, rollup 9.8, imputation 14.3 GiB) against the native path's 2.4 to 3.7 GiB on the same steps are the memory baseline KIT-011, KIT-012 and KIT-017 are judged against; the default path cannot complete a 300-sample journey inside the 16 GiB cap. On the native path the DA step holds the journey's anonymous peak (5.5 to 5.7 GiB: the full protein object plus limma), and normalisation with RUV is the longest step.

## Populations and oracles

`kit/bindings.json` now binds, from control-1, one population per journey step (`journey_<step>`: store, head sha256, entity id, rows) and its oracle (every component, support entry and product fingerprint and the args summary, canonical across runs), 17 steps (the session export leaves no state) plus the three `core` populations. Control-2 reproduced all of them (0 fingerprint mismatches over 18 steps, 0 file mismatches with the canonical file rule). Rows: peptide states 5,846,578 to 5,515,154; `imputed` 11,826,409 (the imputation step writes about twice the rows of its input; recorded for a later ticket); protein states 2,500; `ruv_corrected` and `correlation_filtered` 2,476.

## Three user decisions of 2026-10-10 (science and protocol)

1. **Protocol v1.2, journey ceiling up to 4 h** (attempt 1 showed the data steps scale about 20 x from 8 to 300 samples); the launch record's limit raised to 14,400 s (one constant in `workflowGuiLaunchRecord`), control rebuilt.
2. **RUV sweep: imputation once per sweep** (`prepareRuvCancorInput` / `ruvCancorFromInput`; attempt 3 measured about 970 s per percentage because `ruvCancor` re-imputed the full matrix for each of the 20 percentages; the sweep now takes about 12 min). Exactness test `test-prot-ruv-sweep-imputation.R`.
3. **RUV sweep k validity and RUV-III guards** (attempt 4: the sweep chose 1 % controls with k = 15 while only 3 of those 25 controls were observed in every sample, and RUV-III returned no finite values, which the kernel passed on as an empty table): a percentage whose elbow k exceeds its controls observed in every sample is skipped (`skipped_k_exceeds_complete_controls`); `ruvIII_C_Varying` refuses such a k (`prot_ruv_insufficient_controls`) and an empty result (`prot_ruv_empty_result`). Offline check on the saved 300-sample object ([ruv-sweep-check-300-samples.txt](ruv-sweep-check-300-samples.txt)): percentages 1 to 11 skipped, winner 19 % with k = 44; RUV-III 2,500 x 298. Science note for a later ticket: at 300 samples the elbow k is 44 to 52 for every percentage, every composite score saturates to zero and the winner is decided by the tie-break order.

Plus protocol v1.3 (2026-10-11): journey members tolerate isolated sampling gaps up to 500 ms.

## Harness corrections (member measurements never rerun for a correction)

| # | Member | Cause | Fix |
|---|---|---|---|
| 1 | operation complete | verifier compared page hashes against the oracle's page-count key | `*_page_hashes`; verifier rerun |
| 2 | operation fail | verification record saved from the pre-return variable; head check compared the state id with the input component id | fixed; verifier and head check rerun |
| 3 | journey attempt 1 | pinned browser wait budget tied to the coordinator's elapsed clock (900 s) | per-member driver process (`workflow_kit_member.R`) |
| 4 | journey attempt 2 | browser transport timed out while the tab rendered the 2,500 x 298 protein table | evaluations retried while the renderer is frozen |
| 5 | proactive | pinned 2 h process cap on waits; 10 s completion notifications; fingerprint pass inside the measured window | wait budget = member ceiling; persistent notification log; teardown before fingerprints |
| 6 | journey attempt 5 | table reader's 64-column cap refused the wide protein tables; six process spawns per sample left gaps above 250 ms | column cap = table width; one awk per sample (56 ms mean, 60 ms max) |
| 7 | journey attempt 6 control-2 | between-member file check compared raw paths (timestamped directories, request ids, export bookkeeping) | canonical paths; export bookkeeping and PDFs by path; verifier rerun |

Launch record: journey attempts 1 to 5 are archived under `<proof>/kit/runtime/kit-journey.attempt*` with the causes in `<proof>/kit/kit-journey-notes.txt` (attempt 3's coordinator log was overwritten by the attempt 4 launch; its receipts are intact). Two screen lifecycle members ran outside the ledger at the user's request to prove the corrected driver path; the second completed all 18 steps in 462 s under a 600 s screen ceiling. A stray dry-check process of mine overlapped the first 22 minutes of attempt 1 (elapsed time only).

## What was built

`tools/profiling/workflow_kit.R` (operation cases, `kit-0NN` stages, baseline and production-default journey cases, canonical cross-run fingerprints with memoised record hydration, per-step products, reproducibility and canonical file checks, operation and journey decisions, generalised verifier), `workflow_kit_coordinator.R` (operation-member launch through the installed adapters dispatch; per-member driver process), `workflow_kit_member.R` (new), `workflow_kit_driver.R` (patient driver, notification log, teardown before fingerprints), `run_workflow_kit.sh` (stages, one-awk sampler), `bindings-generator.R` (pair, ceiling, operation lifecycle, population ingestion). `test-workflow-kit-harness.R`: 114 assertions. No CORE-pinned file edited.

## Claim limits

- The production-default member did not complete: its baseline is "guard at imputation" with per-step peaks through the replicate filter; a later comparison against it is completion plus those per-step peaks, not a whole-journey ratio.
- Control-1's sampled anonymous trace has one 280 ms blind window during RUV (protocol v1.3); the whole-slice peak is the kernel's counter.
- The elapsed times were measured on a shared host (a game ran for part of the day); time is reported, never judged.

## Hand-off

KIT-006 binds `full_p14`/`full_primary` operation members through `kit-006` (`bindings$tickets[["006"]]`); KIT-009/010/014/015 bind the `journey_<step>` populations and oracles; KIT-011/012/017 use the baseline pair and the production-default result.
