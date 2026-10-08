# CORE series closeout (CORE-001 to CORE-017)

Closure source: `6d4f7c22`; the closure commit is the commit that adds this file. Its hash is recorded in `tickets/openai/06-CORE-INFRA-PLAN.md` and `tickets/transcriptomics/00-TRANSCRIPTOMICS-PLAN.md`. Every number below comes from a receipt under the proof root `MultiScholaR-dia-proof-20260911-VhC1n5/core/` or an outcome file in this directory tree. No default activation, real-study or production claim is made.

## Outcome per ticket

| Ticket | What it did | Outcome |
|---|---|---|
| CORE-001 | Codec, exact-sum handle and export publication repairs | done (correctness) |
| CORE-002 | State and rollup science-adjacent repairs | done (correctness) |
| CORE-003 | Execution recovery and leak repairs | done (correctness) |
| CORE-004 | CORE harness, bindings and lifecycle proof | done (harness) |
| CORE-005 | File-hash memoisation and writer metadata caching | **reverted**: failed the per-pair peak rule (1.005 / 0.961) in `core-010-control` |
| CORE-006 | Row-group alignment | closed without change: DuckDB rounds row groups to 10,240 |
| CORE-007 | Byte estimates and lazy strings | closed without change: the peak is not in R, and the strings are dictionary-encoded |
| CORE-008 | Scan, annotation and append consolidation | closed without change: 1–3% of the targeted steps |
| CORE-009 | Integer key encoding | closed without change: about 5% of the grouping step measured end to end |
| CORE-010 | Worker cap follows free host memory; adaptive P14 batches | **accepted by user decision** (budget stage `budget_not_binding`) |
| CORE-011 | Engine spill | closed without change: the ordered grouping cannot spill |
| CORE-012 | Summation contract | decision: keep the exact ordered sum (the compensated sum depends on thread count) |
| CORE-013 | Coordination overhead | closed without change: about 15 ms per poll tick, about 0.2 s of syncs per action |
| CORE-014 | Store garbage collection | closed without change; handed to PERF-107 |
| CORE-015 | Dead code, profiles, options, preconditions, design note | done (trimmed scope) |
| CORE-016 | Closeout | this file |
| CORE-017 | Whole-input checks per batch, plus phase timing | stage invalid (host contention) and reverted; **re-applied in CORE-010**; phase timing kept |

## Measured stages (`full_p14`, 24,000,000 precursors, unless stated)

| Stage | Member | Build | Peak GiB | Anonymous peak GiB | Elapsed s | Result |
|---|---|---|---|---|---|---|
| CORE-004 lifecycle (`screen_p14`) | complete / fail | core_004 | 2.74 / 2.71 | — | 34.6 / 27.0 | observed / observed |
| CORE-017 lifecycle re-proof (`screen_p14`) | complete / fail | core_004 | 2.70 / 2.69 | — | 65.9 / 42.3 | observed / observed |
| core-017 (protocol v3) | control-1, control-2 | core_005 | 7.50, 7.88 | 6.48, 6.67 | 144.5, 147.5 | stage `invalid` (control-2 monitor gap 131 ms, host contention) |
| | candidate-1, candidate-2 | core_017 | 7.33, 7.65 | 6.47, 6.49 | 161.6, 159.3 | checks-phase anonymous peak 4.38/4.40 → 3.18/3.20 GiB |
| CORE-010 lifecycle re-proof (`screen_p14`) | complete / fail | core_004 | 2.73 / 2.70 | — | 34.1 / 28.0 | observed / observed |
| core-010-control (16 GiB) | control-1, control-2 | PERF-108 (`jEmE2F`) | 7.94, 8.20 | 6.61, 6.82 | 118.6, 88.1 | CORE-005 **rejected** (per-pair peak 1.005, 0.961) |
| | candidate-1, candidate-2 | core_005 | 7.98, 7.88 | 6.66, 6.58 | 84.8, 90.5 | |
| core-010 budget (8 GiB free at start) | control-1 | PERF-108 | 7.56 | 6.51 | 95.4 | completed: `budget_not_binding` (start drifted to 9.13 GiB free) |
| | candidate-1 | core_010 | 6.02 | 5.15 | 94.3 | 5 batches, cap 5.32 GiB |
| | candidate-2 | core_010 | 5.64 | 4.77 | 93.6 | 6 batches, cap 4.98 GiB |

## Net change from the PERF-108 closeout

- **PERF-108 `full_p14` receipts** (library `jEmE2F`, quiet host): peaks 7.10 / 7.05 GiB, elapsed 87.5 / 87.0 s.
- **With plenty of free memory,** the final build's batch choice is three, identical to PERF-108, and no stage measured it alone at 16 GiB. The series makes no claim of a lower peak under those conditions.
- **With reduced free memory** (one stage, two runs against one PERF-108 control run in the same session): peak 7.56 → 6.02 / 5.64 GiB (−20% / −25%), anonymous 6.51 → 5.15 / 4.77 GiB, wall time 95.4 → 94.3 / 93.6 s. Time is reported and not judged under protocol v5 and later, and PERF-108 measured +46% elapsed for six batches on a quiet host.
- **`full_primary`** (the sample filter) was not measured in this series; no change is claimed for it.

## Which mechanisms mattered

- **Lowered memory:** adaptive complete-sample batching keyed on free host memory (CORE-010), and the per-batch whole-input checks (CORE-017, re-applied), which lowered the checks-phase floor by about 1.2 GiB of anonymous memory.
- **Did not:**
  - integer keys (about 5%);
  - engine spill (cannot spill the ordered lists);
  - row-group alignment, lazy strings, scan consolidation and coordination overhead (premises false or negligible);
  - the compensated sum (rejected for reproducibility).
- **The remaining floor:** about 2.4 GiB before the rollup and about 3.2 GiB anonymous memory in the checks phase. The single largest remaining lever is the exact ordered sum, which was kept deliberately (CORE-012).

## Unfavourable or incomplete results (not omitted)

- **CORE-005** was reverted by a 0.5% per-pair peak breach that is likely noise. The strict per-pair rule has no noise tolerance, which should be decided before it is used again.
- **CORE-017's stage was invalid** under host contention. Its first reading used cache-inclusive peaks and wrongly concluded "about 5%"; the correction is recorded in its outcome.
- **CORE-010's budget stage never reached its declared budget.** The placeholder is sized once and drifted by about 1 GiB, so the control result is uninterpretable. The drift check is deferred to the stage's next reuse.
- **Protocol revisions:** v3 to v7 were made within one session (user decisions, applied prospectively only).

## Closure checks (CORE-016)

- **Installed check** (closure build `browser-smoke.Hz6YwJ`, guarded `gui-bulk-rollup:writer-r-preflight`, `browser-smoke.iUjeFk`): 174 checks passed, 0 failed; all five exact P14 cases exact; 1.79 GB peak, 0 swap.
- **Legacy checkpoint:** a v1 native session archive written on 2026-09-23 (`browser-smoke.lSCYZW`, 765 MB) opens on the closure build with all 15 states readable as their S4 classes, identical to the PERF-108 build.
- **Public API and S4:** `NAMESPACE` is unchanged across the series (`addcccc0^..6d4f7c22`); no `setClass`, `setGeneric`, `setMethod` or `setValidity` line changed. Saved-table and checkpoint formats are unchanged.
- **Source suite and R CMD check:** the full suite and `R CMD check --no-manual --no-tests` were run at the closure source by CORE-015 (`core-015/removed.md`). The suite has no file worse than before CORE-015. The check has no errors; 3 warnings and 3 notes, all pre-existing (it is not "OK").
- **Not run:** the installed `test-workflow-rollup-qualification.R` readback and `test-workflow-saved-installed.R`, which need further guarded members beyond the one ledgered CORE-016 check.
- **Files changed:** `source-reuse.json` lists every `R/` file the series changed, with its commits and tickets. The C provider change (CORE-001) is in `src/`.

## Rollback

Per-ticket `git revert` in reverse order. Accepted references, history and evidence are never deleted. PERF-108 evidence is reproduced by checking out `9df060f0` or earlier.
