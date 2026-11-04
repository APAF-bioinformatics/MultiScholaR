# KIT-081 outcome: the stage-3b measurement tools

Status 2026-11-04: **completed**. Stage 3b. The instruments KIT-080 will measure with are in place and proven end to end on the screen corpora. No package code changed. The members here were recorded and are not counted; they support no memory claim.

- **Commits** on `performance-v2`: protocol v1.6 992ea0ff; instruments dbd5b88b; records 1b0b136f; route correction f1cef413; bindings 13d4bb12 (closure) and 4f2033cf (rerun).
- **Install:** `kit_081` = `core/browser-smoke.6nbYY3` (dbd5b88b).
- **Ticket:** `tickets/kit/KIT-081.md`, revision 1.1.

## What changed

1. **Protocol v1.6** is section 8 of [protocol.md](../protocol.md), frozen on the user's decision of 2026-11-04. Worker steps (intensity, duplicates, normalisation, correlation, DA) are judged by `app_process`. Their whole-slice per-step figures are reported only, and the whole-run peak decides the total-memory claim.
2. **The controller trace.** The sampler reads the controller's cgroup as `controller_anon` at the same cadence. The coordinator reports `controller_peak_anon_bytes` and the per-step `phase_controller_peak_anon_bytes`.
3. **Operation members on a list assay.** The input can be a named list of references, identified by the digest of its sorted element ids and fingerprinted per element.
4. **Suffixed operation keys** `kit-0NN[a-z]`, with a bound ceiling of at most 1,800 s.
5. `kitWorkerPairDecision` applies v1.6's two rules when a ticket binds `worker_steps`.

## Results

| Gate | Evidence | Result |
|---|---|---|
| Controller trace through completion | `kit-lifecycle-lipid`, `kit-lifecycle-metab` | Both `lifecycle_complete_observed`. Peak controller anonymous memory 0.91 and 0.89 GiB; slice 1.86 and 1.85 GiB. Largest sampler gaps 0.09 and 0.07 s. |
| Controller trace through failure | `kit-lifecycle-lipid-fail`, `kit-lifecycle-metab-fail` | Both `lifecycle_failure_observed`. The trace runs until the guard stops the member: 1,384 and 1,241 samples, controller peak 0.67 and 0.66 GiB. Largest gaps 0.09 and 0.06 s. |
| List-assay operation member, suffixed key | `kit-081a` | `accepted`, with its decision computed. Four members verified. On this small corpus `whole` and `partitioned` are equal in memory: slice 1.05 GiB, controller 0.35 GiB, anonymous ratio 1.00, so no reduction is shown (`reduction FALSE`, `not_worse TRUE`). |
| Earlier decisions | `test-workflow-kit-harness.R` | The recorded pair decisions of KIT-017, KIT-036 and KIT-057 replay unchanged, 525 tests pass. KIT-012 has no pair stage JSON to replay, so the ticket's list could not be met for that one. |
| Omics records | [exactness-gates.txt](exactness-gates.txt) | Identical to `kit_079`'s on both omics, as expected with no package change. |
| Regressions | `test-workflow-execution*.R`, four files | 173 pass, 0 fail, 0 error, 16 skip. |
| Package | `R CMD check` | Baseline: 3 WARNINGs, 3 NOTEs, on the tarball of 4f2033cf. |
| Pinned files | `git diff 383567ac HEAD` | No change under `R/`, `NAMESPACE`, `DESCRIPTION` or `man/`. Under `tools/profiling` only the kit's own three files changed; the CORE-pinned runner did not. |

## Failed attempts, kept

- **`kit-081a.route-unbound`:** the first launch refused every member. `operation_case` dropped `route`, so the vehicle compared a list assay with the state's id. Fixed in f1cef413, with a test. This was an instrumentation failure, and the rerun used the correction budget.
- **`kit-lifecycle-metab.sampler-gaps`:** gaps over the frozen limit, clustered in start-up before the import, invalidated this member, so it was rerun once. The rerun was valid with a largest gap of 0.07 s. The cause is not established; the laptop was on battery in low-power mode for every member here, including the valid ones. The 0.5 s and 0.25 s limits were not changed.
- **The harness file in the parallel batch:** `test-workflow-kit-harness.R` left no result line when it ran beside the four execution files. Run alone from the same clean worktree it passes in 20 s. I did not find out why the batch lost its line.

## Limits

- Time is not comparable: every member ran on battery in low-power mode ([power-state.txt](power-state.txt)).
- `kit-081a` proves the instrument, not the benefit of partitioning. KIT-084 and KIT-080 measure that at larger sizes.

## Hand-off

KIT-082 uses the controller trace as its baseline: for metabolomics by step, import 0.51, design 0.53, intensity 0.63, duplicates 0.64, QC finalise 0.65, normalisation 0.57, correlation 0.57, session export 0.61, DA 0.74, report 0.88 GiB.
