# KIT measurement protocol

**Version 1.3 (frozen by KIT-001, 2026-10-09; amended 2026-10-09 before any counted member, user decision: the control is the rebuilt `kit_001` build and KIT-005 adds a production-default baseline member; amended 2026-10-10 by user decision after the first 300-sample journey attempt: the journey member ceiling may be set up to 4 h (14,400 s) instead of 90 min, because the data steps scale about 20 x from 8 to 300 samples; the launch record's deadline limit is raised to match and the control is rebuilt from the `kit_001` source plus that one constant; amended 2026-10-11 by user decision after KIT-005's baseline stage: a journey member's measurement tolerates isolated sampling gaps up to 500 ms, because its whole-slice peak is the kernel's own high-water mark and a gap only blinds the sampled anonymous trace for that window; the number of gaps above 250 ms is recorded with the member; operation and lifecycle-operation members keep the 250 ms rule).** Its SHA-256 is recorded in `kit/bindings.json` and a change to this file is a new version needing a user decision. Carries the CORE protocol v8.1 (`../core/protocol.md`) unless stated.

## 1. Rules carried from CORE

- Memory alone decides; elapsed and CPU time are reported per pair, never judged.
- The measurement works with the host as normally used: no power-profile or application changes are requested; profile and busy CPU are recorded per member, never a refusal.
- A member's measurement is valid with a monitor sampling gap up to 250 ms (the peak is the kernel high-water mark); zero swap; cleanup proven (worker and monitor stopped). A journey member (v1.3) tolerates isolated gaps up to 500 ms; its gaps above 250 ms are counted in the receipt and the sampled anonymous peak is reported with that count.
- Each member runs as a managed application in its own slice: 16 GiB hard cap, 15.5 GiB early stop, 12 GiB `MemAvailable` at start, 2 GiB continuous host floor, eight threads.
- One mechanism per ticket. A unit-test pass is never a performance claim. A valid unfavourable result, an OOM or a timeout is recorded and never retried; one correction rerun of a single member is allowed only for a demonstrated instrumentation or harness failure.
- Builds: one recorded install per candidate on a clean tree; source commit, patch hash, build identity and exact-sum identity recorded in the bindings before any member uses the build.

## 2. Control

The control is the latest accepted KIT build, run in the same session as the candidate. Until KIT-006 accepts, the control is `kit_001`: the `core_018` source (`906c50dc`) plus the KIT journey launch record (`kit-journey-v1`, native DIA sessions in the controller), the RUV-III technical-replicate structure and the removal of the quantile ceiling in the proteomics intensity filters (user decisions 2026-10-09). Its install is recorded in the bindings before any member uses it. A control receipt is never reused across a build or protocol change; controls are re-run.

## 3. Decision rules

- **Not worse:** candidate peak ≤ 1.02 × control peak (the per-pair rule of CORE-005 had no noise tolerance; this one does).
- **Required reduction:** a ticket that exists to lower memory claims it only when candidate ≤ 0.95 × control on the named metric.
- **Exactness:** every component, support entry, args summary and product of the step is hash-identical to the oracle (zero-row safe fingerprints). The oracle for a post-rollup step is produced by KIT-005's counted control journey on `core_018`, which snapshots the store after each step as the bound population and fingerprints its outputs; the second baseline member must reproduce every fingerprint. No in-process full-size oracle runs exist. The `core` populations keep their bound oracles (the P14 page oracle; `oracles.filter_primary` pages, summary, support and design hashes).
- **Single members** (KIT-014, KIT-015, KIT-018) are exactness-and-completion only: oracle hash, peak recorded, no ratio. The memory decision for migrated steps is taken once by the closure baseline pair (KIT-017), per pair as in CORE: the required reduction is judged against the control run in the same session, KIT-005's two baseline members give the noise estimate for both metrics, and no step's anonymous peak may be worse than 1.02 × the control's.
- **Unit evidence gate:** a step is declared only when the tracked `kit-001/step-map.json` (hashed into the bindings) records its independence unit as determined from the kernel's code, with the evidence (function and the lines that do or do not use anything outside the unit), and no open `flag`. A flag is raised only where the code is batch-safe but the intended science might not be, or where the code is ambiguous; it carries a question for the user and waits for the answer. The gate applies to the partitioner proof steps as well (P14, the sample filter, normalisation/RUV, correlation and DA before KIT-004 and KIT-006). A kernel that takes a dataset-wide statistic declares it as a `global_statistics` pre-pass; a partitioned declaration carries `merge` rules for every non-assay component it writes.

## 4. Member kinds

| Kind | Scope | Metrics |
|---|---|---|
| operation member | one registered operation on one bound population in the worker slice | kernel peak (`peak_bytes`), whole-slice anonymous peak, phase peaks, batches, elapsed and CPU (reported) |
| journey member | the KIT-001 driver through the whole DIA journey on the 300-sample synthetic corpus, under a `gui`-mode managed launch with the KIT journey record (native DIA sessions from import), with the per-step settings frozen in the bindings | whole-slice anonymous peak, per-step anonymous peaks (`phase_peak_anon_bytes$<step>`), `resident_end_bytes` (minimum anonymous memory in the last 5 s before teardown, from the coordinator's anonymous trace and a teardown mark), elapsed (reported); ceiling set from the measured step scaling and recorded before launch, never above 4 h (14,400 s; v1.2) |
| baseline pair | the control build (`kit_001` until KIT-006 accepts) versus a candidate, both journey members in one session; used by KIT-005 (baseline, populations and oracles) and KIT-017 (closure) | as journey member, judged per pair |
| production-default baseline | one journey member of the control build launched without a journey record (in-memory `WorkflowState`, as users run it); the memory baseline that KIT-011, KIT-012 and KIT-017 reductions are judged against; no per-step references, so no populations | whole-slice anonymous peak, per-step anonymous peaks, `resident_end_bytes`, elapsed (reported) |
| lifecycle member | completion and forced failure through a changed launch or measurement path, on the screen corpus | observed / not observed |

Selectors: `kit:0NN:(control|candidate):[12]` for operation members, `kit:journey:(control|candidate):[12]` for journey and baseline-pair members, `kit:lifecycle:(complete|fail)`. All KIT selectors live in the KIT harness files; the CORE-pinned files are not edited.

## 5. Populations and oracles

`kit/bindings.json` (KIT-001) carries the `core` populations (`screen_p14`, `full_p14` with the P14 page oracle, `full_primary` with the sample-filter oracle), the PERF-107 synthetic DIA-NN import corpora (8-sample screen, 300-sample full), the KIT harness file hashes and this protocol's hash. KIT-005's counted control journey adds one population per DIA step (store path, `dia-state-head.rds` sha256, input entity id, rows) and the full fingerprint set that is its oracle. The four CORE-pinned harness files are never edited by KIT tickets.

## 6. Run ledger

The series plan's run ledger is the only authority for counted members. No other full-size execution is permitted unless the user amends it.
