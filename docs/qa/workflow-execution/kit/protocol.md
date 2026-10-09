# KIT measurement protocol

**Version 1 (draft, 2026-10-09).** Frozen by KIT-001 before any counted KIT member; its SHA-256 is then recorded in `kit/bindings.json` and a change to this file is a new version needing a user decision. Carries the CORE protocol v8.1 (`../core/protocol.md`) unless stated.

## 1. Rules carried from CORE

- Memory alone decides; elapsed and CPU time are reported per pair, never judged.
- The measurement works with the host as normally used: no power-profile or application changes are requested; profile and busy CPU are recorded per member, never a refusal.
- A member's measurement is valid with a monitor sampling gap up to 250 ms (the peak is the kernel high-water mark); zero swap; cleanup proven (worker and monitor stopped).
- Each member runs as a managed application in its own slice: 16 GiB hard cap, 15.5 GiB early stop, 12 GiB `MemAvailable` at start, 2 GiB continuous host floor, eight threads.
- One mechanism per ticket. A unit-test pass is never a performance claim. A valid unfavourable result, an OOM or a timeout is recorded and never retried; one correction rerun of a single member is allowed only for a demonstrated instrumentation or harness failure.
- Builds: one recorded install per candidate on a clean tree; source commit, patch hash, build identity and exact-sum identity recorded in the bindings before any member uses the build.

## 2. Control

The control is the latest accepted KIT build, run in the same session as the candidate. Until KIT-006 accepts, the control is `core_018` (source `906c50dc`, install `browser-smoke.gQ40C7`). A control receipt is never reused across a build or protocol change; controls are re-run.

## 3. Decision rules

- **Not worse:** candidate peak ≤ 1.02 × control peak (the per-pair rule of CORE-005 had no noise tolerance; this one does).
- **Required reduction:** a ticket that exists to lower memory claims it only when candidate ≤ 0.95 × control on the named metric.
- **Exactness:** every component, support entry, args summary and product of the step is hash-identical to the oracle (zero-row safe fingerprints). The oracle for a post-rollup step is produced by KIT-005's counted control journey on `core_018`, which snapshots the store after each step as the bound population and fingerprints its outputs; the second baseline member must reproduce every fingerprint. No in-process full-size oracle runs exist. The `core` populations keep their bound oracles (the P14 page oracle; `oracles.filter_primary` pages, summary, support and design hashes).
- **Single members** (KIT-014, KIT-015, KIT-018) are exactness-and-completion only: oracle hash, peak recorded, no ratio. The memory decision for migrated steps is taken once by the closure baseline pair (KIT-017), per pair as in CORE: the required reduction is judged against the control run in the same session, KIT-005's two baseline members give the noise estimate for both metrics, and no step's anonymous peak may be worse than 1.02 × the control's.
- **Unit confirmation gate:** a step is declared only when the tracked `kit-001/step-map.json` (hashed into the bindings) records its independence unit as confirmed by the user as science owner, with the decision reference and date. The gate applies to the partitioner proof steps as well (P14, the sample filter, normalisation/RUV, correlation and DA before KIT-004 and KIT-006). A kernel that takes a dataset-wide statistic declares it as a `global_statistics` pre-pass; a partitioned declaration carries `merge` rules for every non-assay component it writes.

## 4. Member kinds

| Kind | Scope | Metrics |
|---|---|---|
| operation member | one registered operation on one bound population in the worker slice | kernel peak (`peak_bytes`), whole-slice anonymous peak, phase peaks, batches, elapsed and CPU (reported) |
| journey member | the KIT-001 driver through the whole DIA journey on the 300-sample synthetic corpus, under a `gui`-mode managed launch, with the journey's state-manager mode and per-step settings frozen in the bindings | whole-slice anonymous peak, per-step anonymous peaks (`phase_peak_anon_bytes$<step>`), `resident_end_bytes` (minimum anonymous memory in the last 5 s before teardown, from the coordinator's anonymous trace and a teardown mark), elapsed (reported); ceiling set from the KIT-001 screen duration scaled to 300 samples and recorded before launch, never above 90 min |
| baseline pair | `core_018` (source `906c50dc`, install `browser-smoke.gQ40C7`) versus a candidate, both journey members in one session; used by KIT-005 (baseline, populations and oracles) and KIT-017 (closure) | as journey member, judged per pair |
| lifecycle member | completion and forced failure through a changed launch or measurement path, on the screen corpus | observed / not observed |

Selectors: `kit:0NN:(control|candidate):[12]` for operation members, `kit:journey:(control|candidate):[12]` for journey and baseline-pair members, `kit:lifecycle:(complete|fail)`. All KIT selectors live in the KIT harness files; the CORE-pinned files are not edited.

## 5. Populations and oracles

`kit/bindings.json` (KIT-001) carries the `core` populations (`screen_p14`, `full_p14` with the P14 page oracle, `full_primary` with the sample-filter oracle), the PERF-107 synthetic DIA-NN import corpora (8-sample screen, 300-sample full), the KIT harness file hashes and this protocol's hash. KIT-005's counted control journey adds one population per DIA step (store path, `dia-state-head.rds` sha256, input entity id, rows) and the full fingerprint set that is its oracle. The four CORE-pinned harness files are never edited by KIT tickets.

## 6. Run ledger

The series plan's run ledger is the only authority for counted members. No other full-size execution is permitted unless the user amends it.
