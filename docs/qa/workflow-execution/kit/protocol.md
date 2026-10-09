# KIT measurement protocol

**Version 1 (draft, 2026-10-09).** Frozen by KIT-005 before any counted KIT member; its SHA-256 is then recorded in `kit/bindings.json` and a change to this file is a new version needing a user decision. Carries the CORE protocol v8.1 (`../core/protocol.md`) unless stated.

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
- **Exactness:** the candidate's output table component pages are hash-identical to the oracle. The oracle is the in-process pre-kit build (`core_018`) run once on the same bound population and fingerprinted into the bindings, or the P14 page oracle for the replacement kernel.
- **Single members** (KIT-014, KIT-015, KIT-016) are exactness-and-completion only: oracle hash, peak recorded, no ratio. The memory decision for migrated steps is taken once by the closure journey pair (KIT-017).
- **Unit confirmation gate:** a step is declared only when `kit-001/step-map.json` marks its independence unit `confirmed` by a decision-record row in the series plan.

## 4. Member kinds

| Kind | Scope | Metrics |
|---|---|---|
| operation member | one registered operation on one bound population in the worker slice | kernel peak (`peak_bytes`), whole-slice anonymous peak, phase peaks, batches, elapsed and CPU (reported) |
| journey member | the KIT-001 driver through the whole DIA journey on the 300-sample synthetic corpus | whole-slice anonymous peak, per-step anonymous peaks (`phase_peak_anon_bytes$<step>`), `resident_end_bytes` (minimum anonymous memory in the last 5 s before teardown), elapsed (reported) |
| lifecycle member | completion and forced failure through a changed launch or measurement path, on the screen fixture | observed / not observed |

Selectors: `kit:0NN:(control|candidate):[12]` for operation members, `kit:journey:(control|candidate):[12]` for journey members, `kit:lifecycle:(complete|fail)`.

## 5. Populations and oracles

`kit/bindings.json` (KIT-001) carries the `core` populations (`screen_p14`, `full_p14` with the P14 page oracle, `full_primary` with the sample-filter oracle), the PERF-107 synthetic 300-sample DIA-NN import corpus, and one population per DIA step snapshotted after that step by the KIT-001 driver (store path, `dia-state-head.rds` sha256, input entity id, rows). KIT-005 records an in-process oracle for every per-step population on `core_018`.

## 6. Run ledger

The series plan's run ledger is the only authority for counted members. No other full-size execution is permitted unless the user amends it.
