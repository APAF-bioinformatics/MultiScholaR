# CORE-013 outcome: closed without change (premise false)

Ticket: `tickets/openai/CORE-013.md` revision 1.1. Decision 2026-10-08 (user): close without code change or counted runs after the premise check below. **No code changed; no guarded run spent.**

## What was measured

`overhead-probe.R` (in this directory) was run in-process with the `core_010` installed library, on the project filesystem (`/home`, LUKS); `/tmp` is tmpfs and would understate `sync`.

| Suspected cost | Measured | Effect on an action |
|---|---|---|
| `systemctl show` in the controller poll (two calls per 50 ms tick) | 7.6 ms per call, about 15 ms per tick | The controller is busy about 23% of one core (of 16) while waiting. The worker runs in its own unit and is not slowed, and other Shiny work waits at most about 15 ms. |
| Durable journal or control write with `sync` (about 30 per action) | 5.1–7.2 ms each over three runs | 0.15–0.21 s per action, against actions of 10–120 s |
| Monitor re-reading `application-active-worker.rds` every 50 ms | well under 1 ms | negligible |

The 100 ms sleep in `pollWorkflowAction` runs once, on the terminal path. The 20 ms readiness loops run during application launch and admission only.

## Why no change

- **Slower polling would make actions finish later.** Polling at 250–500 ms (item 1) delays completion detection by up to that interval, so a candidate would likely finish later, not sooner.
- **Fewer syncs would save little.** Merging journal states and batching control writes (item 3) saves at most about 0.2 s per action, below run-to-run variation, and touches the durable recovery path for no observable gain.
- **Items 2 and 4 are housekeeping** with no measured cost.
- **The stage could not show an effect.** A four-member stage on the sample-filter action, possibly plus a fresh control pair, could not resolve these amounts under protocol v7, where time is reported and not judged.

## Handoff

Nothing to CORE-015 beyond this record.
