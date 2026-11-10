# KIT-112 W5, chain B: install, check, bind, the journeys (M2 and run 7), 2026-11-10

- **Install:** `kit_112` = `browser-smoke.OftnRb`, from a pinned worktree at `5ecc9fd1` (R code as `a53f9060`).
- **`R CMD check --no-tests`:** 3 WARNINGs, 3 NOTEs. The headers equal `kit-111/check.log`'s, and the "checking R code for possible problems" section is byte-identical.
- **Bindings:** `00e758d1`.

| Member | Decision | Exact | Seconds | Dispatches against KIT-111's member | `da` controller peak (GiB) | `da` seconds |
|---|---|---|---|---|---|---|
| 112l | accepted | yes | 165.6 | adds `DA-contrast-output` 1/1; nothing else changed | 0.657 → 0.657 | 14.3 → 19.8 |
| 112m | accepted | yes | 300.5 | adds `DA-contrast-output` 1/1; nothing else changed | 1.018 → 0.983 | 16.0 → 22.6 |
| 112t | accepted | yes | 167.0 | adds `DA-contrast-output` 1/1; nothing else changed | 0.660 → 0.658 | 13.3 → 18.4 |
| 112f | `journey_failed` at `da` (harness decision `invalid`, failing rule `verifier_accepted`) | — | 128.9 | `DA-contrast-output` 1/0; no `da_publication`. Error: "production action failed: … managed action ended: worker_failed". Store head `correlation_filtered`, with no `da_results` | — | — |
| 112a | accepted | yes | 204.1 | identical to 111a | — | — |
| 112b | accepted | yes | 197.0 | identical to 111b | — | — |
| 112c | accepted | yes | 132.6 | identical to 111c | — | — |

Every member is valid: largest sample interval ≤ 0.15 s, swap 0, slice stopped. On the screen members, every other step rose by at most +0.017 GiB, and the slice peak changed by −0.034 to +0.002 GiB.

**M2 (D8's lifecycle test) passes:** the instruments see `DA-contrast-output` through completion (112l/m/t) and through failure (112f). The screen "no worse" rule (the `da` step only, lead decision) holds on all three. The `da` step is 4.5 to 6.6 s slower, the worker's start-up on a small fixture; it is reported, not judged.
