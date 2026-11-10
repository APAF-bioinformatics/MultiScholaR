# KIT-113 W5, chain B: install, check, bindings, journeys (2026-11-10)

- **Install:** `kit_113` = `core/browser-smoke.cmBMvS`, from a pinned worktree at `95600989`. Bindings: `b64bbc64`.
- **Run 5, `R CMD check --no-tests`** ([../check.log](../check.log)): 3 WARNINGs, 3 NOTEs, the same sections as KIT-112's. The "checking R code for possible problems" NOTE is byte-identical to KIT-112's.
- **Run 7 (M1), the journeys:** one attempt each, under the host rule ([../power-state.txt](../power-state.txt)).

| Member | Decision | Exact | Dispatches against 112l/m/t | Elapsed |
|---|---|---|---|---|
| 113l | accepted | yes | identical (added none, removed none, changed none) | 137.2 s |
| 113m | accepted | yes | identical | 236.2 s |
| 113t | accepted | yes | identical | 139.6 s |

No journey dispatches `enrichment`, because none reaches it (D1). The memory figures are in [../memory-outcome.md](../memory-outcome.md).
