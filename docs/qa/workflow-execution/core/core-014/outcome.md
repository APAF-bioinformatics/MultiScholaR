# CORE-014 outcome: closed without change; handed to PERF-107

Ticket: `tickets/openai/CORE-014.md` revision 1.1. Decision 2026-10-08 (user): close without code change or counted runs. **No code changed; no guarded run spent.** The waste is real, but the safe point to remove it lies in the app's session-swap code, which PERF-107 owns.

## What was measured

`store-closure-probe.R` (in this directory, read-only) opens a store and computes the closure that Save archives (`expandWorkflowTableReferences(collectDIANativeReferences(list(snapshot, head)))`). It then classifies every entity outside that closure by descriptor kind and identity.

| Store | Size | Outside the Save closure | What it is |
|---|---|---|---|
| App session `browser-smoke.lSCYZW/session-store` | 1.13 GB | 33 of 50 entities, 0.66 GB (58%) | earlier app working objects (`owner:data_tbl`, `owner:correlation_filtered_obj`, designs, normalised results), the last Save's `native-session-metadata` and `export-metadata`, one superseded rollup table, 18 superseded head snapshots |
| Harness store `wide-native-1` | 12.84 GB | 77 of 147 entities, 4.98 GB (39%) | result tables and products from re-run steps, 34 superseded head snapshots |
| Leftover `.dia-native-import-*` folders beside sessions | 731 MB each | whole folders | previous sessions' stores, left behind after Open |

Superseded head snapshots, the one category nothing else can reference, total about 4 MB.

## Why no change in CORE

- **The app holds some entities outside the closure.** The running app keeps references that the store cannot see: the last Save's metadata (`context$owned_export`, used to resume or re-publish an export) and the working-object roots created by app adoption (`prepareDIAGUIReferences`). A general clean-up keyed on the Save closure could remove data a live session still uses.
- **Restricting the clean-up to provably dead head snapshots** would save almost nothing.
- **Open already starts clean.** It extracts only the archived closure into a new private folder (`openDIANativeSession`). The lasting waste is the replaced folder that each Open leaves behind, plus a live store's growth until the next Save and Open.
- **Items 2–4** (head-publish coalescing, dependency deltas, hard-link copies) had no measured cost behind them and were not pursued.

## Handoff to PERF-107 (via CORE-016)

- **After a successful Open swap,** once the app has released the previous manager and its rollback handler has exited, remove the replaced private session folder (`.dia-native-import-*` or `.dia-gui-adoption-*`).
- **A live-store clean-up** becomes safe only once the app declares the references it holds, so the clean-up's keep-set includes them. Until then, Save and Open remain the clean path.
