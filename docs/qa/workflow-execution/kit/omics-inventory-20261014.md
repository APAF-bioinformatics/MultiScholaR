# Inventory of the kit-relevant idioms in the other omics (2026-10-14)

Recorded before stage 3 is authored (user request after the KIT-022 report: the same grep pass the kit and DIA paths had). Read-only; counts are files, not call sites. Transcriptomics has no module files yet (TRX-005 extracts the notebook kernels first).

| Pattern | Lipidomics (75 files) | Metabolomics (68 files) | What it means for stage 3 |
|---|---|---|---|
| `.GlobalEnv` reads or writes in modules | 11 | 5 | The KIT-013 pattern applies before any declaration: module globals onto the workflow bus; the kernels' own global reads become adapters as in DIA |
| `WorkflowState$new()` as the state owner | 2 | 3 | KIT-030 replaces them with the adapter on the KIT-003 class declarations |
| `saveState(` with a full S4 | 7 | 8 | Every save site becomes a reference commit (`commitDIAReferenceState` through the declaration), as KIT-011 did for DIA |
| `current_s4_object <-` holders | 1 | 1 | Few holders; the KIT-012 rule (references only) is cheap here |
| `renderDT` / `renderDataTable` full-table renders | 10 | 11 | Paged renders from saved tables (KIT-018/KIT-024 pattern) |
| `getDIAActiveManager()` | 0 | 0 | The owner resolution is DIA-named; stage 3 needs the omic-agnostic name (a kit gap: KIT-030 renames or wraps) |
| `checkParamsObjectFunctionSimplify` in kernels | 62 call sites across both | — | Kernel contract (never rewritten); the kit calls kernels explicitly |

Conclusion for the plan: KIT-030 should carry an explicit first item "globals onto the bus and holders to references" (the KIT-013 and KIT-012 patterns) before the state manager switch, and KIT-031's "declarations only" acceptance test needs the owner resolution and the display-table helpers to be omic-agnostic names, which is a kit gap to close in KIT-030 rather than in lipidomics.
