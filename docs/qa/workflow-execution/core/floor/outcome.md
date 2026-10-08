# Shared floor: what the slice holds before any workflow work starts

Read-only diagnostic (2026-10-09), after the CORE closeout. No guarded run, no package change, no counted member. Question: the closeout reports about 2.4 GiB in use before the P14 rollup starts; is it rollup memory or a floor every workflow step pays?

## Answer

It is a floor every step pays, and it is almost entirely package namespaces. Each managed slice runs three R processes (controller, monitor, worker; `workflowLaunchService`), each started as `Rscript -e "MultiScholaR:::...Main(...)"`, so each loads the full package namespace with all 87 imported packages before doing anything. One such process costs about 0.83 GiB of anonymous memory. Three of them are 2.43 GiB, which is the measured pre-rollup floor.

| What one fresh R process holds (anonymous, median of 3, closure build `Hz6YwJ`) | MiB |
|---|---:|
| Bare R | 41 |
| Full `MultiScholaR` namespace (87 imports) | 830 |
| All imports except 11 app-only packages (Glimma, cowplot, ComplexHeatmap, ggraph, Biobase, seqinr, git2r, devtools, circlize, GGally, ggpubr) | 395 |
| A worker-style set (arrow, duckdb, DBI, digest, cli, jsonlite, dplyr, fastduckdb, rlang, tibble, data.table, later) | 124 |

The package's own code is 11 MiB of that; the rest is its import list. Attributed in load order, Glimma (284), cowplot (121), ComplexHeatmap (118), ggraph (75) and Biobase (47) are the largest first-time loads, each pulling in shared Bioconductor and graphics dependencies.

## Evidence from the existing slice traces (50 ms, whole slice, anonymous)

- **Small fixture** (`core-lifecycle-complete`, 8 samples): 0.8 GiB plateau at 6–9 s (first process loaded), a steady ramp to 2.37 GiB by 21 s (the other two processes loading), the rollup itself at 2.4–2.6 GiB, then 1.6–1.68 GiB after the worker exits (controller plus monitor).
- **Full input** (`core-010-control/control-1`, 24,000,000 precursors): the same 0.8 GiB plateau, the same ramp to 2.37–2.43 GiB by 27–32 s before the rollup starts, and the same 1.6–1.68 GiB after the worker exits.
- The floor is therefore independent of input size (within 0.06 GiB between an 8-sample fixture and the full corpus), which is what a fixed per-process cost looks like.

## What it means

- **Workflow-wide:** every action run through the execution layer pays about 2.4 GiB before its own work, on every machine, regardless of input. This is the first series number that applies to the whole workflow rather than one step.
- **Reclaimable in principle:** about 1.2 to 1.4 GiB of it.
  - Moving the 11 app-only packages from Imports to Suggests with `requireNamespace()` at their use sites (about 70 `::` call sites plus 6 `importFrom` lines) would take each process from 0.83 to about 0.4 GiB, so the three-process floor falls from 2.43 to about 1.16 GiB. This also lowers the GUI's own baseline by the same 0.4 GiB.
  - A slim loader for the worker and monitor (controller full, the other two on a worker-style set) would give about 1.05 GiB, but R cannot partially load a namespace, so this needs a separate small package or a sourced bootstrap. More work, less gain than the Suggests move.
- **Not reclaimable here:** the 3.2 GiB checks-phase floor is rollup work on the whole input and is unrelated.
- **Not claimed:** no run has been made with a slimmer namespace, so the reclaim figures are projections from per-process measurements, not measured workflow peaks. The GUI (Shiny) process baseline was not measured; the controller measured here is the harness's managed application, which carries the same namespace.

## Method notes

- `floor-probe.R` launches fresh `Rscript --vanilla` children and reads `RssAnon` from `/proc/self/status` after one `gc()`; `floor-probe.txt` is its recorded output. Each namespace is loaded by its own top-level call, matching how the managed processes are launched. The same loads written inside a `for` loop measured about 130 MiB higher per process; that form was discarded.
- Load-order attribution per package (first measurement, not in the probe) is approximate because shared dependencies are charged to whichever package loads them first. The set-level figures in the table do not have that problem.

## Follow-up

Done as CORE-018 (2026-10-09): stage `core-018-floor` measured the validation-phase floor at 2.38 GiB on the closure build and 1.16 GiB with the eleven packages in Suggests (`core-018/outcome.md`). The original suggestion follows.

## Suggested bounded follow-up (as written before the ticket)

One ticket: move the 11 app-only packages to Suggests behind `requireNamespace()` guards, keep public behaviour identical, and measure the floor on the existing small-fixture lifecycle member (one guarded screen run, about 35 s) against the 2.37 GiB plateau recorded here. Decision rule: accept if the pre-rollup plateau falls by at least 1 GiB with the lifecycle checks unchanged.
