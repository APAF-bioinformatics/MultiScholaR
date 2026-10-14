# KIT-023: the publication phases bind no globals

Status 2026-10-14: **complete** (no counted member; no separate install: the KIT-024 build, KIT-017's candidate, carries this change). Authority: `tickets/kit/00-KIT-PLAN.md` (decision row 2026-10-14, best-practice pass; run ledger KIT-023: none, one install recorded); protocol v1.4. Control `kit_018`; KIT-024's build is KIT-017's candidate.

## What the evidence changed

The ticket was authored to move the peptide and protein figure rendering into the worker so that the application process would bind no globals for the legacy figure writers. Reading the chain before editing showed that the figure writer (`updateProteinFiltering`) takes `project_dirs` as an argument and reads the global only as a fallback, that the controller already passed the directories, and that nothing in the render chain reads `config_list`. The two `rlang::local_bindings(... .env = .GlobalEnv)` calls in the controller's publication phases were dead weight, and moving the rendering (count plots, no heavy namespace) had no measured cost to justify it. Revision 1.1 of the ticket records the withdrawal of the worker-side staging; the repair is the removal.

## What was done

- `publishWorkflowPeptideResult` and `publishWorkflowProteinResult` (`R/utils_workflow_execution.R`) no longer bind `config_list` and `project_dirs` into the global environment; the figure writer receives the directories as before. The file has no `.GlobalEnv` reference left.
- `publishDIAProteinProgress` (`R/utils_workflow_dia_qc_rollup.R`, the in-process GUI path) no longer reads the global `project_dirs`; it builds the directories from the owner's experiment paths. The file has no `.GlobalEnv` reference left.
- Tests: the KIT-010 protein publication test and the KIT-015 publication loop observe the global environment at the figure boundary (`workflowExportBoundary("before_qvalue_figures")`) and assert the publication phase added nothing.
- README: the "Kernel contracts the kit adapts to" entry now records that no global binding exists at publication; the kernels' own fallbacks stay in the kernels and are never reached from the kit.

## Tests (working tree, proof library on `R_LIBS`)

| Suite | Result |
|---|---|
| test-workflow-partitions.R | 301 pass (the KIT-010 publication test observes the global environment at the figure boundary) |
| test-workflow-prot-declarations.R | 117 pass (the KIT-015 publication loop likewise) |
| test-workflow-execution.R / test-perf107-audit.R | 64 (10 skips) / 193 pass |
| `grep -c GlobalEnv` on the two files | 0 and 0 |

## Claim limits

- The multi-session hazard (two sessions publishing at once sharing one global) is removed by construction in the controller; it is asserted by the boundary test, not observed under concurrent sessions.
- The kernels' global fallbacks remain reachable from scripts that call the writers without directories; the kit never does.
- No counted member.
