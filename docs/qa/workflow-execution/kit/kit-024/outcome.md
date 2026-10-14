# KIT-024: kit consistency before closure

Status 2026-10-14: **complete** (no counted member; one install recorded: kit_024 = `core/browser-smoke.70tjXp`, `28075aa1`, bound as KIT-017's candidate; the control stays `kit_018`). Authority: `tickets/kit/00-KIT-PLAN.md` (decision row 2026-10-14, best-practice pass; run ledger KIT-024: none, one install recorded); protocol v1.4. Control `kit_018`; this build is KIT-017's candidate.

## The three items

1. **Every declaration registered** (`R/utils_workflow_operations.R`). The five whole-matrix declarations that run in process only (`normalization-ruv`, `pearson-pairs`, `pearson-selection`, `DA-contrast-output`, `enrichment`) are registry entries with the compute role `in_process` (no handler; the declaration attached; the validator requires exactly that shape), so `workflowOperationNames()` lists every declaration. `assertWorkflowOperationDispatchable` refuses them at `runWorkflowAction` and `beginWorkflowGUIAction` with `workflow_execution_unsupported`; the worker entry point already refused an unknown role. Registry 19 operations.
2. **Server-side paged renders, references from the accepted state.** `diaDisplayTableFilter` (`R/utils_workflow_dia_products.R`) is the table widget's page filter over a saved table: the widget asks for a page and the store serves it (`readDIADisplayTablePage`, `readDIAResultPage`); an in-memory table keeps the widget's own filter. The RUV optimisation grid (`renderProtNormRuvOptimizationTable`), the g:Profiler and clusterProfiler results tables and the DA results table render through it with a visible row count (ordering and searching off: the store serves source order). After an Open, `resolveProtNormRuvOptimizationTable` saves the grid of the accepted normalisation product as a display table on first use, and `restoreProtEnrichResultsFromState` rebuilds the per-contrast tables (and their display tables) from the `enrichment_completed` state's `EnrichmentResults` on first use, with the method inferred from the results (`diaEnrichmentMethodOf`).
3. **Baseline test quarantine** ([quarantined-tests.md](../quarantined-tests.md)). The full suite was run once file by file (343 files); every failing file was rerun at the pre-ticket commit `e972ad02` to separate pre-existing failures from this ticket's. Pre-existing: 36 tests in 25 files, skipped with a reason naming their entry (Q-01 to Q-20: what fails, since when, what would fix it). Repaired instead of quarantined, where the cause was a stale test rather than the code: the DA render direct suite's widget mock (this ticket's server-side arguments), the peptide rollup shared suite's plot mock (the KIT-014 acknowledgement identity), and the saved-QC module-edge test, which read the configuration from the global environment that KIT-013 moved onto the workflow bus (its own `info` argument bug had masked eight expectations). MCI-007.6 passes at the current tree and is not quarantined.

## Tests (working tree, proof library on `R_LIBS`)

| Suite | Result |
|---|---|
| test-workflow-operation-registry.R | 125 pass (new: the five in-process entries, the dispatch refusal, the validator shape) |
| test-workflow-dia-products.R | 123 pass (new: the page filter on store pages and the in-memory fallback; the method inference; the after-Open resolution of the RUV grid and the enrichment tables on a reopened owner) |
| test-prot-05b-norm-module-contracts.R / test-prot-11-enrichment-module-contracts.R | 562 / 1253 pass (widget mocks widened for the server-side arguments) |
| test-module-ci-prot-da.R / test-prot-da-render-direct-shared.R / test-prot-07b-da-handlers-compat.R | 245 / 37 / 55 pass |
| test-workflow-partitions.R / test-workflow-execution.R / test-perf107-audit.R / test-workflow-kit-harness.R | 306 / 64 / 193 / 326 pass |
| Full suite, first pass (before the quarantine was complete) | 343 files: 18 failures and 22 errors in 19 files beyond the known 13 tests; each file rerun at `e972ad02`: one KIT-024 failure (the DA render mock), three stale tests repaired, the rest pre-existing and quarantined |
| Full suite, second pass | 343 files, 29,432 pass, 0 fail, 0 error, 152 skips (every quarantine skip named) |

## Claim limits

- The page filter implements paging only; ordering and searching are disabled on the paged renders (the exports in the widget's buttons cover the visible page; the modules' download handlers keep the full tables).
- The after-Open resolvers run on first use in the session that opened the project; they save display tables (content-addressed) into that project's store.
- The quarantine hides nothing: every skipped test is listed with its cause, and a passing quarantined test is to be un-skipped.
- Five quarantine entries and one repaired test are expectations overtaken by KIT-008, KIT-011, KIT-013 and KIT-015, whose regression checks ran only the files they touched; the quarantine list's closing section says so for the KIT-017 closure, and the whole suite runs once per ticket from here on.
- No counted member.
