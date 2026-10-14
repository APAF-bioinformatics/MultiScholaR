# KIT-018: design checkpoint verification, reports from products, legacy-session adoption

Status 2026-10-14: **implemented** (`6a7304e7`); candidate install and the counted single journey member pending the user's go. Authority: `tickets/kit/00-KIT-PLAN.md` (run ledger KIT-018: 1 single journey member, import to design); protocol v1.3 (single members: exactness and completion only, peak recorded, no ratio). Control: `kit_015` (`core/browser-smoke.hXymGJ`). Predecessor cases: PERF-112, PERF-125, PERF-127, PERF-104 (legacy sessions); audit B32.

## What was verified or built

1. **Design checkpoint default: verified, not changed.** `buildProtDesignStateCheckpoint` commits the design state through `newDIAPeptideCheckpoint` whenever the owner is a reference-states adapter, which the DIA-NN/IQ import installs on Linux (`initializeDIANativeImport`). New test (`test-workflow-native-session.R`): on a real adapter the call returns `raw_data_s4`, the state is a saved reference, the design gained `tech_rep_group`, and the restored object is the `PeptideQuantitativeData` with the imported rows. The in-memory DIA branch is **kept**: it is the path of the Limpa route, of non-Linux hosts and of non-DIA-NN formats (the import leaves the in-memory owner for those; `test-workflow-native-session.R` "native import is the DIA-NN/IQ default and does not migrate Limpa or other formats"). The ticket's "removes any remaining in-memory branch" is therefore not applied; the divergence is recorded here (claim limit as the ticket's item 3 already states for Limpa and non-Linux).
2. **Reports and exports read accepted products; two tables paged.** `resolveDIARuvOptimizationResult(manager)` (`R/utils_workflow_dia_products.R`) reads the RUV optimisation result from the saved normalisation product the current normalised state references (`products$normalization_science`); the study parameter writer (`createWorkflowArgsFromConfig`) and the summary tab's session metadata (`summariseProtSummaryWorkflowMetadata`) take it first, before the module copy or the RDS file. Display tables: `saveDIADisplayTable` saves a bounded data frame as a content-addressed result table (a repeated save is a hit) and `readDIADisplayTablePage` pages it from the store, with an optional equality filter on one column (scanned in 10,000-row reads, rows read by ordinal). The normalisation step saves the RUV optimisation grid (`norm_data$ruv_optimization_table`) and `renderProtNormRuvOptimizationTable` pages from it (first page, 1,000 rows); the enrichment finalisation saves each contrast's g:Profiler and clusterProfiler result frames (`gprofiler_results_table`, `clusterprofiler_results_table`, carried through the contrast switch) and the two results-table renders page from them with the direction filter applied in the store read. The in-memory frames stay for the plots, summaries and downloads (holders are not the resident cost: KIT-012); the renders are what became bounded.
3. **Legacy-session adoption: verified, not changed.** The DA tab's Open calls `ingestDIANativeLegacySession` for a non-owned RDS session (`mod_prot_da_handlers_load_session.R`); it adopts on Linux, DIA, non-Limpa (gate `isDIANativeDefaultEnabled()`, TRUE since KIT-011). The existing adoption test now also asserts the adopted history is a reference history (`usesDIAReferenceStates`; every state a saved reference). Limpa and non-Linux legacy sessions keep the in-memory load (claim limit, audit B32).

## Claim limits (pre-execution audit PE-03, PE-04)

- The rendered report (`RenderReport`, the Rmd) reads the files the workflow published under the project, which are the publication contract; the summary tab's own reads (the parameter writer, the session metadata) take the accepted normalisation product.
- The paged renders show the first page (at most 1,000 rows) of the saved display table; the downloads and plots keep the in-memory frames. After an Open, the table references are absent until the normalisation or enrichment step runs again, and the renders fall back to the in-memory frames (the RUV grid is bounded by the sweep; the enrichment frames are the module's as before).

## Harness: the single journey member (KIT-018)

`tickets$018 = list(kind = "single_journey", build = "kit_018", corpus = "full", stop_after = "design")` expands to one journey case `member` (candidate, no launch record, the production default) whose driver runs the journey from import through `stop_after` only (`workflow_kit_driver.R`); the stop step is part of the case binding; the verifier expects the steps up to it and judges their states against the KIT-005 `journey_<step>` oracles (`journey_import`, `journey_design`); `kitSingleJourneyDecision` accepts on eligibility, completion of those steps, complete populations and an exact verification, rejects on an oracle mismatch named by the verifier, and is invalid otherwise; peak recorded, no ratio, no control. Harness tests: `test-workflow-kit-harness.R` 325 (new: expansion, binding, refusals, decision outcomes).

## Tests (working tree at `6a7304e7`, proof library on `R_LIBS`)

| Suite | Result |
|---|---|
| test-workflow-dia-products.R | 64 pass (new: display table save/page/filter; RUV result from the product) |
| test-workflow-native-session.R | 500 pass, 1 pre-existing failure+error ("native publication recovers after process loss": needs an installed package under load_all, baseline) |
| test-workflow-kit-harness.R | 325 pass |
| test-prot-05b-norm-module-contracts.R | 562 pass |
| test-prot-11-enrichment-module-contracts.R | 1253 pass (three mocks widened for the new seams) |
| test-prot-11a-enrichment-module-characterization.R | 143 pass |
| test-module-ci-prot-enrich.R | 96 pass, MCI-010.2 pre-existing error |
| test-prot-12-summary-module-contracts.R | 446 pass |
| test-general-filemgmt-config-direct-shared / roundtrip / s4 | 33 / 24 / 30 pass |
| test-metab-05-summary-module-characterization.R | 160 pass, 9 skips |
| test-prot-04a-design-module-characterization.R | 230 pass |

## Counted member

Pending: one single journey member (import to design, 300-sample corpus, kit_018 build, exactness against `journey_import` and `journey_design`, completion, peak recorded).
