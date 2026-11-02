# KIT-076 plan

**Work item.** Lipid and metab DA as one saved operation per contrast on a saved per-assay fit, with manifests v2 and
the record labelled by contrast. Exit: the frozen old-vs-new rule holds on every capture, 076l/076m exact against new
oracles, v1 records still replay, one failing contrast leaves the others complete.

## Design decisions

- **Kernel halves (moved code).** `fit<Omic>DaAssayModel` (design, sample intersect, `lmFit`) and
  `apply<Omic>DaContrasts` (term validation, `makeContrasts`, `contrasts.fit`, `eBayes`, `treat`, extraction, q-values;
  takes one or more contrasts). `runTestsContrasts*DA` calls both with all its contrasts, so `runLipidsDA` /
  `runMetabolitesDA` and their all-contrast `fit.eb` are unchanged.
- **Deviation from the ticket's split.** The ticket puts contrast-term validation in the fit half. It moves to the
  contrast half, since item 6 needs an invalid contrast to fail alone and a saved fit must not depend on the contrasts.
- **Also moved out of the orchestrators:** the assay input (expression matrix, annotation) and the per-contrast table
  decoration, so both paths build the same tables.
- **Declarations.** `<prefix>_da_contrast` (whole matrix, one per contrast). It finds or saves one `da_fit` per assay
  (`newWorkflowOmicsDaFitIdentity`, new). `<prefix>_da_output` stays registered, for v1 records; the wrapper no longer
  runs it.
- **Failures.** A contrast failure on an assay is caught in the kernel and recorded under `errors`; store errors are not
  caught. A fit failure drops the assay as before (recorded, not fixed). If every contrast fails the run stops with
  their messages, so the lipid preflight's run-level stop survives now that it runs per contrast.
- **Consumers.** `readWorkflowOmicsDaResults` rebuilds the kernel's list from the manifest. `fit.eb` in that list is a
  coefficient view (`list(coefficients = )`, one column per contrast): Glimma reads only its names and row order, and
  the test fakes already use that shape. The full per-contrast `fit.eb` stays in the saved product. So the exports,
  Glimma, volcanoes, heatmaps, selectors and tables need no edit. A failed contrast adds a warning notification.
- **Record.** v2 entries are one per contrast and assay (`contrast`, `assay`, `entity_id`, `rows`). v1 reading is
  unchanged. Replay accepts either operation and always runs per contrast.

## Scientific proof, frozen before any new code runs

The rule is the ticket's: `~ 0 + group` `identical()`; other designs absolute difference <= 1e-12 with identical calls,
rows and order; anything beyond it is a stop.

**The ticket's named inputs cannot falsify the claim alone.** The screen corpora and e2e lanes each have one contrast,
so one contrast alone is the same computation as v1. Scenarios added on the same data, captured on `kit_075`:

| Scenario | Data | Formula | Contrasts | Rule |
|---|---|---|---|---|
| `screen` | screen head state | `~ 0 + group` | recorded (1) | identical |
| `multi` | same, group x batch (6 levels) | `~ 0 + group` | 3, sharing terms | identical |
| `multi_treat` | same as `multi` | `~ 0 + group` | 3, lfc 0.5 | identical |
| `nonorth` | screen head state + batch | `~ 0 + group + batch` | 2 | <= 1e-12 |
| `e2e` | `lipid_canonical`, `metab_combined`, log2 | `~ 0 + group` | recorded (1) | identical |

Compared: `runLipidsDA` / `runMetabolitesDA` directly (old vs new, whole list); the wrapper's long tables, per-assay
tables and counts; each per-contrast `fit.eb` column against v1's; exported TSV and workbook content.

## Files

Kernels `func_lipid_da_model.R`, `func_metab_da_core_engine.R`, `func_metab_da_orchestration.R`; core
`utils_workflow_omics_declarations.R`; bindings `utils_workflow_{lipid,metab}_declarations.R`; registry maps
`utils_workflow_operations.R`; `utils_workflow_class_declarations.R`; `utils_workflow_analysis_record.R`; finalizers
`mod_lipid_da_analysis_helpers.R`, `mod_metab_da_observer_helpers.R`; recorders only if their reads need it; tests.

## Gates

Stage-3b exactness set against new oracles `lipid76` / `metab76`; proteomics part (record code changes); v1 member
records replay; targeted tests, dependents grep, R CMD check, install. No full suite (next checkpoint is KIT-078).
No counted member.
