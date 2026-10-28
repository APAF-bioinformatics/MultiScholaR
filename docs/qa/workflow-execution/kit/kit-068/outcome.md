# KIT-068: the lipidomics and metabolomics routes on the kit (tooling only)

Status 2026-10-29: **completed**; the three step-map flags were answered by the user on 2026-10-29 as the determinations read (recorded in the maps' rows as `decision`, `open_flags = 0`; bindings regenerated for the maps' hashes). Date 2026-10-29. Authority: the series plan, the stage-3 re-cut review (`../REVIEW-20261028-stage3-recut.md`). Control and build under test: `kit_065` (`core/browser-smoke.InuLC5`, f55bef1b). No package code changed; no install; no counted member.

## What was built

| Item | Where | Commit |
|---|---|---|
| Synthetic corpus generator: `lipid_screen` (two LipidSearch exports LCMS_Pos/LCMS_Neg, 200 lipids + 10 duplicate-adduct rows + 6 standards per assay), `metab_screen` (two MS-DIAL exports LCMS_Pos/GCMS, 300 features + 6 standards per assay); 24 samples = 12 biological x 2 technical replicates, 3 batches, 2 groups, 15 % differential, a batch effect, no missing values; `design.tsv` with `sample, group, batch, tech_rep_group, replicates`; `manifest.json` with per-file sha256 | `tools/profiling/kit_omics_corpus.R`; generated under `<proof>/kit/corpus/lipid-screen` and `metab-screen` by the bindings generator | 5d79cc4f |
| Routes `lipidomics` and `metabolomics` in `kitJourneySteps` (ids `import, design, intensity, duplicates, qc_finalize, normalisation, correlation, session_export, da, report`; kinds `store, omics_qc, omics_norm, da, report`; the normalisation state follows the bound RUV mode); the stage grammar `kit-lifecycle-(lipid|metab)[-fail]` in `kitCases`, `kitMemberSeconds`, `kitSelectorStage`, the coordinator and the runner; the route-aware `kitLifecycleDecision` and the verifier's failure branch for an in-memory member (no store: the head rule recorded as not applicable); `lifecycle_omics` kept by `kitBindings` | `tools/profiling/workflow_kit.R`, `workflow_kit_coordinator.R`, `run_workflow_kit.sh` | 5d79cc4f, 56f1112c, a7d56839 |
| Driver: `runKitOmicsImportActions` (tile, setup, vendor radio, assay names and files through the module's shinyFiles buttons, sanitised names off, process); `runKitDesignActions` generalised on prefix, tab input and the save wait (the DIA default unchanged); `runKitOmicsStep` for every step through the report, each action waiting on the module's own completion text | `tools/profiling/workflow_kit_driver.R` | 5d79cc4f |
| Coordinator: the corpus's `assays` (by assay name) and `vendor`; the fault from the case (`fault_spec`) | `tools/profiling/workflow_kit_coordinator.R` | 5d79cc4f |
| Step maps with line evidence and three flags each | `kit-001/step-map-lipid.json`, `step-map-metab.json` | 5d79cc4f |
| Bindings: corpora `lipid_screen`, `metab_screen`; `lifecycle_omics` (ceiling 1,800 s, build `kit_065`, fault `guard_after_qc_finalize`, normalisation ITSD + cyclic loess + automatic RUV + correlation 0.5); `step_map_lipid(_sha256)`, `step_map_metab(_sha256)` | `bindings-generator.R`, `bindings.json` | 5d79cc4f and the regenerations |
| Test: the routes' step lists and states, the four stages' cases, selectors, ceilings, the decisions on in-memory members | `tests/testthat/test-workflow-kit-harness.R` (471 pass, the two known skips) | 5d79cc4f, 56f1112c |
| Stage summary script extended to lifecycle stages | `kit-053/stage-summary.R` | 5d79cc4f |

The generator's invocation for the regenerations: every bound build as `kit_NNN=<install dir>`, `kit_journey=<proof>/kit/runtime/kit-journey`, `KIT_065_LIFECYCLE=1`, with the kit_065 library on `R_LIBS` (the journey ingestion loads the package). The diff against the previous bindings was the additions, the harness file hashes and the generation dates only.

## The four lifecycle members on kit_065

| Stage | Decision | Steps | Elapsed | Anonymous peak | Notes |
|---|---|---|---|---|---|
| `kit-lifecycle-lipid` | **lifecycle_complete_observed** | 10 of 10 through the report | 226 s | 2.22 GiB | first attempt; verifier passed; `exact = NA` (in-memory) |
| `kit-lifecycle-lipid-fail` | **lifecycle_failure_observed** (attempt 3) | 5 complete, `normalisation` failed under the injected guard | 43 s | | attempts 1 and 2 were harness defects (below) |
| `kit-lifecycle-metab` | **lifecycle_complete_observed** (attempt 2, the user's single rerun) | 10 of 10 through the report; verifier passed | 210 s | 2.19 GiB | attempt 1 completed every step and verified but its sampler measure was invalid (10 gaps above 250 ms, max 1.05 s, inside the RUV sweep; archived); the rerun's max interval 0.07 s, no gaps |
| `kit-lifecycle-metab-fail` | **lifecycle_failure_observed** | 5 complete, `normalisation` failed under the injected guard | 41 s | | first attempt |

Summaries: [kit-lifecycle-lipid-stage.json](kit-lifecycle-lipid-stage.json), [kit-lifecycle-lipid-fail-stage.json](kit-lifecycle-lipid-fail-stage.json), [kit-lifecycle-metab-stage.json](kit-lifecycle-metab-stage.json), [kit-lifecycle-metab-fail-stage.json](kit-lifecycle-metab-fail-stage.json); the archived attempts [kit-lifecycle-lipid-fail.attempt1-no-launch-record-stage.json](kit-lifecycle-lipid-fail.attempt1-no-launch-record-stage.json), [kit-lifecycle-lipid-fail.attempt2-sampler-na-uptime-stage.json](kit-lifecycle-lipid-fail.attempt2-sampler-na-uptime-stage.json) and [kit-lifecycle-metab.attempt1-sampler-gaps-stage.json](kit-lifecycle-metab.attempt1-sampler-gaps-stage.json). The stage directories (and the two renamed attempts) stay under `<proof>/kit/runtime/`.

What the complete members produced (checked in the member projects, synthetic data): the design's cleaned tables carry only the id, annotation and sample columns; ITSD applied (median aggregation), log2 offset 1, cyclic loess, the automatic RUV sweep (lipid: best k 2 at 4 % and 7 %; metab: best k 1 at 3 % and 7 %), the correlation filter kept all 24 samples; retained features lipid 183 / 182, metab 281 / 283; the DA long tables per assay (`de_posmode_*`, `de_negmode_*` / `de_gcms_*`), the publication workbook and the rendered `.docx` report under `results_summary/<omic>/`.

### The two harness corrections (the failure path)

1. **Attempt 1 (`invalid_guard_injection`).** The application permits the guard injection only under a journey launch record with `guard_stop_test` (`workflowGuiLaunchRecord`); the omics case carried none (`record = FALSE`, the production default). The failure member now carries the record; the complete member keeps the production default, which the in-memory modules need nothing from. Commit 56f1112c.
2. **Attempt 2 (`time_guard` at 46 s).** With the injection accepted, the sampler took a sample as the cgroup vanished under the application's guard stop, read no uptime, and the guard comparison `"NA" >= <deadline>` was a string comparison that fired the time guard. The sampler now skips a sample without a numeric uptime. Commit a7d56839. This race could have hit any failure member on a fast stop; the DIA failure member of KIT-065 did not meet it.

Both are harness defects repaired within the ticket's correction rule (a driver or instrumentation defect, repaired and rerun; a module failure is a finding). No module failed in any of the four members.

### The metabolomics complete member's measure (attempt 1) and the single rerun

Attempt 1's decision was `invalid: measurement_valid`: the rule requires every sampler interval at or below 500 ms, and the trace had 10 intervals above 250 ms, two above 500 ms and one of 1.05 s, all inside the 105 s normalisation step (the RUV sweep); the lipid member minutes earlier had none, and the metabolomics code runs no parallel workers, so the cause was host contention during that window. The member's completion and verification stood. Under the ticket's correction rule (driver failures only) it was not rerun by default; the user authorised a single rerun (2026-10-29). The rerun completed every step through the report in 210 s with a maximum sampler interval of 0.07 s and no gaps: `lifecycle_complete_observed`. Attempt 1 is archived (`kit-lifecycle-metab.attempt1-sampler-gaps`), a valid instrumentation failure preserved, not a module finding.

## The step maps' flags (answered 2026-10-29)

Three flags, the same on both maps; the user accepted each as read ("recommended" on all three, 2026-10-29). KIT-031/032 declare these units:

1. **Intensity filter**: row unit after a per-assay quantile pre-pass (one global statistic per assay), with `ceiling()` removed in the declaration (user decision 2026-10-28). The second half of the flag as written, whether the design's cleaned tables carry only the id, annotation and sample columns, is answered by the members' projects: they do (`LipidName, LipidClass, <samples>`; `Alignment ID, Metabolite name, <samples>`), so the quantile sees intensities only.
2. **Internal-standard normalisation**: complete sample with the mean-factor pre-pass (one aggregate per sample, the mean over samples), the standards' row removal in the same operation.
3. **Correlation filter**: whole matrix, the DIA step map's determination for the same statistic.

Resolved before the maps were written: the DA output is the manifest (`lipid.da_results.v1` / `metab.da_results.v1`, the `dia.da_results.v1` pattern), user confirmed 2026-10-29; recorded in the maps' `da` rows without a flag.

Observations recorded in the maps for KIT-031/032 (not flags): the RUV-III replicate structure comes from `ruv_grouping_variable`, not the declared technical replicates (the series' open replicate rule); MS-DIAL alignment ids are unique, so the duplicate step is the identity on that layout; the kit drives the design builder (group, technical replicates from `replicates`, one contrast), not the design import modal, and the corpus's `batch` column is recorded in `design.tsv` only.

## Claim limits

Synthetic corpora; no result on them is a scientific claim. The members prove the instrumentation observes the two routes through completion and failure on the in-memory modules (`exact = NA`, no store, no populations, no head rule). No memory claim.

## Hand-off

KIT-030, 069, 070, 071 take the driver and the corpora for their before/after checks (`run kit-lifecycle-lipid` and `kit-lifecycle-metab` on each install; a store will appear from KIT-070 on, at which point `kitFindStore` and the fingerprint gate need the generalisation KIT-070 names). KIT-031/032 take the step maps once the three flags are answered.
