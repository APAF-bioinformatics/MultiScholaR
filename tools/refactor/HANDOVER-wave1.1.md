# Wave 1.1 Handover

## Goal

Stabilize the first proteomics DA split before further decomposition.

Wave 1 applied the DA split into live `R/`, but several resulting files still
exceed the file-size budget. Wave 1.1 adds characterization tests around the DA
wrapper, DA init handlers, and DA results builders so those files can be split
again safely.

This handover also captures the repo-local `renv` test environment and the
normal verification entrypoints, so a later session can resume without needing
chat history.

## Canonical Context

- Refactor protocol:
  [PLAYBOOK.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/PLAYBOOK.md:1)
- Stabilization queue:
  [GOD_MODULE_STABILIZATION_BACKLOG.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/GOD_MODULE_STABILIZATION_BACKLOG.md:1)
- Current import seam handover:
  [HANDOVER-import-server-seams.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/HANDOVER-import-server-seams.md:1)
- Current normalization seam handover:
  [HANDOVER-norm-server-seams.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/HANDOVER-norm-server-seams.md:1)
- Repo size baseline:
  [AUDIT-file-sizes.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/AUDIT-file-sizes.md:1)
- Wave 1 size baseline:
  [AUDIT-wave1-file-sizes.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/AUDIT-wave1-file-sizes.md:1)
- Filename-coupling audit:
  [AUDIT-filename-coupling.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/AUDIT-filename-coupling.md:1)
- Wave 1 apply handover:
  [HANDOVER-wave1.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/HANDOVER-wave1.md:1)
- `renv` bootstrap script:
  [tools/bootstrap_renv.R](/home/doktersmol/Documents/MultiScholaR/tools/bootstrap_renv.R:1)
- `renv` test runner:
  [tools/test_with_renv.R](/home/doktersmol/Documents/MultiScholaR/tools/test_with_renv.R:1)
- Installed Codex skill:
  [god-module-stabilization](</home/doktersmol/.codex/skills/god-module-stabilization/SKILL.md>)
- Lockfile:
  [renv.lock](/home/doktersmol/Documents/MultiScholaR/renv.lock:1)
- Project env config:
  [.Renviron](/home/doktersmol/Documents/MultiScholaR/.Renviron:1),
  [.Rprofile](/home/doktersmol/Documents/MultiScholaR/.Rprofile:1),
  [.renvignore](/home/doktersmol/Documents/MultiScholaR/.renvignore:1)

## Live Wave 1 State

Wave 1 is already applied to live `R/`.

Wave 1.1 is now applied live for the proteomics DA follow-up slice.

Completed live splits:

- the former [R/func_prot_da_results.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_results.R:1)
  monolith was split into:
  - [R/func_prot_da_results_io.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_results_io.R:1) `529`
  - [R/func_prot_da_results_long_format.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_results_long_format.R:1) `204`
  - [R/func_prot_da_results_summary.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_results_summary.R:1) `359`
  - [R/func_prot_da_results_methods.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_results_methods.R:1) `495`
  - [R/func_prot_da_results.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_results.R:1)
    is now a 7-line breadcrumb stub
- the former [R/mod_prot_da_handlers_analysis.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_da_handlers_analysis.R:1)
  monolith was split into:
  - [R/mod_prot_da_handlers_init.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_da_handlers_init.R:1) `335`
  - [R/mod_prot_da_handlers_load_session.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_da_handlers_load_session.R:1) `261`
  - [R/mod_prot_da_handlers_run_analysis.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_da_handlers_run_analysis.R:1) `516`
  - [R/mod_prot_da_handlers_analysis.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_da_handlers_analysis.R:1)
    is now a 6-line breadcrumb stub
- the former [R/func_prot_da_model.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_model.R:1)
  monolith was split into:
  - [R/func_prot_da_model_wrapper.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_model_wrapper.R:1) `611`
  - [R/func_prot_da_model_volcano_prep.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_model_volcano_prep.R:1) `126`
  - [R/func_prot_da_model_stats.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_model_stats.R:1) `578`
  - [R/func_prot_da_model_methods.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_model_methods.R:1) `596`
  - [R/func_prot_da_model.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_model.R:1)
    is now a 7-line breadcrumb stub
- the former [R/func_prot_da_volcano.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_volcano.R:1)
  monolith was split into:
  - [R/func_prot_da_volcano_glimma.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_volcano_glimma.R:1) `421`
  - [R/func_prot_da_volcano_static.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_volcano_static.R:1) `130`
  - [R/func_prot_da_volcano_write.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_volcano_write.R:1) `190`
  - [R/func_prot_da_volcano_main.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_volcano_main.R:1) `111`
  - [R/func_prot_da_volcano.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_volcano.R:1)
    is now a 7-line breadcrumb stub

Current state of the former oversized Wave 1.1 files:

- no remaining Wave 1.1 DA follow-up file is above the `801-1000` soft-cap band
- the largest remaining files in this slice are now:
  - [R/func_prot_da_model_wrapper.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_model_wrapper.R:1) `611`
  - [R/func_prot_da_model_methods.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_model_methods.R:1) `596`
  - [R/func_prot_da_model_stats.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_model_stats.R:1) `578`
  - [R/func_prot_da_results_io.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_results_io.R:1) `529`
  - [R/mod_prot_da_handlers_run_analysis.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_da_handlers_run_analysis.R:1) `516`

These are acceptable under the current budget, but are candidates for a later
“push acceptable toward ideal” pass if we want to tighten them further.

## Tests Added In Wave 1.1

### 1. DA Compat and Handler Characterization

File:

- [tests/testthat/test-prot-07b-da-handlers-compat.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-prot-07b-da-handlers-compat.R:1)

Coverage added:

- `writeInteractiveVolcanoPlotProteomicsMain()` accepts legacy `de_*` aliases
- clear error for missing DA/DE analysis result inputs
- `da_server_init_handlers()` seeds DA state from workflow state
- contrast normalization from `A-B` to `groupA-groupB`

Extra compatibility note:

- this test now uses a local `with_function_overrides()` helper so it can
  safely override locked bindings when run through `pkgload::load_all()`

### 2. DA Results Characterization

File:

- [tests/testthat/test-prot-07c-da-results-characterization.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-prot-07c-da-results-characterization.R:1)

Coverage added:

- `getTypeOfGrouping()`
- `extractResults()`
- `countStatDaGenesHelper()`
- `getSignificantData()`
- `createDaResultsLongFormat()`

## Completed Import Slice And Next Active Target

Proteomics import stabilization is now complete in live `R/`.

Current live import state:

- the `func_prot_import.R` helper wave is already applied live into:
  - [R/func_prot_import_detection.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_import_detection.R:1)
  - [R/func_prot_import_readers.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_import_readers.R:1)
  - [R/func_prot_import_tmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_import_tmt.R:1)
  - [R/func_prot_import_diann_format.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_import_diann_format.R:1)
  - [R/func_prot_import_organisms.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_import_organisms.R:1)
- wrapper apply wave 1 is live in:
  - [R/mod_prot_import_ui_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_import_ui_helpers.R:1)
  - [R/mod_prot_import_detection_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_import_detection_helpers.R:1)
  - [R/mod_prot_import_path_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_import_path_helpers.R:1)
- wrapper apply wave 2 is live in:
  - [R/mod_prot_import_processing_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_import_processing_helpers.R:1)
  - [R/mod_prot_import_config_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_import_config_helpers.R:1)
  - [R/mod_prot_import_organism_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_import_organism_helpers.R:1)
  - [R/mod_prot_import_orchestration_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_import_orchestration_helpers.R:1)
- wrapper apply wave 3 is live in:
  - [R/mod_prot_import_ui.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_import_ui.R:11)
  - [R/mod_prot_import_server.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_import_server.R:18)
- [R/mod_prot_import.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_import.R:1)
  is now the extractor breadcrumb shell at `64` lines
- [tools/refactor/HANDOVER-import-server-seams.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/HANDOVER-import-server-seams.md:1)
  now serves as the archival import handover rather than an active seam map

Verification already done for the completed import state:

```bash
Rscript tools/test_with_renv.R tests/testthat/test-prot-01-import.R
Rscript tools/test_with_renv.R tests/testthat/test-prot-01b-import-detection-characterization.R
Rscript tools/test_with_renv.R tests/testthat/test-prot-01c-import-module-contracts.R
Rscript tools/test_with_renv.R tests/testthat/test-format-diann.R
Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-import-server-wave2.yml
Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-import-server-wave3.yml
```

Observed result:

- import gate green
- `test-prot-01-import.R`: passed with one expected skip for missing LFS binary
- wave 2 and wave 3 post-apply checks passed
- `mod_prot_import_ui()` and `mod_prot_import_server()` now live in dedicated
  files with the public entry points preserved

The next active god-module target is now proteomics normalization:

- [R/mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1)
- [R/func_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_norm.R:1)
- the first normalization helper apply wave is now live in:
  - [R/mod_prot_norm_qc_support_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_qc_support_helpers.R:1)
  - [R/mod_prot_norm_qc_generation_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_qc_generation_helpers.R:1)
- the second normalization helper apply wave is now live in:
  - [R/mod_prot_norm_workflow_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_workflow_helpers.R:1)
  - [R/mod_prot_norm_ruv_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_ruv_helpers.R:1)
- the first shared correlation-completion seam is now introduced in-place in
  [R/mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:508)
  for correlation result naming, persistence, final QC/filtering updates,
  workflow finalization, and summary text
- the next core correlation seam is now introduced in-place in
  [R/mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:754)
  for correlation input-object resolution, correlation-vector calculation,
  threshold application, and skipped-correlation state routing
- the outer correlation observer seam is now introduced in-place in
  [R/mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:859)
  for apply/skip workflow runners, correlation completion, and shared error
  handling
- the export session seam is now introduced in-place in
  [R/mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1068)
  for export readiness, source-dir resolution, session-data assembly,
  artifact saving, summary writing, and export error handling
- [R/mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1)
  is now `2066` lines and still owns the public wrapper at
  [mod_prot_norm_server()](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1343)

## Verification Done

Both new test files parse.

Targeted sourced-file harness runs passed before the `renv` runner became the
canonical path:

### Compat and Handler Harness

```bash
Rscript -e 'library(testthat); library(shiny); library(tibble); library(dplyr); library(rlang); source("R/allGenerics.R"); source("R/func_general_helpers.R"); source("R/func_pept_s4_objects.R"); source("R/func_prot_s4_objects.R"); source("R/func_lipid_s4_objects.R"); source("R/func_prot_da_volcano.R"); source("R/mod_prot_da_handlers_analysis.R"); testthat::test_file("tests/testthat/test-prot-07b-da-handlers-compat.R")'
```

### Results Harness

```bash
Rscript -e 'library(testthat); library(dplyr); library(tidyr); library(tibble); library(purrr); library(stringr); library(rlang); library(ggplot2); source("R/allGenerics.R"); source("R/func_pept_s4_objects.R"); source("R/func_prot_s4_objects.R"); source("R/func_prot_da_results.R"); testthat::test_file("tests/testthat/test-prot-07c-da-results-characterization.R")'
```

## Repo-Local `renv` State

`renv` is now bootstrapped and synchronized for this repo.

Important files:

- [renv.lock](/home/doktersmol/Documents/MultiScholaR/renv.lock:1)
- [tools/bootstrap_renv.R](/home/doktersmol/Documents/MultiScholaR/tools/bootstrap_renv.R:1)
- [tools/test_with_renv.R](/home/doktersmol/Documents/MultiScholaR/tools/test_with_renv.R:1)

Current behavior:

- `tools/bootstrap_renv.R` hydrates from existing local libraries first, then
  installs remaining CRAN/Bioconductor/remote packages into the project library
  non-transactionally, then snapshots `renv.lock`
- `tools/test_with_renv.R` activates `renv`, runs `pkgload::load_all(..., export_all = TRUE)`,
  then runs either `testthat::test_dir()` or `testthat::test_file()`
- normal Wave 1.1 verification should now go through the `renv` runner, not
  through ad hoc sourced-file harnesses

Exact verification done under `renv`:

```bash
Rscript tools/test_with_renv.R tests/testthat/test-prot-07b-da-handlers-compat.R
Rscript tools/test_with_renv.R tests/testthat/test-prot-07c-da-results-characterization.R
Rscript -e 'library(renv); s <- renv::status(); cat("SYNC=", s$synchronized, "\n", sep="")'
```

Observed result:

- `test-prot-07b-da-handlers-compat.R`: passed
- `test-prot-07c-da-results-characterization.R`: passed
- `renv::status()`: `SYNC=TRUE`
- all `R/*.R` files parse after the full Wave 1.1 DA follow-up split

Bootstrap lessons already folded into tooling:

- the first bootstrap attempt failed to leave a usable library because
  `renv::install()` was transactional and a `V8` failure prevented successful
  packages from being linked
- `tools/bootstrap_renv.R` now uses `transactional = FALSE`
- the first hydrate attempt missed the user library because activated `renv`
  overrides `R_LIBS_USER`; `tools/bootstrap_renv.R` now derives the default
  user library path explicitly

Known runtime noise:

- `timedatectl` / system bus warnings may still appear during test startup in
  this environment; they did not block the targeted tests
- package startup prints the MultiScholaR dependency banner during `load_all()`

## Not Yet Done

- No full-package `devtools::document()` run
- No full package test suite run under `renv`
- No post-Wave-1.1 size-audit refresh has been generated yet

## Exact Next Step

Wave `1.1` itself is complete, the import slice is fully applied live, and the
first two normalization helper waves are now also applied live.

The next exact step is to continue the now-started proteomics normalization
wrapper stabilization:

1. rerun the normalization gate
2. continue from
   [HANDOVER-norm-server-seams.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/HANDOVER-norm-server-seams.md:1)
3. treat the normalization, correlation, and export observers as seam-ready,
   then continue with the reset observer in `mod_prot_norm.R`
4. stop again at a clean seam boundary

## Resume Guardrails

Use these guardrails before continuing with the next god-module target:

- Follow [PLAYBOOK.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/PLAYBOOK.md:1) literally for structural changes:
  - move code by exact source extraction, not by hand-rewriting bodies
  - stage extraction output outside live `R/` first
- Prefer the installed
  [god-module-stabilization](</home/doktersmol/.codex/skills/god-module-stabilization/SKILL.md>)
  skill for survey/stabilize workflow orchestration. It wraps the repo's
  existing `tools/refactor/` scripts and records the current protocol.
- Keep newly introduced helper/function names in the local camelCase style.
- Do not create new live `R/` helper files during exploratory work. If a new helper file is needed, create it in staging first and only promote the reviewed version into `R/`.
- If a targeted test exposes a real runtime bug, fix that bug minimally in live code, then return to the staged extraction workflow for any structural breakup.

When resuming normalization work, start by re-running the current gate:

```bash
Rscript tools/test_with_renv.R tests/testthat/test-prot-05-normalisation.R
Rscript tools/test_with_renv.R tests/testthat/test-prot-05b-norm-module-contracts.R
Rscript tools/test_with_renv.R tests/testthat/test-prot-06-ruv.R
```

Then continue from the current next target:

1. treat proteomics import as complete and archived in
   [HANDOVER-import-server-seams.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/HANDOVER-import-server-seams.md:1)
2. continue proteomics normalization from
   [HANDOVER-norm-server-seams.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/HANDOVER-norm-server-seams.md:1):
   - `mod_prot_norm.R` is classified `high-risk-wrapper` / `needs-seam-introduction`
   - `func_prot_norm.R` is classified `direct-extraction-ready`
   - the first normalization QC helper wave is applied live in
     `mod_prot_norm_qc_support_helpers.R` and
     `mod_prot_norm_qc_generation_helpers.R`
   - the second normalization workflow/RUV helper wave is applied live in
     `mod_prot_norm_workflow_helpers.R` and
     `mod_prot_norm_ruv_helpers.R`
   - the shared, core, and outer correlation shells now delegate through new
     in-file helpers in `mod_prot_norm.R`
   - the export session shell now also delegates through new in-file helpers
     in `mod_prot_norm.R`
   - `test-prot-05b-norm-module-contracts.R` now covers the new helpers directly
   - the next safe target is `observeEvent(input$reset_normalization, ...)`
3. if normalization progress is blocked, do a cleanup pass on the remaining acceptable-but-not-ideal DA files:
   - [R/func_prot_da_results_io.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_results_io.R:1) `529`
   - [R/mod_prot_da_handlers_run_analysis.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_da_handlers_run_analysis.R:1) `516`
   - [R/func_prot_da_model_wrapper.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_model_wrapper.R:1) `611`
   - [R/func_prot_da_model_stats.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_model_stats.R:1) `578`
   - [R/func_prot_da_model_methods.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_da_model_methods.R:1) `596`

## Compact Resume Set

- [tools/refactor/HANDOVER-wave1.1.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/HANDOVER-wave1.1.md:1)
- [tools/refactor/GOD_MODULE_STABILIZATION_BACKLOG.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/GOD_MODULE_STABILIZATION_BACKLOG.md:1)
- [tools/refactor/PLAYBOOK.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/PLAYBOOK.md:1)
- [tools/refactor/HANDOVER-import-server-seams.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/HANDOVER-import-server-seams.md:1)
- [tools/refactor/HANDOVER-norm-server-seams.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/HANDOVER-norm-server-seams.md:1)
- [tools/bootstrap_renv.R](/home/doktersmol/Documents/MultiScholaR/tools/bootstrap_renv.R:1)
- [tools/test_with_renv.R](/home/doktersmol/Documents/MultiScholaR/tools/test_with_renv.R:1)
- [renv.lock](/home/doktersmol/Documents/MultiScholaR/renv.lock:1)
- [.Renviron](/home/doktersmol/Documents/MultiScholaR/.Renviron:1)
- [.renvignore](/home/doktersmol/Documents/MultiScholaR/.renvignore:1)
- [tests/testthat/test-prot-05-normalisation.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-prot-05-normalisation.R:1)
- [tests/testthat/test-prot-05b-norm-module-contracts.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-prot-05b-norm-module-contracts.R:1)
- [tests/testthat/test-prot-06-ruv.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-prot-06-ruv.R:1)
- [tests/testthat/test-prot-01-import.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-prot-01-import.R:1)
- [tests/testthat/test-prot-01c-import-module-contracts.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-prot-01c-import-module-contracts.R:1)
- [god-module-stabilization](</home/doktersmol/.codex/skills/god-module-stabilization/SKILL.md>)

No tests were rerun after the final doc-only updates in this handover.
