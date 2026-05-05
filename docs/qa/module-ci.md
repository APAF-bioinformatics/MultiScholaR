# Module CI Corpus

The module-CI corpus is the per-module layer underneath the browser E2E suite.
E2E tests prove standard workflows compose through the GUI. Module-CI tests prove
that individual modules preserve data integrity across supported parameters,
edge cases, state transitions, and artifact schemas.

## Manifest

The manifest lives at `tests/testdata/module_ci/manifest.json`. Each scenario has
these required fields:

- `scenario_id`
- `ticket_id`
- `item_id`
- `omic`
- `module_family`
- `runtime`
- `ci_lane`
- `fixture`
- `parameters`
- `expected_state`
- `expected_artifacts`
- `required_tests`

Scenario IDs start with their implementation item, for example
`MCI-018.1-lipidomics-import-schema-smoke`. This makes the manifest directly
sliceable at ticket or item boundary.

## Runtime Classes

- `unit-contract`: pure helper/module tests without a browser.
- `module-browser`: isolated browser tests for one module or module family.
- `module-artifact`: generated file and schema fidelity tests.
- `workflow-e2e`: full workflow browser scenarios.
- `release-full`: exhaustive release-gate scenarios.

## Oracles

Shared oracles live in `tests/testthat/helper-module-ci-oracles.R` and validate:

- sample identity;
- feature identity;
- assay provenance;
- design/sample alignment;
- finite numeric outputs;
- duplicate keys;
- non-empty files and table schemas;
- deterministic table checksums for selected input columns;
- expected per-column missingness;
- expected non-finite value sentinels.

## Fixture Packs

`tests/testdata/module_ci/fixture_packs.json` is the fixture catalog added for
MCI-002. It is generated from `tests/testdata/module_ci/generate-fixtures.R`
and is intentionally separate from the workflow E2E manifest so module-level
edge cases do not disturb stable browser fixtures.

The catalog covers every current omic/import form:

- proteomics: DIA, DDA, TMT, and DIA+LIMPA;
- metabolomics: LC and GC;
- lipidomics: LC and GC.

Each omic/import form has all required fixture classes:

- `happy_path`: clean, ordinary module input.
- `missingness`: deterministic missing-value locations and expected row-filter
  candidates.
- `duplicates`: duplicate feature IDs and expected duplicate-key validation.
- `invalid_design`: duplicated/missing design samples and expected design
  validation errors.
- `multi_factor`: group, batch, sex, site, and age design covariates.
- `bad_names`: unsafe sample and feature names that exercise name sanitisation.
- `nonfinite`: `Inf`/`NaN` sentinels that must be rejected before modelling.
- `small_n`: one replicate per group for insufficient-replicate guards.
- `large_enough_for_plots`: plot-scale fixtures labelled nightly-only.
- `multi_assay_mixed`: two-assay payloads for assay provenance checks.

Each pack includes:

- a data TSV under `fixtures/<omic>/`;
- a design TSV under `fixtures/<omic>/`;
- a JSON oracle under `oracles/<omic>/`;
- `ci_lane` and `push_safe` labels;
- source generator metadata;
- expected exported schema columns;
- selected-column MD5 checksums for data and design corruption detection.

Push CI should select `ci_lane = "push"` and `push_safe = true`. Nightly and
release gates may opt into `ci_lane = "nightly"` for larger plotting fixtures.
The test helper rejects non-push-safe packs unless they are explicitly in the
nightly lane.

## Local Runner

The foundation suite is:

```r
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(export_all = TRUE, helpers = FALSE, attach_testthat = FALSE)
testthat::test_dir(
  "tests/testthat",
  filter = "^(module-ci-manifest|module-ci-oracles|module-ci-harness)$",
  reporter = "summary"
)
```

Later MCI tickets add omic-specific scenario files and CI wrappers.

The MCI-002 fixture/oracle suite is:

```r
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(export_all = TRUE, helpers = FALSE, attach_testthat = FALSE)
testthat::test_dir(
  "tests/testthat",
  filter = "^(module-ci-fixture-integrity|module-ci-oracles|e2e-fixture-integrity)$",
  reporter = "summary"
)
```

## CI Partitions

MCI-003 adds local CI entrypoints and the first GitHub Actions workflow for this
corpus:

- `tools/ci/run-module-ci.R` selects module-CI scenarios by `--omic`,
  `--module`, `--runtime`, and `--scenario`, runs the mapped `testthat` files,
  and writes `module-ci-run-manifest.json` into the artifact directory.
- `tools/ci/module-ci-scorecard.R` reads the run manifest and emits
  `module-ci-scorecard.json` plus `module-ci-scorecard.md`. The scorecard counts
  `passed`, `failed`, `skipped`, and `unknown` separately; skipped browser
  scenarios are never counted as passed.
- `.github/workflows/module-ci.yml` defines push/PR jobs, browser smoke, current
  E2E smoke, and nightly/release full gates.

Required push/PR jobs:

- `module-ci unit (foundation)`: foundation manifest/harness checks.
- `module-ci unit (proteomics)`: proteomics fixture/oracle module checks.
- `module-ci unit (metabolomics)`: metabolomics fixture/oracle module checks.
- `module-ci unit (lipidomics)`: lipidomics fixture/oracle module checks.
- `module-ci browser smoke`: browser harness availability and smoke coverage.
- `current E2E smoke`: stable E2E manifest, fixture integrity, and browser
  harness smoke.

Nightly/release gates:

- `all-module-unit-contracts`: every unit-contract module scenario.
- `all-current-e2e-workflows`: every current `e2e-*` test file.
- `alternate-launch-browser-smoke`: launch/test-mode and browser harness paths.
- `report-export-fidelity`: report/export smoke tests and enrichment report path.

Artifact retention is deterministic and limited to CI output directories:

- `tests/testthat/_module_ci_artifacts/**`;
- `tests/testthat/_e2e_artifacts/**`;
- scorecard JSON/Markdown files inside each job artifact directory.

The workflow deliberately does not upload root-level config, ignored local
profiles, renv directories, workbooks, or user-specific files.

Local examples:

```bash
Rscript tools/ci/run-module-ci.R --runtime unit-contract --reporter summary
Rscript tools/ci/run-module-ci.R --omic proteomics --module fixtures --runtime unit-contract --reporter summary
Rscript tools/ci/module-ci-scorecard.R --artifact-dir tests/testthat/_module_ci_artifacts
```

## Proteomics Import Matrix

MCI-004 adds `tests/testthat/test-module-ci-prot-import.R` plus
`tests/testthat/helper-module-ci-prot-import.R`. This suite targets the
proteomics import branch point before design/QC state handoff.

The matrix covers:

- DIA-NN: canonical reports, no-FASTA continuation, duplicate protein IDs,
  zero/NA intensities, unsafe run names and sanitisation, missing required
  columns, extra metadata preservation, and peptide/protein grouping integrity.
- Proteome Discoverer TMT: `Abundance: F...` channel names, simplified
  `Abundance.*` names, duplicate reporter/sample collisions, missing accession
  columns, shared peptide accessions, and all-NA channels.
- MaxQuant LFQ: `LFQ.intensity.*` detection, contaminant/reverse filtering,
  protein-group accessions, missing gene names, duplicate accessions, missing
  intensity/protein columns, and zero-intensity preservation.
- FragPipe LFQ: `MaxLFQ Intensity` columns, regular intensity fallback,
  alternate sample names, missing intensity columns, missing protein IDs, and
  the same downstream LFQ contract as MaxQuant.
- Browser import controls: test mode exposes `search_results_standard` and
  `fasta_file_standard`; production mode keeps shinyFiles controls when
  `shinyFiles` is installed.
- Import state/artifacts: checkpoint capture, workflow type, report-template
  mapping, setup log, state digest, source/results directories, and design-step
  unlock are asserted after valid import.

The import matrix exposed two production hardening changes:

- PD-TMT duplicate abundance channels now fail before pivoting with
  `Duplicate TMT reporter channel/sample names after normalization`.
- MaxQuant imports read `proteinGroups.txt` as explicit tab-delimited files so
  malformed low-column fixtures still reach MaxQuant validation instead of
  failing during delimiter guessing.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic proteomics --module import --runtime unit-contract --reporter summary
```

## Proteomics Design Matrix

MCI-005 adds `tests/testthat/test-module-ci-prot-design.R` plus
`tests/testthat/helper-module-ci-prot-design.R`. This suite targets the sample
assignment and contrast contract between import and every downstream proteomics
module.

The matrix covers:

- Builder payloads for two-group, three-group, unbalanced, multi-factor,
  batch-aware, and interaction-ready designs.
- Removed-sample handling, sample-order independence, sanitized run names,
  suffix-based renames, duplicate rows, missing samples, extra rows, and
  case-drift rejection.
- Contrast append and serialization for valid, reversed, duplicate, self, empty,
  and invalid-term contrasts while preserving the existing `*_vs_*` friendly
  name format and `groupA-groupB` limma contrast format.
- Existing-design import round trips for `design_matrix.tab`, `data_cln.tab`,
  `contrast_strings.tab`, and `config.ini`, including TSV type drift in all-NA
  optional columns.
- Builder artifact persistence for `design_matrix.tab`, `data_cln.tab`,
  `contrast_strings.tab`, `manifest.json`, and `config.ini`.
- Workflow state digest checks that keep design pending until internally
  consistent design/data/contrast artifacts are present.
- Lightweight downstream DA setup validation that confirms sample identity,
  model-matrix contrast terms, required cleaned-data columns, and numeric
  intensity shape before QC or DA consumes the design artifacts.

The design matrix did not require production behavior changes. It codifies the
current proteomics design contract so future refactors preserve:

- every imported sample exactly once unless a rejection scenario is intentional;
- existing builder control IDs such as `prot-design-import-existing`;
- existing artifact filenames and text formats;
- `*_vs_*` friendly names generated from contrast names;
- shared design patterns used by metabolomics and lipidomics modules.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic proteomics --module design --runtime unit-contract --reporter summary
Rscript tools/ci/module-ci-scorecard.R --artifact-dir tests/testthat/_module_ci_artifacts
```
