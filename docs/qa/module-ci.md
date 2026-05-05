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

## Proteomics Peptide QC Matrix

MCI-006 adds `tests/testthat/test-module-ci-prot-peptide-qc.R` plus
`tests/testthat/helper-module-ci-prot-peptide-qc.R`. This suite targets the DIA
peptide-QC contract before peptide-to-protein rollup and downstream protein QC.

The matrix covers:

- Intensity filtering at zero, exact threshold boundary, above-max thresholds,
  missing values, all-zero peptides, and non-finite intensities.
- Q-value filtering for valid thresholds, inclusive `<= threshold` boundary
  behavior, missing q-values, missing output/filter columns, all-pass, and
  all-fail cases.
- Sample and replicate filters for balanced inputs, unbalanced sample counts,
  duplicate peptide/run observations, inclusion-list overrides, and missing
  design rows.
- Deterministic technical-replicate imputation, no-imputation settings,
  all-missing peptides, and minimum-observation guards.
- Plot and summary outputs for density, PCA, RLE, summary counts, and
  degenerate low-variance fixtures.
- State history for `qvalue_filtered`, `precursor_rollup`,
  `intensity_filtered`, `protein_peptide_filtered`, `sample_filtered`,
  `replicate_filtered`, and `imputed`, with a guard that protein QC cannot
  proceed before the peptide `imputed` state exists.

The peptide-QC matrix exposed production hardening changes:

- Peptide intensity/missingness filters now coerce all non-finite values
  (`Inf`, `-Inf`, and `NaN`) to missing before threshold tests.
- Peptide q-value filtering now treats threshold equality as passing, matching
  the UI wording and common FDR cutoff expectation.
- Peptide intensity helpers now resolve S4 slot-backed column parameters before
  tidy evaluation, so module apply steps do not pass literal expressions such as
  `theObject@sample_id` into `dplyr`.
- Peptide sample filtering now honors module-updated object args before legacy
  stack introspection, preventing tiny valid DIA fixtures from being emptied by
  the historical default sample cutoff.
- The downstream protein intensity transition now resolves
  `removeRowsWithMissingValuesPercent` parameters from explicit args or the
  correct S4 `@args` section when the generic is injected through module
  apply-step dependencies. Browser DIA and DIA+LIMPA lanes both advance through
  `protein_intensity_filtered`.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic proteomics --module qc_peptide --runtime unit-contract --reporter summary
```

## Proteomics Protein QC Matrix

MCI-007 adds `tests/testthat/test-module-ci-prot-protein-qc.R` plus
`tests/testthat/helper-module-ci-prot-protein-qc.R`. This suite targets the
protein-QC contract shared by protein-level TMT/LFQ imports and DIA peptide
rollup outputs.

The matrix covers:

- Protein-level TMT and LFQ bypass into `protein_s4_initial`, including invalid
  protein-ID column rejection before state advancement.
- DIA rollup through IQ MaxLFQ and limpa DPC-Quant test-mode fallback,
  including missing optional `limpa`, alias restoration, dropped-sample design
  pruning, report-template metadata, ambiguous/shared peptides, and
  peptide-count guards.
- Accession cleanup for semicolon-delimited protein groups, contaminants,
  canonical isoform choice, duplicate canonical IDs, missing gene names, FASTA
  metadata present, and FASTA metadata absent.
- Protein intensity filtering for missingness thresholds, percentile thresholds,
  all-zero proteins, all-non-finite proteins, and downstream sample/design
  alignment.
- Duplicate protein resolution for mean/median/max-compatible numeric sample
  columns, including sample IDs that do not contain digits.
- Protein replicate filtering for balanced, unbalanced, single-replicate, and
  disagreeing technical-replicate patterns.
- State/artifact checks for `protein_s4_initial`, `protein_s4_created`,
  `protein_accession_cleaned`, `protein_intensity_filtered`,
  `duplicates_removed`, and `protein_replicate_filtered`, plus protein counts,
  QC parameter persistence, exported replicate-filter TSVs, density, PCA, and
  RLE plots.

The protein-QC matrix exposed production hardening changes:

- Protein missingness filtering now coerces every non-finite intensity
  (`Inf`, `-Inf`, and `NaN`) to missing before threshold tests.
- Duplicate protein removal now aggregates all numeric sample columns instead
  of relying on digit-matching column names, preserving valid sample IDs such as
  `Alpha` and `Beta`.
- Protein-level density plotting now dispatches for `ProteinQuantitativeData`
  by reusing the existing protein PCA path and ggplot density renderer.
- Protein PCA now accepts `label_column = NULL`, matching module usage and the
  peptide-QC plotting contract.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic proteomics --module qc_protein --runtime unit-contract --reporter summary
```

## Proteomics Normalization/RUV Matrix

MCI-008 adds `tests/testthat/test-module-ci-prot-norm.R` plus
`tests/testthat/helper-module-ci-prot-norm.R`. This suite targets the
normalization, RUV, correlation-filtering, filtered-session export, and DA reload
contract after protein QC and before differential abundance.

The matrix covers:

- Between-sample normalization for `none`, `cyclicloess`, `quantile`, and
  `scale`, with explicit oracles for already-log values, zeros, negative values,
  `NA`, `Inf`, and `NaN` sentinel handling.
- RUV skip, manual, and automatic branches, including named negative-control
  vectors, too-few-control rejection, invalid imputation-k rejection, optimizer
  weak-effect and plateau boundaries, audit persistence, and R6 state updates.
- RUV-III correction orientation and sample/feature alignment, verifying that
  corrected matrices preserve the expected protein/sample shape and finite
  numeric payload when inputs are finite.
- Correlation filtering for pass-all, fail-one, fail-many, exact-threshold,
  one-sample/no-pair, and missing-sample-correlation cases.
- QC artifact generation for pre-normalization, post-normalization, and
  RUV-corrected PCA, RLE, density, and correlation plots, plus composite layout
  selection and RUV canonical-correlation plot persistence.
- Degenerate constant-feature plot fallbacks for protein PCA, RLE, and density
  so plot generation remains push-safe even when variance is zero.
- Filtered-session export fidelity for R6 current state, complete state map,
  state history, workflow type, report template, normalization method, RUV mode,
  RUV k, correlation threshold, limpa metadata, protein/sample counts, and
  stable filenames such as `filtered_session_data_latest.rds`.
- DA reload smoke for exported skip, manual, and automatic+limpa variants
  without running full DA, confirming that the DA module can restore the S4
  object, R6 state manager, design, contrasts, config, UniProt annotations,
  tab status, and RUV optimization payload.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic proteomics --module normalization --runtime unit-contract --reporter summary
```

## Proteomics Differential Abundance Matrix

MCI-009 adds `tests/testthat/test-module-ci-prot-da.R` plus
`tests/testthat/helper-module-ci-prot-da.R`. This suite targets the DA contract
between filtered-session reload, model fitting, visual inspection, export, and
downstream enrichment/report consumers.

The matrix covers:

- Filtered-session reload for DIA, DIA+limpa, TMT, LFQ MaxQuant, and LFQ
  FragPipe payloads, including workflow type, report template, limpa flag,
  R6 state, design, contrasts, UniProt annotations, RUV audit metadata, and tab
  status restoration.
- Invalid session rejection for missing session files, stale source
  directories, and malformed RDS payloads before DA state or globals are
  mutated.
- Formula and contrast validation for two-group, multi-group, batch-adjusted,
  reversed, raw-only, duplicate, empty, invalid-term, and invalid-formula
  definitions.
- A production validation seam that rejects invalid formulas/contrasts before
  `differentialAbundanceAnalysis()` or result export can run.
- Statistical output schemas for standard limma results, no significant
  proteins, all significant proteins, tied p-values, missing/non-finite
  p-values, and small-n q-value fallback.
- Renderer behavior for result-table significance filters, q-value/logFC
  threshold changes, empty result states, heatmap top-N and clustering options,
  static volcano, and interactive volcano handoff.
- Export fidelity for raw contrast filenames, friendly comparison labels,
  TSV/XLSX long-annot schemas, sample intensity columns, gene annotations,
  UI/report-facing DA parameters, and stable `da_proteins_*_long_annot.*`
  filenames.
- Downstream smoke checks proving accepted DA schemas can be loaded by
  `createDAResultsForEnrichment()`, selected through enrichment-friendly
  contrast labels, and discovered by summary/report DA workbook indexing.

The DA matrix exposed production hardening changes:

- DA filtered-session loading now validates that the payload contains a
  `ProteinQuantitativeData` object, design matrix, contrasts table, and valid
  protein/sample counts before mutating DA or workflow state.
- DA run handling now normalizes contrast tables to raw, full-format, and
  friendly labels, rejects duplicate/empty/malformed contrasts, and validates
  formula terms against the model matrix before statistical execution.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic proteomics --module differential_abundance --runtime unit-contract --reporter summary
```

## Proteomics Enrichment Matrix

MCI-010 adds `tests/testthat/test-module-ci-prot-enrich.R` plus
`tests/testthat/helper-module-ci-prot-enrich.R`. This suite targets the
proteomics enrichment boundary after DA, where backend selection, identifiers,
thresholds, deterministic doubles, browser state, exports, publication workbooks,
and report parameters must remain synchronized.

The matrix covers:

- Backend routing for supported model-organism taxa through `gprofiler2`,
  unsupported numeric taxa through `clusterProfiler`, explicit taxon changes in
  the reactive method state, and invalid taxon rejection before process
  dispatch.
- Identifier handling for UniProt accessions, gene symbols, Ensembl-like IDs,
  duplicate IDs, missing annotations, mixed-case IDs, and no-mappable-ID
  inputs.
- Threshold and display behavior for q-value boundaries, p-value columns,
  source filters, direction filters, top-N table limits, and valid empty-result
  states.
- Deterministic test-mode doubles for gprofiler-like list payloads and
  clusterProfiler-like S4 payloads, including backend-specific pathway TSV
  files and collated all-contrast result state.
- Export fidelity for backend-specific ZIP contents, summary text, selected
  contrast, organism/taxon, thresholds, backend TSV columns, publication
  workbook indexing, and report-facing payload serialization.
- Pre-seeded module smoke that initializes enrichment from existing DA output so
  browser/module checks do not need to rerun full DA for every backend case.

The enrichment matrix exposed production hardening changes:

- Taxon IDs are normalized and validated as positive integer NCBI taxonomy IDs
  before backend routing or process argument construction.
- Enrichment process parameters now reject non-finite cutoffs, out-of-range
  q-value thresholds, and empty correction methods before deterministic or real
  backend execution.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic proteomics --module enrichment --runtime unit-contract --reporter summary
```
