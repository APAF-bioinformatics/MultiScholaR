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

## Proteomics Summary Report Matrix

MCI-011 adds `tests/testthat/test-module-ci-prot-summary.R` plus
`tests/testthat/helper-module-ci-prot-summary.R`. This suite targets the
release-facing summary/report boundary where workflow parameters, publication
artifacts, template selection, rendered reports, and session-state exports must
stay branch-correct.

The matrix covers:

- Workflow-parameter round trips for DIA, DIA+limpa, TMT, MaxQuant LFQ,
  FragPipe LFQ, enrichment-run/no-enrichment branches, RUV skip/run states, and
  altered DA thresholds.
- Publication-copy behavior for DA tables, enrichment pathway TSVs, normalized
  intensity files, QC plots, copied publication figure folders, and explicitly
  missing optional artifacts.
- Report-template selection for `DIANN_report.rmd`,
  `DIANN_limpa_report.rmd`, `TMT_report.rmd`, and `LFQ_report.rmd`, including
  stale-template rejection across branch changes.
- Rendered-report smoke using deterministic render stubs in push-safe module CI.
- Session-state export schema checks for omic type, workflow type, report
  template, DA UI parameters, enrichment UI parameters, and RUV metadata.
- Artifact scorecard coverage for study parameters, publication workbooks,
  copied figures, report file, and session-state export.

The summary matrix exposed production hardening changes:

- Template status now reports DIA, DIA+limpa, TMT, and LFQ templates instead of
  only DIA/TMT.
- Summary report selection now ignores stale explicit templates that do not
  match the current workflow branch.
- Session-state export can serialize workflow/template metadata and parameter
  payloads when workflow data is available.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic proteomics --module summary_report --runtime unit-contract --reporter summary
```

## Metabolomics Import Matrix

MCI-012 adds `tests/testthat/test-module-ci-metab-import.R` plus
`tests/testthat/helper-module-ci-metab-import.R`. This suite targets the
metabolomics import boundary before design, QC, normalization, and DA can
consume assay state.

The matrix covers:

- Custom imports with explicit `Feature.Name` identifiers, annotation-present
  and annotation-absent payloads, alternate sample regex families, extra
  metadata columns, and character-to-numeric sample coercion boundaries.
- MS-DIAL-style imports with expected metadata columns, positive/negative mode
  annotations, duplicate feature IDs, zero intensities, missing intensities, and
  missing identifier rejection.
- Multi-assay imports for `LCMS_Pos+LCMS_Neg`, `GCMS` only, `LCMS_Pos+GCMS`,
  duplicate assay-name rejection, empty secondary assay behavior, and unusual
  assay-name rejection before file artifacts are written.
- Invalid schemas for missing ID columns, no sample columns, non-numeric sample
  values, mismatched sample patterns, and malformed payloads that must fail
  before workflow state advances.
- Browser import controls for standard upload paths in test mode and
  shinyFiles-compatible production control rendering.
- Import artifact fidelity for `data_cln_<assay>.tab`,
  `assay_manifest.txt`, `column_mapping.json`, source upload copies, import
  summary rows, and state-digest import payloads.

The metabolomics import matrix exposed production hardening changes:

- Built workflow payloads are now validated before state mutation: assay names
  must be non-empty, unique, path-safe, and paired with explicit metabolite and
  sample-column mappings.
- Sample intensity columns are accepted as numeric or losslessly
  numeric-coercible strings, then coerced once validation has passed.
- When an experiment source directory is available, import processing writes
  deterministic cleaned-assay, manifest, mapping, summary, and copied-source
  artifacts for downstream audit and report fidelity checks.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic metabolomics --module import --runtime unit-contract --reporter summary
```

## Metabolomics Design Matrix

MCI-013 adds `tests/testthat/test-module-ci-metab-design.R` plus
`tests/testthat/helper-module-ci-metab-design.R`. This suite targets the shared
metabolomics design contract that sits between import and every downstream
multi-assay step.

The matrix covers:

- Builder-created designs for two-group, three-group, unbalanced, batch-aware,
  missing-group, and extra-metadata cases.
- Multi-assay alignment for `LCMS_Pos+LCMS_Neg`, `GCMS`, `LCMS_Pos+GCMS`,
  reordered assay columns, missing assay samples, extra sample-like columns, and
  case-varied sample names.
- Contrast serialization for raw limma terms, friendly names, full
  `friendly=raw` formats, reversed contrasts, duplicate contrasts, invalid
  model terms, and no-contrast preflight behavior.
- Existing-design import from source directories containing
  `design_matrix.tab`, `data_cln_<assay>.tab`, `assay_manifest.txt`,
  `column_mapping.json`, `contrast_strings.tab`, and optional `config.ini`.
- Browser-safe import modal controls for test-mode path input and import
  confirmation.
- Direct DA preflight that proves accepted designs can produce a model matrix
  and that contrast terms are present before analysis dispatch.

The design matrix exposed production hardening changes:

- A shared metabolomics design validation layer now checks one design row per
  expected sample, unique sample rows, assigned groups, replicate assignments,
  unique assay names, mapped sample columns, and sample-like assay columns not
  present in the mapping.
- Contrast tables are normalized to preserve raw contrast expressions,
  friendly labels, and full serialization strings while rejecting empty,
  duplicate, or model-incompatible terms.
- Builder save now runs a non-destructive preflight before mutating workflow
  state or writing artifacts. This keeps corruptive designs from advancing into
  QC/DA while preserving the existing ability to save a design before contrasts
  are finalized.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic metabolomics --module design --runtime unit-contract --reporter summary
```

## Metabolomics QC Matrix

MCI-014 adds `tests/testthat/test-module-ci-metab-qc.R` plus
`tests/testthat/helper-module-ci-metab-qc.R`. This suite targets the
metabolomics QC boundary that mutates multi-assay S4 state before
normalization.

The matrix covers:

- Intensity filtering for threshold boundaries, zeros, missing values,
  all-pass/all-fail branches, and per-assay feature-count preservation.
- Duplicate handling for repeated feature IDs, repeated annotations,
  identical-intensity ties, conflicting intensities, and highest-average
  retention across LC/GC assay lists.
- ITSD detection for absent standards, explicit regex matches, fallback regex
  detection, multi-assay standards, and malformed metadata/search-column
  failures.
- S4 plotting for PCA, Pearson correlations, RLE, density, small-n CV
  fallbacks, constant-feature Pearson fallbacks, and missing-group rejection.
- Finalization state for valid `MetaboliteAssayData`, invalid current state,
  missing design rows, `metab_qc_complete` persistence, QC status completion,
  and normalization handoff history.
- Static browser smoke for the QC orchestrator, intensity, duplicates, ITSD,
  and finalization tabs so control IDs remain wired for Shiny/browser tests.

The QC matrix exposed one production hardening change:

- QC finalization now rejects empty assay payloads and empty design matrices
  before saving `metab_qc_complete`, preventing normalization from receiving a
  syntactically valid but unusable metabolomics S4 object.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic metabolomics --module qc --runtime unit-contract --reporter summary
```

## Metabolomics Normalization/RUV Matrix

MCI-015 adds `tests/testthat/test-module-ci-metab-norm.R` plus
`tests/testthat/helper-module-ci-metab-norm.R`. This suite targets the
metabolomics normalization boundary after QC and before differential abundance,
where multi-assay LC/GC state, normalization parameters, RUV correction,
correlation filtering, and exported sessions must stay coherent.

The matrix covers:

- Log transformation with zeros, missing values, negative values, positive
  offsets, and rejection of invalid offsets while preserving assay names and
  sample/design alignment.
- Between-sample normalization for `none`, `cyclicloess`, `quantile`, and
  `scale`, including finite numeric sanity and explicit persistence of the
  selected method.
- ITSD normalization for regex-based global selection, manual assay-specific
  selections, absent standards, invalid manual standards, invalid aggregation
  methods, ITSD feature-count serialization, and state persistence to
  `metab_itsd_norm`.
- RUV-III control flow for skipped correction, manual optimization, automatic
  optimization, named control vectors, assay-specific best-k payloads, failed
  too-few-control assays, invalid-k rejection, and guarded
  `metab_ruv_corrected` state advancement.
- Correlation filtering for skipped filtering, pass-all, exact threshold
  boundary, fail-one replicate pair, fail-many replicate pairs, small-n/no-pair
  cases, invalid correlation payloads, and explicit handoff state.
- Session export fidelity for `omic_type`, `current_s4_object`, `assay_names`,
  `feature_counts`, normalization/RUV/correlation completion flags, RUV
  optimization results, ITSD selections, QC params, summary text, latest RDS,
  and metadata sidecars.
- DA reload smoke for every accepted export variant: RUV skipped, manual RUV,
  and automatic RUV. Reload verifies S4 restoration, state-manager persistence,
  contrast dropdown payloads, assay dropdown payloads, and formula restoration
  from S4 args.

The normalization matrix is intentionally push-safe and unit-contract scoped.
It avoids browser launch while still exercising the production helper seams used
by the Shiny module:

- `logTransformAssays()`
- `normaliseUntransformedData()`
- `normaliseBetweenSamples()`
- `runMetabNormItsdProgressApplyShell()`
- `runMetabNormBetweenSampleProgressApplyShell()`
- `runMetabNormRuvProgressApplyShell()`
- `runMetabNormRuvOptimizationStep()`
- `runMetabNormRuvCorrectionStep()`
- `dispatchMetabNormApplyCorrelation()`
- `runMetabNormSkipCorrelationObserverEntry()`
- `buildMetabNormExportSessionData()`
- `saveMetabNormExportSessionRdsFiles()`
- `saveMetabNormExportMetadataFiles()`
- `saveMetabNormExportSummaryFile()`
- `runMetabDaLoadSessionObserverShell()`
- `restoreMetabDaLoadedSessionState()`

The suite exposed no required production behavior changes. It codifies the
current metabolomics normalization contract so future refactors preserve:

- LC/GC assay names and sample columns through every normalization branch;
- S4 args and session export fields for ITSD, RUV, normalization, and DA
  formulas;
- non-silent failure behavior for invalid RUV and invalid correlation inputs;
- the `metab_filtered_session_data_latest.rds` payload consumed by DA and
  reports;
- summary and metadata sidecar filenames used by downstream audit/report flows.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic metabolomics --module normalization --runtime unit-contract --reporter summary
```

## Metabolomics DA Matrix

MCI-016 adds `tests/testthat/test-module-ci-metab-da.R` plus
`tests/testthat/helper-module-ci-metab-da.R`. This suite targets the
metabolomics differential-abundance boundary after normalization/correlation
filtering and before summary/report export. It is designed to catch wiring
regressions in session reload, model validation, per-assay analysis,
combined/per-assay rendering, export schemas, and report-facing artifacts.

The matrix covers:

- Filtered-session reload for LC-only, GC-only, and combined LC+GC payloads,
  including state-manager restoration, formula restoration from S4 args,
  contrast dropdowns, assay dropdowns, stale `source_dir` fallback to
  `export_dir`, missing session files, malformed RDS files, and mismatched
  `assay_names` versus S4 assay slots.
- Formula and contrast preflight for two-group, multi-group, batch-aware,
  reversed, duplicate, no-contrast, invalid-term, and invalid-formula cases.
  Invalid DA inputs now fail before the analysis shell can run.
- DA orchestration for significant, non-significant, all-significant,
  no-variance, tied-p-value, and missing-feature-metadata cases while
  preserving feature IDs, assay provenance, raw contrasts, friendly contrasts,
  and per-sample intensity columns.
- Render behavior for combined and per-assay static volcano paths, heatmap
  parameter plumbing, results-table filters, threshold changes, empty outputs,
  max-row truncation, and selector choices for contrast/assay controls.
- Export behavior for TSV and XLSX outputs, stable `de_<mode>_metabolites_*`
  file prefixes, `posmode` and sanitized assay fallbacks such as `gcms`, NumSig
  summary files, contrast columns, friendly names, feature IDs, assay columns,
  and intensity columns.
- Summary/report consumer smoke that reads exported TSV schemas and exercises
  the `copyToResultsSummary()` workbook path used to build
  `DA_results_metabolomics.xlsx`.

The DA matrix exposed three production hardening changes:

- DA analysis input resolution now calls metabolomics design preflight when a
  formula is supplied, so invalid model terms and invalid contrasts are rejected
  before analysis starts.
- Reloaded metabolomics DA sessions now reject non-empty `assay_names` that do
  not match the S4 object's assay slots, preventing stale dropdown state from
  silently diverging from the object being analyzed.
- Metabolomics DA export tables now retain the `assay` column while preserving
  existing filename prefixes and first-six statistical columns.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic metabolomics --module differential_abundance --runtime unit-contract --reporter summary
```

## Metabolomics Summary/Report Matrix

MCI-017 adds `tests/testthat/test-module-ci-metab-summary.R` plus
`tests/testthat/helper-module-ci-metab-summary.R`. This suite targets the
release-facing metabolomics summary/report boundary where workflow parameters,
publication artifacts, deterministic report rendering, and session-state
exports must preserve LC, GC, and mixed LC+GC provenance.

The matrix covers:

- Parameter roundtrip for LC-only, GC-only, and combined LC+GC workflows,
  including normalization method, ITSD settings, RUV mode/result payloads, DA
  q-value/logFC thresholds, formula strings, report template, and assay layout.
- Publication-copy behavior for DA workbooks, enrichment workbooks, QC figures,
  DA figure directories, NumSig tables, normalized/RUV TSV and RDS outputs, and
  mixed-assay DA artifacts.
- Required versus optional artifact handling: missing design, contrast, or
  study-parameter artifacts now fail explicitly, while missing optional plots
  are reported as warnings without blocking copy completion.
- Deterministic report rendering using push-safe render stubs for
  `metabolomics_report.rmd`, plus missing-output behavior.
- Study/session parameter payload assertions for assay names, feature counts,
  sample counts, friendly and raw contrasts, selected normalization/RUV/ITSD/DA
  parameters, workflow type, and report template.
- Reloadable summary session-state exports and scorecard entries for study
  parameters, DA/enrichment publication workbooks, normalized outputs, copied
  figures, rendered report, and session RDS.

The summary/report matrix exposed three production hardening changes:

- Metabolomics session-state export now includes workflow/report metadata and a
  parameter payload when workflow data is available, while preserving the prior
  minimal export shape when workflow data is absent.
- Publication copy now classifies required copy failures before setting the
  summary module to copied, so missing design/contrast/study-parameter artifacts
  cannot be silently treated as success.
- Metabolomics normalized TSV outputs are copied to `Publication_tables`
  alongside the existing normalized/RUV RDS outputs when present.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic metabolomics --module summary_report --runtime unit-contract --reporter summary
```

## Lipidomics Import Matrix

MCI-018 adds `tests/testthat/test-module-ci-lipid-import.R` plus
`tests/testthat/helper-module-ci-lipid-import.R`. This suite targets the
lipidomics import module as a format-sensitive, multi-assay contract instead of
only checking the canonical E2E upload lane.

Coverage includes:

- LipidSearch import for LCMS_Pos, LCMS_Neg, GCMS-named, missing `LipidName`,
  missing `LipidClass`, duplicate lipid IDs, extra metadata, zero intensities,
  and non-numeric sample payloads.
- MS-DIAL import for expected columns, missing annotation, ion-mode metadata,
  duplicate features, zero/missing intensities, and exclusion of numeric
  annotation columns such as total scores from sample-column selection.
- Custom/vendor override behavior for explicit mappings, case-insensitive
  custom lookup, unsupported detected formats falling back through the default
  reader, and mismatch between selected and detected formats.
- Multi-assay routing for absent assay2, empty assay2, duplicate assay names,
  unusual assay names, and assay2 reader parity with the active format.
- Invalid input short-circuiting for missing ID columns, missing sample columns,
  malformed/insufficient files, no numeric sample data, and mismatched sample
  naming before `workflow_data` or downstream tab state advances.
- Source artifact and state digest fidelity for `data_cln_*.tab`,
  `assay_manifest.txt`, `column_mapping.json`,
  `lipidomics_import_summary.tsv`, import source copies, workflow type, and
  design unlock.

The import matrix exposed three production hardening changes:

- Lipid import processing now validates assay names, ID columns, sample-column
  presence, and numeric sample payloads before committing `workflow_data`.
- Successful imports now write source artifacts directly from the import module,
  matching the downstream design import contract and recording the artifact
  digest in `processing_log$setup_import`.
- The Shiny process observer now passes `experiment_paths` and source file paths
  through to the processing shell, so browser and unit-contract paths exercise
  the same artifact-writing route.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic lipidomics --module import --runtime unit-contract --reporter summary
```

## Lipidomics Design Matrix

MCI-019 adds `tests/testthat/test-module-ci-lipid-design.R` plus
`tests/testthat/helper-module-ci-lipid-design.R`. This suite targets the
lipidomics design boundary between import and QC/DA for every loaded lipid
assay.

Coverage includes:

- Builder payloads for balanced two-group, unbalanced, three-group, batch-aware,
  metadata-rich, and missing-group designs.
- Multi-assay alignment across LCMS_Pos, LCMS_Neg, GCMS, and combined assay
  layouts, including reordered, missing, extra, and case-varied sample columns.
- Contrast serialization for raw, friendly, reversed, duplicate, invalid, empty,
  and design-save-without-final-contrast cases while preserving
  `*_vs_*` friendly names and `groupA-groupB` contrast expressions.
- Existing-design import from source directories with `design_matrix.tab`,
  `assay_manifest.txt`, `data_cln_*.tab`, `contrast_strings.tab`,
  `column_mapping.json`, and browser-safe import controls.
- Artifact fidelity for design, contrast, assay manifest, column mapping,
  manifest, config, and per-assay cleaned data files.
- DA preflight checks that reject invalid formulas, invalid contrast terms, and
  sample/design drift before QC triggers, S4 state persistence, or tab status
  advancement.

The design matrix exposed two production hardening changes:

- Lipid design builder saves now validate design/sample alignment and DA contrast
  feasibility before mutating workflow data, writing artifacts, creating
  `lipid_raw_data_s4`, or advancing QC.
- Existing lipid design imports now validate the imported design, assay manifest,
  column mapping, contrast table, and formula locally before committing workflow
  state or assigning globals.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic lipidomics --module design --runtime unit-contract --reporter summary
```

## Lipidomics QC Matrix

MCI-020 adds `tests/testthat/test-module-ci-lipid-qc.R` plus
`tests/testthat/helper-module-ci-lipid-qc.R`. This suite targets the lipidomics
QC boundary after design and before normalization, where intensity filtering,
duplicate resolution, ITSD evaluation, plots, and final QC state must stay
multi-assay safe.

Coverage includes:

- Intensity filtering for percentile/proportion thresholds, zeros, missing
  values, non-finite values, all-pass/all-fail branches, and per-assay feature
  count preservation.
- Duplicate handling for repeated lipid IDs, repeated lipid annotations/classes,
  tied intensities, conflicting intensities, and state-save statistics from the
  duplicate-resolution module helper.
- ITSD metrics for global regexes, assay-specific matching via annotation
  columns, absent standards, fallback pattern matching, multiple matches, and
  invalid regex input.
- S4 plotting smoke for PCA, Pearson, RLE, density, small-n CV fallbacks,
  constant-feature Pearson fallbacks, and missing-group rejection.
- Finalization state for valid `LipidomicsAssayData`, invalid current state,
  empty design, empty assay payloads, `lipid_qc_complete` persistence, QC status
  completion, and normalization unlock.
- Static browser smoke for the lipid QC orchestrator, intensity, duplicates,
  ITSD, and finalization tabs so control IDs remain wired for Shiny/browser
  tests.

The QC matrix exposed two production hardening changes:

- Lipid intensity filtering now rejects non-finite sample intensities before
  threshold calculations, while still allowing ordinary missing values.
- Lipid QC finalization now validates assay payloads, design rows, and
  assay/design alignment before saving `lipid_qc_complete`; successful
  finalization also unlocks the normalization tab when it was disabled/locked.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic lipidomics --module qc --runtime unit-contract --reporter summary
```

## Lipidomics Normalization/RUV Matrix

MCI-021 adds `tests/testthat/test-module-ci-lipid-norm.R` plus
`tests/testthat/helper-module-ci-lipid-norm.R`. This suite targets the
lipidomics normalization boundary after QC and before DA, where ITSD
normalization, log transforms, between-sample normalization, RUV, correlation
filtering, export, and DA reload all have to preserve multi-assay provenance.

Coverage includes:

- Normalization methods `none`, `cyclicloess`, `quantile`, and `scale`, plus
  log offsets, zeros, negative values, missing values, constant features,
  per-assay shifts, and sample/design alignment.
- ITSD normalization for regex-selected standards, manual per-assay standards,
  absent standards, invalid manual IDs, invalid aggregation, ITSD removal, and
  serialized per-assay ITSD feature/count metadata.
- RUV skip, manual, and automatic branches, including named control vectors,
  per-assay `k` payloads, optimizer traces, failed assay optimization,
  invalid-k correction failure, and no `lipid_ruv_corrected` state save on RUV
  failure.
- Correlation filtering for skipped filtering, pass-all, exact threshold
  boundary, fail-one, fail-many, small-n/empty correlation payloads, invalid
  payloads, tab-status completion, and explicit filtered-state handoff.
- Session export for skip/manual/automatic variants, validating `omic_type`,
  current S4 class, assay names, feature counts, normalization flags, RUV
  payloads, ITSD payloads, metadata sidecars, summary text, and
  `lipid_filtered_session_data_latest.rds` compatibility.
- Lipid DA reload smoke for every accepted normalized session variant, including
  state-manager restoration, formula restoration, contrast dropdowns, assay
  dropdowns, and table assay choices.

The normalization matrix exposed one production hardening change:

- Lipid normalized-session export now counts samples by intersecting assay
  columns with the design matrix sample IDs, instead of subtracting two metadata
  columns from the assay width. This preserves correct sample counts when lipid
  assays include additional metadata such as `LipidClass`.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic lipidomics --module normalization --runtime unit-contract --reporter summary
```

## Lipidomics DA Matrix

MCI-022 adds `tests/testthat/test-module-ci-lipid-da.R` plus
`tests/testthat/helper-module-ci-lipid-da.R`. This suite targets the lipidomics
DA boundary after normalized-session export and before summary/report
consumption, where filtered sessions, design formulas, contrasts, multi-assay
results, visualizations, and publication artifacts must remain wired together.

Coverage includes:

- Filtered-session reload for single-assay, dual-assay, and LCMS-plus-GCMS
  layouts, including stale source fallback, missing file notifications,
  malformed RDS rejection, state-manager restoration, formula restoration,
  contrast dropdowns, assay dropdowns, and explicit assay-name/S4 mismatch
  rejection for lipidomics sessions.
- Formula and contrast preflight for two-group, multi-group, batch-aware,
  reversed, duplicate-friendly-name, empty, invalid-term, and invalid-formula
  cases. Invalid inputs are rejected before the DA engine mutates state or marks
  analysis complete.
- DA orchestration across LCMS_Pos, LCMS_Neg, and GCMS with deterministic edge
  cases for mixed significance, no significant features, all significant
  features, no variance, tied p-values, and missing annotation metadata.
- Long-table schema fidelity for lipid IDs, lipid names, assay provenance,
  raw contrasts, friendly contrasts, numerator/denominator labels, q-values,
  significance calls, and per-sample intensity columns.
- Render contracts for combined and per-assay volcano inputs, heatmap state and
  cluster storage, summary text, significance filters, empty result tables, max
  row handling, and all DA selector updates.
- Export/report fidelity for per-assay TSV/XLSX files, publication volcano
  outputs, heatmap no-significant fallbacks, NumSigDE tables, mode-specific file
  prefixes, and report-consumer schema smoke.

The DA matrix exposed production hardening changes:

- Lipid DA now runs a reusable preflight before analysis, validating
  assay/design alignment, formula terms, and non-empty contrast tables before
  invoking the DA engine.
- Lipid filtered-session reload now rejects lipidomics sessions whose serialized
  `assay_names` disagree with the loaded S4 object's `lipid_data` assays.
- Lipid DA results now fall back to the lipid ID as `lipid_name` when an assay
  has no annotation column, preserving export/report schemas for metadata-poor
  imports.
- Lipid DA table rendering now returns an empty state after filters remove all
  rows, instead of rendering a blank DataTable.
- Lipid DA success finalization now preserves selector state and disk-write
  status when output directories are intentionally absent, instead of returning
  early from the disk-write branch.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic lipidomics --module differential_abundance --runtime unit-contract --reporter summary
```

## Lipidomics Summary/Report Matrix

MCI-023 adds `tests/testthat/test-module-ci-lipid-summary.R` plus
`tests/testthat/helper-module-ci-lipid-summary.R`. This suite targets the final
lipidomics user-visible fidelity surface: summary parameter capture,
publication copying, report rendering, session export, and artifact scorecards
after DA has generated assay-specific outputs.

Coverage includes:

- Parameter roundtrip for single LCMS_Pos, LCMS_Pos+LCMS_Neg, and
  LCMS_Pos+GCMS layouts, with `none`, `scale`, and `quantile` normalization,
  RUV skip/manual/automatic branches, ITSD skip/manual/regex branches, q-value
  thresholds, raw contrast expressions, and friendly contrast names.
- Study-parameter fidelity for assay names, lipid counts per assay, sample
  count, lipid metadata columns, ITSD settings, log-transform settings,
  between-sample normalization, RUV-III status, per-assay RUV `k`/control
  counts, DA thresholds, and saved integration S4 `@args`.
- Publication-copy fidelity for combined DA workbooks, per-assay DA inputs,
  RUV/normalized TSV/RDS outputs, composite QC figures, stage-specific QC
  plots, ITSD figures, legacy volcano/heatmap directories, lipid heatmaps,
  lipid NumSigDE tables, Study_report design/contrast/study-parameter copies,
  and copied normalized table schemas.
- Missing-artifact behavior for optional report assets versus required study
  artifacts. Optional plot/directory gaps now surface as warning-status copies,
  while missing contrasts, design matrix, or study parameters fail explicitly
  before the summary state marks files as copied.
- Report rendering with deterministic push-safe stubs, preservation of
  `lipidomics_report.rmd`, download-handler filename fidelity, successful
  non-empty HTML output, and missing-render-output failure classification.
- Session-state export for reloadable RDS payloads containing omic type,
  workflow type, report template, assay names, feature counts, sample count,
  contrast count, parameter payloads, UI parameters, RUV optimization payloads,
  and saved S4 args.
- Scorecard coverage for study parameters, DA workbook, RUV-normalized TSV/RDS,
  QC figures, NumSigDE table, session state, and rendered report.

The summary/report matrix exposed production hardening changes:

- Lipid summary export now enriches session-state RDS files with workflow
  metadata, assay names, feature counts, sample counts, contrast counts,
  parameter payloads, RUV optimization results, and S4 args instead of exporting
  only basic UI state.
- Lipid publication copying now classifies required and optional copy failures:
  required study artifacts stop the workflow explicitly, while optional plot or
  directory gaps remain visible as warnings without blocking report generation.
- Lipidomics-specific publication copying now preserves stage QC plots, ITSD
  figures, heatmap folders, NumSigDE tables, and normalized/RUV assay outputs
  under the existing results-summary directory conventions.
- Push-safe observer shells now exercise save, publication copy, report render,
  and session export behavior without a live browser session, making the
  summary/report workflow runnable from module CI on every push.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic lipidomics --module summary_report --runtime unit-contract --reporter summary
```

## Cross-Module Integrity Sentinels

MCI-024 adds `tests/testthat/test-module-ci-sentinels.R` plus
`tests/testthat/helper-module-ci-sentinels.R`. This suite provides reusable
sentinels for every omic/module matrix so per-module-valid outputs cannot drift
silently before the next module consumes them.

Coverage includes:

- Sample identity digests that preserve order by default and fail unexplained
  drops, additions, reorders, or renames. Legitimate transformations must be
  declared through `expected_dropped` or named `expected_renamed` policies.
- Feature identity sentinels for proteomics, metabolomics, and lipidomics
  fixture packs, with reusable extraction from tables and S4 objects. Expected
  duplicate-resolution drops or ID renames can be declared explicitly.
- Assay provenance sentinels for metabolomics and lipidomics S4 objects,
  session exports, DA tables, and report parameter payloads, ensuring assay
  names remain aligned through multi-assay workflows.
- Parameter fidelity sentinels that compare live workflow state with serialized
  payloads for workflow type, report template, design formula, normalization,
  RUV, ITSD, DA thresholds, and enrichment backend choices.
- Shared-state collision sentinels for one-session multi-omic runs, including
  `project_dirs`, `config_list`, selected omic state, and path isolation across
  proteomics, metabolomics, and lipidomics.
- Artifact schema sentinels for TSV, XLSX, RDS, and HTML/report files. These
  validate non-empty outputs, required columns, row counts, RDS classes, and
  report text/template fingerprints.

The sentinel layer is additive: it lives in test helpers and does not require
production modules to expose new user-facing APIs. Future MCI tickets can adopt
the helpers at each module boundary by creating a digest at module entry,
asserting it at module exit, and declaring expected drops/renames in oracle
sidecars when a transformation is intentional.

Runner example:

```bash
Rscript tools/ci/run-module-ci.R --omic all --module cross_module_integrity --runtime unit-contract --reporter summary
```

## CI Enforcement And Release Gates

MCI-025 upgrades `.github/workflows/module-ci.yml` from an entrypoint smoke
workflow into the enforced module-CI/release-gate corpus. The workflow remains
additive until branch protection is configured, but the job names are stable so
they can be made required checks when the branch policy is enabled.

Push and pull-request checks:

- `module-ci unit (foundation)` runs manifest/harness foundation checks.
- `module-ci unit (cross-module-integrity)` runs the MCI-024 corruption
  sentinels and fixture-integrity bridge.
- `module-ci unit (proteomics)` runs all proteomics unit-contract module
  scenarios: import, design, peptide QC, protein QC, normalization/RUV, DA,
  enrichment, and summary/report.
- `module-ci unit (metabolomics)` runs all metabolomics unit-contract module
  scenarios: import, design, QC, normalization/RUV/ITSD, DA, and
  summary/report.
- `module-ci unit (lipidomics)` runs all lipidomics unit-contract module
  scenarios: import, design, QC, normalization/RUV/ITSD, DA, and
  summary/report.
- `module-ci browser smoke (...)` runs the shared browser harness plus one
  representative module surface for each current omic.
- `representative E2E (...)` runs one full GUI workflow per omic:
  proteomics DIA, metabolomics LC/GC, and lipidomics canonical LCMS.
- `E2E fixture and harness smoke` keeps the manifest, fixture, and browser
  harness contracts hot on every push/PR.

Nightly and release checks:

- `full gate (all-module-unit-contracts)` runs every unit-contract module-CI
  scenario in the manifest.
- `full gate (all-current-e2e-workflows)` runs every current `e2e-*` workflow.
- `full gate (alternate-launch-browser-smoke)` covers `run_app`, app-dir, and
  test-mode launch surfaces.
- `full gate (report-export-fidelity)` checks report/export smoke paths and
  enrichment report fidelity.
- `full gate (cross-omic-coexistence)` runs cross-omic E2E/shared-state checks
  and the MCI-024 sentinels together.
- `release candidate promotion gate` runs for `release/**` branches, `v*` tags,
  and manual `gate=release|all` dispatches. It fails unless the full gate matrix
  succeeds and writes a promotion summary artifact.

Artifacts are deliberately limited to CI output directories:

- module-CI artifacts under `tests/testthat/_module_ci_artifacts/**`;
- E2E/browser artifacts under `tests/testthat/_e2e_artifacts/**`;
- `module-ci-run-manifest.json`;
- `module-ci-scorecard.json`;
- `module-ci-scorecard.md`;
- `release-gate-summary.md`.

The workflow does not upload root-level ignored or local-development paths such
as `renv/`, `renv.lock`, `.Rprofile`, `.ticket-config.json`, `dev/`, or
`Workbooks/`.

Scorecards now expose both machine-readable and human-readable triage fields:
scenario ID, ticket ID, omic, module, runtime, CI lane, status/result, failure
reason, artifact count, and artifact paths. The Markdown scorecard is appended
to the GitHub Actions step summary for quick triage; the JSON scorecard remains
the stable machine artifact.

Local release-gate reproduction:

```bash
Rscript tools/ci/run-module-ci.R --runtime unit-contract --reporter summary
Rscript tools/ci/run-module-ci.R --omic proteomics --runtime unit-contract --reporter summary
Rscript tools/ci/run-module-ci.R --omic metabolomics --runtime unit-contract --reporter summary
Rscript tools/ci/run-module-ci.R --omic lipidomics --runtime unit-contract --reporter summary
Rscript tools/ci/module-ci-scorecard.R --artifact-dir tests/testthat/_module_ci_artifacts
```

Failure triage should start with the scorecard row, then the run manifest, then
the uploaded artifact directory named by that row. Promotion to a final release
tag should wait for the release-candidate gate to pass on the release branch or
manual release dispatch; any follow-up fix can be committed normally and the
gate rerun before tagging.

## Impact-Aware CI Routing

MCI-026 through MCI-033 replace the fixed push/PR module matrix with a
fail-closed impact router. The release/nightly/full gates remain complete
corpus runs; only routine push and pull-request checks are dynamically narrowed
when the changed paths are safely understood.

The routing contract has three owned pieces:

- `tools/ci/impact-map.json` is the versioned source ownership map. Each rule
  has a stable ID, path regexes, owner metadata, a reason string, an impact
  level, module-CI rows, optional browser rows, optional E2E lanes, and
  escalation flags for cross-omic or report/export checks.
- `tools/ci/detect-impact.R` resolves explicit `--changed-files` input or
  `--base-ref/--head-ref` git diffs into stable JSON plus a Markdown summary.
  It validates rule shape before trusting the map and fails closed to broad
  impact on resolver errors, malformed regexes, all-zero push bases, unknown
  risky sources, unknown CI/test files, and invalid maps.
- `.github/workflows/module-ci.yml` runs `detect-impact` on push, pull request,
  and manual `gate=push|all`, then fans out impacted module, browser, E2E,
  cross-omic, and report/export jobs from the resolver matrices.

The resolver JSON is intentionally suitable for GitHub Actions without a
translation layer. Required fields include `schema_version`, `resolver_version`,
`map_schema_version`, `impact_level`, `changed_files`, `matched_rules`,
`no_op_files`, `suppressed_no_op_files`, `unmatched_files`, `module_matrix`,
`browser_matrix`, `e2e_matrix`, `run_cross_omic`, `run_report_export`,
`impact_summary`, and `errors`. Matrix rows include artifact directories so the
impact scorecard can be joined to module-CI and E2E run manifests.

Impact levels are conservative:

- `none` means only safe no-op paths matched, such as docs-only changes with no
  risky file in the same change set.
- `targeted` means module-owned source changes selected only the relevant omic,
  module family, and E2E lanes.
- `omic` is reserved for omic-wide changes that are broader than one module but
  do not require all omics.
- `broad` means shared state, file management, report/export, CI, tests,
  fixtures, unknown source paths, resolver errors, or mixed risky changes
  selected all representative safety gates.
- `full` is available for future always-full release-sensitive surfaces; full
  release/nightly coverage is still driven by the existing full-gate workflow
  triggers.

Current source ownership policy:

- `R/mod_lipid_*.R` and `R/func_lipid_*.R` route to lipidomics module families
  and the lipidomics canonical E2E lane.
- `R/mod_metab_*.R` and `R/func_metab_*.R` route to metabolomics module
  families and the LC/GC/combined metabolomics E2E lane set.
- `R/mod_prot_*.R` and `R/func_prot_*.R` route to proteomics module families
  and the appropriate DIA, DIA+LIMPA, TMT/LFQ/FragPipe, enrichment, or report
  E2E lane.
- Shared workflow state, R6/S4 state, file management, config, general design,
  reporting, enrichment dispatch, plotting, QC helpers, fixtures, tests, and CI
  code escalate broad so all omics and representative E2Es remain protected.
- Unknown tracked `R/*.R`, unknown test files, unknown CI files, delete/rename
  ambiguity, or resolver failures cannot produce `impact_level=none`.

Local routing examples:

```bash
Rscript tools/ci/detect-impact.R --changed-files R/mod_lipid_norm_server.R --output /tmp/mci-impact-targeted-lipid.json --summary /tmp/mci-impact-targeted-lipid.md
Rscript tools/ci/detect-impact.R --changed-files R/mod_metab_da_server.R --output /tmp/mci-impact-targeted-metab.json --summary /tmp/mci-impact-targeted-metab.md
Rscript tools/ci/detect-impact.R --changed-files R/func_prot_limpa_qc_helpers.R --output /tmp/mci-impact-targeted-prot-limpa.json --summary /tmp/mci-impact-targeted-prot-limpa.md
Rscript tools/ci/detect-impact.R --changed-files R/utils_workflow_state.R --output /tmp/mci-impact-broad-state.json --summary /tmp/mci-impact-broad-state.md
Rscript tools/ci/detect-impact.R --changed-files docs/qa/module-ci.md --output /tmp/mci-impact-docs.json --summary /tmp/mci-impact-docs.md
Rscript tools/ci/detect-impact.R --changed-files R/new_unmapped_surface.R --output /tmp/mci-impact-unknown-source.json --summary /tmp/mci-impact-unknown-source.md
```

Local impacted job reproduction:

```bash
Rscript tools/ci/run-module-ci.R --omic lipidomics --module normalization --runtime unit-contract --reporter summary --artifact-dir tests/testthat/_module_ci_artifacts/impact/module/lipidomics-normalization
Rscript tools/ci/run-e2e-ci.R --lane lipid_canonical --filter '^e2e-lipidomics-canonical$' --reporter summary --artifact-dir tests/testthat/_e2e_artifacts/impact/lipidomics-canonical
Rscript tools/ci/run-module-ci.R --runtime unit-contract --reporter summary
Rscript tools/ci/module-ci-scorecard.R --artifact-dir tests/testthat/_module_ci_artifacts
```

The E2E manifest at `tests/testdata/e2e/manifest.json` now declares routing
metadata for every current lane: `test_filter`, `module_families`,
`touchpoints`, `critical_shared`, `report_export`, and `cross_omic_packs`.
`tools/ci/run-e2e-ci.R` accepts one or more `--lane` values, one or more
`--filter` values, writes `e2e-run-manifest.json`, and respects the existing
direct `testthat` invocation model for full `filter = "^e2e-"` nightly runs.
The runner also supports `--dry-run true` for router self-tests only; CI does
not set dry-run for impacted E2E jobs.

Impact scorecards and summaries should be read in this order:

- Start with `impact-routing` artifacts and `impact.md`. Confirm the impact
  level, changed files, matched rules, no-op files, suppressed no-op files,
  unmatched files, selected module rows, selected browser rows, selected E2Es,
  resolver errors, and local reproduction commands.
- For each selected module row, inspect the corresponding
  `module-ci-impact-*` artifact, `module-ci-run-manifest.json`, and
  `module-ci-scorecard.md`.
- For each selected E2E row, inspect the `e2e-impact-*` artifact and
  `e2e-run-manifest.json`.
- If cross-omic or report/export flags are true, inspect the dedicated
  `e2e-impact-cross-omic` or `e2e-impact-report-export` artifacts.

Failure triage rules:

- False narrow routing: add or broaden an impact-map rule before trusting the
  push result. Reproduce with `detect-impact.R`, verify selected matrices, then
  add a focused `ci-impact-*` test.
- False broad routing: keep the broad result until a narrower rule is justified
  by source ownership and E2E touchpoints. Broad is noisy but safe; narrow can
  hide regressions.
- Resolver error: treat as broad until the error is fixed. The summary should
  include the resolver error and broad fallback reason.
- Missing map ownership: any new source family must add a rule, reason, owner,
  module consequences, E2E consequences, and source-map coverage test in the
  same change.
- Downstream module/E2E failure: debug from the module or E2E scorecard first;
  only change the router if the lane was incorrectly selected.

Maintenance rules:

- New modules must add module-CI manifest scenarios, impact-map source rules,
  resolver tests, and docs in the same change.
- New omics must define module families, fixture packs, E2E lanes, cross-omic
  state expectations, and branch-protection guidance before enabling targeted
  routing.
- New E2E lanes must declare `test_filter`, module families, touchpoints,
  report/export coverage, cross-omic packs, and fixture paths in the E2E
  manifest.
- New shared helpers should default to broad impact until ownership is proven
  narrower by tests and reviewer consensus.
- Routing docs, impact maps, CI scripts, helper tests, workflow YAML, fixtures,
  and manifests are self-protecting paths and should trigger broad checks.

Branch protection guidance:

- Require `detect impact` and `impact routing summary` for push/PR protection
  once dynamic routing is enabled on GitHub.
- Treat `module-ci impacted (*)`, `module-ci browser impacted (*)`, and
  `e2e impacted (*)` as matrix jobs whose concrete names may vary by change;
  the stable summary job is the branch-protection anchor.
- Keep `release candidate promotion gate` required for release branches and
  version tags. It remains blocked on full gates, not targeted routing.
- Use manual `workflow_dispatch gate=all` when touching unfamiliar shared
  infrastructure, moving source files across omics, changing fixture schemas, or
  preparing a final release.

Residual risk is path-based ownership. The router cannot infer semantic
coupling from code at runtime; it depends on maintained regex ownership,
complete E2E touchpoint metadata, and fail-closed defaults. When modules move or
helpers become shared, update the map first and accept a temporary broad route
until targeted coverage is proven.
