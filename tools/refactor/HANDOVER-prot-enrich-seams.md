# Proteomics Enrichment Seam Map

## Goal

Stabilize `mod_prot_enrich.R` one bounded checkpoint at a time while keeping
`mod_prot_enrich_server()` as the public wrapper identity.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/mod_prot_enrich.R`
- classification: `complete`
- next step: `Archive this target as stabilized; no further god-module checkpoint is queued for mod_prot_enrich.R.`
- latest completed checkpoint: `prot_enrich_run_dependency_annotation_wave29_apply`
- latest completed checkpoint summary:
  applied the reviewed exact-source run-dependency and annotation helper wave
  from `tools/refactor/manifest-prot-enrich-server-wave29.yml` into
  `R/mod_prot_enrich_server_helpers.R`, removed the moved helper cluster from
  `R/mod_prot_enrich.R`, reran the post-apply checker, and reran the focused
  enrichment gate green with `1243` passing expectations.
- current recovery state:
  - live `R/mod_prot_enrich_server_helpers.R` now carries the accumulated
    exact-source helper waves 1 through 25 plus the live wave-27
    organism-filter helper cluster
  - live `R/mod_prot_enrich_server.R` now carries the exact-source
    `mod_prot_enrich_server()` entrypoint extracted from
    `R/mod_prot_enrich.R`
  - `tools/refactor/manifest-prot-enrich-server-wave26.yml` remains the
    reviewed exact-source manifest for the applied server-entrypoint wave
  - `tools/refactor/manifest-prot-enrich-server-wave27.yml` now remains the
    reviewed exact-source manifest for the applied organism-filter helper wave
    for
    `resolveProtEnrichOrganismMapping()`,
    `applyProtEnrichOrganismFilter()`, and
    `persistProtEnrichOrganismFilterMetadata()`
  - `tools/refactor/manifest-prot-enrich-server-wave28.yml` now remains the
    reviewed exact-source manifest for the applied context/contrast helper wave
    for
    `resolveProtEnrichCurrentS4Object()`,
    `buildProtEnrichContrastChoices()`,
    `resolveProtEnrichRawContrastName()`,
    `createProtEnrichRawContrastNameReactive()`,
    `resolveProtEnrichSelectedContrastResults()`, and
    `resolveProtEnrichSelectedDaResults()`
  - `tools/refactor/manifest-prot-enrich-server-wave29.yml` now remains the
    reviewed exact-source manifest for the applied run-dependency and
    annotation helper wave for
    `resolveProtEnrichRunDependencies()`,
    `resolveProtEnrichOutputDirectories()`,
    `resolveProtEnrichUniprotAnnotations()`, and
    `resolveProtEnrichAnnotationMatching()`
  - `tools/refactor/staging/prot-enrich-server-wave26/R/mod_prot_enrich_server.R`
    remains the reviewed staged exact-source source for the live apply
  - `tools/refactor/staging/prot-enrich-server-wave26/collate-prot-enrich-server-wave26.txt`
    remains the reviewed helper/server ordering artifact for the applied wave
  - `tools/refactor/staging/prot-enrich-server-wave27/R/mod_prot_enrich_server_helpers.R`
    now holds the reviewed staged exact-source organism-filter helper wave
    extracted from live `R/mod_prot_enrich.R`
  - `tools/refactor/staging/prot-enrich-server-wave27/collate-prot-enrich-server-wave27.txt`
    now records the staged helper target ordering used for the applied wave-27
    organism-filter helper checkpoint
  - `tools/refactor/staging/prot-enrich-server-wave28/R/mod_prot_enrich_server_helpers.R`
    now remains the reviewed staged exact-source context/contrast helper wave
    extracted from live `R/mod_prot_enrich.R` and used for the live wave-28
    apply
  - `tools/refactor/staging/prot-enrich-server-wave28/collate-prot-enrich-server-wave28.txt`
    now records the staged helper target ordering used for the live wave-28
    context/contrast helper apply
  - `tools/refactor/staging/prot-enrich-server-wave29/R/mod_prot_enrich_server_helpers.R`
    now holds the reviewed staged exact-source run-dependency and annotation
    helper wave extracted from live `R/mod_prot_enrich.R`
  - `tools/refactor/staging/prot-enrich-server-wave29/collate-prot-enrich-server-wave29.txt`
    now records the staged helper target ordering for the reviewed wave-29
    extraction checkpoint
  - live `R/mod_prot_enrich.R` no longer contains
    `prepareProtEnrichAnalysisBodySetup()`,
    `runProtEnrichAnalysisBody()`, or `mod_prot_enrich_server()`
  - live `R/mod_prot_enrich.R` no longer contains
    `resolveProtEnrichOrganismMapping()`,
    `applyProtEnrichOrganismFilter()`, or
    `persistProtEnrichOrganismFilterMetadata()`
  - live `R/mod_prot_enrich.R` no longer contains
    `resolveProtEnrichCurrentS4Object()`,
    `buildProtEnrichContrastChoices()`,
    `resolveProtEnrichRawContrastName()`,
    `createProtEnrichRawContrastNameReactive()`,
    `resolveProtEnrichSelectedContrastResults()`, or
    `resolveProtEnrichSelectedDaResults()`
  - live `R/mod_prot_enrich.R` no longer contains
    `resolveProtEnrichRunDependencies()`,
    `resolveProtEnrichOutputDirectories()`,
    `resolveProtEnrichUniprotAnnotations()`, or
    `resolveProtEnrichAnnotationMatching()`
  - live `R/mod_prot_enrich.R` now delegates the
    display/status output registration cluster via
    `setupProtEnrichDisplayStatusOutputBootstrap()` while preserving
    `mod_prot_enrich_server()` as the public wrapper identity
  - live `R/mod_prot_enrich.R` now delegates the
    results/summary output registration cluster via
    `setupProtEnrichResultsSummaryOutputBootstrap()` while preserving
    `mod_prot_enrich_server()` as the public wrapper identity
  - live `R/mod_prot_enrich.R` now delegates the
    remaining observer registration cluster via
    `setupProtEnrichObserverRegistrationBootstrap()` while preserving
    `mod_prot_enrich_server()` as the public wrapper identity
  - live `R/mod_prot_enrich.R` now delegates the
    remaining run/results/plot/download registration shell via
    `setupProtEnrichRunOutputDownloadBootstrap()` while preserving
    `mod_prot_enrich_server()` as the public wrapper identity
  - live `R/mod_prot_enrich_server_helpers.R` now also defines
    `setupProtEnrichDisplayStatusOutputBootstrap()` for the wrapper shell
    alongside the existing display/status setup seams
  - live `R/mod_prot_enrich_server_helpers.R` now also defines
    `setupProtEnrichResultsSummaryOutputBootstrap()` for the wrapper shell
    alongside the existing per-output setup helpers
  - live `R/mod_prot_enrich_server_helpers.R` now also defines
    `setupProtEnrichObserverRegistrationBootstrap()` for the wrapper shell
    alongside the existing observer setup seams
  - live `R/mod_prot_enrich_server_helpers.R` now also defines
    `setupProtEnrichRunOutputDownloadBootstrap()` for the wrapper shell
    alongside the existing run observer, results/summary output, plot
    output, and results download setup seams
  - live `R/mod_prot_enrich_server_helpers.R` now also contains the live
    wave-28 context/contrast helper cluster for
    `resolveProtEnrichCurrentS4Object()`,
    `buildProtEnrichContrastChoices()`,
    `resolveProtEnrichRawContrastName()`,
    `createProtEnrichRawContrastNameReactive()`,
    `resolveProtEnrichSelectedContrastResults()`, and
    `resolveProtEnrichSelectedDaResults()`
  - live `R/mod_prot_enrich_server_helpers.R` now also contains the live
    wave-29 run-dependency and annotation helper cluster for
    `resolveProtEnrichRunDependencies()`,
    `resolveProtEnrichOutputDirectories()`,
    `resolveProtEnrichUniprotAnnotations()`, and
    `resolveProtEnrichAnnotationMatching()`
  - `DESCRIPTION` `Collate:` now preserves the live helper/server ordering as
    `mod_prot_enrich_server_helpers.R`, `mod_prot_enrich_server.R`,
    `mod_prot_enrich.R`
  - `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-prot-enrich-server-wave26.yml`
    passes after the live wave-26 apply
  - `Rscript tools/refactor/verify_refactor.R --manifest tools/refactor/manifest-prot-enrich-server-wave27.yml`
    passes for the reviewed wave-27 organism-filter helper manifest
  - `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-prot-enrich-server-wave27.yml`
    passes after the live wave-27 organism-filter helper apply
  - `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-prot-enrich-server-wave28.yml`
    passes after the live wave-28 context/contrast helper apply
  - `Rscript tools/refactor/verify_refactor.R --manifest tools/refactor/manifest-prot-enrich-server-wave29.yml`
    passes for the reviewed wave-29 run-dependency and annotation helper
    manifest
  - `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-prot-enrich-server-wave29.yml`
    passes after the live wave-29 run-dependency and annotation helper apply
  - `Rscript tools/test_with_renv.R tests/testthat/test-prot-11-enrichment-module-contracts.R`
    passes again with `1243` expectations after the live wave-29 apply
  - current target classification is now `complete` at `418` lines with `1`
    top-level function, max function length `351`, no `moduleServer()`
    calls, and no direct `observeEvent()` or `render*()` calls left in the
    target file
- next bounded stop point:
  target complete; `mod_prot_enrich.R` is now a bounded UI-only file in the
  ideal size band, so no additional stabilization checkpoint is queued for
  this handover

## Existing Safety Net

- `Rscript tools/test_with_renv.R tests/testthat/test-prot-11-enrichment-module-contracts.R`

## Notes

- No loop-control artifacts were changed in this iteration.
- Root cause of the wave-21 block was tooling, not extraction logic:
  the old live apply path called `extract_blocks.R --write-targets` in a way
  that rewrote shared helper targets from only the current manifest entries.
- `tools/refactor/extract_blocks.R` now has a preserve-existing-targets mode
  for live apply.
- `tools/refactor/apply_wave.py` now runs transactionally:
  it backs up touched sources/targets, preserves pre-existing target symbols,
  and restores the repo files automatically if the apply or post-apply checks
  fail.
- The installed skill `scripts/apply_wave.py` was synced to the hardened repo
  version so future loop iterations use the safer apply path.
- The staged wave-12 helper file now lives at
  `tools/refactor/staging/prot-enrich-server-wave12/R/mod_prot_enrich_server_helpers.R`.
- The staged wave-13 helper file now lives at
  `tools/refactor/staging/prot-enrich-server-wave13/R/mod_prot_enrich_server_helpers.R`
  and was the reviewed exact-source source for the live wave-13 apply.
- The staged wave-14 helper file now lives at
  `tools/refactor/staging/prot-enrich-server-wave14/R/mod_prot_enrich_server_helpers.R`
  and was the reviewed exact-source source for the live wave-14 apply.
- The staged wave-15 helper file now lives at
  `tools/refactor/staging/prot-enrich-server-wave15/R/mod_prot_enrich_server_helpers.R`
  and was the reviewed exact-source source for the live wave-15 apply.
- The staged wave-16 helper file now lives at
  `tools/refactor/staging/prot-enrich-server-wave16/R/mod_prot_enrich_server_helpers.R`
  and was the reviewed exact-source source for the live wave-16 apply.
- The staged wave-17 helper file now lives at
  `tools/refactor/staging/prot-enrich-server-wave17/R/mod_prot_enrich_server_helpers.R`
  and was the reviewed exact-source source for the live wave-17 apply.
- The staged wave-17 collate artifact now lives at
  `tools/refactor/staging/prot-enrich-server-wave17/collate-prot-enrich-server-wave17.txt`.
- The staged wave-18 helper file now lives at
  `tools/refactor/staging/prot-enrich-server-wave18/R/mod_prot_enrich_server_helpers.R`
  and was the reviewed exact-source source for the live wave-18 apply.
- The staged wave-18 collate artifact now lives at
  `tools/refactor/staging/prot-enrich-server-wave18/collate-prot-enrich-server-wave18.txt`.
- The staged wave-19 helper file now lives at
  `tools/refactor/staging/prot-enrich-server-wave19/R/mod_prot_enrich_server_helpers.R`.
- The staged wave-19 collate artifact now lives at
  `tools/refactor/staging/prot-enrich-server-wave19/collate-prot-enrich-server-wave19.txt`.
- The staged wave-20 helper file now lives at
  `tools/refactor/staging/prot-enrich-server-wave20/R/mod_prot_enrich_server_helpers.R`.
- The staged wave-20 collate artifact now lives at
  `tools/refactor/staging/prot-enrich-server-wave20/collate-prot-enrich-server-wave20.txt`.
- The staged wave-21 helper file now lives at
  `tools/refactor/staging/prot-enrich-server-wave21/R/mod_prot_enrich_server_helpers.R`.
- The staged wave-21 collate artifact now lives at
  `tools/refactor/staging/prot-enrich-server-wave21/collate-prot-enrich-server-wave21.txt`.
- The staged wave-22 helper file now lives at
  `tools/refactor/staging/prot-enrich-server-wave22/R/mod_prot_enrich_server_helpers.R`.
- The staged wave-22 collate artifact now lives at
  `tools/refactor/staging/prot-enrich-server-wave22/collate-prot-enrich-server-wave22.txt`.
- The staged wave-23 helper file now lives at
  `tools/refactor/staging/prot-enrich-server-wave23/R/mod_prot_enrich_server_helpers.R`.
- The staged wave-23 collate artifact now lives at
  `tools/refactor/staging/prot-enrich-server-wave23/collate-prot-enrich-server-wave23.txt`.
- The staged wave-24 helper file now lives at
  `tools/refactor/staging/prot-enrich-server-wave24/R/mod_prot_enrich_server_helpers.R`.
- The staged wave-24 collate artifact now lives at
  `tools/refactor/staging/prot-enrich-server-wave24/collate-prot-enrich-server-wave24.txt`.
- The staged wave-25 helper file now lives at
  `tools/refactor/staging/prot-enrich-server-wave25/R/mod_prot_enrich_server_helpers.R`
  and was the reviewed exact-source source for the live wave-25 apply.
- The staged wave-26 server entrypoint file now lives at
  `tools/refactor/staging/prot-enrich-server-wave26/R/mod_prot_enrich_server.R`.
- The staged wave-26 collate artifact now lives at
  `tools/refactor/staging/prot-enrich-server-wave26/collate-prot-enrich-server-wave26.txt`.
- The staged wave-27 organism-filter helper file now lives at
  `tools/refactor/staging/prot-enrich-server-wave27/R/mod_prot_enrich_server_helpers.R`.
- The staged wave-27 collate artifact now lives at
  `tools/refactor/staging/prot-enrich-server-wave27/collate-prot-enrich-server-wave27.txt`.
- The staged wave-28 context/contrast helper file now lives at
  `tools/refactor/staging/prot-enrich-server-wave28/R/mod_prot_enrich_server_helpers.R`.
- The staged wave-28 collate artifact now lives at
  `tools/refactor/staging/prot-enrich-server-wave28/collate-prot-enrich-server-wave28.txt`.
- The staged wave-29 run-dependency and annotation helper file now lives at
  `tools/refactor/staging/prot-enrich-server-wave29/R/mod_prot_enrich_server_helpers.R`.
- The staged wave-29 collate artifact now lives at
  `tools/refactor/staging/prot-enrich-server-wave29/collate-prot-enrich-server-wave29.txt`.
- `R/mod_prot_enrich_server_helpers.R` now contains the exact-source helper
  accumulation from waves 1 through 25 after the live wave-25 apply.
- `R/mod_prot_enrich_server_helpers.R` now also contains the live wave-27
  organism-filter helper cluster for `resolveProtEnrichOrganismMapping()`,
  `applyProtEnrichOrganismFilter()`, and
  `persistProtEnrichOrganismFilterMetadata()`.
- Live `R/mod_prot_enrich.R` no longer contains the wave-24 analysis-body
  setup helper `prepareProtEnrichAnalysisBodySetup()`.
- Live `R/mod_prot_enrich.R` no longer contains the wave-25 analysis-body
  wrapper helper `runProtEnrichAnalysisBody()`.
- Live `R/mod_prot_enrich.R` no longer contains the wave-27 organism-filter
  helper cluster for `resolveProtEnrichOrganismMapping()`,
  `applyProtEnrichOrganismFilter()`, and
  `persistProtEnrichOrganismFilterMetadata()`.
- Live `R/mod_prot_enrich.R` no longer contains the wave-16 render-helper
  cluster for `renderProtEnrichGprofilerResultsTable()`,
  `renderProtEnrichGprofilerPlot()`, and
  `renderProtEnrichClusterProfilerPlot()`.
- Live `R/mod_prot_enrich.R` no longer contains the wave-18 download-helper
  cluster for `buildProtEnrichResultsDownloadFilename()` and
  `writeProtEnrichResultsDownloadArchive()`.
- Live `R/mod_prot_enrich.R` no longer contains the wave-19 analysis-body
  persistence-helper cluster for `buildProtEnrichAnalysisResultsPayload()`,
  `propagateProtEnrichResultsArgs()`, and
  `propagateProtEnrichUiParams()`.
- Live `R/mod_prot_enrich.R` no longer contains the wave-20 analysis-body
  completion-helper cluster for `updateProtEnrichStateManagerUiParams()`,
  `saveProtEnrichCompletedState()`, and
  `completeProtEnrichTabStatus()`.
- Live `R/mod_prot_enrich.R` no longer contains the wave-21 analysis-body
  completion-feedback cluster for `completeProtEnrichProgress()`,
  `notifyProtEnrichCompletion()`, `notifyProtEnrichAnalysisError()`,
  `logProtEnrichAnalysisError()`, `messageProtEnrichAnalysisError()`, and
  `reportProtEnrichAnalysisError()`.
- Live `R/mod_prot_enrich.R` no longer contains the wave-22 analysis-body
  persistence/finalization cluster for
  `captureProtEnrichPostProcessResults()`,
  `persistProtEnrichAnalysisResults()`, and
  `finalizeProtEnrichAnalysisBodyResults()`.
- Live `R/mod_prot_enrich.R` no longer contains the wave-23 analysis-body
  execution/results cluster for
  `resolveProtEnrichAnalysisInputColumns()`,
  `buildProtEnrichProcessEnrichmentsArgs()`,
  `prepareProtEnrichProcessExecution()`,
  `executeProtEnrichProcessEnrichments()`, and
  `buildProtEnrichAllContrastResults()`.
- Live `R/mod_prot_enrich_server_helpers.R` now carries the accumulated
  exact-source helper waves 1 through 23 after the live wave-23 apply;
  `DESCRIPTION` `Collate:` remained unchanged.
- `DESCRIPTION` `Collate:` remained unchanged because
  `R/mod_prot_enrich_server_helpers.R` was already loaded before
  `R/mod_prot_enrich.R`.
