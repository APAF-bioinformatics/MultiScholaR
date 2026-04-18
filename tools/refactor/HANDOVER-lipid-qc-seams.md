# Lipid QC Seam Map

## Goal

Document the active safe stabilization seams for
[func_lipid_qc.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_qc.R:1)
while keeping the live lipid-QC entry points behaviorally stable.

## Current Position In The Flow

- April 13, 2026 classification for
  [func_lipid_qc.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_qc.R:1)
  is `review` and `direct-extraction-ready`.
- The first bounded helper wave is now applied live via
  [tools/refactor/manifest-lipid-qc-wave1.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-qc-wave1.yml:1)
  into
  [R/func_lipid_qc_filtering_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_qc_filtering_helpers.R:1),
  covering:
  - `calculateLipidFilteringAssayMetrics()`
  - `prepareLipidFilteringContext()`
  - `finalizeLipidFilteringStep()`
  - `resolveLipidFilteringPlotSaveDir()`
  - `saveLipidFilteringPlots()`
  - `returnLipidFilteringPlots()`
- Live
  [updateLipidFiltering()](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_qc.R:307)
  now delegates to those helpers without changing the outer wrapper contract.
- The live collate artifact now exists at
  [tools/refactor/collate-lipid-qc-wave1.txt](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/collate-lipid-qc-wave1.txt:1).
- The second bounded helper wave is now applied live via
  [tools/refactor/manifest-lipid-qc-wave2.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-qc-wave2.yml:1)
  into
  [R/func_lipid_qc_progress_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_qc_progress_helpers.R:1),
  covering:
  - `getFilteringProgressLipidomics()`
  - `updateFilteringProgressLipidomics()`
  - `countUniqueLipids()`
  - `countLipidsPerSample()`
  - `calculateLipidMissingness()`
  - `calculateLipidSumIntensityPerSample()`
  - `calculateTotalUniqueLipidsAcrossAssays()`
- The live collate artifact for the second wave now exists at
  [tools/refactor/collate-lipid-qc-wave2.txt](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/collate-lipid-qc-wave2.txt:1).
- The third bounded helper wave is now applied live via
  [tools/refactor/manifest-lipid-qc-wave3.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-qc-wave3.yml:1)
  into
  [R/func_lipid_qc_reporting_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_qc_reporting_helpers.R:1),
  covering:
  - `calculateLipidCVs()`
  - `getLipidInternalStandardMetrics()`
  - `generateLipidFilteringPlots()`
- The live collate artifact for the third wave now exists at
  [tools/refactor/collate-lipid-qc-wave3.txt](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/collate-lipid-qc-wave3.txt:1).
- The fourth bounded helper wave is now applied live via
  [tools/refactor/manifest-lipid-qc-wave4.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-qc-wave4.yml:1)
  into
  [R/func_lipid_qc_support_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_qc_support_helpers.R:1),
  covering:
  - `resolveLipidDuplicateFeaturesByIntensity()`
  - `calculateLipidPairCorrelation()`
- The live collate artifact for the fourth wave now exists at
  [tools/refactor/collate-lipid-qc-wave4.txt](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/collate-lipid-qc-wave4.txt:1).
- `DESCRIPTION` `Collate:` now includes
  `func_lipid_qc_filtering_helpers.R`
  and
  `func_lipid_qc_progress_helpers.R`
  and
  `func_lipid_qc_reporting_helpers.R`
  and
  `func_lipid_qc_support_helpers.R`
  after `func_lipid_qc.R`.

## Existing Safety Net

- Focused characterization gate:
  - [tests/testthat/test-lipid-01-qc-filtering-helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tests/testthat/test-lipid-01-qc-filtering-helpers.R:1)
- Replay command:
  - `Rscript -e "source('R/allGenerics.R'); source('R/func_general_s4_objects.R'); source('R/func_lipid_s4_objects.R'); source('R/func_lipid_import.R'); source('R/func_lipid_qc.R'); source('R/func_lipid_qc_progress_helpers.R'); source('R/func_lipid_qc_reporting_helpers.R'); source('R/func_lipid_qc_filtering_helpers.R'); source('R/func_lipid_qc_support_helpers.R'); testthat::test_file('tests/testthat/test-lipid-01-qc-filtering-helpers.R', stop_on_failure = TRUE)"`
- Gate coverage now characterizes:
  - `getFilteringProgressLipidomics()`
  - `updateFilteringProgressLipidomics()`
  - `countUniqueLipids()`
  - `countLipidsPerSample()`
  - `calculateLipidMissingness()`
  - `calculateLipidSumIntensityPerSample()`
  - `calculateTotalUniqueLipidsAcrossAssays()`
  - `prepareLipidFilteringContext()`
  - `calculateLipidFilteringAssayMetrics()`
  - `calculateLipidCVs()` via the assay-metrics characterization shell
  - `getLipidInternalStandardMetrics()` via the assay-metrics characterization shell
  - `finalizeLipidFilteringStep()`
  - `generateLipidFilteringPlots()` via the finalize-and-plot characterization shell
  - `resolveLipidFilteringPlotSaveDir()`
  - `saveLipidFilteringPlots()`
  - `returnLipidFilteringPlots()`
  - `resolveLipidDuplicateFeaturesByIntensity()`
  - `calculateLipidPairCorrelation()`

The source-based gate is used in this lane because
[tools/test_with_renv.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/test_with_renv.R:1)
cannot run here without a checked-in `renv/activate.R`, and full `load_all()`
currently stops on absent optional imports in the sandbox.

## Next Safe Checkpoint

- No further seam work is required for this backlog target.
- If later cleanup is desired, treat the remaining
  [R/func_lipid_qc.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_qc.R:1)
  wrapper/method shells as optional follow-up rather than an active blocker.

## Safe Compact Checkpoint

It is safe to compact here.

## Notes

- No prior target handover existed for this lane.
- The active stop point is now a live wave-4 apply for the remaining
  duplicate-resolution and pair-correlation helpers, under the same passing
  focused characterization gate.
- The post-apply checker and the focused source-based lipid-QC helper gate both
  reran green after the live wave-4 apply.
- Live
  [R/func_lipid_qc.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_qc.R:1)
  now measures `474` lines while
  [R/func_lipid_qc_filtering_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_qc_filtering_helpers.R:1)
  measures `398` lines and
  [R/func_lipid_qc_progress_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_qc_progress_helpers.R:1)
  measures `296` lines and
  [R/func_lipid_qc_reporting_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_qc_reporting_helpers.R:1)
  measures `731` lines and
  [R/func_lipid_qc_support_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_qc_support_helpers.R:1)
  measures `159` lines.
- `func_lipid_qc.R` is now back inside the ideal size band for this backlog and
  no longer blocks the lipidomics QC lane.
