# General File Management Seam Map

## Goal

Freeze bounded helper surfaces for `func_general_filemgmt.R` while exact-source
waves drain the file without changing public wrapper behavior.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R`
- classification: `direct-extraction-ready`
- checkpoint reached: `general_filemgmt_na_duplicate_live_seam`
- next step: `None required for this target; R/func_general_filemgmt.R is now within the ideal band at 495 lines with only setupDirectories() and loadDependencies() remaining live.`

## Existing Safety Net

- `tests/testthat/test-general-filemgmt-path-contracts.R`
- `tests/testthat/test-general-filemgmt-config-contracts.R`
- `tests/testthat/test-general-filemgmt-config-roundtrip-contracts.R`
- `tests/testthat/test-general-filemgmt-export-contracts.R`
- `tests/testthat/test-general-filemgmt-utility-contracts.R`
- `tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R`
- `tests/testthat/test-general-filemgmt-s4-contracts.R`
- `tests/testthat/test-prot-12-summary-module-contracts.R`

## Notes

Manual target bucket 12 for shared file-management stabilization.

- This target previously had no dedicated handover; this file now records the
  first safe stop point for `R/func_general_filemgmt.R`.
- April 18, 2026 stabilize-mode characterization checkpoint added
  [tests/testthat/test-general-filemgmt-path-contracts.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-general-filemgmt-path-contracts.R:1)
  to freeze the initial paths/directories surface before any staged wave.
- The new focused gate exercises:
  - `setupDirectories()`
  - `setupAndShowDirectories()`
  - `getProjectPaths()`
  - `createDirectoryIfNotExists()`
  - `createDirIfNotExists()`
- The gate passed through:
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
- The first staged extraction wave is now authored in
  [tools/refactor/manifest-general-filemgmt-paths-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-paths-wave1.yml:1)
  and verified cleanly against the live source tree.
- The staged wave materialized
  [func_general_filemgmt_path_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-paths-wave1/R/func_general_filemgmt_path_helpers.R:1)
  with:
  - `getProjectPaths()`
  - `createDirectoryIfNotExists()`
  - `createDirIfNotExists()`
- The staged collate preview is recorded in
  [collate-general-filemgmt-paths-wave1.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-paths-wave1/collate-general-filemgmt-paths-wave1.txt:1).
- April 18, 2026 wave 1 was applied live from
  [tools/refactor/manifest-general-filemgmt-paths-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-paths-wave1.yml:1)
  into:
  - [R/func_general_filemgmt_path_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_path_helpers.R:1)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:697)
- The live apply removed these top-level helpers from
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:697)
  and materialized them in
  [R/func_general_filemgmt_path_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_path_helpers.R:1):
  - `getProjectPaths()`
  - `createDirectoryIfNotExists()`
  - `createDirIfNotExists()`
- The live collate artifact is recorded in
  [tools/refactor/collate-general-filemgmt-paths-wave1.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-filemgmt-paths-wave1.txt:1),
  and `DESCRIPTION` now collates
  `R/func_general_filemgmt_path_helpers.R` immediately before
  `R/func_general_filemgmt.R`.
- The focused gate reran green after the live apply through:
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
- April 18, 2026 the second bounded directory wave was staged from
  [tools/refactor/manifest-general-filemgmt-directory-wave2.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-directory-wave2.yml:1)
  into
  [tools/refactor/staging/general-filemgmt-directory-wave2/R/func_general_filemgmt_directory_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-directory-wave2/R/func_general_filemgmt_directory_helpers.R:1)
  for:
  - `createOutputDir()`
- April 18, 2026 wave 2 was applied live into:
  - [R/func_general_filemgmt_directory_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_directory_helpers.R:1)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:742)
- The live apply removed `createOutputDir()` from
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:742)
  and materialized it in
  [R/func_general_filemgmt_directory_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_directory_helpers.R:1).
- The live collate artifact is recorded in
  [tools/refactor/collate-general-filemgmt-directory-wave2.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-filemgmt-directory-wave2.txt:1),
  and `DESCRIPTION` now collates
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  and then `R/func_general_filemgmt.R`.
- The focused gate reran green again after the live apply through:
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
- April 18, 2026 the third bounded directory wave was staged from
  [tools/refactor/manifest-general-filemgmt-directory-wave3.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-directory-wave3.yml:1)
  into
  [tools/refactor/staging/general-filemgmt-directory-wave3/R/func_general_filemgmt_directory_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-directory-wave3/R/func_general_filemgmt_directory_helpers.R:1)
  for:
  - `setupAndShowDirectories()`
- April 18, 2026 wave 3 was applied live into:
  - [R/func_general_filemgmt_directory_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_directory_helpers.R:33)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1434)
- The live apply removed `setupAndShowDirectories()` from
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1434)
  and materialized it in
  [R/func_general_filemgmt_directory_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_directory_helpers.R:33).
- The live collate artifact is recorded in
  [tools/refactor/collate-general-filemgmt-directory-wave3.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-filemgmt-directory-wave3.txt:1),
  with the existing `DESCRIPTION` collate order unchanged:
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  and then `R/func_general_filemgmt.R`.
- The focused gate reran green again after the live apply through:
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
- `setupDirectories()` remains live in
  [func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:465),
  while `setupAndShowDirectories()` now lives in
  [func_general_filemgmt_directory_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_directory_helpers.R:33).
  This checkpoint intentionally stops before any follow-up seam that would
  split the larger multi-omic `setupDirectories()` body or widen into config,
  results IO, or Rmd helpers.
- April 18, 2026 one bounded live seam split the
  `setupDirectories()` input parsing and omic-type validation block into
  [parseSetupDirectoriesOmicTypes()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:190)
  while keeping
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:465)
  as the public wrapper and structural work surface.
- The focused gate reran green after the seam through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
  with `46` passing expectations.
- April 18, 2026 one bounded live seam split the `setupDirectories()`
  per-omic configuration switch into
  [getSetupDirectoriesOmicConfig()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:233)
  while keeping
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:465)
  as the public wrapper and structural work surface.
- The focused gate reran green again after the seam through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
  with `46` passing expectations.
- Post-seam classification stayed `review` / `direct-extraction-ready`; the
  file now reports `85` top-level functions with a maximum helper span of
  `328` lines.
- April 18, 2026 one bounded live seam split the `setupDirectories()`
  path-list assembly and named-return-path block into
  [buildSetupDirectoriesPathList()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:335)
  while keeping
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:465)
  as the public wrapper and structural work surface.
- The focused gate reran green again after the seam through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
  with `46` passing expectations.
- Post-seam classification stayed `review` / `direct-extraction-ready`; the
  file now reports `86` top-level functions with a maximum helper span of
  `259` lines.
- April 18, 2026 the fourth bounded `setupDirectories()` helper wave was
  authored in
  [tools/refactor/manifest-general-filemgmt-setupdirectories-wave4.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-setupdirectories-wave4.yml:1)
  and staged cleanly into
  [tools/refactor/staging/general-filemgmt-setupdirectories-wave4/R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-setupdirectories-wave4/R/func_general_filemgmt_setupdirectories_helpers.R:1)
  for:
  - `parseSetupDirectoriesOmicTypes()`
  - `getSetupDirectoriesOmicConfig()`
  - `buildSetupDirectoriesPathList()`
- The staged collate preview is recorded in
  [tools/refactor/staging/general-filemgmt-setupdirectories-wave4/collate-general-filemgmt-setupdirectories-wave4.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-setupdirectories-wave4/collate-general-filemgmt-setupdirectories-wave4.txt:1)
  with the intended order:
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`.
- The focused gate reran green after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
  with `46` passing expectations.
- Post-stage classification of the live target stayed `review` /
  `direct-extraction-ready`; the file still reports `86` top-level functions
  with a maximum helper span of `259` lines.
- April 18, 2026 the fourth bounded `setupDirectories()` helper wave was
  applied live from
  [tools/refactor/manifest-general-filemgmt-setupdirectories-wave4.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-setupdirectories-wave4.yml:1)
  into:
  - [R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_setupdirectories_helpers.R:1)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:225)
- The live apply removed these top-level helpers from
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:225)
  and materialized them in
  [R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_setupdirectories_helpers.R:1):
  - `parseSetupDirectoriesOmicTypes()`
  - `getSetupDirectoriesOmicConfig()`
  - `buildSetupDirectoriesPathList()`
- The live collate artifact is recorded in
  [tools/refactor/collate-general-filemgmt-setupdirectories-wave4.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-filemgmt-setupdirectories-wave4.txt:1),
  and `DESCRIPTION` now collates
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`,
  and then `R/func_general_filemgmt.R`.
- The focused gate reran green after the live apply through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
  with `46` passing expectations.
- Post-apply classification of the live target stayed `review` /
  `direct-extraction-ready`; the file now reports `83` top-level functions
  with a maximum helper span of `259` lines.
- April 18, 2026 one bounded live seam split the `setupDirectories()`
  existing-directory overwrite/reuse decision block into
  [handleSetupDirectoriesExistingDirs()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:198)
  while keeping
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:280)
  as the public wrapper and structural work surface.
- The focused gate was extended with a non-interactive reuse-existing
  characterization in
  [tests/testthat/test-general-filemgmt-path-contracts.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-general-filemgmt-path-contracts.R:56)
  and reran green through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
  with `51` passing expectations.
- Post-seam classification stayed `review` / `direct-extraction-ready`; the
  file now reports `84` top-level functions with a maximum helper span of
  `259` lines.
- This stop point keeps the directory-materialization and script-copy body live
  inside
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:280),
  so the next bounded checkpoint should stage one exact-source helper wave for
  `handleSetupDirectoriesExistingDirs()` into
  [R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_setupdirectories_helpers.R:1)
  without widening into results IO or Rmd helpers.
- April 18, 2026 the fifth bounded `setupDirectories()` helper wave was
  authored in
  [tools/refactor/manifest-general-filemgmt-setupdirectories-wave5.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-setupdirectories-wave5.yml:1)
  and staged cleanly into
  [tools/refactor/staging/general-filemgmt-setupdirectories-wave5/R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-setupdirectories-wave5/R/func_general_filemgmt_setupdirectories_helpers.R:1)
  for:
  - `handleSetupDirectoriesExistingDirs()`
- The staged collate preview is recorded in
  [tools/refactor/staging/general-filemgmt-setupdirectories-wave5/collate-general-filemgmt-setupdirectories-wave5.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-setupdirectories-wave5/collate-general-filemgmt-setupdirectories-wave5.txt:1)
  with the intended order:
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`.
- The focused gate reran green after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
  with `51` passing expectations.
- Post-stage classification of the live target stayed `review` /
  `direct-extraction-ready`; the file still reports `84` top-level functions
  with a maximum helper span of `259` lines.
- This stop point keeps live `R/` sources unchanged for the new seam, so the
  next bounded checkpoint should review and, if approved, apply
  `manifest-general-filemgmt-setupdirectories-wave5.yml` into
  [R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_setupdirectories_helpers.R:1)
  while keeping the remaining directory-materialization and script-copy body
  live inside
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:280).
- April 18, 2026 the reviewed `setupDirectories()` wave5 helper slice was
  applied live into:
  - [R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_setupdirectories_helpers.R:244)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:196)
- The live apply moved `handleSetupDirectoriesExistingDirs()` out of
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:226)
  while keeping `setupDirectories()` as the public wrapper and primary work
  surface.
- The emitted collate artifact now exists at
  [tools/refactor/collate-general-filemgmt-setupdirectories-wave5.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-filemgmt-setupdirectories-wave5.txt:1).
- The focused gate reran green after the live apply through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
  with `51` passing expectations.
- Post-apply classification of the live target stayed `review` /
  `direct-extraction-ready`; the file now reports `83` top-level functions
  with a maximum helper span of `259` lines.
- April 18, 2026 one bounded live seam split the duplicated
  directory-materialization and conditional script-copy branch inside
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:311)
  into
  [materializeSetupDirectoriesStructure()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:199)
  while keeping `setupDirectories()` as the public wrapper and primary work
  surface.
- The focused gate reran green after the live seam through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
  with `51` passing expectations.
- Post-seam classification of the live target stayed `review` /
  `direct-extraction-ready`; the file now reports `84` top-level functions
  with a maximum helper span of `259` lines.
- The next bounded checkpoint should stage one exact-source helper wave for
  [materializeSetupDirectoriesStructure()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:199)
  into
  [R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_setupdirectories_helpers.R:1)
  while keeping path-list construction, return assembly, and print flow live
  inside
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:389)
  and not widening into results IO or Rmd helpers.
- April 18, 2026 the sixth bounded `setupDirectories()` helper wave was
  authored in
  [tools/refactor/manifest-general-filemgmt-setupdirectories-wave6.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-setupdirectories-wave6.yml:1)
  and staged cleanly into
  [tools/refactor/staging/general-filemgmt-setupdirectories-wave6/R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-setupdirectories-wave6/R/func_general_filemgmt_setupdirectories_helpers.R:1)
  for:
  - `materializeSetupDirectoriesStructure()`
- The staged collate preview is recorded in
  [tools/refactor/staging/general-filemgmt-setupdirectories-wave6/collate-general-filemgmt-setupdirectories-wave6.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-setupdirectories-wave6/collate-general-filemgmt-setupdirectories-wave6.txt:1)
  with the intended order:
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`.
- The staged helper target parsed cleanly after extraction.
- The focused gate reran green after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
  with `51` passing expectations.
- Post-stage classification of the live target stayed `review` /
  `direct-extraction-ready`; the file still reports `84` top-level functions
  with a maximum helper span of `259` lines.
- This stop point keeps live `R/` sources unchanged for the new seam, so the
  next bounded checkpoint should review and, if approved, apply
  `manifest-general-filemgmt-setupdirectories-wave6.yml` into
  [R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_setupdirectories_helpers.R:1)
  while keeping path-list construction, return assembly, and print flow live
  inside
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:389)
  and not widening into results IO or Rmd helpers.
- April 18, 2026 the reviewed `setupDirectories()` wave6 helper slice was
  applied live via
  [tools/refactor/manifest-general-filemgmt-setupdirectories-wave6.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-setupdirectories-wave6.yml:1)
  into:
  - [R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_setupdirectories_helpers.R:302)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:227)
  for:
  - `materializeSetupDirectoriesStructure()`
- The live collate preview is recorded in
  [collate-general-filemgmt-setupdirectories-wave6.txt](/home/doktersmol/Documents/MultiScholaR/collate-general-filemgmt-setupdirectories-wave6.txt:1)
  with the applied order:
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`.
- The post-apply checker passed through
  `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-general-filemgmt-setupdirectories-wave6.yml`.
- The focused gate reran green after apply through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
  with `51` passing expectations.
- Post-apply classification of the live target stayed `review` /
  `direct-extraction-ready`; the file now reports `83` top-level functions
  with a maximum helper span of `259` lines.
- The next bounded checkpoint should stage one exact-source helper wave for
  the final print/report block inside
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:227)
  while keeping `all_created_paths` return assembly and per-omic
  orchestration live in `R/func_general_filemgmt.R` and not widening into
  results IO or Rmd helpers.
- April 18, 2026 one bounded live seam split the final print/report block
  inside
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:284)
  into
  [printSetupDirectoriesSummary()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:200)
  while keeping `setupDirectories()` as the public wrapper and primary work
  surface.
- The focused gate reran green after the live seam through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
  with `51` passing expectations.
- Post-seam classification of the live target stayed `review` /
  `direct-extraction-ready`; the file now reports `84` top-level functions
  with a maximum helper span of `259` lines.
- This seam makes the final print/report surface symbol-addressable, so the
  next bounded checkpoint should stage one exact-source helper wave for
  [printSetupDirectoriesSummary()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:200)
  into
  [R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_setupdirectories_helpers.R:1)
  while keeping `all_created_paths` return assembly and per-omic
  orchestration live in
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:284)
  and not widening into results IO or Rmd helpers.
- April 18, 2026 the seventh bounded `setupDirectories()` helper wave was
  staged from
  [tools/refactor/manifest-general-filemgmt-setupdirectories-wave7.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-setupdirectories-wave7.yml:1)
  into
  [tools/refactor/staging/general-filemgmt-setupdirectories-wave7/R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-setupdirectories-wave7/R/func_general_filemgmt_setupdirectories_helpers.R:1)
  for:
  - `printSetupDirectoriesSummary()`
- The staged collate preview is recorded in
  [tools/refactor/staging/general-filemgmt-setupdirectories-wave7/collate-general-filemgmt-setupdirectories-wave7.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-setupdirectories-wave7/collate-general-filemgmt-setupdirectories-wave7.txt:1)
  with the intended order:
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`.
- The staged helper target parsed cleanly after extraction.
- The focused gate reran green again after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
  with `51` passing expectations.
- April 18, 2026 the reviewed `setupDirectories()` wave7 helper slice was
  applied live via
  [tools/refactor/manifest-general-filemgmt-setupdirectories-wave7.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-setupdirectories-wave7.yml:1)
  into:
  - [R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_setupdirectories_helpers.R:387)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:228)
  for:
  - `printSetupDirectoriesSummary()`
- The live apply moved `printSetupDirectoriesSummary()` out of
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:228)
  while keeping `setupDirectories()` as the public wrapper and primary work
  surface.
- The emitted collate artifact now exists at
  [tools/refactor/collate-general-filemgmt-setupdirectories-wave7.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-filemgmt-setupdirectories-wave7.txt:1).
- The post-apply checker passed through
  `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-general-filemgmt-setupdirectories-wave7.yml`.
- The focused gate reran green after apply through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
  with `51` passing expectations.
- Post-apply classification of the live target stayed `review` /
  `direct-extraction-ready`; the file now reports `83` top-level functions
  with a maximum helper span of `259` lines.
- The next bounded checkpoint should introduce one live seam for the
  per-omic path-definition and timestamped-directory initialization block
  inside
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:228),
  extracting only the `publication_graphs_dir_base` / `qc_dir_base` /
  `current_omic_paths_def` assembly and timestamped `dir.create()` call into
  a new top-level helper while keeping `materializeSetupDirectoriesStructure()`,
  `buildSetupDirectoriesPathList()`, return assembly, and print delegation
  live in `setupDirectories()` and not widening into results IO or Rmd
  helpers.
- April 18, 2026 one bounded live seam split the per-omic path-definition and
  timestamped-directory initialization block inside
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:269)
  into
  [initializeSetupDirectoriesOmicPaths()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:201)
  while keeping `setupDirectories()` as the public wrapper and primary work
  surface.
- The focused gate reran green after the live seam through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
  with `51` passing expectations.
- Post-seam classification of the live target stayed `review` /
  `direct-extraction-ready`; the file now reports `84` top-level functions
  with a maximum helper span of `259` lines.
- This seam makes the per-omic path-definition block symbol-addressable, so
  the next bounded checkpoint should stage one exact-source helper wave for
  [initializeSetupDirectoriesOmicPaths()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:201)
  into
  [R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_setupdirectories_helpers.R:1)
  while keeping `materializeSetupDirectoriesStructure()`,
  `buildSetupDirectoriesPathList()`, return assembly, and print delegation
  live in
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:269)
  and not widening into results IO or Rmd helpers.
- April 18, 2026 the eighth bounded `setupDirectories()` helper wave was
  staged from
  [tools/refactor/manifest-general-filemgmt-setupdirectories-wave8.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-setupdirectories-wave8.yml:1)
  into
  [tools/refactor/staging/general-filemgmt-setupdirectories-wave8/R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-setupdirectories-wave8/R/func_general_filemgmt_setupdirectories_helpers.R:1)
  for:
  - `initializeSetupDirectoriesOmicPaths()`
- The staged collate preview is recorded in
  [tools/refactor/staging/general-filemgmt-setupdirectories-wave8/collate-general-filemgmt-setupdirectories-wave8.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-setupdirectories-wave8/collate-general-filemgmt-setupdirectories-wave8.txt:1)
  with the intended order:
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`.
- The staged helper slice parsed cleanly after extraction through
  `Rscript -e 'parse(file = "tools/refactor/staging/general-filemgmt-setupdirectories-wave8/R/func_general_filemgmt_setupdirectories_helpers.R")'`.
- The focused gate reran green after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
  with `51` passing expectations.
- Live classification of the target remained `review` /
  `direct-extraction-ready`; [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  still reports `84` top-level functions with a `259`-line maximum helper
  span.
- The next bounded checkpoint should review and apply wave8 via
  [tools/refactor/manifest-general-filemgmt-setupdirectories-wave8.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-setupdirectories-wave8.yml:1)
  into:
  - [R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_setupdirectories_helpers.R:1)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:201)
  for `initializeSetupDirectoriesOmicPaths()` while keeping
  `materializeSetupDirectoriesStructure()`,
  `buildSetupDirectoriesPathList()`, return assembly, and print delegation
  live in
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:269)
  and not widening into results IO or Rmd helpers.
- April 18, 2026 the reviewed `setupDirectories()` wave8 helper slice was
  applied live via
  [tools/refactor/manifest-general-filemgmt-setupdirectories-wave8.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-setupdirectories-wave8.yml:1)
  into:
  - [R/func_general_filemgmt_setupdirectories_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_setupdirectories_helpers.R:442)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:266)
  for:
  - `initializeSetupDirectoriesOmicPaths()`
- The live apply moved `initializeSetupDirectoriesOmicPaths()` out of
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:266)
  while keeping `setupDirectories()` as the public wrapper and primary work
  surface.
- The emitted collate artifact now exists at
  [tools/refactor/collate-general-filemgmt-setupdirectories-wave8.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-filemgmt-setupdirectories-wave8.txt:1).
- The post-apply checker passed through
  `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-general-filemgmt-setupdirectories-wave8.yml`.
- The focused gate reran green after apply through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-path-contracts.R`
  with `51` passing expectations.
- Post-apply classification of the live target stayed `review` /
  `direct-extraction-ready`; the file now reports `83` top-level functions
  with a maximum helper span of `259` lines.
- The next bounded checkpoint should add one focused characterization gate for
  [readConfigFile()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:564)
  and
  [readConfigFileSection()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:832)
  before drafting any config-helper extraction wave, because the current
  `test-general-filemgmt-path-contracts.R` surface only freezes the
  `setupDirectories()` / path-contract behavior.
- April 18, 2026 one bounded config-reader characterization checkpoint added
  [tests/testthat/test-general-filemgmt-config-contracts.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-general-filemgmt-config-contracts.R:1)
  to freeze the `readConfigFile()` /
  `readConfigFileSection()` surface before any config-helper staging work.
- The focused config gate exercises:
  - `readConfigFile()`
  - `readConfigFileSection()`
- The live source now matches the documented call contract by accepting the
  `file_type` parameter in
  [readConfigFile()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:564),
  which keeps
  [readConfigFileSection()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:832)
  runnable without widening the checkpoint into extraction work.
- The focused config gate passed through:
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-config-contracts.R`
  with `35` passing expectations.
- The next bounded checkpoint should stage one config-helper wave for
  [readConfigFile()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:564)
  and
  [readConfigFileSection()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:832)
  into `R/func_general_filemgmt_config_helpers.R` while keeping config
  round-trip writers, results IO, and Rmd sourcing helpers live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 the first bounded config-helper wave was staged from
  [tools/refactor/manifest-general-filemgmt-config-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-config-wave1.yml:1)
  into
  [tools/refactor/staging/general-filemgmt-config-wave1/R/func_general_filemgmt_config_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-config-wave1/R/func_general_filemgmt_config_helpers.R:1)
  for:
  - `readConfigFile()`
  - `readConfigFileSection()`
- The staged collate preview is recorded in
  [tools/refactor/staging/general-filemgmt-config-wave1/collate-general-filemgmt-config-wave1.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-config-wave1/collate-general-filemgmt-config-wave1.txt:1)
  with the intended order:
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`,
  `R/func_general_filemgmt_config_helpers.R`.
- The staged helper target parsed cleanly after extraction through
  `Rscript -e 'parse(file = "tools/refactor/staging/general-filemgmt-config-wave1/R/func_general_filemgmt_config_helpers.R")'`.
- The focused config gate reran green after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-config-contracts.R`
  with `35` passing expectations.
- Live classification of the target stayed `review` /
  `direct-extraction-ready`; [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  still reports `83` top-level functions with a `259`-line maximum helper
  span.
- This stop point keeps live `R/` sources unchanged for the new seam, so the
  next bounded checkpoint should review and, if approved, apply
  [tools/refactor/manifest-general-filemgmt-config-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-config-wave1.yml:1)
  into:
  - [R/func_general_filemgmt_config_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_config_helpers.R:1)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:564)
  for `readConfigFile()` and `readConfigFileSection()` while keeping config
  round-trip writers, results IO, and Rmd sourcing helpers live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 the reviewed config wave1 helper slice was applied live via
  [tools/refactor/manifest-general-filemgmt-config-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-config-wave1.yml:1)
  into:
  - [R/func_general_filemgmt_config_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_config_helpers.R:1)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:556)
  for:
  - `readConfigFile()`
  - `readConfigFileSection()`
- The live apply moved the config readers out of
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  while keeping config round-trip writers, results IO, and Rmd sourcing live
  in the primary wrapper file.
- The emitted collate artifact now exists at
  [tools/refactor/collate-general-filemgmt-config-wave1.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-filemgmt-config-wave1.txt:1)
  and `DESCRIPTION` now collates:
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`,
  `R/func_general_filemgmt_config_helpers.R`,
  and then `R/func_general_filemgmt.R`.
- The post-apply checker passed through
  `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-general-filemgmt-config-wave1.yml`.
- The focused config gate reran green after apply through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-config-contracts.R`
  with `35` passing expectations.
- Post-apply classification of the live target is now `direct-extraction-ready`;
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  reports `4983` lines, `81` top-level functions, and a `210`-line maximum
  helper span.
- The next bounded checkpoint should add one config round-trip
  characterization gate for:
  - [formatConfigList()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:778)
  - [updateConfigParameter()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:2093)
  - [createStudyParametersFile()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:2553)
  - [createWorkflowArgsFromConfig()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:3203)
  before drafting any second config-helper wave.
- April 18, 2026 one bounded config round-trip characterization checkpoint
  added
  [tests/testthat/test-general-filemgmt-config-roundtrip-contracts.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-general-filemgmt-config-roundtrip-contracts.R:1)
  to freeze the `formatConfigList()` / `updateConfigParameter()` /
  `createStudyParametersFile()` / `createWorkflowArgsFromConfig()` surface
  before staging any second config-helper wave.
- The focused round-trip gate exercises:
  - `formatConfigList()`
  - `updateConfigParameter()`
  - `createStudyParametersFile()`
  - `createWorkflowArgsFromConfig()`
  - env-to-`study_parameters.txt` propagation for updated config values.
- The focused round-trip gate passed through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-config-roundtrip-contracts.R`
  with `24` passing expectations.
- The existing config-reader gate reran green again through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-config-contracts.R`
  with `35` passing expectations.
- Live `R/` sources stayed unchanged during this checkpoint; post-checkpoint
  classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  still reports `4983` lines, `81` top-level functions, and a `210`-line
  maximum helper span.
- The next bounded checkpoint should stage one second config-helper wave for:
  - [formatConfigList()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:778)
  - [updateConfigParameter()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:2093)
  - [createStudyParametersFile()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:2553)
  - [createWorkflowArgsFromConfig()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:3203)
  into
  [R/func_general_filemgmt_config_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_config_helpers.R:1)
  while keeping results IO, report/project export, and Rmd sourcing helpers
  live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 the second bounded config-helper wave was authored in
  [tools/refactor/manifest-general-filemgmt-config-wave2.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-config-wave2.yml:1)
  and staged cleanly into
  [tools/refactor/staging/general-filemgmt-config-wave2/R/func_general_filemgmt_config_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-config-wave2/R/func_general_filemgmt_config_helpers.R:1)
  for:
  - `formatConfigList()`
  - `updateConfigParameter()`
  - `createStudyParametersFile()`
  - `createWorkflowArgsFromConfig()`
- The staged collate preview is recorded in
  [tools/refactor/staging/general-filemgmt-config-wave2/collate-general-filemgmt-config-wave2.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-config-wave2/collate-general-filemgmt-config-wave2.txt:1)
  with the intended order:
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`,
  `R/func_general_filemgmt_config_helpers.R`.
- The staged helper target parsed cleanly after extraction.
- The focused round-trip gate reran green after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-config-roundtrip-contracts.R`
  with `24` passing expectations.
- The existing config-reader gate reran green again through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-config-contracts.R`
  with `35` passing expectations.
- Live `R/` sources stayed unchanged during staging; post-stage classification
  remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  still reports `4983` lines, `81` top-level functions, and a `210`-line
  maximum helper span.
- April 18, 2026 the approved second bounded config-helper wave was applied
  live from
  [tools/refactor/manifest-general-filemgmt-config-wave2.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-config-wave2.yml:1)
  into:
  - [R/func_general_filemgmt_config_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_config_helpers.R:1)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
- The live apply removed these top-level helpers from
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  and materialized them in
  [R/func_general_filemgmt_config_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_config_helpers.R:300):
  - `formatConfigList()`
  - `updateConfigParameter()`
  - `createStudyParametersFile()`
  - `createWorkflowArgsFromConfig()`
- The live collate artifact is recorded in
  [tools/refactor/collate-general-filemgmt-config-wave2.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-filemgmt-config-wave2.txt:1),
  and `DESCRIPTION` now collates
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`,
  `R/func_general_filemgmt_config_helpers.R`,
  and then `R/func_general_filemgmt.R`.
- The focused round-trip gate reran green after the live apply through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-config-roundtrip-contracts.R`
  with `24` passing expectations.
- The existing config-reader gate reran green again through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-config-contracts.R`
  with `35` passing expectations.
- Post-apply classification stayed `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  now reports `3258` lines, `55` top-level functions, and a `210`-line
  maximum helper span.
- Existing summary-module coverage already exercises
  [copyToResultsSummary()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:795)
  and
  [RenderReport()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1725)
  through
  [tests/testthat/test-prot-12-summary-module-contracts.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-prot-12-summary-module-contracts.R:1),
  so the next stop point can stay in staged-wave mode.
- The next bounded checkpoint should draft and verify one report-helper staging
  wave for:
  - [copyToResultsSummary()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:795)
  - [downloadReportTemplate()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1670)
  - [RenderReport()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1725)
  into `R/func_general_filemgmt_report_helpers.R` while keeping
  [sourceRmdFileSimple()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:339),
  [sourceRmdFile()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:355),
  [saveListOfPdfs()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:310),
  [savePlot()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:480),
  [save_plot()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:543),
  and
  [write_results()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:551)
  live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 the first bounded report-helper wave was authored in
  [tools/refactor/manifest-general-filemgmt-report-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-report-wave1.yml:1)
  and staged cleanly into
  [tools/refactor/staging/general-filemgmt-report-wave1/R/func_general_filemgmt_report_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-report-wave1/R/func_general_filemgmt_report_helpers.R:1)
  for:
  - `copyToResultsSummary()`
  - `downloadReportTemplate()`
  - `RenderReport()`
- The staged collate preview is recorded in
  [tools/refactor/staging/general-filemgmt-report-wave1/collate-general-filemgmt-report-wave1.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-report-wave1/collate-general-filemgmt-report-wave1.txt:1)
  with the intended order:
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`,
  `R/func_general_filemgmt_config_helpers.R`,
  `R/func_general_filemgmt_report_helpers.R`.
- The staged helper target parsed cleanly after extraction.
- The focused summary-module gate reran green after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-12-summary-module-contracts.R`
  with `436` passing expectations.
- Live `R/` sources stayed unchanged during staging; post-stage classification
  remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  still reports `3258` lines, `55` top-level functions, and a `210`-line
  maximum helper span.
- The next bounded checkpoint should review and, if approved, apply
  [tools/refactor/manifest-general-filemgmt-report-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-report-wave1.yml:1)
  into:
  - [R/func_general_filemgmt_report_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_report_helpers.R:1)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:795)
  for `copyToResultsSummary()`, `downloadReportTemplate()`, and
  `RenderReport()` while keeping
  [sourceRmdFileSimple()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:339),
  [sourceRmdFile()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:355),
  [saveListOfPdfs()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:310),
  [savePlot()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:480),
  [save_plot()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:543),
  and
  [write_results()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:551)
  live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 the reviewed report-helper wave1 slice was applied live from
  [tools/refactor/manifest-general-filemgmt-report-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-report-wave1.yml:1)
  into:
  - [R/func_general_filemgmt_report_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_report_helpers.R:1)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
- The live apply materialized the first bounded `report handling` helper slice
  out of `R/func_general_filemgmt.R` for:
  - [copyToResultsSummary()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_report_helpers.R:16)
  - [downloadReportTemplate()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_report_helpers.R:731)
  - [RenderReport()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_report_helpers.R:786)
- The emitted collate artifact now exists at
  [tools/refactor/collate-general-filemgmt-report-wave1.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-filemgmt-report-wave1.txt:1),
  and `DESCRIPTION` now collates
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`,
  `R/func_general_filemgmt_config_helpers.R`,
  `R/func_general_filemgmt_report_helpers.R`,
  and then `R/func_general_filemgmt.R`.
- The focused summary-module gate reran green after the live apply through
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-12-summary-module-contracts.R`
  with `436` passing expectations.
- Post-apply classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  now reports `2285` lines, `39` top-level functions, and a `210`-line
  maximum helper span.
- The next bounded checkpoint should draft and verify one follow-up
  report-helper staging wave for:
  - [saveListOfPdfs()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:310)
  - [sourceRmdFileSimple()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:339)
  - [sourceRmdFile()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:355)
  into `R/func_general_filemgmt_report_helpers.R` while keeping
  [savePlot()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:480),
  [save_plot()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:543),
  and
  [write_results()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:551)
  live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 the follow-up report-helper wave is now recorded in
  [tools/refactor/manifest-general-filemgmt-report-wave2.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-report-wave2.yml:1)
  and staged cleanly into
  [tools/refactor/staging/general-filemgmt-report-wave2/R/func_general_filemgmt_report_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-report-wave2/R/func_general_filemgmt_report_helpers.R:1)
  for:
  - `saveListOfPdfs()`
  - `sourceRmdFileSimple()`
  - `sourceRmdFile()`
- The staged collate preview is recorded in
  [tools/refactor/staging/general-filemgmt-report-wave2/collate-general-filemgmt-report-wave2.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-report-wave2/collate-general-filemgmt-report-wave2.txt:1)
  with the intended order:
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`,
  `R/func_general_filemgmt_config_helpers.R`,
  `R/func_general_filemgmt_report_helpers.R`.
- The staged helper target parsed cleanly after extraction, and the focused
  summary-module gate reran green after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-12-summary-module-contracts.R`
  with `436` passing expectations.
- Live `R/` sources stayed unchanged during staging; post-stage classification
  remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  still reports `2285` lines, `39` top-level functions, and a `210`-line
  maximum helper span.
- The next bounded stop point is a reviewed apply of wave2 via
  [tools/refactor/manifest-general-filemgmt-report-wave2.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-report-wave2.yml:1)
  into:
  - [R/func_general_filemgmt_report_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_report_helpers.R:1)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  for `saveListOfPdfs()`, `sourceRmdFileSimple()`, and `sourceRmdFile()`
  while keeping
  [savePlot()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:480),
  [save_plot()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:543),
  and
  [write_results()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:551)
  live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 the reviewed report-helper wave2 slice was applied live via
  [tools/refactor/manifest-general-filemgmt-report-wave2.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-report-wave2.yml:1)
  into:
  - [R/func_general_filemgmt_report_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_report_helpers.R:1)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  for:
  - `saveListOfPdfs()`
  - `sourceRmdFileSimple()`
  - `sourceRmdFile()`
- The live apply materialized the second bounded `report handling` helper
  slice out of `R/func_general_filemgmt.R` for:
  - [saveListOfPdfs()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_report_helpers.R:981)
  - [sourceRmdFileSimple()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_report_helpers.R:1003)
  - [sourceRmdFile()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_report_helpers.R:1018)
- The emitted collate artifact now exists at
  [tools/refactor/collate-general-filemgmt-report-wave2.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-filemgmt-report-wave2.txt:1),
  and `DESCRIPTION` `Collate:` remains:
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`,
  `R/func_general_filemgmt_config_helpers.R`,
  `R/func_general_filemgmt_report_helpers.R`,
  and then `R/func_general_filemgmt.R`.
- The post-apply checker passed through:
  `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-general-filemgmt-report-wave2.yml`.
- The focused summary-module gate reran green after the live apply through
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-12-summary-module-contracts.R`
  with `436` passing expectations.
- Post-apply classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  now reports `2231` lines, `36` top-level functions, and a `210`-line
  maximum helper span.
- The next bounded stop point is one export-helper characterization checkpoint
  for:
  - [savePlot()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:426)
  - [save_plot()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:489)
  - [write_results()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:497)
  before drafting any third report/results-IO wave, because the current
  focused summary-module gate does not directly freeze that file-writer
  surface.
- Keep
  [testRequiredFiles()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:324),
  `testRequiredFilesWarning()`, `testRequiredArguments()`, `parseType()`,
  `parseString()`, `parseNumeric()`, `parseCharacter()`, `parseLogical()`,
  `parseToList()`, `isArgumentDefined()`, and
  [loadDependencies()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:536)
  live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 one bounded export-helper characterization checkpoint added
  [tests/testthat/test-general-filemgmt-export-contracts.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-general-filemgmt-export-contracts.R:1)
  to freeze:
  - [savePlot()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:426)
  - [save_plot()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:489)
  - [write_results()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:497)
- The characterization gate reran green through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-export-contracts.R`
  with `26` passing expectations.
- The existing summary-module gate reran green again through
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-12-summary-module-contracts.R`
  with `436` passing expectations.
- Live `R/` sources stayed unchanged during this checkpoint; post-checkpoint
  classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  still reports `2231` lines, `36` top-level functions, and a `210`-line
  maximum helper span.
- This stop point keeps
  [savePlot()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:426),
  [save_plot()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:489),
  and
  [write_results()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:497)
  live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1),
  so the next bounded checkpoint should stage one exact-source export-helper
  wave for those symbols into
  [R/func_general_filemgmt_report_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_report_helpers.R:1)
  without widening into validation helpers or `loadDependencies()`.
- April 18, 2026 one bounded export-helper staging wave is now recorded in
  [tools/refactor/manifest-general-filemgmt-report-wave3.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-report-wave3.yml:1)
  and staged cleanly into
  [tools/refactor/staging/general-filemgmt-report-wave3/R/func_general_filemgmt_report_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-report-wave3/R/func_general_filemgmt_report_helpers.R:1)
  for:
  - `savePlot()`
  - `save_plot()`
  - `write_results()`
- The staged collate preview is recorded in
  [tools/refactor/staging/general-filemgmt-report-wave3/collate-general-filemgmt-report-wave3.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-report-wave3/collate-general-filemgmt-report-wave3.txt:1)
  with the intended order:
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`,
  `R/func_general_filemgmt_config_helpers.R`,
  `R/func_general_filemgmt_report_helpers.R`.
- The staged export-helper slice parsed cleanly after extraction and now
  materializes:
  - [savePlot()](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-report-wave3/R/func_general_filemgmt_report_helpers.R:20)
  - [save_plot()](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-report-wave3/R/func_general_filemgmt_report_helpers.R:82)
  - [write_results()](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-report-wave3/R/func_general_filemgmt_report_helpers.R:89)
- The focused export-helper gate reran green after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-export-contracts.R`
  with `26` passing expectations.
- Live `R/` sources stayed unchanged during staging; post-stage
  classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  still reports `2231` lines, `36` top-level functions, and a `210`-line
  maximum helper span.
- The next bounded stop point is a reviewed apply of wave3 via
  [tools/refactor/manifest-general-filemgmt-report-wave3.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-report-wave3.yml:1)
  into:
  - [R/func_general_filemgmt_report_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_report_helpers.R:1)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  for `savePlot()`, `save_plot()`, and `write_results()` while keeping
  [testRequiredFiles()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:324),
  `testRequiredFilesWarning()`, `testRequiredArguments()`, `parseType()`,
  `parseString()`, `parseList()`, `isArgumentDefined()`, and
  [loadDependencies()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:536)
  live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 the reviewed report wave3 helper slice was applied live via
  [tools/refactor/manifest-general-filemgmt-report-wave3.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-report-wave3.yml:1)
  into:
  - [R/func_general_filemgmt_report_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_report_helpers.R:1)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  for:
  - [savePlot()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_report_helpers.R:1053)
  - [save_plot()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_report_helpers.R:1115)
  - [write_results()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_report_helpers.R:1122)
- The live apply moved the bounded export-helper/file-writer surface into
  [R/func_general_filemgmt_report_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_report_helpers.R:1)
  while keeping
  [testRequiredFiles()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:324),
  [testRequiredFilesWarning()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:337),
  [testRequiredArguments()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:349),
  [parseType()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:363),
  [parseString()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:375),
  [parseList()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:389),
  [isArgumentDefined()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:402),
  and
  [loadDependencies()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:447)
  live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- The emitted collate artifact now exists at
  [tools/refactor/collate-general-filemgmt-report-wave3.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-filemgmt-report-wave3.txt:1)
  with the applied order:
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`,
  `R/func_general_filemgmt_config_helpers.R`,
  `R/func_general_filemgmt_report_helpers.R`.
- The post-apply checker passed through:
  `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-general-filemgmt-report-wave3.yml`.
- The focused export-helper gate reran green after apply through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-export-contracts.R`
  with `26` passing expectations.
- Post-apply classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  now reports `2142` lines, `33` top-level functions, and a `210`-line
  maximum helper span.
- The next bounded stop point is one utility-helper characterization
  checkpoint for:
  - [testRequiredFiles()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:324)
  - [testRequiredFilesWarning()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:337)
  - [testRequiredArguments()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:349)
  - [parseType()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:363)
  - [parseString()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:375)
  - [parseList()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:389)
  - [isArgumentDefined()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:402)
  before drafting any shared-helper extraction wave into
  [R/func_general_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_helpers.R:1),
  while keeping
  [loadDependencies()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:447)
  and the downstream S4/NA-validation helpers live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 the bounded utility-helper characterization checkpoint added
  [tests/testthat/test-general-filemgmt-utility-contracts.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-general-filemgmt-utility-contracts.R:1)
  to freeze:
  - [testRequiredFiles()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:324)
  - [testRequiredFilesWarning()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:337)
  - [testRequiredArguments()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:349)
  - [parseType()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:363)
  - [parseString()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:375)
  - [parseList()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:389)
  - [isArgumentDefined()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:402)
- The focused utility-helper gate reran green through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-utility-contracts.R`
  with `15` passing expectations.
- The characterization intentionally freezes the current caller-visible
  behavior before any shared-helper extraction wave:
  - `testRequiredFiles()` and `testRequiredArguments()` log once per missing
    entry and call `q()` once per missing entry.
  - `testRequiredFilesWarning()` logs once per missing file without quitting.
  - `parseType()`, `parseString()`, and `parseList()` currently return the
    original list without caller-visible mutation.
  - `isArgumentDefined()` currently treats a named `NULL` entry as defined.
- Live `R/` sources stayed unchanged during characterization; post-checkpoint
  classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  still reports `2142` lines, `33` top-level functions, and a `210`-line
  maximum helper span.
- The next bounded stop point is to draft and stage one shared-helper wave
  into [R/func_general_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_helpers.R:1)
  for:
  - [testRequiredFiles()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:324)
  - [testRequiredFilesWarning()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:337)
  - [testRequiredArguments()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:349)
  - [parseType()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:363)
  - [parseString()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:375)
  - [parseList()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:389)
  - [isArgumentDefined()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:402)
  while keeping
  [loadDependencies()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:447)
  and the downstream S4/NA-validation helpers live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 the first bounded shared-helper wave was authored in
  [tools/refactor/manifest-general-filemgmt-utility-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-utility-wave1.yml:1)
  and staged cleanly into
  [tools/refactor/staging/general-filemgmt-utility-wave1/R/func_general_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-utility-wave1/R/func_general_helpers.R:1)
  for:
  - [testRequiredFiles()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:324)
  - [testRequiredFilesWarning()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:337)
  - [testRequiredArguments()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:349)
  - [parseType()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:363)
  - [parseString()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:375)
  - [parseList()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:389)
  - [isArgumentDefined()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:402)
- The staged collate preview is recorded in
  [tools/refactor/staging/general-filemgmt-utility-wave1/collate-general-filemgmt-utility-wave1.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-utility-wave1/collate-general-filemgmt-utility-wave1.txt:1)
  with the targeted order:
  `R/func_general_helpers.R`.
- The focused utility-helper gate reran green after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-utility-contracts.R`
  with `15` passing expectations.
- Live `R/` sources stayed unchanged during staging; post-stage
  classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  still reports `2142` lines, `33` top-level functions, and a `210`-line
  maximum helper span.
- The next bounded stop point is a reviewed apply of
  [tools/refactor/manifest-general-filemgmt-utility-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-utility-wave1.yml:1)
  into:
  - [R/func_general_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_helpers.R:1)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  for the staged validation and parser helpers while keeping
  [loadDependencies()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:447)
  and the downstream S4/NA-validation helpers live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 the reviewed apply attempt for
  [tools/refactor/manifest-general-filemgmt-utility-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-utility-wave1.yml:1)
  was blocked cleanly because
  `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-general-filemgmt-utility-wave1.yml`
  reported duplicate symbols already shared between
  [R/func_general_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_helpers.R:1003)
  and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:599):
  - `extract_experiment`
  - `updateMissingValueParameters`
  - `chooseBestProteinAccession_s3`
- The blocked apply ran through the skill wrapper
  `python3 .../god-module-stabilization/scripts/apply_wave.py`
  and restored repo files from backup, so live `R/` sources stayed unchanged.
- April 18, 2026 a replacement dedicated utility-helper wave was authored in
  [tools/refactor/manifest-general-filemgmt-utility-wave2.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-utility-wave2.yml:1)
  and staged cleanly into
  [tools/refactor/staging/general-filemgmt-utility-wave2/R/func_general_filemgmt_utility_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-utility-wave2/R/func_general_filemgmt_utility_helpers.R:1)
  for:
  - [testRequiredFiles()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:324)
  - [testRequiredFilesWarning()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:337)
  - [testRequiredArguments()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:349)
  - [parseType()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:363)
  - [parseString()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:375)
  - [parseList()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:389)
  - [isArgumentDefined()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:402)
- The staged collate preview is recorded in
  [tools/refactor/staging/general-filemgmt-utility-wave2/collate-general-filemgmt-utility-wave2.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-utility-wave2/collate-general-filemgmt-utility-wave2.txt:1)
  with the targeted order:
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`,
  `R/func_general_filemgmt_config_helpers.R`,
  `R/func_general_filemgmt_report_helpers.R`,
  `R/func_general_filemgmt_utility_helpers.R`.
- The focused utility-helper gate reran green after restaging through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-utility-contracts.R`
  with `15` passing expectations.
- Live `R/` sources stayed unchanged during the replacement staging; post-stage
  classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  still reports `2142` lines, `33` top-level functions, and a `210`-line
  maximum helper span.
- April 18, 2026 the reviewed utility wave2 helper slice was applied live via
  [tools/refactor/manifest-general-filemgmt-utility-wave2.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-utility-wave2.yml:1)
  into:
  - [R/func_general_filemgmt_utility_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_utility_helpers.R:1)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  for:
  - `testRequiredFiles()`
  - `testRequiredFilesWarning()`
  - `testRequiredArguments()`
  - `parseType()`
  - `parseString()`
  - `parseList()`
  - `isArgumentDefined()`
- The live apply moved the bounded validation and parser helper surface into
  [R/func_general_filemgmt_utility_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_utility_helpers.R:1)
  while keeping
  [loadDependencies()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:374),
  [extract_experiment()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:526),
  `updateMissingValueParameters()`, `chooseBestProteinAccession_s3()`, and the
  downstream S4/NA-validation helpers live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- The emitted collate artifact now exists at
  [tools/refactor/collate-general-filemgmt-utility-wave2.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-filemgmt-utility-wave2.txt:1)
  with the applied order:
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`,
  `R/func_general_filemgmt_config_helpers.R`,
  `R/func_general_filemgmt_report_helpers.R`,
  `R/func_general_filemgmt_utility_helpers.R`.
- `DESCRIPTION` now collates
  `R/func_general_filemgmt_utility_helpers.R` immediately before
  `R/func_general_filemgmt.R`.
- The post-apply checker passed through
  `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-general-filemgmt-utility-wave2.yml`.
- The focused utility-helper gate reran green after apply through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-utility-contracts.R`
  with `15` passing expectations.
- Post-apply classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  now reports `2069` lines, `26` top-level functions, and a `210`-line
  maximum helper span.
- The next bounded stop point is one duplicate-resolution characterization
  checkpoint for
  [extract_experiment()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:526)
  before any live dedupe between
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  and
  [R/func_general_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_helpers.R:1003),
  while keeping
  [loadDependencies()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:374),
  `updateMissingValueParameters()`, `chooseBestProteinAccession_s3()`, and the
  downstream S4/NA-validation helpers live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 one bounded duplicate-resolution characterization checkpoint
  added
  [tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R:1)
  to freeze the current `extract_experiment()` duplicate surface before any
  live dedupe between
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:526)
  and
  [R/func_general_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_helpers.R:1034).
- The focused duplicate-resolution gate exercises:
  - `DESCRIPTION` collate ordering for `func_general_filemgmt.R` before
    `func_general_helpers.R`
  - structural parity between the two live `extract_experiment()` definitions
  - package-level behavior parity for range/start/end/error/warning paths
- The focused duplicate-resolution gate passed through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R`
  with `15` passing expectations.
- Post-checkpoint classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  still reports `2069` lines, `26` top-level functions, and a `210`-line
  maximum helper span.
- The next bounded stop point is one live duplicate-resolution seam for
  [extract_experiment()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:526),
  removing the duplicate implementation from
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  in favor of the later-loading canonical helper in
  [R/func_general_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_helpers.R:1034),
  while keeping
  [loadDependencies()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:374),
  `updateMissingValueParameters()`, `chooseBestProteinAccession_s3()`, and the
  downstream S4/NA-validation helpers live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 one bounded live duplicate-resolution seam removed the
  duplicate `extract_experiment()` implementation from
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1),
  leaving the canonical package-level helper in
  [R/func_general_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_helpers.R:1034).
- The focused duplicate-resolution gate now exercises:
  - `DESCRIPTION` collate ordering that keeps `func_general_filemgmt.R` ahead
    of `func_general_helpers.R`
  - absence of a live top-level `extract_experiment()` definition in
    `R/func_general_filemgmt.R`
  - package-level range/start/end/error/warning parity against the canonical
    helper in `R/func_general_helpers.R`
- The focused duplicate-resolution gate reran green after the live seam through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R`
  with `12` passing expectations.
- Post-seam classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  now reports `2013` lines, `25` top-level functions, and a `210`-line
  maximum helper span.
- The next bounded stop point is one duplicate-resolution characterization
  checkpoint for
  [updateMissingValueParameters()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:544)
  before any live dedupe between
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  and
  [R/func_general_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_helpers.R:1236),
  while keeping
  [loadDependencies()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:374),
  [chooseBestProteinAccession_s3()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:716),
  and the downstream S4/NA-validation helpers live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 one bounded duplicate-resolution characterization checkpoint
  extended
  [tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R:1)
  to freeze the current `updateMissingValueParameters()` duplicate surface
  before any live dedupe between
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:544)
  and
  [R/func_general_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_helpers.R:1236).
- The focused duplicate-resolution gate now also exercises:
  - both live top-level `updateMissingValueParameters()` definitions
  - canonical helper signature coverage for `function_name` and
    `grouping_variable`
  - package-level parity with the canonical helper for the default
    `removeRowsWithMissingValuesPercent` path and the helper-only
    `peptideIntensityFiltering` path
- The focused duplicate-resolution gate reran green after the checkpoint
  through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R`
  with `26` passing expectations.
- Post-checkpoint classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  still reports `2013` lines, `25` top-level functions, and a `210`-line
  maximum helper span.
- April 18, 2026 one bounded live duplicate-resolution seam removed the
  duplicate `updateMissingValueParameters()` implementation from
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1),
  leaving the canonical package-level helper in
  [R/func_general_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_helpers.R:1236).
- The focused duplicate-resolution gate now exercises:
  - `DESCRIPTION` collate ordering that keeps `func_general_filemgmt.R` ahead
    of `func_general_helpers.R`
  - absence of live top-level `extract_experiment()` and
    `updateMissingValueParameters()` definitions in
    `R/func_general_filemgmt.R`
  - package-level `extract_experiment()` parity against the canonical helper
    in `R/func_general_helpers.R`
  - package-level `updateMissingValueParameters()` signature and
    default/helper-only parity against the canonical helper in
    `R/func_general_helpers.R`
- The focused duplicate-resolution gate reran green after the live seam through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R`
  with `24` passing expectations.
- Post-seam classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  now reports `1868` lines, `24` top-level functions, and a `210`-line
  maximum helper span.
- The next bounded stop point is one duplicate-resolution characterization
  checkpoint for
  [chooseBestProteinAccession_s3()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:571)
  before any live dedupe between
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  and
  [R/func_general_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_helpers.R:1397),
  while keeping
  [loadDependencies()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:374),
  [updateRuvParameters()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:512),
  and the downstream S4/NA-validation helpers live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 one bounded duplicate-resolution characterization checkpoint
  extended
  [tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R:1)
  to freeze the current `updateRuvParameters()` duplicate surface before any
  live dedupe between
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:512)
  and
  [R/func_prot_norm_optimization_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_norm_optimization_helpers.R:289).
- The focused duplicate-resolution gate now also exercises:
  - `DESCRIPTION` collate ordering that keeps `func_general_filemgmt.R` ahead
    of `func_prot_norm_optimization_helpers.R`
  - both live top-level `updateRuvParameters()` definitions
  - package-level parity with both duplicates for the current populated-control
    and empty-control mutation paths
- The focused duplicate-resolution gate reran green after the checkpoint
  through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R`
  with `54` passing expectations.
- Post-checkpoint classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  still reports `1494` lines, `10` top-level functions, and a `210`-line
  maximum helper span.
- April 18, 2026 one bounded live duplicate-resolution seam removed the
  duplicate `updateRuvParameters()` implementation from
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1),
  leaving the canonical package-level helper in
  [R/func_prot_norm_optimization_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_norm_optimization_helpers.R:289).
- The focused duplicate-resolution gate now exercises:
  - `DESCRIPTION` collate ordering that keeps `func_general_filemgmt.R` ahead
    of `func_prot_norm_optimization_helpers.R`
  - absence of a live top-level `updateRuvParameters()` definition in
    `R/func_general_filemgmt.R`
  - package-level parity against the canonical normalization helper for the
    current populated-control and empty-control mutation paths
- The focused duplicate-resolution gate reran green after the live seam through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R`
  with `51` passing expectations.
- Post-seam classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  now reports `1477` lines, `9` top-level functions, and a `210`-line maximum
  helper span.
- April 18, 2026 one bounded S4-helper characterization checkpoint added
  [tests/testthat/test-general-filemgmt-s4-contracts.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-general-filemgmt-s4-contracts.R:1)
  to freeze the public `createWorkflowArgsFromConfig()` serialization surface
  for metabolomics and lipidomics S4 parameter extraction before any live
  apply of the S4-helper wave.
- The focused S4 gate now exercises:
  - [extractMetabS4Params()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:532)
    through multi-assay `MetaboliteAssayData` workflow-args output
  - [extractLipidS4Params()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:765)
    through `LipidomicsAssayData` workflow-args output
- The focused S4 gate reran green after the checkpoint through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-s4-contracts.R`
  with `30` passing expectations.
- April 18, 2026 the first bounded S4-helper wave was staged from
  [tools/refactor/manifest-general-filemgmt-s4-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-s4-wave1.yml:1)
  into
  [tools/refactor/staging/general-filemgmt-s4-wave1/R/func_general_filemgmt_s4_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-s4-wave1/R/func_general_filemgmt_s4_helpers.R:1)
  for:
  - `extractMetabS4Params()`
  - `extractLipidS4Params()`
- The staged collate preview is recorded in
  [tools/refactor/staging/general-filemgmt-s4-wave1/collate-general-filemgmt-s4-wave1.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-filemgmt-s4-wave1/collate-general-filemgmt-s4-wave1.txt:1)
  and places `R/func_general_filemgmt_s4_helpers.R` before
  `R/func_general_filemgmt_config_helpers.R`.
- Live classification of the target remained `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  still reports `1477` lines, `9` top-level functions, and a `210`-line
  maximum helper span because staging left live `R/` unchanged.
- The next bounded stop point is to review and, if approved, apply
  [tools/refactor/manifest-general-filemgmt-s4-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-s4-wave1.yml:1)
  into:
  - [R/func_general_filemgmt_s4_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_s4_helpers.R:1)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  while keeping
  [loadDependencies()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:374)
  and the downstream NA-validation helpers live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 the reviewed S4 wave1 slice was applied live from
  [tools/refactor/manifest-general-filemgmt-s4-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-filemgmt-s4-wave1.yml:1)
  into:
  - [R/func_general_filemgmt_s4_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_s4_helpers.R:1)
  - [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
- The live apply materialized the first bounded S4 helper slice out of
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  for:
  - [extractMetabS4Params()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_s4_helpers.R:19)
  - [extractLipidS4Params()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt_s4_helpers.R:248)
- The emitted collate artifact now exists at
  [tools/refactor/collate-general-filemgmt-s4-wave1.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-filemgmt-s4-wave1.txt:1),
  and `DESCRIPTION` now collates
  `R/func_general_filemgmt_path_helpers.R`,
  `R/func_general_filemgmt_directory_helpers.R`,
  `R/func_general_filemgmt_setupdirectories_helpers.R`,
  `R/func_general_filemgmt_s4_helpers.R`,
  `R/func_general_filemgmt_config_helpers.R`,
  `R/func_general_filemgmt_report_helpers.R`,
  `R/func_general_filemgmt_utility_helpers.R`,
  and then `R/func_general_filemgmt.R`.
- The focused S4 gate reran green after the live apply through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-s4-contracts.R`
  with `30` passing expectations.
- Post-apply classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  now reports `1021` lines, `7` top-level functions, and a `122`-line maximum
  helper span.
- The next bounded stop point is one duplicate-resolution characterization
  checkpoint for:
  - [checkPeptideNAPercentages()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:538)
  - [validatePostImputationData()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:665)
  - [getProteinNARecommendations()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:719)
  - [checkProteinNAPercentages()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:813)
  - [validatePostImputationProteinData()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:948)
  against the canonical package-level helpers in:
  - [R/func_prot_qc_peptide_support.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_qc_peptide_support.R:21)
  - [R/func_peptide_qc_imputation.R](/home/doktersmol/Documents/MultiScholaR/R/func_peptide_qc_imputation.R:143)
  - [R/func_prot_qc_reporting_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_qc_reporting_helpers.R:21)
  while keeping
  [loadDependencies()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:374)
  live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 one bounded duplicate-resolution characterization checkpoint
  extended
  [tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R:1)
  to freeze the remaining NA-validation duplicate surface before any live
  dedupe between
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:538)
  and the later-loading canonical helpers in:
  - [R/func_prot_qc_peptide_support.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_qc_peptide_support.R:21)
  - [R/func_peptide_qc_imputation.R](/home/doktersmol/Documents/MultiScholaR/R/func_peptide_qc_imputation.R:143)
  - [R/func_prot_qc_reporting_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_qc_reporting_helpers.R:21)
  for:
  - [checkPeptideNAPercentages()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:538)
  - [validatePostImputationData()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:665)
  - [getProteinNARecommendations()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:719)
  - [checkProteinNAPercentages()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:813)
  - [validatePostImputationProteinData()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:948)
- The focused duplicate-resolution gate now also exercises:
  - `DESCRIPTION` collate ordering that keeps `func_general_filemgmt.R` ahead
    of `func_prot_qc_peptide_support.R`,
    `func_peptide_qc_imputation.R`, and
    `func_prot_qc_reporting_helpers.R`
  - both live top-level NA-analysis, recommendation, and validation
    definitions in `R/func_general_filemgmt.R`
  - package-level parity with the canonical helpers for the current peptide
    analysis/validation path and protein analysis/recommendation/validation
    path
- The focused duplicate-resolution gate reran green after the checkpoint
  through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R`
  with `98` passing expectations.
- Post-checkpoint classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  still reports `1021` lines, `7` top-level functions, and a `122`-line
  maximum helper span.
- The next bounded stop point is one live duplicate-resolution seam removing
  the five remaining NA-validation duplicate implementations from
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  in favor of the later-loading canonical helpers in:
  - [R/func_prot_qc_peptide_support.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_qc_peptide_support.R:21)
  - [R/func_peptide_qc_imputation.R](/home/doktersmol/Documents/MultiScholaR/R/func_peptide_qc_imputation.R:143)
  - [R/func_prot_qc_reporting_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_qc_reporting_helpers.R:21)
  for:
  - [checkPeptideNAPercentages()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:538)
  - [validatePostImputationData()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:665)
  - [getProteinNARecommendations()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:719)
  - [checkProteinNAPercentages()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:813)
  - [validatePostImputationProteinData()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:948)
  while keeping
  [loadDependencies()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:374)
  live in
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1).
- April 18, 2026 one bounded live duplicate-resolution seam removed the five
  remaining NA-validation duplicate implementations from
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  in favor of the later-loading canonical helpers in:
  - [R/func_prot_qc_peptide_support.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_qc_peptide_support.R:21)
  - [R/func_peptide_qc_imputation.R](/home/doktersmol/Documents/MultiScholaR/R/func_peptide_qc_imputation.R:143)
  - [R/func_prot_qc_reporting_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_qc_reporting_helpers.R:21)
  for:
  - `checkPeptideNAPercentages()`
  - `validatePostImputationData()`
  - `getProteinNARecommendations()`
  - `checkProteinNAPercentages()`
  - `validatePostImputationProteinData()`
  while keeping
  [loadDependencies()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:374)
  live.
- The focused duplicate-resolution gate now exercises:
  - `DESCRIPTION` collate ordering that keeps `func_general_filemgmt.R` ahead
    of `func_prot_qc_peptide_support.R`,
    `func_peptide_qc_imputation.R`, and
    `func_prot_qc_reporting_helpers.R`.
  - Absence of live top-level `checkPeptideNAPercentages()`,
    `validatePostImputationData()`, `getProteinNARecommendations()`,
    `checkProteinNAPercentages()`, and
    `validatePostImputationProteinData()` definitions in
    `R/func_general_filemgmt.R`.
  - Package-level parity against the canonical peptide and protein NA helper
    implementations on the current analysis, recommendation, and validation
    paths.
- The focused duplicate-resolution gate reran green after the live seam
  through
  `Rscript tools/test_with_renv.R tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R`
  with `88` passing expectations.
- Post-seam classification remains `direct-extraction-ready`, and
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1)
  now reports `495` lines, `2` top-level functions, and a `122`-line maximum
  helper span.
- This completes the bounded stabilization target for
  [R/func_general_filemgmt.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:1);
  only
  [setupDirectories()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:229)
  and
  [loadDependencies()](/home/doktersmol/Documents/MultiScholaR/R/func_general_filemgmt.R:374)
  remain live, so this handover is now archival context for bucket 12.
