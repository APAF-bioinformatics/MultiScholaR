# Lipid Import Seam Map

## Goal

Apply bounded exact-source helper waves for
[func_lipid_import.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import.R:1)
while keeping the live lipid-import detection and mapping behavior stable.

## Current Position In The Flow

- April 14, 2026 classification for
  [func_lipid_import.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import.R:1)
  is `direct-extraction-ready`.
- The first bounded helper wave is now applied live via
  [tools/refactor/manifest-lipid-import-wave1.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-import-wave1.yml:1)
  into
  [R/func_lipid_import_detection.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import_detection.R:1).
- The applied wave covers:
  - `detectLipidomicsFormat()`
  - `findLipidMatchingColumn()`
  - `getLipidomicsColumnDefaults()`
  - `validateColumnMapping()`
  - `validateLipidColumnMapping()`
- The live collate artifact now exists at
  [tools/refactor/collate-lipid-import-wave1.txt](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/collate-lipid-import-wave1.txt:1).
- `DESCRIPTION` now collates
  [R/func_lipid_import_detection.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import_detection.R:1)
  immediately after
  [R/func_lipid_import.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import.R:1).
- The second bounded helper wave is now applied live via
  [tools/refactor/manifest-lipid-import-wave2.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-import-wave2.yml:1)
  into
  [R/func_lipid_import_readers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import_readers.R:1).
- The applied wave covers:
  - `importMSDIALData()`
  - `importLipidMSDIALData()`
  - `importLipidSearchData()`
- The live collate artifact now also exists at
  [tools/refactor/collate-lipid-import-wave2.txt](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/collate-lipid-import-wave2.txt:1).
- `DESCRIPTION` now collates
  [R/func_lipid_import_readers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import_readers.R:1)
  immediately after
  [R/func_lipid_import_detection.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import_detection.R:1).
- The third bounded helper wave is now applied live via
  [tools/refactor/manifest-lipid-import-wave3.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-import-wave3.yml:1)
  into
  [R/func_lipid_import_core.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import_core.R:1).
- The applied wave covers:
  - `createLipidomicsAssayData()`
  - `getLipidQuantData()`
- The live collate artifact now also exists at
  [tools/refactor/collate-lipid-import-wave3.txt](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/collate-lipid-import-wave3.txt:1).
- `DESCRIPTION` now collates
  [R/func_lipid_import_core.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import_core.R:1)
  immediately after
  [R/func_lipid_import_readers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import_readers.R:1).
- Live
  [R/func_lipid_import.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import.R:1)
  now measures `21` LOC as a breadcrumb shell, and
  [R/func_lipid_import_detection.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import_detection.R:1)
  measures `367` LOC, while
  [R/func_lipid_import_readers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import_readers.R:1)
  measures `260` LOC, while
  [R/func_lipid_import_core.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import_core.R:1)
  measures `132` LOC.

## Existing Safety Net

- Focused source-based characterization gate:
  [tests/testthat/test-lipid-00-import-detection-helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tests/testthat/test-lipid-00-import-detection-helpers.R:1)
- Replay command:
  - `Rscript -e "options(device = function(...) grDevices::pdf(file = tempfile(fileext = '.pdf'))); source('R/func_general_helpers.R'); source('R/func_lipid_import.R'); source('R/func_lipid_import_detection.R'); source('R/func_lipid_import_readers.R'); source('R/func_lipid_import_core.R'); testthat::test_file('tests/testthat/test-lipid-00-import-detection-helpers.R', stop_on_failure = TRUE)"`
- Post-apply wave checkers:
  - `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-lipid-import-wave1.yml`
  - `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-lipid-import-wave2.yml`
  - `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-lipid-import-wave3.yml`

## Notes

- No prior target handover existed for this lane.
- The legacy shell in
  [R/func_lipid_import.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import.R:1)
  is now intentionally reduced to a breadcrumb stub after the third helper wave.
- Lipid-import helper stabilization is complete for this backlog target; the
  remaining filename-coupled dev scripts noted in the backlog stay optional
  follow-up work and are not gating this lane.
- `Rscript tools/test_with_renv.R ...` is not currently usable in this worktree
  because `renv/activate.R` is absent, so the focused gate remains source-based
  for this lane.
