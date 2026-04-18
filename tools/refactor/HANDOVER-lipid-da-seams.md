# Lipid DA Seam Map

## Goal

Apply the bounded exact-source lipid-DA helper waves for
[func_lipid_da.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_da.R:1)
while keeping live lipid-DA behavior stable.

## Current Position In The Flow

- April 14, 2026 classification for
  [func_lipid_da.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_da.R:1)
  is `direct-extraction-ready`.
- The first bounded helper wave is now applied live via
  [tools/refactor/manifest-lipid-da-wave1.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-da-wave1.yml:1)
  into:
  - [R/func_lipid_da_model.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_da_model.R:1)
  - [R/func_lipid_da_results.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_da_results.R:1)
- The applied wave covers:
  - `runTestsContrastsLipidDA()`
  - `runLipidsDA()`
  - `createLipidDaResultsLongFormat()`
- The live collate artifact now exists at
  [tools/refactor/collate-lipid-da-wave1.txt](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/collate-lipid-da-wave1.txt:1).
- The second bounded helper wave is now applied live via
  [tools/refactor/manifest-lipid-da-wave2.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-da-wave2.yml:1)
  into:
  - [R/func_lipid_da_plotting.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_da_plotting.R:1)
  - [R/func_lipid_da_export.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_da_export.R:1)
- The applied second wave covers:
  - `generateLipidDAVolcanoPlotGlimma()`
  - `generateLipidDAHeatmap()`
  - `generateLipidDAVolcanoStatic()`
  - `outputLipidDaResultsAllContrasts()`
- The live collate artifact for the second wave now exists at
  [tools/refactor/collate-lipid-da-wave2.txt](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/collate-lipid-da-wave2.txt:1).
- `DESCRIPTION` now collates the four lipid-DA helper files after
  [func_lipid_da.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_da.R:1).
- Live
  [func_lipid_da.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_da.R:1)
  is reduced to `155` LOC after the second exact-source removal, while
  [R/func_lipid_da_plotting.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_da_plotting.R:1)
  is `773` LOC and
  [R/func_lipid_da_export.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_da_export.R:1)
  is `415` LOC.

## Existing Safety Net

- Focused source-based characterization gate:
  - [tests/testthat/test-lipid-02-da-core-helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tests/testthat/test-lipid-02-da-core-helpers.R:1)
- Replay command:
  - `Rscript -e "source('R/allGenerics.R'); source('R/func_general_s4_objects.R'); source('R/func_lipid_s4_objects.R'); source('R/func_lipid_da.R'); source('R/func_lipid_da_model.R'); source('R/func_lipid_da_results.R'); source('R/func_lipid_da_plotting.R'); source('R/func_lipid_da_export.R'); source('R/func_lipid_import.R'); testthat::test_file('tests/testthat/test-lipid-02-da-core-helpers.R', stop_on_failure = TRUE)"`
- Gate coverage now characterizes:
  - `runTestsContrastsLipidDA()`
  - `runLipidsDA()`
  - `createLipidDaResultsLongFormat()`
- The focused source gate reran green after the first live apply and again
  after the second live apply.

## Notes

- No prior target handover existed for this lane.
- April 14, 2026 lane-control hardening note:
  the overnight issue was not a lipid DA code failure. After the earlier
  `func_lipid_qc` commit, generic backlog fallback drifted the worktree into a
  proteomics bucket. The loop and supervisor now persist lane
  `scopePrefixes` (`R/func_lipid_`, `R/mod_lipid_`) so unattended fallback
  selection stays inside lipid families only.
- This checkpoint intentionally leaves `getCountsTable()` and
  `getLipidQuantData()` in live
  [func_lipid_da.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_da.R:1).
- `getLipidQuantData()` is currently a duplicate-symbol surface between
  [func_lipid_da.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_da.R:123)
  and
  [func_lipid_import.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import.R:156);
  it stays out of the second live apply as an explicit follow-up compatibility
  decision, but it no longer blocks this backlog target because the public
  wrapper file is back inside the ideal size band.
- The focused source gate avoids a sandbox dependency on `limma` by
  characterizing `runLipidsDA()` with a stubbed
  `runTestsContrastsLipidDA()` return surface, while the validation-path test
  still exercises the live undefined-level guard in
  `runTestsContrastsLipidDA()`.
- Lipid-DA helper stabilization is complete for this backlog target.
- Any later work on `getLipidQuantData()` ownership should be handled as a
  separate compatibility cleanup shared with lipid import/QC consumers, not as
  another required stabilization checkpoint for
  [func_lipid_da.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_da.R:1).
- The lipid worktree lane no longer stops after this target:
  - `stabilization-loop.py queue` was used to seed the remaining scoped lipid
    workflow targets into the same lane state
  - the active unattended follow-on target is now
    [func_lipid_import.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_import.R:1)
    under run `manual-func-lipid-import-iter-016`
  - the scoped queue now continues through the remaining `R/func_lipid_*` and
    `R/mod_lipid_*` surfaces instead of truthfully parking as `completed`
    after only the initial manual seed set
