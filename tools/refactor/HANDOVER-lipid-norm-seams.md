# Lipid Norm Seam Map

## Goal

Document the completed direct-helper stabilization checkpoints for
[func_lipid_norm.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_norm.R:1)
after the second live helper apply, while keeping the wrapper-bound config
helper in the separate `mod_lipid_norm.R` lane.

## Current Position In The Flow

- April 14, 2026 classification for
  [func_lipid_norm.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_norm.R:1)
  is `direct-extraction-ready`.
- The first bounded helper wave is recorded via
  [tools/refactor/manifest-lipid-norm-wave1.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-norm-wave1.yml:1).
- That reviewed wave is now applied live from
  [func_lipid_norm.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_norm.R:1)
  into
  [R/func_lipid_norm_ruv_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_norm_ruv_helpers.R:1),
  covering:
  - `runLipidPerAssayRuvOptimization()`
  - `extractLipidBestKPerAssay()`
  - `extractLipidCtrlPerAssay()`
  - `buildLipidCombinedRuvTable()`
- Manifest verification, live extraction, and post-apply checks all passed for
  the wave.
- The new live helper file currently measures `283` lines.
- The live collate artifact was emitted at
  [tools/refactor/collate-lipid-norm-wave1.txt](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/collate-lipid-norm-wave1.txt:1).
- `DESCRIPTION` `Collate:` now includes
  `func_lipid_norm_ruv_helpers.R`.
- The second bounded helper wave is now recorded via
  [tools/refactor/manifest-lipid-norm-wave2.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-norm-wave2.yml:1)
  and is now applied live into
  [R/func_lipid_norm_support_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_norm_support_helpers.R:1),
  covering:
  - `buildLipidItsdSelectionTable()`
  - `generateLipidQcPlots()`
- Manifest verification, live extraction, and post-apply checks all passed for
  the second wave.
- The new live helper file currently measures `321` lines.
- The live collate artifact for that second wave was emitted at
  [tools/refactor/collate-lipid-norm-wave2.txt](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/collate-lipid-norm-wave2.txt:1).
- `DESCRIPTION` `Collate:` now includes
  `func_lipid_norm_ruv_helpers.R` and `func_lipid_norm_support_helpers.R`.
- Live
  [func_lipid_norm.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_norm.R:1)
  is reduced to `87` lines and now only retains `buildLipidNormConfig()`.

## Existing Safety Net

- Focused gate:
  [tests/testthat/test-lipid_norm_exclusion.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tests/testthat/test-lipid_norm_exclusion.R:1)
- Replay command:
  - `Rscript -e "testthat::test_file('tests/testthat/test-lipid_norm_exclusion.R', stop_on_failure = TRUE)"`
- Post-apply replay:
  - `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-lipid-norm-wave2.yml`

## Next Safe Checkpoint

- Treat the direct-helper stabilization lane for
  [func_lipid_norm.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_norm.R:1)
  as complete.
- Keep `buildLipidNormConfig()` deferred to the separate
  [mod_lipid_norm.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_norm.R:1)
  wrapper lane because it is Shiny-input-bound and is not yet reused outside
  that server path.
- If lipid-normalization stabilization resumes, start from
  [mod_lipid_norm.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_norm.R:1)
  instead of reopening helper extraction in
  [func_lipid_norm.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/func_lipid_norm.R:1).

## Notes

- No prior target handover existed for this lane.
- This checkpoint stops after one clean live apply boundary with the focused
  lipid-normalization gate green.
- The staged wave-2 artifacts remain available for audit in
  [tools/refactor/staging/wave2_lipidomics_norm_support_helpers/R/func_lipid_norm_support_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/staging/wave2_lipidomics_norm_support_helpers/R/func_lipid_norm_support_helpers.R:1).
