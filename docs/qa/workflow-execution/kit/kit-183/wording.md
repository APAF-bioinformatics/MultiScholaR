# KIT-183 wording list (VAL-41: implemented as proposed; the owner edits before release)

Tab (`R/utils_integration_gui.R`): `enrich_test`, `enrich_test_camera`, `enrich_test_gsea`, `enrich_test_camera_pr`
(not shown under O1), `enrich_test_camera_note`, `enrich_test_gsea_note`, `enrich_test_camera_pr_note` (not shown
under O1), `enrich_values_ruv_note`, `enrich_caution`, `enrich_factor_strength`, `enrich_test_few_samples`,
`enrich_test_counts`, `enrich_direction_oriented` (added: the first sentence of `enrich_direction`, for camera
results), `enrich_direction_camera`, `enrich_seed` (relabelled).

Record and Methods (`R/utils_integration_record.R`): `col_enrich_score` (relabelled), `col_enrich_strength` (added:
the summary sheet's column for the factor-strength sentence), `enrich_string_test_note`, `enrich_ranking_value_camera`,
`enrich_ranking_value_camera_pr`, `methods_test_name_camera`, `methods_test_name_camera_pr`, `methods_test_name_gsea`,
`methods_test_settings_camera`, `methods_test_settings_camera_pr`, `methods_test_settings_gsea`,
`methods_enrich_camera`, `methods_enrich_camera_pr`, `methods_enrich_limits_tests`.

Refusals (`R/utils_integration_enrichment.R`, 8.4): `integration_enrichment_too_few_complete`,
`integration_enrichment_too_few_samples_for_test`, `integration_enrichment_test_not_offered`.

Numbers written in from the simulation (results.md): R_cam 6.6 %, I_cam 5.8 to 7.4 %, S_cam 12.8 %, R_gsea_0.3 43.0 %,
I_gsea_0.3 41.9 to 44.1 %, R_cpr_0.3 29.9 %, R_cam_6 5.0 %.

Replaced after the implementation audit (AUDIT-KIT-183-impl I01, I02, I10; VAL-85): `enrich_test_camera_note` (the
range 0.4 to 6.6 % and the 10.1 % case), `methods_enrich_limits_tests` (the same, and the pre-ranked figures for
metabolites and genes), `enrich_factor_strength` (the pre-ranked test's 60.1 %), `enrich_test_gsea_note` (13.5 % for
metabolites to 43.0 % for genes), `enrich_test_camera_pr_note` (not offered in this release), `enrich_test_few_samples`
(the quantity named; up to 9.9 % for the most variable features), `enrich_test_counts` ("ranked features with an
identifier"); NEWS.
