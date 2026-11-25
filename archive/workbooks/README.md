# Archived R Markdown workbooks (not installed, not supported)

These eight workbooks (five DIA proteomics, two metabolomics, one multi-omics integration) were moved here from
`inst/workbooks/` on 2026-11-25 by the owner's decision (methods ledger VAL-209; ticket KIT-234). They are kept for a
later version, which will bring them back for educational use. They are not part of the built or installed package
(`^archive$` in `.Rbuildignore`), and the application no longer copies them into projects.

State when archived:

- Each of the seven proteomics and metabolomics workbooks stopped with an error at one or more steps against this
  version (KIT-231 G0, run-tests 45511537888, observed chunk by chunk, not whole workbooks; two DIA limpa design
  chunks were not observed). The records are on branch `origin/kit-231`, `docs/qa/workflow-execution/kit/kit-231/`.
- The integration workbook rendered end to end when last run (KIT-143, runs 44864211198 and 44875316282) and carries
  KIT-184's paragraph on adjusting for the removed RUV factors. The functions it calls (`integrationFindInputs()` to
  `integrationReplayProject()`) remain exported.

The repair inventory for bringing them back: tickets KIT-230, KIT-231 and KIT-232 (parked), KIT-231's gate records,
and the assertions KIT-234 removed with the workbooks: KIT-203 S11 and KIT-173 T10's loop (the five DIA workbooks pass
`significance_rule = "bh_whole_family"` and `go_propagation = "is_a_part_of_regulates"` to `processEnrichments()`),
KIT-143's static test and KIT-154's workbook-code checks (`test-integration-workbook.R`), the render of the
integration workbook (`test-integration-workbook-mofa.R`) and KIT-184's four phrases
(`test-integration-removed-factors.R`). To bring them back: `git mv archive/workbooks inst/workbooks`.
