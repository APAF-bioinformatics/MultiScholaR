# KIT-172: declared changes visible at an offset above 0 (written before the control request, section 7(c) step 2)

Changed step per lane: `norm` for `lipid_canonical`, `metab_lc`, `metab_gc`, `metab_combined` (the log step runs in the
`norm` step of every lane). `<omic>` is `lipidomics` or `metabolomics`, the module's output prefix in the frozen files.

| Capture key | Steps | Expected value on the head with `log_offset = 1` supplied |
|---|---|---|
| `outputs/<omic>_workflow-norm-log2_recorded_offset` | before `norm` (`import_refusal` to `norm_refusal`) | `"<NULL>"` (absent on the base) |
| same | from `norm` on | `<ui> <sha256>` of the rendered `shiny::helpText(sprintf(omicGuiLog2Texts()[["recorded_offset"]], "1"))`, computed in the head child and written to `head1/expected-recorded-offset.txt` |
| `outputs/<omic>_workflow-de-da_summary_stats` (the DA summary, `output$da_summary_stats`, `R/mod_omic_da.R:233`) | before `da` | the base's value (`"<req>"`) |
| same | from `da` on | the base's text plus the line `Not tested (no value in a compared group): 0` (the lane fixtures hold no missing or zero sample cell; the script checks the fixture files before the passes) |

Every other key must be identical to the base's at every step; no F1b fix changes a capture on these lanes (fixture
intensities 2.4e3 and above, no missing or zero cell: the QC counts read the same, and the volcano's sentence appears
only when an assay holds missing cells). Not captured as content (`<present>` only): the report, its copied template,
the study parameters file, the session export's summary text and the analysis record, whose wording changes are
therefore not visible to the control.

Classes in `gui-entry-control.tsv`: `frozen-base` rows `integration_artifacts` or `metab_fixture_4` (lane-range rule:
metabolomics lanes from `norm` on); `base-head1` rows `declared_change (KIT-172)`; `head1-head0` rows from `norm` on
`plain_log2 (KIT-172)`; `base-head0` rows the class of the same key's `base-head1` row, if any, plus
`plain_log2 (KIT-172)`. Any row classed `unexplained` stops the ticket.
