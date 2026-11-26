# KIT-268 outcome notes (G0 stage)

Step-0 mismatch for the ticket writer: ticket 2.1 cites `R/mod_prot_da_handlers_run_analysis.R:6-8` for "the kernel in process before
the publication request"; at 4d588ba3 lines 6-8 are `update_contrast_choices`, and that statement is the comment at lines 4-5. All
other cited lines were re-read and hold.

Re-freeze (run 45715855490): see runs.tsv; the oracle is not committed; the captured file is refreeze/run/out/lipid_canonical.json for the decision, outside tests/.

Correction pass (coordinator ruling, 2026-11-26), step (2): `test-omic-enrich-module.R`, the six head-only failures
("(b) in process: refusals", fields `values` of `pre_da/tab_before_da`, `run_before_da`, `da`, `open`, `run_no_contrast`,
`run_taxon_abc`). The failing element is `lines(actual[[11]])`; in the frozen capture (`omicEnrichW1Frozen`, e.g. test lines
846 to 855) element 11 of `values` is `stringdb_plot = "<text> b54fa6cb..."`, which is the hash the failure expects. Its renderer is
`renderProtEnrichStringDbPlot` (R/mod_omic_enrich_panels.R:596), a static plotly "STRING-DB coming soon". It is NOT the proteomics
q-value notice. Nothing of this ticket touches it. The file was 1011/59 on every record since checkpoint-20261126, so the head's
1005/65 is unexplained: STOPPED per the ruling ("if it is anything else, stop"), no fix and no run sent.

Prediction slip, test-kit236-report-sentences.R: head 313 passed / 0 failed; predicted 307 / 0. Derivation: base side (the tests-only
commit) 307 passed / 6 failed = 289 existing + 18 passing T4 expectations, 6 failing; on head the 6 failing now pass: 307 + 6 = 313.

Correction pass, step (a) (local search of the committed records, no run): the head hash of element 11 in R2, `<text> a3f2703c...`, appears
in no earlier record; the same element (`stringdb_plot`) showed a different hash, `<text> 62479f92...`, in kit-201 r3's compare-failures (on
the base side only, the "refusals table_removed" steps), while the frozen value is `b54fa6cb...`. So that element has varied between
runs before. By reading: the renderer (R/mod_omic_enrich_panels.R:596) is `plotly::plot_ly() |> add_text(...)` with a constant text;
whether the widget's JSON carries anything that varies between processes (plotly gives each `plot_ly()` call generated ids) was not read
in plotly's source in this session and stays inferred. R2's 1011/59 to 1005/65 stays on record as observed and unexplained.
