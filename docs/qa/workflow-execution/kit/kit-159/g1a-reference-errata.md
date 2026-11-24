# G1(r): the KIT-172 key defect (coordinator's ruling on re-audit 3, R3.B3, 2026-11-24)

**Ruling.** For the 18 groups below, G1(r)'s own rows are adopted. KIT-172's committed
`kit-172/k156-gui-entry-failures.tsv` (SHA-256 0429bb99...35ce) is not edited: it is a merged
ticket's evidence with its own copy and hash. The exact pairs are committed as
`harness/msrval/inst/extdata/g1r-key-errata.tsv`. `valAdoptReferenceRows` accepts a listed pair
only where lane, step, class and the earlier key all match. Any other disagreement still stops.

**The defect.** KIT-172's reader kept waldo's line window in one key. The source line is
`lines(actual$outputs$metabolomics_workflow-norm-correlation_filter_summary)[2:8] vs ...`. Here
`)` closes `lines(` and `[2:8]` is the displayed window; neither is part of the output id.

**The 18 pairs.** Each lane covers the steps `correlation`, `session_export`, `da_load`, `da`,
`summary_refusal`, `summary_save`, `summary_copy`, `summary_report` and `summary_export`. All 18
are class `metab_fixture_4`.

| lane | KIT-172's key | G1(r)'s key |
|---|---|---|
| `metab_gc` (9 steps) | `outputs/metabolomics_workflow-norm-correlation_filter_summary)[2:8]` | `outputs/metabolomics_workflow-norm-correlation_filter_summary` |
| `metab_combined` (9 steps) | `outputs/metabolomics_workflow-norm-correlation_filter_summary)[2:11]` | `outputs/metabolomics_workflow-norm-correlation_filter_summary` |

**Everything else agrees.** With the window removed, all 18 groups agree on lane, step, class,
fields, the full key set, oracle and source hash (re-audit 3, R3.B2).

**The oracle check.** The id without the window appears 20 times in each of
`tests/testdata/gui-entry/metab_gc.json`, `metab_combined.json` and `metab_lc.json`. It never
appears with `)[`.

**The truncation flags.** KIT-172 flags 16 summary groups as `names_list_truncated` = `yes`:
`summary_save`, `summary_copy`, `summary_report` and `summary_export` of `lipid_canonical`,
`metab_lc`, `metab_gc` and `metab_combined`. G1(r) flags them `no`, correctly: CI's waldo prints
every differing key on its own line. The adoption rule allows this; it is not one of the 18.

**Who read the keys.** KIT-172's own control never read them. `kit-172/gui-entry-control.R` uses
only `reference$step` per lane, so its counts and conclusions do not depend on the key text.

**Effect on verdicts.** None. The (lane, step) pairs, the classes and the per-(step, class)
counts are unchanged: `metab_fixture_4` gives 10 per metabolomics lane.
