# KIT-226 outcome (ticket r1.2; owner's decisions VAL-196)

State: R1 and R2 done. Stopped, as ordered, **before** the control C, any re-freeze commit and the shared full
suite (KIT-215, KIT-216 and KIT-220 move stored captures in the same files; the coordinator orders the re-dumps).

## Commits (branch `kit-226`, from `origin/v0.9-staging` cf8af245)

| Commit | What |
|---|---|
| 09d92168 | tests only: new `tests/testthat/test-kit226-gprofiler-columns.R` (T1 to T5, G1 to G3); the three declared expectations `test-prot-11-enrichment-module-contracts.R:1709`, `:1711`, `:1712`; the declared fourth `test-enrichment-significance.R:497` (archive line, OD-5) |
| 77ab881f | fix: `R/utils_workflow_enrichment.R` (builder, `@noRd` block), `R/func_general_enrichment_engine_helpers.R` (`gprofilerAdjustmentMethod`, figure labels, published TSV column), `R/mod_omic_enrich_panels.R:596` (archive line); `NEWS.md`; ledger `VAL-<merge>` |

## Observed

**R1, run 45504002943, at 09d92168 (no `R/` change).** `test-kit226-gprofiler-columns.R` 19 passed, 19 failed. Each
failure is at the ticket's stated base value: T1 `Count` = 1392 1780 892 131 418 (twice); T2 `pvalue`, `qvalue`
present; T3 `p_adjust_method` NULL (for `g_SCS`, `fdr`, absent, mixed and empty); T4 the same on the restore path;
T5 the archive TSV has `pvalue`/`qvalue`, `Count` = `term_size`, no `p_adjust_method`; G2 the published TSV has no
`p_adjust_method`, axis "-log10(p-value)", tooltip "P-value: ". G1, G3 and G2's data/mapping guards passed.
`test-enrichment-significance.R` failed only at `:497` (summary line "Q-value Cutoff:").
`test-prot-11-enrichment-module-contracts.R`: **0 expectations run**; the file is quarantined by a whole-file skip at
line 1 (Q-163, KIT-107 W3). The three declared expectation changes are therefore inert on base and head; T1 to T3
cover the same behaviour. Stop rule 1 not triggered (no failing-first row passed on base).

**R2, run 45506344960, 77ab881f against 09d92168, 23 files, 3 shards per side, check on both** (`r2/compare.tsv`):
- fixed: `test-kit226-gprofiler-columns.R` (38 passed, 0 failed), `test-enrichment-significance.R` (191/0);
- same: 18 files (all other R2 files, `test-workflow-class-declarations.R` added because it calls the plot builders);
- new failures, every one at a stored reference awaiting its re-freeze (by case, step and field, from `failures.txt`):
  - `test-omic-enrich-module.R` W1: `gprofiler` run values/bus/display/returned, filter_up/down/all values, download
    values and download; `organism_filter` no_mapping/run and mapped/run values/bus/display/returned; `restore` run
    values/bus/display/returned, restore values and restore; `kernel_gprofiler` run values (both plot leaves `:1842`,
    `:1847` and the table leaf), bus, display, returned and files (only the two pathway TSVs, OD-3). Download shows
    exactly the two predicted leaves (`enrichment_analysis_summary.txt`, `gprofiler2_results.tsv`). Two further
    tests, (R3, A107-04) and (R3, A107-08), compare their run with the same frozen `display` leaves
    (`kernel_gprofiler` and `gprofiler` run display); no further stored value;
  - `test-omics-enrich-core.R`: the core-side gprofiler and kernel_gprofiler run steps and the restore capture, which
    compare with the same W1 frozen leaves;
  - `test-gui-entry-proteomics.R` lane prot_dia: exactly the six `gprofiler_results_table` leaves (steps enrichment,
    summary_refusal, summary_save, summary_copy, summary_report, summary_export); the four base failures (integration
    files `.tsv` against `.parquet`, pre-existing) unchanged on head.
- Every clusterProfiler and kernel_clusterprofiler W1 failure, and the KIT-113 W1 failures, are identical on base and
  head (pre-existing).
- R CMD check: base and head 0 ERROR, 3 WARNING, 3 NOTE, findings identical.

## Inferred

- The fields that moved are those section 5 predicts (W1 `values`, `bus`, `display`, `returned`, `download`,
  `restore`; under OD-3 `files`); no `history`, `states`, `product`, `notes` or `services` field of a g:Profiler case
  moved. Leaf-by-leaf equality with section 5 (stop rule 2) is control C's to show; this run shows fields and the
  first differing leaves only.
- Stop rule 3: G1, G2 (data and mappings) and G3 pass on head: not triggered.

## Acceptance (ticket section 4 to 6)

| Item | Status | Evidence |
|---|---|---|
| T1 to T5 fail first, pass on head | met | R1 45504002943; R2 45506344960 |
| G1, G2 (data, mappings), G3 unchanged | met | pass on both sides |
| Three declared expectations `test-prot-11…:1709/:1711/:1712` | not demonstrable | file quarantined (Q-163); edits made, inert |
| `test-enrichment-significance.R:497` | met | fails on base, passes on head |
| R2 files: no verdict change outside stored leaves | met (by field) | `r2/compare.tsv` |
| R CMD check unchanged | met | `r2/check-compare.txt` |
| Control C, re-freeze, shared R3 | pending | held for the coordinator |

## Stored references awaiting re-freeze (not edited)

`test-omic-enrich-module.R`: the 38 leaves of section 5 (i), the two plot leaves (`:1842`, `:1847`), the summary leaf
(`:1259`), and the two `kernel_gprofiler` run `files` leaves of the pathway TSVs (OD-3; VAL-196 (4) does not name
them: coordinator to confirm). `tests/testdata/gui-entry/prot_dia.json`: 6 leaves. `test-omics-enrich-core.R`
reads the W1 frozen values, so it follows that re-freeze with no edit of its own.
`tests/testthat/helper-module-ci-prot-enrich.R:167-170` keeps its own fixture with the old mapping (left as is).

## Core-hours

R1 0.92, R2 11.13: 12.05 against 17.0 expected for R1 and R2 (ticket 23.0 with C), ceiling 304.67.

## Not checked

- Whether today's live g:Profiler service states `significance_threshold_method`, and its spelling under `fdr` (no
  service called; only the fixture strings are tested). If it does not, live tables read "not recorded".
- Leaf-by-leaf equality of the moved values with section 5 (control C).
- The three declared expectations (quarantined file).
- The shared full suite.

## For the change report

- **What was wrong.** On the g:Profiler route the enrichment results table copied g:Profiler's `p_value`, which the
  service has already adjusted for multiple testing (g:SCS by default), into three columns named `pvalue`, `p.adjust`
  and `qvalue`, and put the term's total size into `Count`. The in-app table, its Copy/CSV/Excel buttons, the tab's
  'Download results' archive (`gprofiler2_results.tsv`) and a reopened project showed an adjusted value as a nominal
  p-value and as a q-value, and a term's size as the number of the user's genes in it. The published per-comparison
  TSVs and the workbook said nowhere that `p_value` is adjusted. The figure (app, published PNG/HTML, report) labelled
  the adjusted value "-log10(p-value)" / "P-value"; the archive summary called the cutoff "Q-value Cutoff". Every
  g:Profiler enrichment user since the legacy pipeline. No term, flag or order was wrong.
- **Code.** Base: `R/utils_workflow_enrichment.R:154-161` (cf8af245; legacy `R/mod_prot_enrich.R:1450-1456` at
  744a2b27), `R/func_general_enrichment_engine_helpers.R:200`, `:228`, `:398-403`, `R/mod_omic_enrich_panels.R:596`.
  Fix 77ab881f: the builder (`buildProtEnrichAllContrastResults`), `gprofilerAdjustmentMethod`, the figure, the
  published TSV writer (`generate_enrichment_plots`) and the archive line.
- **What changed in results.** Tables: `pvalue` and `qvalue` removed, `p.adjust` = the service's value, new
  `p_adjust_method` (the answer's own statement, e.g. `g_SCS`, or "not recorded"), `Count` = `intersection_size`.
  Published TSV/xlsx: new `p_adjust_method` column. Figure: axis and hover text say "adjusted p-value"; points
  unchanged. Projects saved before the fix show the corrected table on Open; files already downloaded keep the old
  labels.
- **Evidence.** R1 45504002943, R2 45506344960 (deterministic g:Profiler double; no live service).
