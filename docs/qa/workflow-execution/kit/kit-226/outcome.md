# KIT-226 outcome (ticket r1.2; owner's decisions VAL-196, wording corrected by VAL-208)

State: R1, R2, the merge of staging 16112917 (KIT-215 in), control C and the re-freeze are done. The merge check rides
a merge train and the shared full suite is pending (no own full suite, coordinator).

## Commits (branch `kit-226`, from `origin/v0.9-staging` cf8af245)

| Commit | What |
|---|---|
| 09d92168 | tests only: new `tests/testthat/test-kit226-gprofiler-columns.R` (T1 to T5, G1 to G3); the three declared expectations `test-prot-11-enrichment-module-contracts.R:1709`, `:1711`, `:1712`; the declared fourth `test-enrichment-significance.R:497` (archive line, OD-5) |
| (merge) | staging 16112917 merged; ledger only conflicted (staging's rows kept, this row after VAL-229) |
| 257349af | control C script `redump.R` (records only) |
| da3ca8d0 | records: prediction before control C |
| 1e2641be | re-freeze: exactly the 49 leaves control C showed (VAL-196 with VAL-208) |
| 77ab881f | fix: `R/utils_workflow_enrichment.R` (builder, `@noRd` block), `R/func_general_enrichment_engine_helpers.R` (`gprofilerAdjustmentMethod`, figure labels, published TSV column), `R/mod_omic_enrich_panels.R:596` (archive line); `NEWS.md`; ledger `VAL-257` |

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
- same: 18 files (all other R2 files, `test-workflow-class-declarations.R` added because it calls the plot builders); 15 of them informative (`test-prot-11-enrichment-module-contracts.R` and `test-prot-11a-enrichment-module-characterization.R` run no test, quarantines Q-163 and Q-164; `test-e2e-proteomics-enrichment.R` skips both its tests without a browser);
- new failures, every one at a stored reference awaiting its re-freeze (by case, step and field, from `r2/failures.txt`; base listing `r2/failures-base.txt`):
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
| Control C (stop rule 2) | met | 45587822866, 45587827135, 45587831399; `control/leaves.tsv` |
| Re-freeze of exactly the changed leaves | met | 1e2641be (49 lines) |
| Merge check, shared R3 | pending | rides a merge train |

## Stored references re-frozen (1e2641be; authority VAL-196 with VAL-208)

`test-omic-enrich-module.R`: the 38 leaves of section 5 (i), the two plot leaves (`:1842`, `:1847`), the summary leaf
(`:1259`), and the two `kernel_gprofiler` run `files` leaves of the pathway TSVs (`:1923`, `:1926`; named by VAL-208).
`tests/testdata/gui-entry/prot_dia.json`: 6 leaves. `test-enrichment-significance.R:497` (named by VAL-208) was
changed in the tests-first commit. `test-omics-enrich-core.R` reads the W1 frozen values, so it follows that
re-freeze with no edit of its own.
`tests/testthat/helper-module-ci-prot-enrich.R:167-170` keeps its own fixture with the old mapping (left as is).

## Attribution and re-freeze go leaf by leaf (audit K226-I04, I07)

The fields are not disjoint in two places, so every difference is attributed, and every re-freeze edits, by leaf:
- W1 `omicEnrichW1Frozen[["kernel_gprofiler"]]$run$files` holds KIT-220's NumSig `.tab` leaf (`:1914`) and this
  ticket's two pathway TSV leaves (`:1923`, `:1926`); the W1 test compares field by field (`:3669`).
- `tests/testdata/gui-entry/prot_dia.json` is compared step by step (`helper-gui-entry.R:414`): KIT-216's leaves
  (`head`, `history`, `states`, `outputs$…da_status`) and this ticket's `outputs$…gprofiler_results_table` fail in one
  expectation in the same six steps.

Order (coordinator): KIT-215 merges first; control C, the re-freeze and the merge-check run follow on the
coordinator's word. The three stored values VAL-196 (4) does not name (`:1923`, `:1926`,
`test-enrichment-significance.R:497`) are allowed through the control (coordinator's ruling (a); recorded as VAL-208,
which corrects VAL-196's wording: the re-freeze authority is VAL-196 with VAL-208).

**Control C must show exactly these 49 changed leaves** (base = staging's head when C runs; head = `kit-226` merged
onto it; base dumps equal the stored value at each; two head dumps in separate processes identical):
- W1 `test-omic-enrich-module.R`, 43 leaves:
  - `gprofiler` 12: run `values[["gprofiler_results_table"]]`, `bus`, `display` ×2, `returned` ×2
    (`all_enrichment_results`, `gprofiler_results`); `filter_up`, `filter_down`, `filter_all` `values`; `download`
    `values`; `download` `gprofiler2_results.tsv` and `enrichment_analysis_summary.txt` (`:1259`, `:1260`);
  - `organism_filter` 12: sessions `no_mapping` and `mapped`, the six run leaves each;
  - `restore` 9: the six run leaves; `restore` step `values`; `restore` ×2;
  - `kernel_gprofiler` 10: run `values` table and both plot leaves (`:1842`, `:1847`); `bus`; `display` ×2;
    `returned` ×2; `files` ×2 (`:1923`, `:1926`).
- `prot_dia.json`, 6 leaves: `outputs$proteomics_workflow-enrichment_analysis-gprofiler_results_table` in steps
  enrichment, summary_refusal, summary_save, summary_copy, summary_report, summary_export.
- Unchanged: `current_contrast_results`, every `history`, `states`, `product`, `notes` leaf, `omicEnrichW113Frozen`,
  `omicEnrichW1FrozenUi`. Pre-existing differences (recorded in `kit-226/control/leaves.tsv`, left as they are):
  `prot_dia`'s four integration `.tsv`/`.parquet` steps; the clusterProfiler, kernel_clusterprofiler and KIT-113 W1
  failures.
- R2's head values C should match: table `<text> 258d7618…`; bus `dcc3ee2f940bd0fa`; display
  `10x23 8b713e2fe8c8dfa5`, `<saved> 10 542f2d5f5bec9ac4`; returned `37aa971ecab5066f`, `8b713e2fe8c8dfa5`; summary
  `c5b5570d…`, archive TSV `71e8ee7c…`; kernel plots `e29b944a…` ×2, table `ef47d5f9…`, bus `e262161c5371ef44`,
  display `10x22 cb4a9aea9d03ca0d`, `<saved> 10 b3c9a33ab15f1fd0`, returned `f49c6bf6d7f341f5`, `cb4a9aea9d03ca0d`,
  pathway TSVs `a7a539d1…` ×2.
- The merge-check run then holds R2's 23 files on the merge with staging plus any merged partner's new test file,
  `check: true`: the new file 38/0, the three moved files back to staging's pre-existing failure lists, check
  identical to base; then the shared R3 (ticket section 6).

## Before control C: does KIT-215's merge move the prediction? (stated before sending)

Staging 16112917 (merge train 3, KIT-215 = VAL-227) merged into `kit-226` (merge commit; only the ledger conflicted,
resolved by keeping staging's rows and putting this ticket's `VAL-257` row after VAL-229; no code or test
conflict). Read, not run:
- `git diff cf8af245 16112917` does not touch `tests/testthat/test-omic-enrich-module.R`,
  `tests/testdata/gui-entry/prot_dia.json`, `test-omics-enrich-core.R`, `helper-gui-entry.R`,
  `test-enrichment-significance.R`, `R/func_general_enrichment_engine_helpers.R` or `R/utils_workflow_dia_enrichment.R`:
  no stored value this control reads was re-frozen by the train.
- KIT-215's edits to `R/utils_workflow_enrichment.R` and `R/mod_omic_enrich_panels.R` add the not-tested and no-list
  fields, notices and summary lines only when a list was not tested or a direction had no list; its own control moved
  no stored capture, so the W1 g:Profiler cases and the prot_dia lane have none.
- Prediction therefore unchanged: the same 49 leaves, with R2's head values. Any other leaf that differs base/head is
  stop rule 2.

## Control C (coordinator's word; observed)

Script `redump.R` (KIT-220's form: each test file in its own process with `MULTISCHOLAR_KIT107_W1_DUMP` and
`MULTISCHOLAR_GUI_ENTRY_RECORD_DIR`); files `test-omic-enrich-module.R` and `test-gui-entry-proteomics.R`. Base
3a48351b = staging 16112917 plus the script only (run 45587822866); head da3ca8d0 = `kit-226` merged with 16112917
plus the script, twice in separate jobs (45587827135, 45587831399). Compared leaf by leaf by `compare_dumps.py`
(text parse of the dumps and JSON; no R), table `control/leaves.tsv` (100 differing leaves):
- **Changed base/head: exactly 49 leaves**, the 43 W1 leaves and 6 `prot_dia.json` leaves listed above; at each the
  base dump equals the stored value, the two heads are equal, and the head value is R2's (e.g. table
  `<text> 258d7618…`, kernel plots `e29b944a…` ×2, pathway TSVs `a7a539d1…` ×2, summary `c5b5570d…`, archive TSV
  `71e8ee7c…`).
- **No other leaf changed**: 51 pre-existing differences (stored ≠ base = head): 27 W1 (the clusterProfiler and
  kernel_clusterprofiler cases: status text, notes, product parameters, display, files, returned, services, states,
  plots; one clusterprofiler download summary) and 24 `prot_dia` (the four summary steps' integration `.tsv` against
  `.parquet` artifacts). Left as they are.
- The two head W1 dumps are byte-identical; the `omicEnrichW1FrozenUi` dump is identical on all three runs.
- Test counts inside the script (both files): base 1263 expectations, 23 failed; heads 1263, 55 failed each (the
  32 extra are the not-yet-re-frozen leaves).
- Stop rule 2 not triggered. The re-freeze 1e2641be edits exactly those 49 lines; re-running the comparison with the
  re-frozen files against head 1 leaves only the 51 pre-existing differences.

## Core-hours

R1 0.92, R2 11.13, C 7.54 (113.0 job-minutes): 19.59 against 23.0 expected, ceiling 304.67.

## Not checked

- A live g:Profiler answer was not observed: whether today's service states `significance_threshold_method`, and
  whether it spells the default `gSCS` or `g_SCS` and BH `fdr` or otherwise (no service called; the stated values
  in T3 are fixture strings). By the coordinator's ruling (d) there is no fallback to the requested method: the
  column shows what the answer carries, and "not recorded" otherwise.
- The enrichment workbook carrying the column (inferred from its reader,
  `R/func_general_filemgmt_results_summary_helpers.R:689-698`); the Copy/CSV/Excel exports; the application in a
  browser; the published files and report of a run made before the fix (not rebuilt; change report).
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
  Fix 77ab881f: `R/utils_workflow_enrichment.R:145`, `:159` (`p_adjust_method` per direction), `:168-174`
  (`p.adjust` once, `Count` from `intersection_size`; `pvalue`/`qvalue` no longer made), roxygen `:104-115`;
  `R/func_general_enrichment_engine_helpers.R:200`, `:228` (labels), `:233-244` (`gprofilerAdjustmentMethod`),
  `:414-415` (published TSV); `R/mod_omic_enrich_panels.R:596` (archive line). Legacy figure labels: 744a2b27
  `R/func_general_enrichment.R:3003`, `:3029`.
- **What changed in results.** Tables: `pvalue` and `qvalue` removed, `p.adjust` = the service's value, new
  `p_adjust_method` (the answer's own statement, e.g. `g_SCS`, or "not recorded"), `Count` = `intersection_size`.
  Published TSV/xlsx: new `p_adjust_method` column. Figure: axis and hover text say "adjusted p-value"; points
  unchanged. On Open, a project saved before the fix shows the corrected table and figure, and a 'Download results'
  archive made after the fix is corrected (both are rebuilt from the stored g:Profiler answer; a project that stored
  its figure, saved before KIT-119, keeps the old axis label in the tab). The published per-comparison TSVs, the
  workbook, the PNG/HTML figures and the report of an enrichment run made before the fix keep the old labels until
  the enrichment is run again (Summary copies the files the run wrote); files already downloaded keep them too.
- **Evidence.** R1 45504002943, R2 45506344960 (deterministic g:Profiler double; no live service). Limits: the
  package's g:Profiler doubles state no correction, so every tested table reads "not recorded" and the stated values
  `g_SCS`/`fdr` are fixture strings (T3); whether today's service states the method (and spells the default `gSCS`
  or `g_SCS`) was not observed; the stored references are re-frozen by control C; the shared full suite is pending.
