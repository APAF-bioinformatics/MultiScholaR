# KIT-172 outcome: plain log2 for lipidomics and metabolomics (state: stopped at the GUI-entry control's mechanism rule)

Ticket revision 1.6. Branch `kit-172-plain-log2` from `v0.9-staging` a60c4e9a. Every run on `run-tests` (R 4.6.1,
reference BLAS, ubuntu-22.04); run ids, commits and job minutes in `runs.tsv`.

| Commit | Content |
|---|---|
| c565b8f5 | tests and scripts only (`R/` equal to base): T1a to T10, fixture F172, the contract rewrites of 7(a), the scripts |
| d8ee5388 | plain log2 (method, generic, declared operation, tab, module, actions, export, study parameters, reports, workbooks, DA summary count, man, NEWS); F1b fixes not yet applied |
| 8c02d61a | F1b fixes (first form) |
| 6e4ef83b | T9 base-saved projects (`tests/testdata/kit-172/base-projects.tar.gz`, sidecar `.sha`), R1 evidence, T8 volcano test on a stand-in builder |
| dd5e406c | correction pass after R2: F1b counts read the object's log flag; pre-log sites restored; tests and scripts corrected |
| 9931beea | control script correction (the control's one correction pass) |
| 5250032c | 7(b) re-dump pasted; R2 and R3 evidence; `test-metab-qc-shared.R` stubs take `...` |
| (records) | this file, `runs.tsv`, control mechanism evidence, the ledger row; the head1 wrapper fix prepared, not run |

## Observed

**Failing first (R1, c565b8f5, R/ = base).** Every new test red on base at the stated value: T1a (1) to (3) fail for
every configuration (72 per class); T2 logFC 0.126-type value; T3, T4's offset-0 case and the 7(a) rewrites error with
"`offset` must be a single positive numeric value"; T5 offset 1 everywhere, `norm-log_offset` rendered, no
`logOffset`; T6 the difference is 0, so the rule log2(1 + 1/x) fails by 1.42e-4; T7's not-tested count absent; T8
reads; T9 fixture absent and `omicNormRecordedLog2Offset` absent; T10 no log2 value below 0 (log2(x + 1) > 0).
`r1/tests.tsv`, `r1/failures.txt`.

**Middle commit (R1b, d8ee5388).** T7's QC detected and missing counts and T8 red before the F1b fixes (detected 220
against 465; missingness 12.2 against 10.4); T7's steps, the T9 pre-M3 case and T10 green. `r1/r1b-*`.

**Fix (R2 6e4ef83b, R3 dd5e406c).** Green at R2: T1a (both classes, every configuration and factor; no
`KIT172-T1A-TP` line, so no excess in t or p), T2, T4, T5, T6, T7, T9 (reopen, replay mismatches none, recorded-offset
sentence; new project offset 0), T10 (both routes, both classes, and its three negative cases), the 7(a) rewrites,
and 12 dependent files. Red at R2 and corrected: T3's forced error (a column of one-element lists was coerced without
error; now a column of functions), `test-omic-gui-config.R` (pins of offset 1 and of the DA summary fields: rewritten,
class (a)), `test-metab-01-...`, `test-lipid-01-...` (their raw-data fixtures hold zeros: the first F1b form counted a
raw zero as detected; the counts now read the object's log flag, and on intensities a zero stays not detected, so
both files pass unchanged). R3 green on those files; `test-metab-qc-shared.R` red because its stubs of the two QC
helpers lacked the new argument (stubs take `...` in 5250032c; not re-run).

**T1b and T11 (R2 scripts, measured, no pass rule).** `r2/t1b-summary.tsv`, `r2/t1b-seeds.tsv` (840 rows, no
error), `r2/t11.tsv`. First seed 8.8 minutes; the 20 lipid seeds and the metab seed 83.1 minutes in all; seed
172101 on metab gave the lipid tables exactly.

| T1b row (lipid, 20 seeds) | seeds with a call change at q < 0.05 | seeds changed (call, k or share) | largest change in t | rule of 6.4 |
|---|---|---|---|---|
| `scale` | 20 | 20 | 7.46 | own ticket (D6) |
| legacy RUV `automatic` | 20 | 20 | 210.6 | own ticket (D8) |
| legacy RUV `manual` | 20 | 20 | 10.03 | own ticket (D8) |
| `groups_all_features` (gate on) | 0 | 0 | 8.9e-13 | measured |
| `first_pass_nulls` (gate on) | 0 | 0 | 2.4e-12 | measured |

T11 (default route, offset 1 against offset 0): T7's concentration fixture: largest change of a log2 value 6.55, of a
logFC 29.7, of t 1,680; 170 calls change (165 at offset 1, 16 at offset 0); legacy RUV k 7 to 2, control share 17 % to
1 %. F172: 6.90, 196.7, 8.0; 37 calls change (0 to 37); k 1 to 2, share 12 % to 1 %. Lipid and metab identical. These
sizes include the legacy RUV route's own unit dependence (T1b).

**7(b) re-dump (R3 script).** `r3/redump-diff.txt`: the lipid and metab tab HTML digests; `norm-log_offset` removed
from the inputs; `norm-log2_recorded_offset` added to the outputs; and, beyond the ticket's list, the module
scenarios' frozen captures, which run at the rendered defaults and so at plain log2: the outputs digest (new output),
the step log (start line without the offset, the values-set-to-missing line), states, provenance, normData and data
digests, the export's "Log2 offset: 0" and its rds and object digests. Not re-run after pasting.

**7(c) control (C1, C2).** Two mechanism failures, nothing committed: C1 the script stopped before any pass; C2
(correction) checked the fixtures, compiled both trees, ran `lipid_canonical` on base (child SHA a60c4e9a) and on the
head at the default (1.9 minutes each), and the head-at-offset-1 child failed: its wrapper read `omicsUiInputDefaults`
from the package namespace, but it is a test helper (`tests/testthat/helper-workflow-omics-residency.R:67`). Stop
rule of 7(c): stopped for the owner. The one-line fix (`packaged <- env$omicsUiInputDefaults`) is on the branch, not
run. `control-mechanism/`.

**Reference rows.** `k156-gui-entry-failures.tsv`: 34 rows (lipid 4, metab 30) read from run 44953187101's preserved
`failures.txt` (SHA-256 1b918062...bde), in G1(r)'s row form (lane, step, class, fields, keys, source, SHA-256);
metabolomics summary steps classed `metab_fixture_4` by the lane-range rule; testthat truncated the integration
artifact names in the summary steps, flagged in `names_list_truncated`.

## F1b (consumers of post-log values)

| Site | Class | Done |
|---|---|---|
| filtering figures' detected and missing counts (`countLipidsPerSample`, `calculateLipidMissingness`, `countMetabolitesPerSample`, `calculateMissingness`) | changes a displayed number on log2 states (workbooks call them after the log step) | fixed: they take `logged` from the object's `args$log_transformed` (T7, T8) |
| store missing percent (`omicsAssayMissingPercent`), finalise tab statistics (`omicQcAssayStatsFrame`) | read only states before the log step (QC tab) | left (first form changed them; restored) |
| CV helpers (`Intensity > 0`), `log2(pmax(x, 1))` QC plots (4 sites) | display only, on intensities in the app; in the workbooks' log-state figures CV of log2 values is unit-dependent under either offset | listed, not fixed (not one line) |
| Glimma volcano's fill of missing cells with 0 | display only | the fill stays; the page says where missing cells are drawn when the assay has any (T8 on a stand-in builder; the real widget is not built on T7's fixture, before and after: "arguments imply differing number of rows: 1, 0", recorded) |
| `scale` (`normalizeMedianAbsValues`), RUVIIIC's zero fill | measured by T1b | warning shipped (D6); own tickets by the rule |
| workbooks' `updateMetaboliteFiltering` on the log state | display | covered by the first row |
| DA not-tested count | displayed number | added to the DA summary (count of rows with no fold change; see divergences) |

## F2 (report only, handed to KIT-177)

`R/func_prot_s4_missingness.R:88, 382`; `R/utils_workflow_dia_state.R:1410, 1435` (TMT and label-free design step,
`log2(.x + 1)`); `R/func_pept_s4_missingness.R:193, 204`; `R/func_pept_s4_norm_methods.R:724`. No new site found.
F1 re-searched at the head: every reader of `log_transform_offset` and `input$log_offset` is in section 8's table.
F3: answered (no regeneration).

## Divergences from the ticket text (for the audit)

- The DA summary's not-tested count counts result rows with no fold change in the selected contrast and assay
  (limma's untested rows), not `findOmicsRuvUnobservedGroups` itself: the tab holds only the long table, and adding a
  stored field would move stored products. Under `~ 0 + group` the two agree (T7 checks it against the test's own count).
- T1b and T11 ran in the fix commit's request (R2), not in the failing-first request: they measure head code.
- The re-dump also moved the module scenarios' frozen captures (above).
- T1b's new designs: `groups_all_features` and `first_pass_nulls` (k rule `num_sv_be`); `technical_replicates`
  cannot run on F172 (no technical replicates), `supplied_controls` needs supplied controls.

## Acceptance criteria

| Criterion | State |
|---|---|
| failing first, then green (T1a to T10) | met for T1a, T2, T3 (after its test fix), T4, T5, T6, T7, T8, T9, T10 (R1, R1b, R2, R3) |
| T1a numerical findings | none (no excess) |
| T1b table and rule, to the user | delivered here; rule: own tickets for `scale` and legacy RUV |
| F1b list with tests | above |
| `scale` warning (D6) | shipped, T10 |
| F2, F1, F3 | above |
| base-saved project reopens and replays | met (T9, R2) |
| control of 7(c) and re-freeze | **not met**: stopped at the mechanism rule after its correction pass |
| full suite base against head with check | not run |
| gate-off script; rendered-panel diff | not run; the panel diff is in `r3/redump-diff.txt` (HTML digests and input/output lists) |
| bundle lane | not run |
| NEWS, man, reports, texts, wording list | done (wording in `omicGuiLog2Texts()`, `R/utils_omic_gui_config.R`) |
| ledger row | written (`VAL-<merge>`) |

## Not checked

- The GUI-entry control and re-freeze; the four lipid and metabolomics GUI-entry oracles are not re-frozen.
- The full suite and R CMD check base against head; the bundle lane; the gate-off script.
- `test-omic-norm-module.R` after the pasted re-dump; `test-metab-qc-shared.R` after its stub fix.
- Whether the recorded-offset output's `<ui>` capture equals the value computed in the head child (head1 never ran).
