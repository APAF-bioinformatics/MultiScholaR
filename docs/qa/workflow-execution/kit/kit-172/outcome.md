# KIT-172 outcome: plain log2 for lipidomics and metabolomics

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

## The GUI-entry control (7(c)), after the owner's rule VAL-54

C1 (45136470598), C2 (45137945090) and C3 (45141094821) were instrument failures (no comparison produced): C1 the
`compareLane` extraction; C2 the head-at-offset-1 wrapper looked up `omicsUiInputDefaults` in the package namespace (it is
a test helper); C3 the same wrapper stopped on a UI with no inputs (`names()` NULL). Fixes: 9931beea, e485208b (also:
the expected recorded-offset capture read through `renderUI` and the lanes' own output reader), e623d6e2. No expected
value, threshold, declared comparison or stop rule changed. Evidence: `control-mechanism/`.

C4 (45142934736, e623d6e2) produced the comparison (`control-run4/`: `gui-entry-control.tsv`, `gui-entry-digests.tsv`,
`declared-keys.tsv`, the three passes' records; runner CPU AMD EPYC 7763; nothing withheld). Result:

- Step 1 (base against the frozen file): **stop.** The metabolomics lanes' differing steps equal the reference list
  (10 each, `metab_fixture_4`). `lipid_canonical` differs at 6 steps, the list has 4: `da` and `summary_refusal` also
  differ, both only in the `artifacts` key `<store entries>`. Observed: at those steps the base's and the frozen file's
  `<store entries>` hold the same paths with the same counts; only their order differs (the frozen file has
  `lipid-da-fit-<ID>` before `lipid_between_sample_normalisation-<ID>`, the base after the `lipid_*` entries).
  Inferred, not tested: the oracles were written inside `test_that`, where testthat sets the C collation (`-` sorts
  before `_`), and the control's children run the lanes outside `test_that` under the runner's locale, so `table()`
  orders the paths differently wherever such names occur. If so, the re-frozen files written by the head-at-default
  child carry the same locale ordering and would not match a later run under testthat.
- Step 2 (head with `log_offset = 1` against base): conforms. Differences only on the two declared keys, at every
  step; each observed value equals its expected value (`declared-keys.tsv`: 136 of 136 equal; the recorded-offset
  sentence's capture `<ui> 916f8a8c...` as computed).
- Step 3 (head at the default against head at offset 1): conforms. No difference before `norm`; every difference
  from `norm` on is `plain_log2 (KIT-172)`.
- Under 7(c) a differing step list stops the ticket. The re-frozen files are not committed. Not re-run (a result).

**C5 (45148367369, 8fb3c215), the one run VAL-55 allows.** The only change: each child process sets
`LC_COLLATE` to C at its start, as testthat does inside `test_that` (testthat 3.3.2, `local_reproducible_output`:
`withr::local_collate("C")`, read in testthat's source at tag v3.3.2). Result: **pass, steps 1 to 3.** Step 1: the
base's differing steps equal the reference list lane by lane (`lipid_canonical` 4, `integration_artifacts`;
`metab_lc`, `metab_gc`, `metab_combined` 10 each, `metab_fixture_4`), which shows C4's two extra steps were the key
order of `<store entries>` under the children's locale. Step 2: differences only on the two declared keys; 160 of 160
declared values equal their expected value. Step 3: nothing differs before `norm`; every difference from `norm` on is
`plain_log2 (KIT-172)`. Nothing withheld; `scripts.tsv` ok. The four re-frozen oracles (source SHA 8fb3c215, written
under C collation by `guiEntryExpectLane`) are committed in a4f60900. The control record is `gui-entry-control.tsv`
(first line: runner CPU AMD EPYC 9V45; changed step `norm` for every lane), with `gui-entry-digests.tsv`,
`declared-keys.tsv`, `gui-entry-control-status.tsv` and the passes' records in `control-run5/`.

What the re-frozen files absorb besides this ticket's change: on `lipid_canonical` the integration-input artifact
names and the `<store entries>` counts of the summary steps (KIT-147 to KIT-153); on the three metabolomics lanes
`metab_fixture_4` from `norm` on (with the integration artifact names in their summary steps). They are
intensity-scale (fixture intensities 2.4e3 and above) and say nothing about values below 1.

## Corrections to the records

- The commit message of 5250032c cites "run 45136468961 family" for R3: the R3 run is **45136474750**, and the commit
  message of dd5e406c cites "run 45125620735 family" for R2: R2 is **45123630549** (`runs.tsv` is correct).
- Two rows ran past their own ceilings: R2's four test jobs used 68.2 job-minutes = 4.55 core-hours against the
  targeted row's ceiling of 4 (expected 2); the T1b and T11 scripts job ran 103.3 job-minutes = 6.89 core-hours, inside
  the T1b row's ceiling of 7.5 but above its expected 4.4, because the parallel part took 74 minutes against the
  first seed's projection of 44 (four seed processes sharing the runner's four cores with T11 done).
- The 7(b) re-dump moved, besides the listed interface values, the module scenarios' frozen captures
  (`omicNormFrozen[["lipid"|"metab"]][["skip"|"manual"|"automatic"]]`): they drive the server with the tab's rendered
  defaults, which no longer carry `norm-log_offset`, so they run at offset 0. Moved: the `outputs` digest at every
  step (the new output), the step log (the start line, the values-set-to-missing line), the `states`, `provenance`,
  `normdata` and per-slot `state` and `data` digests from the run on, and in the `skip` scenarios' export the summary
  line "- Log2 offset: 0" with the `rds_values` and `current_s4_object` digests; and the KIT-110 W3 test that compares
  the registered outputs with `omicNormFrozenUi$outputs`. Added to the ticket's 7(b) with a History note.

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
| control of 7(c) and re-freeze | met by C5 (VAL-55): steps 1 to 3 pass; oracles re-frozen (a4f60900); C4's stop recorded |
| full suite base against head with check | not run |
| gate-off script; rendered-panel diff | not run; the panel diff is in `r3/redump-diff.txt` (HTML digests and input/output lists) |
| bundle lane | not run |
| NEWS, man, reports, texts, wording list | done (wording in `omicGuiLog2Texts()`, `R/utils_omic_gui_config.R`) |
| ledger row | written (`VAL-<merge>`) |

## Not checked

- The full suite and R CMD check base against head; the bundle lane; the gate-off script.
- `test-omic-norm-module.R` after the pasted re-dump; `test-metab-qc-shared.R` after its stub fix.
