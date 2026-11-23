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
1 %. F172: 6.90, 196.7, 8.0; 37 calls change (0 to 37); k 1 to 2, share 12 % to 1 %. Lipid and metab identical. Both fixtures have no true difference between the groups (`kit172Values(..., planted = FALSE)`; the concentration fixture is a monotone rescaling of the same draw), so every call at either offset is a false discovery: 165 and 16 on T7's fixture, 0 and 37 on F172 (one draw each). How much comes from the log step and how much from the legacy RUV route is not separated (inferred: mostly the RUV route, whose k and controls move). At plain log2, T1b's `legacy_automatic` row moves fold changes by up to 1,067 log2 units where cells are missing.

On the four GUI-entry lanes (intensities 2.4e3 and above) no call changes (their DA summaries are identical at offset 1 and at plain log2, control run 45148367369); the change at the log step is at most about 6e-4 log2 (a bound from the smallest value, not measured).

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

- B1 used 15.55 core-hours against the bundle row's ceiling of 12; R1, R1b and R2 together used 16.41 against the targeted and failing-first rows' ceilings of 11.5. The total is 111.17 core-hours (1,667.5 job-minutes), inside VAL-68's 'about 112'.
- R5 (45167057104) was beyond the allowance (VAL-68: one diagnosis run; one correction pass), sent without the owner's word; recorded for the owner (VAL-72 records it as a departure). D2 (45225156715) likewise. Neither changed an expected value.
- The 7(b) extension to the module scenarios' captures is ratified by the owner (VAL-72).

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

## Full suite, check, gate-off and bundle lane

**Full suite (F1, 45151310949, head a4f60900 against base a60c4e9a, 14 shards per side, `sva` in both).** 452
files: 441 `same`; 2 `fixed` (`test-gui-entry-lipidomics.R` 4 to 0 and `test-gui-entry-metabolomics.R` 30 to 0
failing expectations: the re-frozen oracles pass under testthat); 4 `only-on-head` (the KIT-172 files, all green); 2
`count-changed` (the rewritten log-transform contract files); 3 `new-failure`:
- `test-metab-03-norm-support-characterization.R`: the `buildNormConfig` pin of offset 1, class (a); rewritten to 0
  (dfb2ad9f), `same` in R4.
- `test-omics-module-globals.R`: the new `<<-` in the two log methods; the counts are kept in an environment
  (dfb2ad9f, 55ff3063), green in R5.
- `test-omic-da-module.R`: 22 failing expectations on both sides become 51 on head. The frozen DA-module captures of
  `lipid_main`, `metab_main`, `metab_no_friendly`, `lipid_managed`, `metab_managed` and the table-filter test (c)
  hold the DA summary text, which now ends with "Not tested (no value in a compared group): N" (the declared change
  of 6.3 item 1); no number in them moves otherwise. Not re-recorded (no exception names this file): a row for the
  plan's section 11.
`test-omic-norm-module.R` (the pasted re-dump) and `test-metab-qc-shared.R` (the stub fix): 639 and 19 passed on head.
R CMD check: base and head 0 ERROR, 3 WARNING, 3 NOTE; one NOTE line only on head, a record path of this ticket over
100 bytes, renamed (dfb2ad9f); R4's check: findings identical on both sides, none only on head. Gate-off identity:
the head's gate-off pages identical to the base's, renders deterministic, the integration module the unavailable
panel (`full/gate-off-verdict.txt`).

Not re-run after F1: the full suite, the GUI-entry lanes and the bundle at dfb2ad9f and 55ff3063. Their code change
(the count tally in an environment in place of `<<-`) does not change any value; R4 and R5 covered the affected files.

**Class (a) rewrites (old and new expectation):**

| File | Old expectation | New expectation |
|---|---|---|
| `test-lipid-09b-s4-log-transform-assays-characterization.R` ("rejects non-positive offsets") | offset 0 refused with "single positive numeric value" | offset 0 accepted; offset -1 refused, class `omics_log2_offset_invalid` |
| `test-metab-02r-log-transform-assays-characterization.R` | redefined generic `offset = 1`; offset 0 refused; source-text test of the clamp and `args` lines | generic `offset = 0`; offset 0 accepted, -1 refused by class; behavioural test: clamp kept at offset 1, cells of 0 or below missing at offset 0, `args` record the offset |
| `test-metab-s4-norm-core-shared.R` | redefined generic `offset = 1`; offset 0 refused | generic `offset = 0`; -1 refused by class |
| `test-module-ci-lipid-norm.R`, `test-module-ci-metab-norm.R` (MCI-021.1, MCI-015.1) | offset 0 refused | -1 refused by class |
| `test-omic-gui-config.R` | `log2$offset` 1; DA `summary` without `not_tested` | 0; `not_tested = TRUE` for lipid and metab |
| `test-metab-03-norm-support-characterization.R` | `buildNormConfig` default offset 1 | 0 |
| `test-metab-qc-shared.R` | stubs of the two QC helpers without `...` | stubs take `...` (the new `logged` argument) |

**Plan section 11 rows (references this ticket moves and does not re-record):**

| File | Oracle or line | Change |
|---|---|---|
| `tests/testthat/test-omic-da-module.R` | the frozen captures (`values` fields) of the lipid and metab DA scenarios `lipid_main`, `metab_main`, `metab_no_friendly`, `lipid_managed`, `metab_managed`, and test (c) | 29 more failing expectations on runners (22 to 51); the DA summary text gains one line, "Not tested (no value in a compared group): N" |
| `tests/testthat/test-e2e-lipidomics-canonical.R`, `test-e2e-metabolomics-lc-gc.R` | the browser lanes | measured on the bundle lane (B1): 88/0 and 190/0 on ubuntu-22.04 and windows-2025, as on the laptop; no row |

**Bundle lane (B1, 45151293709, at a4f60900; reference run 44998302733).** Red, at the same step as the reference
run on every platform. Linux (ubuntu-22.04, ubuntu-24.04): `lipid_canonical`, `metab_lc`, `metab_gc` and
`metab_combined` have no difference from the re-frozen oracles. Windows (windows-2025): the four lanes show the reference run's differences except its integration-artifact and `metab_normalized` ones, which are gone: `session_summary` (2 steps each) and `norm_log` at the same 8 steps (`session_export` to `summary_export`), which the reference run lists as discrete differences (oracle text against a Windows digest) and B1 as digest differences (both sides now digests, the log being longer than 1,000 characters). No state, table or artifact differs.
`omic_worker_da` 49 failing expectations against the reference's 23 (the DA summary's new line, as on runners).

## Declared before sending: the Windows norm_log diagnosis (VAL-68)

- Purpose: read the text behind the Windows-only `norm_log` digest difference of the four lanes (B1). Diagnosis
  only: no package code, oracle or ticket test changes.
- Binding: run branch `run-k172-win` = the ticket head 55ff3063 plus instrumentation only: the suite stage of
  `run-win-k147` (13f755bb: `distributable.yml` stage `suite`, `packaging/verify/run-suite*.R`) and one diagnosis
  file, `tests/testthat/test-k172-normlog-diag.R`, that wraps the lanes' capture function to keep the `norm_log`
  text whole and returns it in failure messages (SUITEFAIL lines of the job log), with the digest the capture would
  write and the re-frozen oracle's. Four lanes under testthat (C collation, as the suite runs files), then
  `lipid_canonical` and `metab_gc` again under the runner's own locale (`LC_ALL ""`), as the bundle's lane step runs
  lanes outside testthat.
- Request: `distributable.yml` dispatched on `run-k172-win`, stage `suite`, platform `windows`, files
  `test-k172-normlog-diag.R`. Linux not run (a separate job; not free).
- Ceiling: one Windows suite job, about 60 minutes (bundle build about 40, lanes about 20), about 4 core-hours; total
  limit about 112 (VAL-68). One mechanics rerun at most if no text comes back.
- Decides: the cause as observed, the smallest fix, and whether any number changes.

## The Windows norm_log diagnosis (VAL-68): observed

D1 (45219858480) returned the four lanes' `norm_log` text at step `norm` on Windows (the bundle's R, installed
package, under testthat and again under the runner's own locale): it equals the re-frozen oracle's text exactly
(866, 866, 823 and 868 characters). The re-frozen oracles hold this log as text up to `correlation` and as a digest
from `session_export` on (the log passes 1,000 characters there); B1's eight differing steps per lane are exactly the
digest steps. D1 had captured the wrong step, so D2 (45225156715, the one rerun) returned the text at
`session_export` for `lipid_canonical` and `metab_gc`.

Observed in D2: the last line, "Exported comprehensive session data to: <path>" (`R/mod_omic_norm.R:629`), carries the
session file's path unmasked, `D:/a/_temp/st/1/RtmpKsbQ1W/gui-entry-<lane>-<HEX>/project/...`, where Linux writes
`<ROOT>/project/...`. The text is ASCII, without carriage returns; every other line equals Linux's. The capture's mask
(`guiEntryMask`, `tests/testthat/helper-gui-entry.R:51-58`) replaces the lane's root and `tempdir()` as
`normalizePath()` gives them, which on Windows is the backslash form, so the forward-slash form of the same path in
this line is left as it is. The digest of the log therefore depends on the Windows runner's temporary path.

Observed in the reference run's Windows log (job 137134081006): `norm_log` differs at 8 steps on each of the four lanes, kind discrete, the runner's value a `<text>` digest against the oracle's text, consistent with the same unmasked path pushing the base's 934 to 939 characters past 1,000. This ticket's line changed the kind of the difference, not its presence. That the reference run's 8 steps are `session_export` to `summary_export` is inferred (the job log keeps no step names).

Numbers: none change. The difference is in the text of a log output only; every state, table and artifact of the four
lanes is equal on Windows at these steps (B1 lists only `norm_log` and the known `session_summary` family).

Smallest fix (not made; diagnosis only, and it is a test helper): in `guiEntryMask`, also replace the forward-slash
forms of the roots and of `tempdir()` (`normalizePath(x, winslash = "/", mustWork = FALSE)`), so the capture masks
the path on Windows as on Linux. It moves no oracle (Linux text already masked). Alternatively the log line could
name the file relative to the project, which would change the line's text and so the oracles.

Diagnosis cost: D1 1.66 and D2 1.47 core-hours; total 111.2 (VAL-68 limit about 112). D2 was beyond the allowance (VAL-68: one diagnosis run; one correction pass), sent without the owner's word; recorded for the owner (VAL-72 records it as a departure).

## Budget

108.0 core-hours used: R1 3.90, R1b 1.04, R2 11.47, R3 2.47, C1 0.58, C2 0.97, C3 0.83, C4 1.46, C5 1.21, F1 64.46,
R4 3.46, R5 0.64, B1 15.55. This is above the owner's limit of 106 by 2.0. The full suite and the bundle lane were sent
together at 23.9 core-hours used, against the ticket's expected 30.4 and 11.2; they took 64.46 (14 shards per side,
33 jobs, longest 80 minutes) and 15.55 (the Windows job 122.6 minutes against the reference's 79). The ticket's
ceiling of 189.5 is not passed.

## F1b (consumers of post-log values)

| Site | Class | Done |
|---|---|---|
| filtering figures' detected and missing counts (`countLipidsPerSample`, `calculateLipidMissingness`, `countMetabolitesPerSample`, `calculateMissingness`) | changes a displayed number on log2 states (workbooks call them after the log step) | fixed: they take `logged` from the object's `args$log_transformed` (T7, T8) |
| store missing percent (`omicsAssayMissingPercent`), finalise tab statistics (`omicQcAssayStatsFrame`) | display only; read the head state, so on a reopened normalised project a log2 value of exactly 0 (a measured 1) is shown as missing | listed, not fixed (plan section 13 candidate) |
| filtering figures on a replayed offset-1 state | display only (workbooks, pre-M3 imports with zeros): a clamped zero, written as log2(1) = 0, now counts as detected | listed; the exact rule would read `logged` as `log_transformed` and offset 0 |
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
| full suite base against head with check | run (F1); every difference classed above; two corrections verified by R4 and R5 |
| gate-off script; rendered-panel diff | gate-off identity TRUE (F1); the panel diff is in `r3/redump-diff.txt` (HTML digests and input/output lists) |
| bundle lane | run (B1): Linux clean for the four lanes; Windows: `norm_log` differs from `session_export` on, diagnosed (D1, D2) as the capture's path mask missing the forward-slash Windows temp path in the export line; no number changes; fix in the test helper proposed, not made |
| NEWS, man, reports, texts, wording list | done; the wording list is below |
| ledger row | written (`VAL-<merge>`) |

## Wording list (VAL-41; the owner edits before release)

| String | File |
|---|---|
| the five texts of `omicGuiLog2Texts()`: log2 section sentence, recorded-offset sentence, values-set-to-missing log line, volcano sentence, `scale` warning | `R/utils_omic_gui_config.R` (`omicGuiLog2Texts`) |
| report Methods sentence ("Intensities were log2-transformed; a value of zero was read as not detected ...") and "... intensities were transformed as log2(x + {offset_val})." | `inst/reports/lipidomics/lipidomics_report.rmd`, `inst/reports/metabolomics/metabolomics_report.rmd` (log transformation text) |
| "* Transform: log2" | `R/func_general_filemgmt_study_parameter_writer.R` (both log-transformation sections) |
| "Not tested (no value in a compared group): %d" | `R/mod_omic_da_panels.R` (`omicDaSummaryText`) |
| "Log2 offset (0 = plain log2)" | `R/utils_workflow_operations.R` (`lipid_log2`, `metab_log2` labels) |
| "Applying log2 transformation" (without the offset) | `R/utils_workflow_omics_actions.R` (`omicsLog2StartText`) |
| the two workbook comments (plain log2; the `scale` sentence) | `inst/workbooks/metabolomics/metabolomics_workflow_starter.rmd`, `metabolomics_workflow_experienced.rmd` |

"non-positive, non-finite" in the step-log line is vocabulary the UI audience rule asks to avoid; the owner's wording.

## Not checked

- The report's log2 sentence and the study-parameter line are not tested; no report was rendered.

- Whether the base's Windows log in the bundle's lane step carried the same unmasked path (only the suite environment was observed).
- The full suite, GUI-entry lanes and bundle at the final commits dfb2ad9f and 55ff3063 (code change without value change; R4, R5 only).
