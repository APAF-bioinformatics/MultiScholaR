# Merge train 15 of 2026-11-26: KIT-190, KIT-204, KIT-217

**This train is STACKED on merge train 14.** It was branched from `origin/merge-train-20261126g` at its tip
`e69b411ea1d4a61abe5d20cfe5cc30e9106d78f5` (train 14's declaration commit; train 14's code head d1018f6b: KIT-248),
which is stacked on train 13 (KIT-254, f8c871f3), which is stacked on train 12 (landed with staging bbbfdc1a). Staging
was bbbfdc1a (trains 10, 11 and 12 landed; ledger through VAL-297). Trains 13 and 14 were not landed and train 14's
run (45676428221) was out when this train was built. **Train 15 lands only if trains 13 and 14 land with their code
unchanged.** If either is rebuilt or its code changes, or if train 14's run is not green, this train is rebuilt and its
check repeated; this train's evidence stands on e69b411e (base) and the head named below.

## Train

Base: `e69b411ea1d4a61abe5d20cfe5cc30e9106d78f5` (train 14's tip, fetched before building).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-190 | `origin/kit-190-ruv-control-selection` acc97cb0 | 588c960c | `VAL-<train15-1>` (the ticket), `VAL-<train15-1b>` (KIT-094 snapshot, binding-steps), `VAL-<train15-1c>` (declarations snapshot); `VAL-<train15-1d>` reserved for a further coordinator row, not written here | AUDIT-KIT-190-impl-20261126, closure check (e) |
| 2 | KIT-204 | `origin/kit-204` 7e947ce2 | ab5ba116 | `VAL-<train15-2>` | AUDIT-KIT-204-impl-20261125, closure check |
| 3 | KIT-217 | `origin/kit-217` 9d2f9ed0 | 4308c5ae | `VAL-<train15-3>` | AUDIT-KIT-217-impl-20261125, closure check |

Head (bound by the run): `4308c5aee` (full hash in the request). This file, the predictions table and the records
commit after the run are under `docs/qa/` only.

### Conflicts met and their resolutions

No conflict under `R/`, `tests/` or `inst/` except the one authorised resolution in the replay-marker list.

- **KIT-190** (conflict: `docs/methods-decision-ledger.md` only). The branch carries staging bbbfdc1a, so the merge also
  brings staging's landed rows VAL-283 to VAL-297 and the base's rows are train 12's placeholders
  `VAL-<train12-1>`, `VAL-<train12-2>` (landed as VAL-295 and VAL-297, with VAL-296 the correction row), then train
  13's `VAL-<train13-1>`, `VAL-<train13-2>`, train 14's `VAL-<train14-1>`. Resolved as trains 12 and 14 resolved the same
  situation: staging's landed rows VAL-283 to VAL-297 (VAL-295 and VAL-297 replace the placeholder copies, so the KIT-249
  and KIT-236 rows are not listed twice; the landed rows are the placeholders' text with the landing details), then the
  base's `VAL-<train13-1>`, `VAL-<train13-2>`, `VAL-<train14-1>` unchanged, then KIT-190's three rows renamed
  `VAL-<train15-1>`, `VAL-<train15-1b>`, `VAL-<train15-1c>`. This differs from a literal "keep both" only in not
  duplicating the two landed rows. `NEWS.md`, `R/utils_workflow_analysis_record.R` (the marker list, KIT-190's two
  rows after `treat_moderation`, with the roxygen line), `R/utils_workflow_omics_declarations.R`
  (`omicsRuvOperation`'s `fixed`), `R/utils_workflow_omics_actions.R` and `R/utils_workflow_saved_exports.R` merged
  automatically.
- **KIT-204** (conflict: ledger only; `NEWS.md` and `R/utils_workflow_dia_da.R` automatic). Base rows kept, then
  KIT-204's row as `VAL-<train15-2>`.
- **KIT-217** (conflicts: ledger; `R/utils_workflow_analysis_record.R`). `R/utils_workflow_omics_declarations.R` and
  `R/utils_workflow_operations.R` merged automatically (the adjacent insertions are all present, the base's first).
  - Ledger: the base's rows, then KIT-217's row as `VAL-<train15-3>`.
  - Marker list (authorised resolution (b)). Conflict hunk, base side (KIT-190's two rows closing the list with `))`)
    against KIT-217's one row closing with `))`: kept all, KIT-190's `ruv_k_bound` row now ends `),` and KIT-217's
    row follows it and closes the list once:

    ```r
        ruv_k_bound = list(option = "MultiScholaR.ruv_k_bound", current = "sample_size", earlier = "none"),
        # KIT-217 (VAL-173 (2); owner 2026-11-25): ... ("group_counts"); records made before did not count missing
        # values ("not_counted").
        intensity_filter_missing = list(option = "MultiScholaR.intensity_filter_missing", current = "group_counts",
            earlier = "not_counted"))
    ```

    The resulting order is `missing_value_comparison`, `enrichment_background`, `intensity_filter_samples`,
    `duplicate_resolution_samples`, `ruv_feature_levels`, `treat_moderation`, `ruv_control_grouping`, `ruv_k_bound`,
    `intensity_filter_missing`. The roxygen line (auto-merged) names KIT-190's two rows and KIT-217's. The functions that
    follow the list (`workflowTreatModerationCurrent`, `workflowTreatModerationRecord`, `workflowDaReplayRuleOptions`,
    `workflowRecordIntensityFilterEarlier`, `workflowReplayRule`) are all kept, each with its own closing brace.
    Braces checked by reading the region and by counting `{ }` and `( )` outside strings and comments: the merged
    file has the same imbalance as the base, KIT-190's and KIT-217's versions of that file (all 1 and 2, a counting
    artefact of the regex literals, identical on every side); `R/utils_workflow_omics_declarations.R`: the merged file
    equals KIT-217's count (-1 and 3), KIT-217's block having been added to the base's count of -2 and 2 by itself, and
    KIT-190's `fixed = list(ruv_feature_levels = "restored", ruv_control_grouping = "kept", ruv_k_bound =
    "sample_size")` is intact.
- `DESCRIPTION` and `NAMESPACE` unchanged.

### Path lengths

Checked by script before sending: 234 paths are new since the base; none exceeds 100 bytes with the `MultiScholaR/`
prefix; the longest is exactly 100 (`docs/qa/workflow-execution/kit/kit-190/bound-lc/re-selftest/selftest-effect_shifted.tsv`).

## Declaration (written before the run is sent)

**Purpose.** The three members' merge checks on the merged tree, once, on the train, stacked on train 14.

**Scope.** 166 files, deduplicated, base and head, R CMD check on both sides:
- KIT-190: the 106 files of its R2 (`k190-r2-audit-687ba6f3`, run 45651851022) plus the staging and train files its
  register (e) names: `test-kit257-norm-own-files.R`, `test-workflow-saved-exports.R`, `test-kit236-report-sentences.R`,
  `test-kit235-lipid-metab-report-parser.R`, `test-kit227-report-ruv-sentences.R`, `test-kit254-omics-correlation-skip.R`,
  `test-kit248-correlation-skip-record.R`, `test-omics-norm-state-guards.R`, `test-workflow-operation-registry.R`
  (the new file is `test-ruv-control-anova-grouping.R`; the brief's `test-kit190-*.R` is that file);
- KIT-204: the 37 files of its R2d (`kit204-r2d-d792847f`) plus `test-kit249-prot-da-qvalue-fallback.R` and
  `test-kit257-norm-own-files.R`;
- KIT-217: the 57 files of its run 45666403751 (`kit217-mt-13742d43`);
- overlaps: 6 files between KIT-190 and KIT-204, 30 between KIT-190 and KIT-217, 6 between KIT-204 and KIT-217.

**Run settings.** `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` and `MULTISCHOLAR_KIT105_W3_DUMP=1` on both sides;
`extra_packages` limpa, mixOmics (train 14's list) **and sva**: the new KIT-190 file (T7 (v) and a second test) and
four other files in scope `skip_if_not_installed("sva")`, and the declared 139/0 holds only with sva (KIT-190's R2
listed it); this is a deliberate addition to train 14's list, the same on both sides, for the same reason train 14 added
limpa and mixOmics (a skip would be an instrument defect, not a result). Reference BLAS; 8 shards per side; per-file kill
limits 120 minutes for the four GUI-entry files and `test-prot-norm-core.R`, 90 for `test-workflow-saved-da.R` (not
thresholds), 45 default.

**Ceiling.** 2,200 job-minutes (146.7 core-hours = 2,200 x 4 / 60): KIT-190's R2 ceiling 700 (its request
`k190-r2-audit-687ba6f3`: 106 files x 4 shards per side, check) + KIT-204's R2d ceiling 600 (`kit204-r2d-d792847f`) +
KIT-217's merge-train comparison ceiling 900 (`kit217-mt-13742d43`, run 45666403751) = 700 + 600 + 900. The overlap
(42 of the files shared) makes the union cheaper than the sum.

**Base values.** Per file in `merge-train-20261126h-predictions.tsv`: the declared head value of the latest train whose
scope holds the file (train 14, then 13, then the run tables of train 12 and 11), else the base column of the members' own
runs (KIT-204 R2d, KIT-217's run, KIT-190's R2, in that order of recency), else `same` with the count read from the run.
Four files are on no record and read from the run: `test-module-ci-trx-qc.R`, `test-workflow-omics-owner.R`,
`test-workflow-saved-exports.R`, `test-kit172-log2-rule.R`. Counts from older bases may read differently in files that
need `sva`, `limpa` or `mixOmics` (this run installs all three); the verdict is base against head in this run.

### Predictions (the three registers, combined per shared file)

- `test-ruv-control-anova-grouping.R` (KIT-190): only on head, 139/0/0/0 (0 skipped needs sva).
- `test-kit204-prot-da-not-tested.R` (KIT-204): only on head, 81/0/0/0.
- `test-kit217-intensity-filter-missing.R` (KIT-217): only on head, 344/0/0/0.
- `test-workflow-omics-ruv-cancor-once.R` (KIT-190): on THIS base, which lacks KIT-190's test edit, base fails at `:73`
  (14 passed, 2 failed) and head 16/0 (verdict `fixed`). Relative to a base that held the edited test it would be `same`
  16/0.
- `test-omic-norm-module.R` (KIT-190): `same`, 639/0 on both sides: the base test holds the old values and the base
  code produces them; the head test holds the 12 re-dumped provenance values (the 12 one-line hunks of its branch) and
  the head code produces them. On a base that held the old values but head's code, the 12 would fail; that pairing does
  not occur here. Risk: a move of those 12 values by KIT-254 or KIT-248 in the base would show as failures on head.
- `test-prot-norm-core.R` (KIT-190; trains 11 to 14 in the base): 1548/306 on both sides, `message-changed`, only the
  KIT-190 labels (the `lipid_ruv` and `metab_ruv` states); train 14's 12 labels are on both sides. Any other label is a
  finding.
- `test-workflow-omics-binding-steps.R` (KIT-190 and KIT-217): base 58/1, head 57/2, `message-changed`. KIT-190 moves
  messages only (2 parameters, 4 outcomes, 58/1 on both); KIT-217 moves 58/1 to 57/2 (the filter step's description and
  parameters in `records`, two more `report_sources` names). The count is KIT-217's.
- `test-workflow-omics-declarations.R` (KIT-190 and KIT-217): base 45/23, head 43/25, `message-changed`. Each ticket alone
  predicts 45/23 to 43/25; both newly failing expectations are the same two (`current[[omic]]$values`, test line 51, one
  per omic: KIT-190 by the RUV step's five fields, KIT-217 by the filter step's value gaining two entries), so the
  combination is 43/25, not 41/27. If it reads 41/27, the two tickets' new failures are distinct expectations: a
  failure of this combination, neither ticket's declared count.
- `test-workflow-operation-registry.R` (KIT-217; KIT-248 in the base): base 1099/2, head 1107/2 (+4 in the KIT-039
  labels loop and +4 in the KIT-027 defaults loop), `count-changed`; the same two `integration_mofa_model` failures on
  both sides.
- `test-omic-da-module.R` (KIT-204): base 1289/80, head 1288/81, `message-changed`: one more failing expectation
  (`prot_main` `volcano_unticked` values) and 12 message changes.
- `test-metab-qc-shared.R` (KIT-217): 20/0 to 23/0, `count-changed`.
- `test-omic-qc-module.R` (KIT-217): 160/0 to 156/4, `new-failure` (test (b)'s four QC-page digests).
- `test-gui-entry-proteomics.R` (KIT-217): 168/25 to 92/101, `new-failure`, `states` only beyond base's own failures.
- `test-gui-entry-lipidomics.R` 35/0 to 19/16 and `test-gui-entry-metabolomics.R` 105/0 to 57/48 (KIT-217),
  `new-failure`, `states` only, from `qc_intensity` to `summary_export` (16 steps per lane). KIT-190 predicts these
  files `same` alone and adds no move.
- `test-gui-entry-transcriptomics.R` `same` 16/13; `test-kit214-sample-columns.R` `same` 176/0.
- `test-workflow-figure-paths.R`: `same`, except the known flake of staging (KIT-204's C01): a failure at `:99`,
  `:100` or `:471` whose raw message quotes a random temporary folder whose name contains "259" is that flake, not a
  result; it is recorded as such if it occurs.
- Every other file `same` at the base's counts in the predictions table.
- R CMD check: findings identical on both sides, nothing only on head; the portable-file-names NOTE the same count on
  both sides.

### Decision and attempts

- Every prediction holds: TRAIN GREEN, landing only after trains 13 and 14 land with their code unchanged.
- Stop rule: any move OUTSIDE the union of the three registers' declared moves (the list above) is a finding: it is
  reported with its files and counts, not fixed and not re-run here. A move of a declared kind with a different count is
  recorded as a miss.
- One run; one re-send only for an instrument defect with no result judged. An HTTP 500/502 at dispatch with no run
  created: started once by hand with the README's `gh workflow run` command, recorded.

## Result

Run 45681065939 (request `merge-train-20261126h-4308c5ae`, run-tests commit 0e351e1c, dispatched normally), head
4308c5ae against train 14's tip e69b411e, sent once; conclusion success; about 947.9 job-minutes of the 2,200 ceiling
(tests 918.7, check 28.3, plan and results 0.8; summed from the API's job start and end times). Tables in
`merge-train-20261126h/` (`compare.tsv`, `tests.tsv`, `tests-base.tsv`, `check-compare.txt`, `check-summary.txt`,
`meta.txt`, `compare-failures.txt.gz`).

Observed, 166 files: 152 `same`, 14 not `same`. Compared against the declaration, file by file:

| file | base p/f/e/s | head p/f/e/s | verdict | prediction |
|---|---|---|---|---|
| test-ruv-control-anova-grouping.R | absent | 139/0/0/0 | only-on-head | KIT-190, held |
| test-kit204-prot-da-not-tested.R | absent | 81/0/0/0 | only-on-head | KIT-204, held |
| test-kit217-intensity-filter-missing.R | absent | 344/0/0/0 | only-on-head | KIT-217, held |
| test-workflow-omics-ruv-cancor-once.R | 16/0/0/0 | 16/0/0/0 | same | **MISS**: declared base 14/2 failing at `:73`, verdict `fixed`; base passes, so `same` (the register's own statement, "the row-order rule passes on staging", was right; the brief's was not) |
| test-omic-norm-module.R | 639/0/0/0 | 639/0/0/0 | same | KIT-190, held (the 12 re-dumped values pass; train 14's base values untouched) |
| test-prot-norm-core.R | 1548/306/0/0 | 1548/306/0/0 | message-changed | KIT-190, held: 2 failures change message, the lipid and metab RUV declaration test; 303 identical |
| test-workflow-omics-binding-steps.R | 58/1/0/0 | 57/2/0/0 | message-changed | KIT-190 + KIT-217, held |
| test-workflow-omics-declarations.R | 45/23/0/0 | 43/25/0/0 | message-changed | combined prediction 43/25 held (not 41/27: the two tickets' new failures are the same two expectations) |
| test-workflow-operation-registry.R | 1099/2/0/0 | 1107/2/0/0 | count-changed | KIT-217, held; the same two failures, identical on both sides |
| test-omic-da-module.R | 1289/80/0/0 | 1288/81/0/0 | message-changed | KIT-204, held (one more failing expectation, 12 message changes in a 13-block difference) |
| test-metab-qc-shared.R | 20/0/0/0 | 23/0/0/0 | count-changed | KIT-217, held |
| test-omic-qc-module.R | 160/0/0/0 | 156/4/0/0 | new-failure | KIT-217, held |
| test-gui-entry-proteomics.R | 168/25/0/0 | 92/101/0/0 | new-failure | KIT-217, held |
| test-gui-entry-lipidomics.R | 35/0/0/0 | 19/16/0/0 | new-failure | KIT-217, held |
| test-gui-entry-metabolomics.R | 105/0/0/0 | 57/48/0/0 | new-failure | KIT-217, held |
| **test-kit254-omics-correlation-skip.R** | 379/0/0/0 | 376/3/0/0 | new-failure | **UNDECLARED MOVE (finding)**: predicted `same` |

`test-workflow-figure-paths.R`: 0 failures on both sides (the "259" flake did not occur). Every other `same` file has the
base counts of the predictions table (0 differences by script over `compare.tsv`); the four files declared "read from
the run" read `test-module-ci-trx-qc.R` 38/0/1/0, `test-workflow-omics-owner.R` 78/1/0/0,
`test-kit172-log2-rule.R` 212/0/0/0, `test-workflow-saved-exports.R` (see `compare.tsv`), all `same`. No skip in the new
files (sva present). `test-kit214-sample-columns.R` 176/0 and `test-gui-entry-transcriptomics.R` 16/13 `same`.

**The finding: `test-kit254-omics-correlation-skip.R`** (KIT-254 in the base, KIT-190 in the head). Three failures, all
`lipid_ruv` (the E7 tests with RUV): the run's sample list after the correlation filter at `t_hi` reads
`A_1 A_2 A_3 A_4 B_2 B_3 B_4` where the test expects `A_2 A_3 A_4 A_5 A_6 B_1 B_2 B_3 B_4` (precondition "Apply(t_hi)
lipid_ruv", and the re-apply check). Inferred, not shown: KIT-190 chooses other RUV control features, so the
RUV-corrected values on KIT-254's lipid fixture change and the correlation filter keeps other samples; the expectations
were written against the earlier control rule. The metabolomics RUV case passes. Neither ticket's register lists this
file as moving (KIT-190's register named it as a scope file and predicted `same`; it was absent from its R2, which ran on
staging without KIT-254). Not fixed and not re-run here (the stop rule). Someone with the authority over the test has
to decide: re-pin its expected samples under KIT-190's rule or order the merges differently.

**R CMD check:** 0 ERRORs, 3 WARNINGs, 3 NOTEs on both sides; all 6 findings identical, none only on head or only on
base (held).

Verdicts:
- KIT-204: holds exactly.
- KIT-217: holds exactly (344/0; 168/25 -> 92/101; lipid 16 of 35 and metab 48 of 105; test (b)'s four; registry +8).
- KIT-190: holds except (1) the cancor-once base count (a miss in the declaration, favourable: base 16/0 not 14/2) and
  (2) the undeclared `test-kit254-omics-correlation-skip.R` failures, an interaction with train 13 (KIT-254).
- Check: identical.

Train verdict: NOT GREEN as declared, on those two items; every other prediction held. The coordinator decides; the
train is not rebuilt or re-run here. It lands, if at all, only after trains 13 and 14 land with their code unchanged.

Not checked: trains 13 and 14's own results (train 14's run 45676428221 was out when this one was sent; its result was
not read); any file outside the 166; the full suite; the digest values (read from the run).

## After the focused audit: declaration of resolution (a) (written before any run of it)

Source: `tickets/kit/audits/AUDIT-TRAIN15-kit254-e7-20261126.md` (verdict: no product defect; KIT-190's decided control
rule, VAL-110 and VAL-126, moved the `lipid_ruv` expected sample lists of
`tests/testthat/test-kit254-omics-correlation-skip.R`) and the coordinator's ruling of 2026-11-26: resolution (a), a
DECLARED TEST-ONLY edit of that one file. Corrections to the result section above, from the audit: the "7 against 9"
lists are a waldo diff window; the full lists are 9 samples on the head against 11 expected, of 12 (at 0.960024 the head
removes A_5, A_6 and B_1 and keeps A_1); and the file has no metabolomics RUV case ("the metabolomics RUV case passes" is
wrong: `metab` runs without RUV). The cancor-once miss (`test-workflow-omics-ruv-cancor-once.R` 16/0 `same`, not 14/2
`fixed`) was the coordinator's brief error, not a failure of the code.

**Step 1: the G0c run (instrument).** `docs/qa/workflow-execution/kit/kit-254/g0b/g0b-ruv-thresholds.R`, UNCHANGED, as a
scripts-only request: `package_commit` 4308c5aeed00b5551e88be9983032f0c7bacbaf5, `base_commit` `none`, no tests, no
extra packages, no env, `args` and `out` `kit254-g0c`, `out_pattern` `\.(txt|tsv)$`, timeout 25 minutes, ceiling 30
job-minutes (G0b 45628705535 took 12). Purpose: the head's `t_hi`, `t_lo` and the C3 and II0 removed sets of
`lipid_ruv` by G0's rule (`k254Choose`), and the head's best-partner table. Predictions (the audit's, not read from the
run): the script runs unchanged; `t_hi` = 0.95 with C3 removing A_5, A_6, B_1 (SHOWN bound in the audit); `t_lo` and II0
are not predicted, but G0's rule gives a `t_lo` that removes a strict subset of C3's set. Stop rules: if the script does
not run unchanged on this head, that is an instrument question: stop and tell the coordinator, do not edit the script;
if G0's own rule finds no `t_lo` removing a strict subset (its stop rule), stop and tell the coordinator.

**Step 2: the test-only edit.** From that run's output only: `kit254T$lipid_ruv` (`t_hi`, `t_lo`),
`kit254HiRuns$lipid_ruv = setdiff(kit254All, <C3 removed>)`, `kit254LoRuns$lipid_ruv = setdiff(kit254All, <II0
removed>)` and the comments around them; `t_lo` is re-derived, not re-pinned to the head's observed list. Its own
tests-only commit, with the G0c table under `kit-254/g0c/` and a dated addendum in `kit-254/outcome.md`. Prediction: the
file 379 passed, 0 failed on the edited head.

**Step 3: one test-file run.** `test-kit254-omics-correlation-skip.R` alone, base e69b411e and the edited head, both
sides, `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` and `MULTISCHOLAR_KIT105_W3_DUMP=1`, limpa, mixOmics, sva, reference BLAS,
ceiling 60 job-minutes. Prediction: base 379/0, head 379/0, `same`. Any other result is recorded, then stop; no second
edit.

**The other 165 verdicts and R CMD check stand** from run 45681065939 (the edit touches one test file; `kit254T`,
`kit254HiRuns` and `kit254LoRuns` are referenced nowhere else, the audit's argument; each file runs in its own process
and the check runs with `--no-tests`). Trains 16, 17 and 18 are not rebuilt or touched.

## After the focused audit: result

**G0c (run 45704341824, request `kit254-g0c-4308c5ae`, run-tests commit after 0e351e1c; 13.5 job-minutes of 30).** The
G0b script ran unchanged on 4308c5ae (exit 0, 91 s). Tables in `kit-254/g0c/`. G0's rule: `t_hi` = 0.8, `t_lo` = 0.7,
stop rule 2 false. Scan (removed samples): 0, 0.3, 0.5, 0.7 none; 0.8 B_1; 0.9 B_1; 0.95 A_5, A_6, B_1; 0.99 all 12.
C3 Apply(0.8) removes B_1; II0 Apply(0.7) removes none. The audit's SHOWN bound (0.95 removes A_5, A_6, B_1) is confirmed;
the script's rule picks the lowest threshold that removes a sample, 0.8, and a `t_lo` that removes none. (The audit
guessed `t_hi` 0.95; that was a guess it labelled as such, and it did not hold.)

**The edit (commit 3865e5cf, tests-only):** `kit254T$lipid_ruv` = c(hi = 0.8, lo = 0.7); `kit254HiRuns$lipid_ruv` =
`setdiff(kit254All, "B_1")`; `kit254LoRuns$lipid_ruv` = `kit254All` (unchanged value); comments cite run 45704341824.
New tip of the train branch for landing: the records commit after 3865e5cf (3865e5cf is the code head).

**File run (run 45707353015, request `kit254-e7-train15-3865e5cf`; 39.7 job-minutes of 60):** head 3865e5cf against base
e69b411e, `test-kit254-omics-correlation-skip.R` alone, env, limpa, mixOmics, sva: base 379/0/0/0, head 379/0/0/0,
verdict `same`, no failing expectation. Prediction held. Tables in `e7/`.

**Final verdict table:** the 165 other files and R CMD check stand from run 45681065939 (the edit touches one test file;
its constants are referenced nowhere else; each file runs in its own process; the check runs with `--no-tests`;
AUDIT-TRAIN15-kit254-e7-20261126, section (4)). `test-kit254-omics-correlation-skip.R` is now `same` 379/0. Remaining
deviation: the cancor-once declaration (the coordinator's brief error; 16/0 `same`). Train 15 verdict: GREEN after
resolution (a). Job-minutes in total: 947.9 + 13.5 + 39.7 = 1,001.1.
