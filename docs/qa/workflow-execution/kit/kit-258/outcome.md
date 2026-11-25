# KIT-258 outcome

Ticket `tickets/kit/KIT-258.md` r1.2 (Re-check 1 R1-01 to R1-04 applied). Branch `kit-258` from `v0.9-staging` at
**63143aa3** (KIT-214 merged, VAL-256). The owner's questions OD-1 (does an intensity of exactly 0 count as measured)
and OD-2 (the strict ordering) were OPEN during R1 and the first S1. **Answered by the owner 2026-11-26 (VAL-284):
OD-1 (b), a zero is NOT measured, for the count and for the mean, as at import; OD-2 (a), the strict ordering as
written.** The rest of this file (E5, E6, the final fix, the second control S1f, R2) follows those answers.

Labels: OBSERVED = read from a run's returned tables; PREDICTED = written before the run was read.

## Commits

| Commit | What |
|---|---|
| ddbedc3f | tests only: new `test-kit258-duplicate-completeness.R`; declared edits of `test-kit214-sample-columns.R` (`k214Earlier`, `k214OracleKept`), `test-metab-01m-...` and `test-metab-02u-...` source pins. Nothing under `R/` |
| 81fdcf7c | fix as far as VAL-254 ("measured" = not NA), pending OD-1 and OD-2; S1 ran on it |
| ec849390 | merge of `v0.9-staging` 12701142: one conflict, the replay-marker list (every staging row in staging's order, then `duplicate_resolution_rank`); NEWS and the ledger merged without conflict |
| e2842edd | tests for VAL-284: E5 (zeros), E6 (the 12-sample example); declared edits of `k214OracleKept` (a zero is not measured) and the 01m pins. Nothing under `R/` |
| 77a719b8 | final fix: `duplicateResolutionMeasured` / `duplicateResolutionMeasuredMean` (neither NA nor 0) and `keepDuplicateResolutionRow` at the four sites; the exported kernels' help; NEWS; ledger row VAL-<merge> |
| (R2 base) 5dfe2bc8, branch `kit-258-r2-base` | staging 12701142 with the four final test files of e2842edd; nothing under `R/` |

## R1 (tests-only commit ddbedc3f, base none; 1 job, ceiling 45 job-minutes)

Purpose: failing-first. Binds ddbedc3f (R/ = staging 63143aa3). Decides: whether the fix may be written against these
tests (every failing row at its stated base value, every guard passing).

PREDICTED (written before the run was read):
- `test-kit258-duplicate-completeness.R`:
  - E1, E2, E3 (2 omics x 2 variants = 12 tests): each fails on the kept row only, at KIT-214's value: E1 row 1 'A',
    E2 row 2 'A', E3 row 2 'Q'; the one-row-per-id expectation passes.
  - E4 kernels (1 test): fails on id D for both kernels ('A' for 'B'); id E passes. E4 S4 (2 tests): fail on D ('A').
  - T1 (b) (2 tests): `kept` 'A', `replayed_kept` 'A', the rank marker NULL in the record and the replayed record;
    mismatches empty and the samples marker pass.
  - T2 (1 test): the rank marker NULL for lipid and metab; worker status and samples marker pass.
  - G1, G3 (8 tests), G2 (2 tests), T1 (a) (2 tests): pass.
- `test-kit214-sample-columns.R`: passes in full (its G1 objects hold only equally measured duplicates).
- `test-metab-01m-...`: the source-pin test fails on its three new pins (`n_measured`, `keepDuplicateResolutionRow()`,
  the `select`); every other expectation passes.
- `test-metab-02u-...`: the source-pin test fails on its one new pin (`keepDuplicateResolutionRow()`); every other
  expectation passes.

OBSERVED (run 45617808183, 14.5 job-minutes; `r1/tests.tsv`, `r1/failures.txt`): **every prediction held.**
- `test-kit258-duplicate-completeness.R` 76 passed, 26 failed, 0 errors: E1 to E3 (12 tests) fail only on the kept
  row, at E1 "1" "A", E2 "2" "A", E3 "2" "Q"; E4 kernels "A" for both kernels, E4 S4 "A" in both omics; T1 (b) `kept`
  and `replayed_kept` "A" and the rank marker NULL in the record and the replayed record, both omics; T2 the rank
  marker NULL for lipid and metab. G1, G2, G3, T1 (a) and every other expectation pass.
- `test-kit214-sample-columns.R` 176 passed, 0 failed (with the declared edits).
- `test-metab-01m-...` 15 passed, 3 failed: exactly the three new pins. `test-metab-02u-...` 16 passed, 1 failed:
  the one new pin.

Verdict: R1 met (A1). The fix may be written against these tests.

## Fix as far as the owner's decided rule (pending OD-1 and OD-2)

The ticket runs S1 on the fix commit (section 10: "base = staging after KIT-214, head = the fix"), so the fix is written
only as far as VAL-254 decided: at the four sites the most measured (not `NA`, today's definition, unchanged) values
among the averaged columns first, then the highest mean of the measured values, then the first in the table
(`slice_max(with_ties = FALSE)`, unchanged); one internal helper `keepDuplicateResolutionRow` (`R/func_omics_qc_shared.R`)
reads the replay marker `duplicate_resolution_rank` (`measured_then_mean` / `mean`), appended to
`workflowReplayRuleMarkers`; the step's `fixed` gains it. Not written yet (wait for OD-1 and OD-2): E5, E6, NEWS, the
ledger row, the help, log and exported-function texts, `man/`.

## S1 (control on the fix commit, base 63143aa3; 1 scripts job, ceiling 60 job-minutes)

Script `dup-rank-control.R` (KIT-214's control re-pointed as Re-check 1 R1-02 specifies: inputs and variants kept,
classes by completeness on the base side, criteria A4 (a) to (e)). Decides A4 and stop rules 1 and 2.

PREDICTED (written before the run was read):
- (a), (b): every row equal; non-vacuous for lipid and metabolite; every equally measured group keeps its row.
- Classes `e` (unequally measured duplicates): e1_lipid, e1_metab, lipid_msdial_dup_import; e1's kept row does not
  change (row B is both more measured and higher).
- (c), (d1), (d2): every row equal.
- (e): exactly one id changes its kept row, in both variants: lipid_msdial_dup_import / LCMS_Pos / P001, base row 2
  ('TG 52:3', 3 of 4 measured, mean 177.333333) to head row 1 ('PC 34:1', 4 of 4, mean 150.5). Any other changed id is
  stop rule 1.

OBSERVED (run 45620444058, head 81fdcf7c, base 63143aa3; 16.9 job-minutes, script 287 s, both sides exit 0;
`s1/compare.tsv`, `s1/changed-rows.tsv`, `s1/classes.tsv`): **every prediction held; no stop rule.**
- (a) 18/18 and (b) 20/20 rows equal base against head; non-vacuous for both omics; 28/28 equally measured groups keep
  their row.
- Classes `e`: exactly e1_lipid, e1_metab, lipid_msdial_dup_import (as predicted).
- (c) head = KIT-258 oracle 44/44, base = KIT-214 oracle 44/44, head whole = partitioned 22/22.
- (d1) 44/44, (d2) 44/44.
- (e) exactly one id, both variants: lipid_msdial_dup_import / LCMS_Pos / P001, base row 2 'TG 52:3' (3 measured,
  mean 177.333333) to head row 1 'PC 34:1' (4 measured, mean 150.5).

Verdict: A4 met for the rule as decided (measured = not NA). It is re-judged if OD-1 (b) or OD-2 (b) changes the rule.

## S1f (control on the final fix, VAL-284; combined with R2, scripts job; ceiling 60 job-minutes)

The rule changed from what S1 measured (zeros), so the control runs again on the final fix (base 63143aa3 as before;
out `s1f/`). Script changes: "measured" = neither NA nor 0 for the classes, the groups and the KIT-258 oracle; the
KIT-214 oracle keeps the mean of the non-missing values (zeros averaged); the groups table adds that mean.

PREDICTED (written before the run was read):
- (a), (b), (c), (d1), (d2): every row equal; (a), (b) and the group check non-vacuous.
- Classes `e`: e1_lipid, e1_metab, lipid_msdial_dup_import, **and now lipid_msdial_dup_code** (built in code, its
  row 2 holds WT_1 = 0, so it is measured in 3 of 4 against row 1's 4 of 4). Every other input is `a` or `b`: no other
  duplicated row of any input holds a zero (module-CI duplicates complete and non-zero,
  `helper-module-ci-lipid-qc.R:41-53`, `helper-module-ci-metab-qc.R:39-50`; imported tables have their zeros turned to
  missing by the import).
- (e) exactly one id changes its kept row, in both variants, the same as S1: lipid_msdial_dup_import / LCMS_Pos /
  P001, base row 2 'TG 52:3' (3 measured, mean 177.333333) to head row 1 'PC 34:1' (4 measured, mean 150.5).
  lipid_msdial_dup_code does NOT change: base (KIT-214) keeps row 1 (150.5 against 133, the zero averaged), head keeps
  row 1 (4 measured against 3). Why: the zero rule changes counts only where a duplicated row holds a zero, and in
  these fixtures that is the code-built table, where both rules already agree.

OBSERVED (run 45676413310, scripts job 10.2 job-minutes, script 149 s, exit 0; `s1f/`): **every prediction held; no
stop rule.** (a) 18/18, (b) 18/18, groups 26/26 (non-vacuous); classes `e` exactly e1_lipid, e1_metab,
lipid_msdial_dup_code, lipid_msdial_dup_import; (c) 44/44, 44/44, 22/22; (d1) 44/44; (d2) 44/44; (e) exactly
lipid_msdial_dup_import / LCMS_Pos / P001, base row 2 'TG 52:3' (3 measured, mean 177.333333) to head row 1 'PC 34:1'
(4 measured, mean 150.5), both variants. lipid_msdial_dup_code keeps row 1 on both sides.

## R2 (fix 77a719b8 line against tests-only base 5dfe2bc8; 3 shards per side, check; ceiling 600 job-minutes)

Scope: the ticket's grep recounted at the branch point gives 27 files plus the new one (staging added
`test-kit201-treat-moderation.R`, `test-kit224-ruv-feature-levels.R`, `test-kit252-managed-failure-text.R`, which
name `workflowReplayRuleMarkers`), plus the four declared: `test-workflow-omics-binding-steps.R`,
`test-workflow-omics-residency.R`, `test-gui-entry-lipidomics.R`, `test-gui-entry-metabolomics.R`: **32 files**.

PREDICTED (written before the run was read):
- `test-kit258-duplicate-completeness.R`: `fixed` (base fails E1 to E3, E5, E6 on the step (20 tests), the scripted
  E4, E5, E6 (3 kernel tests, 6 S4 tests), T1 (b) (2), T2 (1); head passes every test).
- `test-metab-01m-...` and `test-metab-02u-...`: `fixed` (base fails the new pins only: 4 and 1).
- `test-workflow-omics-binding-steps.R` and `test-workflow-omics-declarations.R`: `message-changed`, counts unchanged
  (the laptop snapshots, already failing on base: the duplicate step records one more fixed parameter, the registry's
  `fixed` one more entry).
- Every other file `same`, including `test-kit214-sample-columns.R` (passes on both).
- Check: findings identical on both sides.
- Stop rule 3 on any other verdict.

OBSERVED (run 45676413310, head e337e72b, base 5dfe2bc8; 238.1 job-minutes in all with S1f and both checks, ceiling
660; `r2/compare.tsv`, `r2/check-compare.txt`):
- `test-kit258-duplicate-completeness.R` `fixed` (base 100 passed, 42 failed; head 142 passed, 0 failed).
- `test-metab-01m-...` `fixed` (base 4 failed: the new pins; head 19/0); `test-metab-02u-...` `fixed` (base 1; head 17/0).
- `test-workflow-omics-binding-steps.R` `message-changed`, counts unchanged (58/1 both sides; `r2/binding-steps.txt`:
  the duplicate step's recorded parameters, as predicted).
- **One verdict not as predicted, declared here (stop rule 3):** `test-workflow-omics-declarations.R:45` is `same`
  (45/23 both sides), not `message-changed`. Its registry expectation (the KIT-073 frozen dump,
  `fixtures/kit073-omics-declarations-frozen.rds`) fails on both sides with the identical text
  "`actual$lipidomics$lipid_duplicate_resolution$declaration$fixed` is a list / `expected...` is NULL": the base
  already has KIT-214's one-entry `fixed`, the head a two-entry one, and the diff printer names only the type, so the
  message cannot show the added entry. The snapshot does differ, as predicted; the comparator cannot see it. Its other
  22 failures (pre-existing last-digit and journey differences of the laptop snapshot) are identical. Counts are
  unchanged. Nothing else is affected; no further run.
- Every other file `same`, including `test-kit214-sample-columns.R` (176/0 both). Pre-existing failures identical on
  both sides: `test-kit224-ruv-feature-levels.R` (1), `test-workflow-operation-registry.R` (2).
- Check: 0 errors, 3 warnings, 3 notes on both sides; findings identical (6), none only on head.

## R3 (the ticket's one correction pass: test E7, tests only): NOT RUN

- Purpose: the mean half of VAL-284 (1) had no falsifying test (audit I02). E7: id V, row A (0, 300, 300, NA) against
  row B (NA, NA, 250, 250), both measured in 2 of 4; A's mean of its measured values 300 (200 with the zero averaged),
  B's 250; the decided rule keeps A, KIT-214's rule keeps B.
- Not run. The base side needs a tests-only commit on `5dfe2bc8` pushed to a base branch; that push was refused by the
  permission system on 2026-11-26 and is parked for the owner. No commit or run for R3 is recorded in `runs.tsv`.
- E7's failure under KIT-214's rule is shown by reading only until the owner allows the base push (KIT-214 averages
  A's non-missing values to 200, below B's 250, and keeps B).
- E7's first execution is the merge train's run, head side only. Prediction there: `test-kit258-duplicate-completeness.R`
  only on head, 162 passed, 0 failed (142 + 20: E7 adds seven tests, four on the declared step with three
  expectations each, one on the two kernels with four, two on the S4 methods with two each). Any other count, or a
  failure in E7: stop and report; no code change.

## Acceptance criteria

| # | Verdict | Evidence |
|---|---|---|
| A1 | met | R1 (45617808183): every failing-first row at KIT-214's value, guards and `test-kit214-sample-columns.R` passing; E5, E6 fail on R2's base side (`r2/tests-base.tsv`, 42 failing expectations in the new file) |
| A2 | met | R2 head: the new file 142/0; the declared edits (section 6 items 1 to 4) are the only edits of existing tests |
| A3 | met, with one declared verdict | R2: `fixed` for the new file, 01m, 02u; binding-steps `message-changed` counts unchanged; declarations `same` (explained above); 27 others `same`; check identical |
| A4 | met | S1 (45620444058, rule before VAL-284) and S1f (45676413310, final rule): every criterion equal, exactly the predicted change |
| A5 | met | NEWS (Scientific changes, after KIT-214's entry); ledger row VAL-<merge> after VAL-294; `wording.md`; this section and "For the change report" |

## Wording, U1

- `wording.md`: the user-facing texts as drafts for the owner (VAL-41).
- U1 (Import Existing Design does not apply the import's zero-to-missing rule, for every later step: intensity filter,
  log2, missingness, duplicates) is **routed, not fixed here**: plan section 13 candidate row (coordinator's ruling
  2026-11-25), parked row 38 (3). This ticket makes only the duplicate step treat a zero as not measured.

## Not checked

- No R on the laptop; every result above comes from a runner.
- Real data: how often duplicated identifiers occur, and how often their rows differ in completeness (VAL-254).
- The Import Existing Design route end to end with zeros (U1, routed): only the duplicate step's handling is tested
  (E5, built in code).
- The user-facing help and log texts in the application: drafts only (`wording.md`).
- `test-workflow-omics-declarations.R`'s registry difference is inferred from the code and the failure text, not shown
  by the comparator.

## For the change report

**Main part (legacy, `main`).** On `main` (and before KIT-214), lipidomics and metabolomics duplicate resolution kept,
for an identifier listed more than once in an assay, the row with the highest mean over every numeric column of the
table, identifier and annotation columns included (the step's columns, every numeric column: `main:R/mod_lipid_qc_duplicates.R:261` and
`main:R/mod_metab_qc_duplicates.R:261`, passed to the kernel at `:274`; the kernel's mean and keep:
`main:R/func_lipid_qc.R:267`, `:275` and `main:R/func_metab_qc.R:267`, `:275`; the same code at
`main:R/func_lipid_s4_objects.R:4634` and `main:R/func_metab_s4_objects.R:4642`). On `main` the import kept an
intensity of exactly 0 as a value (it counted zeros as missing only in its summary,
`main:R/func_lipid_import.R:822`), so a zero entered that mean and lowered the row's average. In 0.9 a zero is read
as missing at import (KIT-118 M3) and, since KIT-258, is not measured by the duplicate rule on every route. Who was affected: users whose table repeats
an identifier (custom imports, Import Existing Design, objects built in code; MS-DIAL default exports have unique
identifiers). What changed in results: such an identifier can now keep a different row, with its own name and
annotation, and every later step (internal-standard normalisation, RUV, the moderated test, the BH family) sees that
row.

**Appendix (0.9 line).** KIT-214 (VAL-173 (1)) restricted the mean to the design's samples but kept the highest mean
alone, which prefers a row seen in 2 of 6 samples averaging 210,000 over one seen in 6 of 6 averaging 150,000.
KIT-258 (VAL-254, VAL-284) keeps the row with the most measured values first (a value is measured when neither missing
nor exactly 0, for the count and for the mean, as at import), then the highest mean of its measured values, then the
first row in the table. The ordering is strict: one more measured sample always wins, so a complete noise-level row
(12 of 12 at 2,000) is kept over a real feature missing in one sample (11 of 12, mean 590,909); the owner accepted
this consequence (VAL-284 (2); test E6). A duplicated internal standard is resolved by the same rule in the
application, and its kept row feeds the ITSD normalisation; the scripted S4 methods keep and rename duplicated
standards instead (unchanged).

- Code: base `63143aa3` `R/func_omics_qc_shared.R:155`, `R/func_lipid_qc_support_helpers.R:59`,
  `R/func_lipid_s4_duplicate_helpers.R:305`, `R/func_metab_s4_qc_methods.R:363` (`slice_max` of the mean); fix
  `77a719b8` `R/func_omics_qc_shared.R:105-138` (`duplicateResolutionMeasured`, `duplicateResolutionMeasuredMean`,
  `keepDuplicateResolutionRow`), called at `R/func_omics_qc_shared.R:205`, `R/func_lipid_qc_support_helpers.R:66`,
  `R/func_lipid_s4_duplicate_helpers.R:311`, `R/func_metab_s4_qc_methods.R:369`; marker in
  `R/utils_workflow_analysis_record.R` (`duplicate_resolution_rank`), fixed by the step at
  `R/utils_workflow_omics_declarations.R:321`.
- Replay: records made before KIT-214 (no marker) replay every numeric column and the mean; records made under KIT-214
  (`duplicate_resolution_samples = design`) the design's samples and the mean; records made now add
  `duplicate_resolution_rank = measured_then_mean`.
- Evidence: R1 45617808183; S1 45620444058 and S1f 45676413310 (15 fixtures, both variants, base against head: exactly
  one identifier changed, P001 in the MS-DIAL lipid fixture through the import, 'TG 52:3' (3 of 4, mean 177.33) to
  'PC 34:1' (4 of 4, mean 150.5)); R2 45676413310 (32 files, check). Limits: constructed and fixture data only; no
  real data; the frequency of such duplicates is unknown.
