# KIT-213 outcome: study_parameters.txt labels the store-scoped inputs (VAL-169)

Ticket r1.2 (`tickets/kit/KIT-213.md`), register `audits/AUDIT-KIT-213-plan-20261125.md` (re-check 2: may start, no
repair outstanding, so step 0 changed no ticket text). Branch `kit-213` from `v0.9-staging` e305fab3. No R was run on
the laptop; every run below is a GitHub-hosted `run-tests` run (`runs.tsv`).

## Commits

| Commit | What |
|---|---|
| c0369d5f, 8583c127, 00583734, b0bb6aca | G0 script `g0.R` and its instrument repairs (docs/qa only; nothing under `R/`) |
| 67631c03 | tests-only: `tests/testthat/test-kit213-study-parameters-labels.R`, control C's `control-w1.R` |
| 196b0e5c | the fix (`R/func_general_filemgmt_study_parameter_writer.R`), NEWS, ledger row `VAL-341` |
| a68fd69c | control C script repair (one `test_file` per side) before its only run |
| ed4e227f | records (this directory), ledger row names control C's run |
| fa335eb6 | the re-freeze (VAL-169 (1), the owner's exception by file; ledger row `VAL-341`): three `study_parameters` leaves of `omicSummaryW1Frozen` (lines 2700, 2786, 3227) |
| 36a8eccc | merge of staging e6791275 (KIT-167, ledger VAL-178 to VAL-188); ledger conflict resolved by placing KIT-213's row after VAL-188 |
| 6cff8354 | merge of staging 4d588ba3; NEWS and ledger: staging's entries first, then KIT-213's; `R/`, `tests/` merged without conflict; against staging the test file differs only at the three re-frozen leaves (now lines 2532, 2618, 3059 at 6cff8354; 2533, 2619, 3060 on the merge with staging 176f9d3f) |

Branch `run-k213-g0c` (9695502a) holds 00583734 plus b0bb6aca's `g0.R` only, so that G0's last instrument rerun ran
on the base's package code.

## What was observed

**G0** (run 45471271185, journeys part; `g0.tsv`, `g0/`). On the base's package code each omic's journey ran twice in
fresh project roots and saved: lipidomics and metabolomics (KIT-207 G1's journey), transcriptomics (RUV automatic: k 2
chosen, the RUVIII export accepted, so E3 needs no deviation), proteomics (KIT-107 F-16's fixture with the correlation
filter, the DA, its publication and the enrichment; every save `done`, study parameters written).
- Lines holding a store entity id (no line held the store's incarnation): lipidomics and metabolomics line 99
  `pairs [input]` (block `[7. <prefix>_correlation_filtered: <prefix>_pearson_selection]`); transcriptomics lines 34
  `trx_ruv_grid.logcpm`, 42 `logcpm`, 43 `grid` (block `[5. trx_ruvs: trx_ruvs]`); proteomics line 100 `pairs`
  (`[3. correlation_filtered: pearson-selection]`) and line 112 `results` (`[4. enrichment_completed: enrichment]`,
  one element `groupB-groupA`). Exactly the L lines; no other line.
- Raw differences between the two projects before any mask (`g0/g0-diff.tsv`): the save's `Timestamp:` (line 8, every
  omic), proteomics' `  timestamp = ` (line 51) and its `Timestamp:` twin (line 92), and the L lines. Nothing else;
  no line holds the project root (root-mask count 0 on every line). **Stop rule 5 not fired.**
- L2 and L5 reached through reference states without a package change. **Stop rule 4 not fired.**
- L4 (`trx_ruviii_export`'s `logcpm`, `grid`) cannot occur: the export records only `ruviii_m_source`
  (`acceptTrxRuviiiExport`, R/utils_workflow_trx_declarations.R:668-669, the only code that records the transcriptomics
  RUVIII export, has stored only that parameter since KIT-095), and G0 shows the export block printing
  `ruviii_m_source` only. The analysis-level label code (and its "computed in this analysis" wording) is unreachable on
  every current route (implementation audit I01).

**G0 W1 part** (instrument reruns 1 and 3, runs 45473660748 and 45479163619; `g0/g0-w1-run*.tsv`). Cases whose saved
study parameters hold an L line: `trx_main` (writes at `save_before_export` and `save`) and `trx_managed` (`save`),
each with the three trx L lines. No proteomics case (their owners commit `correlation_filtered` directly, no
pearson-selection), no lipidomics or metabolomics case (6 writes, 0 L lines), `trx_technical` writes none.

**R1** (run 45473753944, tests-only 67631c03): 28 passed, 27 failed. E1 to E4 fail at their base values: the
differing masked lines are exactly the L lines (lipid and metab 93; trx 28, 36, 37; prot 92, 104, the `  timestamp = `
line masked), the label lines absent; G5 fails on every project; E4's helper absent. G3, G4 pass. **Stop rule 1 not
fired.**

**R2 + control C** (run 45481567724, fix a68fd69c against 67631c03; `r2/`, `control/`):
- 38 files: 35 `same`; the new file `fixed` (28/27 to 56/0); `test-omic-summary-module.R` 1675/28 to 1672/31 and
  `test-omics-summary-core.R` 676/25 to 673/28, `message-changed` only by the three listed `study_parameters` leaves
  (`r2/compare-failures.txt`). Check identical on both sides (0 errors, 3 warnings, 3 notes; nothing only on head).
- G2: the 367 per-line digests of E1 to E4 (lines without an L parameter, after E1's masks) are identical base and
  head (`r2/kit213-lines-*.txt`). G3 digest 5fcb377bb5c4d6be both sides; G4's `fasta_sha256` line identical.
- On head every pair of projects is equal after E1's masks, the exact label lines are present, G5 finds no
  store-bound line.
- Control C (`control/verdict.txt`, `leaves.tsv`): the masked texts differ base against head only at L lines (12
  lines over 4 writes); head against head identical; (i) the base dump equals the stored value at all 3 changed
  leaves; (ii) the two head dumps are byte-identical; (iii) the changed leaves are only the `study_parameters` leaves
  of `trx_main`/`save_before_export`, `trx_main`/`save`, `trx_managed`/`save`. Seven pre-existing leaves (base dump
  differs from stored: trx `files`, `integration` leaves) are recorded and left untouched. **Stop rule 3 not fired.**

**Re-freeze** (fa335eb6): exactly those three lines, to the head dump's values (`47 lines 387b8047e011288e`,
`50 lines bbadf5a52d70875f` twice); `git show fa335eb6` is 3 insertions, 3 deletions.

**R3a** (run 45490237295, full suite on the merge 36a8eccc with staging e6791275, no check; `r3a/compare.tsv`): 477
files `same`, the new file only on head 56/0; both W1 files `same` (1675/28 and 676/25, the laptop-reference failures
on both sides); no timeout or crash. The watcher that should have read it stopped on GitHub API rate limits, so its
result was read only on 2026-11-26.

**R3** (the coordinator's run on the merge 6cff8354 with staging 4d588ba3, with check): predictions committed in
`r3-predictions.md` before it was sent. Run 45700894419 (`r3/`): 511 files `same`, the new file only on head 56/0;
`test-omic-summary-module.R` 1675/28 and `test-omics-summary-core.R` 676/25 `same`; check identical on both sides
(0 errors, 3 warnings, 3 notes; nothing only on head); no timeout or crash. Predictions 1 to 4 met; prediction 5
(cost about 84 core-hours) missed: 1,423.87 job-minutes = 94.92 core-hours, inside the R3 row's ceiling (R3a + R3 =
176.22 of 186.67).

## Re-freeze by file (VAL-169 (1); ledger row `VAL-341`)

| File | Object | Leaf | Stored (base) | Re-frozen (head dump, run 45481567724) |
|---|---|---|---|---|
| `tests/testthat/test-omic-summary-module.R` | `omicSummaryW1Frozen` | `trx_main` / `save_before_export` / `study_parameters` | `47 lines 96b3fdbde9b1bb30` | `47 lines 387b8047e011288e` |
| same | same | `trx_main` / `save` / `study_parameters` | `50 lines 6309448be78f095c` | `50 lines bbadf5a52d70875f` |
| same | same | `trx_managed` / `save` / `study_parameters` | `50 lines 6309448be78f095c` | `50 lines bbadf5a52d70875f` |

Nothing else in that file or any other stored reference was changed. `test-omics-summary-core.R` asserts the same
object and is not edited. KIT-212's re-freeze (VAL-172, VAL-186; now on staging) touched other leaves of the same object
(Copy and Report leaves; `prot_owner`, `prot_owner_variants`): disjoint from these three.

## Acceptance criteria

| # | Verdict | Evidence |
|---|---|---|
| A1 | met | `g0.tsv` (store-bound, masked and L lines per omic and project), `g0/g0-diff.tsv`, `g0/g0-w1-run2.tsv`, `g0/g0-w1-run4.tsv` |
| A2 | met | R1 45473753944 (failing at stated values; G3, G4 pass); R2 45481567724 (56/0 on head) |
| A3 | met | control C in 45481567724; re-freeze commit fa335eb6 edits exactly the three changed leaves |
| A4 | met | R2 `r2/compare.tsv`: only the new file and the two W1 files at the listed leaves; check identical |
| A5 | met | R3 45700894419 on the merge with staging 4d588ba3 (`r3/compare.tsv`, `r3/check-compare.txt`); also R3a 45490237295 on e6791275 |
| A6 | met | NEWS `### Fixes` entry with VAL-169 (2)'s limit; ledger row `VAL-341`; no version, NAMESPACE or man/ change (`git diff e6791275 36a8eccc -- DESCRIPTION NAMESPACE man` empty) |

## Cost

G0 and its three instrument reruns 50.91 job-minutes (3.39 core-hours; G row 6.00, instrument row 18.00); R1 13.15
(0.88; 4.00); R2 with control C 248.65 (16.58; R2 36.00 + C 4.00); R3a 1,219.45 (81.30) and R3 1,423.87 (94.92;
together 176.22 of the R3 row's 186.67). Total 2,956.03 job-minutes = 197.07 core-hours against the ticket's expected
103.5 and ceiling 481.33. The excess over the expectation is the second full suite, sent on the coordinator's
instruction after staging moved by 553 commits.

## Inferred, not executed

- Two projects on the same build and platform now write the same file apart from the run-bound lines G0 listed:
  executed for one journey per omic on one runner image, not across machines or builds.
- The analysis-level label (L4) and the `NA` source-state wording ("in this analysis") are unreachable today (I01);
  E4 executes only the `results` `NA` element through the helper directly.

## Deviations from the ticket text

- G0 ran from its own commit before the tests-only commit (its printed step numbers are the E3/E4 literals), not from
  the tests-only commit; `R/` was the base's in every G0 run.
- G0's W1 part needed all three VAL-54 instrument reruns (`runs.tsv`): an `env` name clash in `evalq`; testthat's
  `teardown_env()` outside a test file; `desc` cannot select the computed W1 test names. The last ran the 13 lipid and
  metab cases in one `test_file` and attributes writes by source directory, not by case. No expected value changed.
- Control C runs one `test_file` per side with the listed cases (the file's other tests run too) instead of one per
  case; writes are compared in order.
- The helper's signature is `studyParameterRecordValueText(parameter, step_order)`: the ticket's third argument
  `operation` is unused by the labels of section 7 and was left out.
- When the step order is unknown (`NA`), the label reads "... computed in this analysis" (the ticket says "the label
  without 'in step <n>'").

## Not checked

- L4 (`trx_ruviii_export` analysis lines): unreachable on every current route (I01), so not executable.
- Windows; other runner images; two different builds.
- The other five report templates' parsers (as the ticket).
- KIT-215 and KIT-212 are on staging 4d588ba3; R3 (45700894419) is the run on the tree holding KIT-212's re-freeze and
  KIT-215's new file, which discharges the ticket's section 6 and 8 rules for them (R3 ran
  `test-kit215-not-tested-lists.R`, `same`). KIT-216, KIT-217 and KIT-263 are not on 4d588ba3. By
  `git merge-tree` only: kit-213 (6cff8354) merges against origin/merge-train-20261126k (KIT-263) and
  origin/merge-train-20261126h (KIT-217) with a conflict in `docs/methods-decision-ledger.md` only; the writer and
  `test-omic-summary-module.R` auto-merge. Not executed.

## For the change report

**New to 0.9.** The analysis-record section of study_parameters.txt exists only on the 0.9 line (KIT-033; `main` has
none), so no released version printed these lines; only 0.9 development builds, including the v0.9.0 staging tag,
were affected (implementation audit I05).

**What was wrong.** In `study_parameters.txt`, the analysis-record section printed four kinds of parameter that
point to an earlier result of the same analysis (the correlation filter's sample pairs in lipidomics, metabolomics
and proteomics; the transcriptomics TMM logCPM and RUV choice; the proteomics enrichment's DA tables) as a 64-digit
hexadecimal value. That value is the result's identifier inside the session's private store, a SHA-256 over a
descriptor that includes the store's random identity, drawn whenever the application creates a store (every session
that imports data). It looked like a checksum of the result but was not one, and two users running the same analysis
on the same data with the same version got different files. Affected: every lipidomics, metabolomics and proteomics
run with the correlation filter, every transcriptomics run with RUV, every proteomics enrichment run. No number, no
result and no other file was affected; nothing reads these values back (replay, reports).

**Code.** Base e305fab3: `R/func_general_filemgmt_study_parameter_writer.R:991-1003` printed
`workflowRecordValueText(parameter$value)` for every record parameter. Fix 196b0e5c: the same lines call
`studyParameterRecordValueText(parameter, order)` (`:990-1006`), new helper at `:1052-1080`.

**What changed in results.** One text field per affected line, for example `pairs [input] = the sample-pair
correlations computed in step 7 of this analysis`, `grid [input] = the RUV choice made in step 5 of this analysis`,
`results [input] = groupB-groupA=the differential abundance results of contrast groupB-groupA in this analysis`.
Nothing else in the file, and no other file.

**Limit (VAL-169 (2)), stated in NEWS.** `analysis_record.json` still names each result by its store identifier, and
its `products[].identity_sha256` values are computed over those identifiers: two projects running the same analysis
write different records, and a record replays only inside a copy of the project that made it. Three writers make a
file named `study_parameters.txt`: `createWorkflowArgsFromConfig` (changed here), `createStudyParametersFile`
(exported, no record block, unchanged) and the integration's `integrationStudyParameterLines`, whose "Checksum of the
model's identity" is project-specific too (unchanged; hand-off H2).

**Evidence.** G0 45471271185 (two projects per omic, raw differences); R1 45473753944; R2 and control C 45481567724;
R3a 45490237295; R3 45700894419. Simulation settings: none (deterministic fixtures: F172 seed 172001, MCI-032's trx corpus, MCI-009's
proteomics object); limits as in "Inferred" above.
