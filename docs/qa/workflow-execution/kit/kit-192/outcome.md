# KIT-192 outcome: transcriptomics DA export without a legacy RUV grid

Ticket `tickets/kit/KIT-192.md` revision 1.2 (C01 of the r1.1 closure check applied at step 0). Decisions: VAL-87 (the
ticket), VAL-102 (D1 (a) nothing written in place of the grid summary; D2 (a) the two labels stay out; D3 the budget).
Branch `kit-192-trx-da-export` from `v0.9-staging` at 51959495baf68b549fb773954b19c002049d7132 (the branch point and
every base below). Nothing was run on the laptop; every run is a `run-tests` request (R 4.6.1, KIT-134 snapshot library,
`blas: reference`, `ubuntu-22.04`, `extra_packages: ["sva"]`, `MULTISCHOLAR_TRX_MANIFEST_NOTE=1`).

## Commits

| Commit | Content |
|---|---|
| 31c9dd64 | Tests only: T1 to T6 appended to `tests/testthat/test-trx-da-without-ruv.R`, T7b in `tests/testthat/test-module-ci-trx-summary.R`, S1 `export-routes-control.R` |
| e8f78773 | Tests only: T7b's printed paragraph matched with `perl = TRUE` (R1 instrument defect, below) |
| 816173df | The fix (`R/utils_workflow_trx_declarations.R`, `computeWorkflowTrxDaPublicationAction`: the grid summary only when `grid[["ruv_design", exact = TRUE]]` is NULL, the RUV publication's test); NEWS; ledger row `VAL-<merge>` |
| 0ef0b641 | Records of R1, R1b, R2, S1 (R3 binds this commit) |
| (records) | R3 and this final outcome |

## Re-read at the branch point (step 0)

- `git diff --name-only 147adde8 51959495` moves only `docs/methods-decision-ledger.md`; every cited line of the ticket
  holds. Defect lines read at 51959495: `R/utils_workflow_trx_declarations.R:910-912` (the unguarded summary write),
  `:827-833` (the RUV publication's guard), `:812-819` (`trxPublicationResult`), `R/func_trx_publication.R:249`.
- Scope search (section 9) at 51959495: 29 test files and 2 helpers match; less the four proteomics-only files
  (`test-prot-05b-norm-module-contracts.R`, `test-prot-05d-norm-module-characterization.R`, `test-prot-norm-core.R`,
  `test-omic-norm-module.R`), plus the three lanes (`test-module-ci-trx-app.R`, `test-e2e-transcriptomics-screen.R`,
  `test-gui-entry-transcriptomics.R`): the 28 files the ticket lists.

## A7: the other omics (read at 51959495; no code change)

- Proteomics: `produceDIADaPublication` (`R/utils_workflow_dia_da.R:584`) calls `outputDaResultsAllContrasts` and
  `assertDIADaPublicationFiles` (`:605`); `R/func_prot_da_results_methods.R` has 399 lines, none matching
  `ruv|grid|separation|best_k`.
- Lipidomics and metabolomics: `computeWorkflowOmicsDaPublicationAction` (`R/utils_workflow_omics_declarations.R:1714`)
  calls `binding$da_export`, which is `outputLipidDaResultsAllContrasts` (`R/utils_workflow_lipid_declarations.R:108`)
  and `outputMetabDaResultsAllContrasts` (`R/utils_workflow_metab_declarations.R:99`); `R/func_lipid_da_export.R` and
  `R/func_metab_da_results_output.R` have no match for `ruv|separation|best|round(`.
- The other readers of `separation_score` under `R/` are normalisation and sweep code (proteomics, peptide, lipid,
  metabolite RUV optimisation, records) and `trxCcaPanels` on the legacy branch of the RUV publication; none is a DA
  publication.

## Observed

**R1** (45292234512, tests-only 31c9dd64, `r1-tests.tsv`, `r1-failures.txt`). Each failing-first test failed at an
assertion, at its stated value:
- T1, T3: the tail's third element `"[warning] Warning: Could not save results to disk: non-numeric argument to
  mathematical function"` against the success message; the no-warning assertion.
- T2: the success message absent among the `message` notifications; the warning present.
- T4: both of the above in process and managed (4 failures).
- T6: `display$saved` FALSE, `display$message` "non-numeric argument to mathematical function", `trxDaExport` success
  FALSE.
- T5 passed (legacy grid: tail, files, summary present with its first line and the grid line).
- Every file assertion and every "summary absent" assertion passed on base (as stated).
- T7b **errored** before any statement was compared: `regexpr("No correction for unwanted variation.{0,1200}", text)`
  "invalid regular expression, reason 'Invalid contents of {}'" (the default TRE engine refuses a repeat count above
  255). The expression only builds the paragraph printed in a failure's `info`.
- `MCI-032.1` failed twice at the laptop-recorded `trx_edge_12` digests (RUVIII newY, the integration TSV), as in
  KIT-186's R1 (runner against laptop references, KIT-156); not this ticket's.

**Instrument defect (VAL-54).** The T7b error is in the test's mechanics, with no statement judged. Fixed in e8f78773
by `perl = TRUE` on that one `regexpr` (PCRE allows the count); the statements, `squash`, the comparison and the stop
rule are unchanged; the ceiling is unchanged. One instrument run used of three.

**R1b** (45294352775, e8f78773, `r1b-tests.tsv`, `r1b-failures.txt`). T7b passed on base: all five statements found
after `squash`, the seed an integer, both forbidden phrases absent (9 expectations). Only the two `MCI-032.1` runner
reference failures. T7b's stop rule (a W3 defect) did not fire.

**R2** (45296626686, fix 816173df against 51959495; `r2-compare.tsv`, `r2-compare-failures.txt`,
`r2-check-compare.txt`, `r2-tests.tsv`; S1 in `s1/`).
- Scope: 26 files `same`; `test-trx-da-without-ruv.R` `count-changed` (T1 to T6 passing on head: 5, 6, 6, 13, 6, 10);
  `test-module-ci-trx-summary.R` `count-changed` (T7b, 9). No `new-failure`, no `message-changed`, no `compare-error`.
  217 failure entries on each side, identical (runner against laptop references).
- Check: base and head 0 ERRORs, 3 WARNINGs, 3 NOTEs; 6 findings identical, none only on head.
- S1 (`s1/export-routes.tsv`, `s1/export-routes-compare.tsv`, `s1/console.txt`): every route ran on both sides; 60 of
  60 rows hold their criterion.
  - Routes 1 to 3 (legacy grid): every field identical, including the grid summary's digest with its timestamp line
    removed (present on both sides; routes 1 and 2 the same digest, route 3 its own) and the 12-file listing.
  - Routes 4 to 7 (none in process, none managed, `groups_all_features` k 2 in process, `first_pass_nulls` k 1
    managed): (i) base last notification the save warning, head last the success message, on each; (ii) the rest
    identical; (iii) errors, head, k (0, 0, 2, 1), DA table digest, the 8-file listing and the `.csv`/`.tsv` digests
    identical; the grid summary absent on both sides.

**R3** (45301883983, 0ef0b641 against 51959495, full suite, 20 shards per side, no check;
`r3-full-suite-compare.tsv`, `r3-full-suite-compare-failures.txt`). 455 files: 453 `same`; `count-changed`
`test-trx-da-without-ruv.R` and `test-module-ci-trx-summary.R` (the new tests passing); no `new-failure`,
`message-changed`, `compare-error`, `timeout` or crash. 899 failure entries on each side, identical (runner against
laptop references, KIT-156).

No stored reference moved.

## After R3: staging gained code

Between R2 and the end of R3, `v0.9-staging` moved from 51959495 by KIT-187's merge (proteomics `R/` files, three
proteomics report templates, one `man/` page, three proteomics test files, its records), `NEWS.md` and the ledger. Of these,
only `NEWS.md` and the ledger are files of 7.3 (text); no code file of 7.3 changed, nor `tests/testthat/helper-module-ci-trx-tabs.R`. By the merger's rule (6.4): R2b (R2's tests and
check on the merge) is due; S1 is not repeated; R3b is not due if the only hand-resolved conflicts are in `NEWS.md` and
the ledger and R2b shows 6.2's verdicts. The merge and R2b were then done at the coordinator's request.

**Merge and R2b.** After the text repairs of the implementation audit (0dae776b: I02, I03), `v0.9-staging` at
e962775d was merged into the branch with a true merge (9d264708). The only conflict was the ledger: staging's VAL-104
to VAL-114 kept in order, KIT-192's row after them (still `VAL-<merge>`). NEWS merged cleanly, the three entries in
place under the one unreleased heading. No code file of 7.3 nor `helper-module-ci-trx-tabs.R` differs on staging
since 51959495, so S1 was not repeated. R2b (45311259029, 9d264708 against e962775d; `r2b-compare.tsv`,
`r2b-compare-failures.txt`, `r2b-check-compare.txt`): 26 `same`, `count-changed` for `test-trx-da-without-ruv.R`
and `test-module-ci-trx-summary.R`, nothing else; 217 failure entries identical on both sides; check findings identical
(6), none only on head. R3b is not due (rule 6.4: only the ledger hand-resolved, R2b on 6.2's verdicts).

## Budget

1,270.1 job-minutes (tests, check and scripts jobs; the orchestration jobs, about 2.8 min, not counted): R1 12.8, R1b 13.7 (the VAL-54 instrument run, inside the ceiling), R2 234.2, R3 1,009.4 = 84.67 core-hours against 83.7 expected and 456.4 ceiling. R2b (conditional, outside that ceiling): 191.7 job-minutes = 12.78 core-hours against 15.5 expected and 36.0 ceiling; longest job 34.6 min. Longest jobs: R2 tests 35.5 min (ceiling 60), check 11.8 and 14.9 (90), S1 15.2 (60); R3 53.8 (70). No overrun.

## Acceptance criteria

| # | Verdict | Evidence |
|---|---|---|
| A1 | met | R1 (T1-T4, T6 at stated values; T5 passes), R1b (T7b passes on base after the instrument fix), R2 head (all pass) |
| A2 | met | `s1/export-routes-compare.tsv`: 60/60 |
| A3 | met | `r2-compare.tsv` and, on the merge, `r2b-compare.tsv`: the verdicts of 6.2 exactly |
| A4 | met | `r2-check-compare.txt`, `r2b-check-compare.txt`: none only on head |
| A5 | met | `r3-full-suite-compare.tsv`: 453 `same`, the two predicted `count-changed`; no other verdict |
| A6 | met | NEWS: this entry under the unreleased heading, KIT-186's sentence deleted; ledger row `VAL-<merge>`; `DESCRIPTION`, `packaging/electron/package.json`, lock files, `man/`, `NAMESPACE` unchanged (`git diff --stat 51959495`) |
| A7 | met | section above |

## Inferred (not observed)

- That the managed k = 0 case behaves as `none` managed (same guard, same handler; A11 of the plan audit).
- That a gate-off user sees nothing change (routes 1 to 3 show it for the three legacy routes driven).

## Not checked

- Figure bytes (PNG, PDF) base against head; the figure code is above the changed line.
- Windows (no platform code touched); real data (none used; no number changes).
- F1 (a legacy-style state without a RUVs product reaching the export), F3, F4 (earlier-run outputs left in place):
  recorded in the ticket, not run.
- The runner's Pandoc version (T7b matched with `squash`, so either way).
- Shown by reading only, not by a test or S1: the gate on with the legacy grid chosen; the design `supplied_controls` and the k rule `sweep_bounded`; an export after a replay; the managed worker as a separate process (the tests run its entry point in this process, with the service layer doubled). Each reads the same product field, `ruv_design`, which only `trxChooseRuvDesign` sets (`R/utils_workflow_trx_declarations.R:476-480`).

## Merge record (coordinator)

Merged into `v0.9-staging` on 2026-11-24; the ledger row written `VAL-<merge>` above is VAL-122 (numbered at merge).
