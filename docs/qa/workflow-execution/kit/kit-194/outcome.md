# KIT-194 outcome: lipidomics, metabolomics and transcriptomics normalisation state guards

Ticket revision 1.2 (B01 to B05 applied at step 0; D1 to D4 answered, VAL-105). Branch
`kit-194-omics-normalisation-guards` from `v0.9-staging` 447ba678. Every run on `run-tests` (R 4.6.1, reference BLAS,
ubuntu-22.04); ids, commits and job minutes in `runs.tsv`. **State: stopped after R1 at acceptance criterion A3** (one
regression-scope file changed; owner decision needed, below). R2 (the full suite) has not run.

| Commit | Content |
|---|---|
| 5c727c10 | tests first: `tests/testthat/test-omics-norm-state-guards.R` (F1, T1 to T10, the classes block) and `fj-confirm.R`; `R/` equal to base |
| c846908d | tests first, three instrument repairs after G0 (feasibility.md); `R/` equal to base |
| 09f05774 | the fix: `R/mod_omic_norm.R`, `R/utils_workflow_trx_norm.R`, `R/utils_workflow_omics_actions.R`, `R/mod_omic_qc_tabs.R`, `R/mod_trx_qc_filter.R`, `R/mod_omic_da.R`; NEWS; ledger row `VAL-<merge>` |

## Observed

**G0 (run 45304509096, base).** F-j reproduced by all four executed routes on both omics: after a run, a Revert to
`<prefix>_qc_complete` then Apply (J1) or Skip (J2), Apply then Revert (J3) and Apply then Reset (J4), the DA ran with
the tab's status "complete" on a table equal to `<prefix>_qc_complete` (maximum 1,373,458.9; the normalised maximum
20.38), 514 result rows each (`fj-confirm.tsv`). No stop rule fired. F1 held in every case on both paths; F2 (trx second
run: an error, `W_2` kept) and F3 (lipid/metab second run with ITSD by pattern: no error, normalised twice) recorded
(`feasibility.md`).

**R1 (run 45307575111), fix 09f05774 against tests-only c846908d.**
- The new file: base 358 passed, 280 failed, **0 errors**; head 647 passed, 0 failed (`fixed`). Base failures per test:
  `r1-new-file-base-failures.tsv`. T9, T10, the classes' no-op cases and every "both commits" case passed on both.
- Effect (`effect.tsv`, synthetic F172, a mechanism check, not an effect size on real data): a second run's normalised
  maximum 4.33 on base against 20.11 on head (equal to the first run's; head re-ran from `<prefix>_qc_complete`, in
  process and managed); J1's DA input maximum 1,373,458.9 on base, refused on head.
- Regression scope: 44 of 45 other files `same`; known base failures identical (test-omic-da-module 53,
  test-omic-summary-module 28, test-gui-entry-transcriptomics 13, and others, identical messages).
  **`test-metab-03d-norm-module-characterization.R`: `new-failure`, 2 expectations** in "mod_metab_norm_server preserves
  export prereq and success behavior": the export's files are not written on head (`r1-compare.tsv`).
- R CMD check: 0 ERROR, 3 WARNING, 3 NOTE on each side; findings only on head: none (`r1-check-compare.txt`).

## Inferred (read in the code, not run)

The new failure's harness commits a plain object as `metab_norm_complete` and presses "Run Normalization Pipeline" in a
fresh tab (`tests/testthat/test-metab-03d-norm-module-characterization.R:320-352`). By 7.3 (R3) and VAL-105 (1), a run
pressed at a normalisation state with no remembered start is refused ("this tab does not know which table the earlier
normalisation started from"); the run therefore does not happen, Skip finds no completed run, and the export is
refused as before the run. This is the behaviour T3 (d) asks for; the ticket's 6.2 expected the file to stay `same`
because it read the harness as "one run per harness", not its start state. Two ways forward, both outside what the
ticket authorises, so the work stops here (one correction pass is left):
1. Edit the existing test's start state to `metab_qc_complete` (two lines, 320 and 327): its export expectations are kept
   and run on a table the run may normalise. Not in 7.8's file list.
2. Narrow R3: refuse at a normalisation state only when the head has normalised lineage or was adopted; a state with a
   normalisation's name but no such lineage (made only by scripts or fixtures) runs on the head as R1. F1 shows every
   state a run commits carries lineage, so T3 (d) still refuses. A change of 7.3's decided rule.

## Acceptance criteria

| # | Status | Evidence |
|---|---|---|
| A0 | met | ticket r1.2: decisions answered (VAL-105), D1 (a) text removed, closure section 18 |
| A1 | met | `feasibility.md`, `fj-confirm.tsv`: G0 all routes reproduced; F1 passes; F2, F3 recorded |
| A2 | met | R1: base 280 failures, 0 errors, head 0 failures (`fixed`); T9, T10 pass on both |
| A3 | **not met** | `test-metab-03d-norm-module-characterization.R` `new-failure` (above); all other files `same` |
| A4 | met | `r1-check-compare.txt`: nothing only on head |
| A5 | not run | R2 waits for the decision on A3 |
| A6 | met | `effect.tsv` |
| A7 | met on the branch | NEWS under the unreleased heading; ledger row `VAL-<merge>`; `DESCRIPTION`, `NAMESPACE`, `man/`, `packaging/`, lock files unchanged |
| A8 | met | this file, `runs.tsv` |

## Divergences recorded

- 7.3: the remembered start is found as the active-history state whose current reference has the start's entity id,
  not by its save label; the same condition (name in the history, entity unchanged), without depending on the label.
- 7.5: the A01 derivation check on the DA fallback applies when the owner holds a `<prefix>_normalized`; an owner with
  none (a fixture's single state, `test-omic-da-module.R`'s no-contrasts and preflight cases) is judged by the name list,
  as proteomics. Without this those regression-scope cases would change (read before R1; R1 shows them `same`).
- B02, B03: the managed `failed` callbacks show a notification only for the QC guard's refusal; a worker's own failure is
  already reported by the step, so it is not shown twice (the sub-tabs' and trx's visible behaviour on other failures
  unchanged).
- The trx run's refusal is also written to the tab's log (as lipid/metab's).

## Not checked

- The real worker's busy refusal and start-time QC check under a press that actually waits (stubs only; the worker
  journey presses at idle and stayed `same`).
- R2, the full suite base against head.
- Enrichment downstream of an affected DA; the report text of an affected record; F-m, F-p, F-q, F-r (candidates).
- Lipid's recorded log2 argument (D4 (b2)) and the state names older builds wrote.
