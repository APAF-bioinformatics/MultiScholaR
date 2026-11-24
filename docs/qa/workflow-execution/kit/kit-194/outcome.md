# KIT-194 outcome: lipidomics, metabolomics and transcriptomics normalisation state guards

Ticket revision 1.3 (r1.2: B01 to B05 applied at step 0, D1 to D4 answered, VAL-105; r1.3: the correction pass after
R1 stopped at A3, coordinator's ruling option 1). Branch `kit-194-omics-normalisation-guards` from `v0.9-staging`
447ba678, with `v0.9-staging` 19997323 merged in (8290ddb8). Every run on `run-tests` (R 4.6.1, reference BLAS,
ubuntu-22.04); ids, commits and job minutes in `runs.tsv`. **State: done; ready for the implementation audit.**

| Commit | Content |
|---|---|
| 5c727c10 | tests first: `tests/testthat/test-omics-norm-state-guards.R` (F1, T1 to T10, the classes block) and `fj-confirm.R`; `R/` equal to base |
| c846908d | tests first, three instrument repairs after G0 (feasibility.md); `R/` equal to base |
| 09f05774 | the fix: `R/mod_omic_norm.R`, `R/utils_workflow_trx_norm.R`, `R/utils_workflow_omics_actions.R`, `R/mod_omic_qc_tabs.R`, `R/mod_trx_qc_filter.R`, `R/mod_omic_da.R`; NEWS; ledger row `VAL-<merge>` |
| cc3adec4 | records of G0 and R1 (the A3 stop) |
| 43fffa9c | the correction pass: `test-metab-03d-norm-module-characterization.R` lines 320 and 327 only (start state `metab_norm_complete` to `metab_qc_complete`; no expectation changed); the NEWS and ledger sentence on lineage-free normalisation-named states |
| 8290ddb8 | true merge of `v0.9-staging` 19997323 (carries KIT-187); ledger conflict resolved with staging's rows VAL-106 to VAL-121 kept and this ticket's `VAL-<merge>` after them; NEWS merged without conflict; no file under `R/` this ticket edits changed on staging |

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

**The A3 stop (read in the code).** The new failure's harness committed a plain object as `metab_norm_complete` and
pressed "Run Normalization Pipeline" in a fresh tab (`test-metab-03d-norm-module-characterization.R:320-352`); by 7.3's R3
branch and VAL-105 (1) that press is refused, so no run, no Skip and no export followed. The ticket's 6.2 had misread
the harness's start state. Coordinator's ruling (option 1, as the ticket's one correction pass): the decided rule stands;
the fixture's start state becomes `metab_qc_complete` (two lines; every expectation unchanged); ticket r1.3 records it
(6.2, 7.8 row 6b, 12.2).

**R1c, the correction pass (run 45313963471), merged head 8290ddb8 against `v0.9-staging` 19997323 plus the ticket's
tests (branch `run-k194-r1c-base`, bfd0bae1; `R/` equal to staging).** New file base 358 passed, 280 failed, 0 errors;
head 647 passed, 0 failed (`fixed`). All 45 other files `same`; `test-metab-03d-norm-module-characterization.R` 13
passed on both. R CMD check: findings identical, none only on head (`r1c-check-compare.txt`). Effect lines on head as R1
(second run maximum 20.11, started from `<prefix>_qc_complete`).

**R2 (run 45319721518), full suite, 8290ddb8 against 19997323.** 457 files `same`, 1 `only-on-head`
(`test-omics-norm-state-guards.R`, 647 passed, 0 failed); no new failure, no changed message; 35 files fail identically on base and head (920 failed expectations and 13 errors on each side) (known laptop-recorded references); no timeout or crash; longest job 50.5 minutes (`r2-compare.tsv`).

**After the implementation audit (register `AUDIT-KIT-194-impl-20261124.md`).** Text repairs I01, I03, I04, I05, I10
(b9079e3d, no code); `v0.9-staging` 0b926295 merged (4295b497; ledger conflict only, staging's rows to VAL-134 kept and
this ticket's `VAL-<merge>` after them; NEWS merged without conflict; none of this ticket's `R/` files changed on
staging). **R1b (run 45328973335)**, 4295b497 against staging 0b926295 plus the ticket's tests (branch
`run-k194-r1b-base`, 4e8bc055): the new file `fixed` (as R1c), the 46 other files `same` (`test-module-ci-trx-summary.R`,
added for KIT-192, 82 passed and 3 failed identically on both sides), check findings identical (`r1b-compare.tsv`,
`r1b-check-compare.txt`). R2 is not repeated (6.3: only NEWS and the ledger were hand-resolved).

**Fixtures building a lineage-free normalisation-named start (surveyed for the ruling; none changed).** Each passes
because it does not press Run (or QC) at that state:
- `test-lipid-12a-norm-module-characterization.R:118` (owner at `lipid_normalized`): the lipid twin; its live blocks
  start the module without pressing Run; the blocks that pressed Run (`:640` on) are quarantined (Q-96).
- `test-lipid-12-norm-module-contracts.R:2605`, `:2662`, `:2738` (`lipid_normalized`, `lipid_ruv_corrected`): Skip,
  correlation and Reset contracts on mocked handlers; the file is skipped (0 tests run on base and head).
- `helper-module-ci-lipid-da.R:165`, `helper-module-ci-metab-da.R:161`, `test-metab-02aa-da-display-filter-characterization.R:2084`,
  `test-omic-da-module.R` lipid/metab fixtures (`omicsTestOwner(..., "<prefix>_correlation_filtered")`): DA tab only; the
  fallback cases are accepted by the name list (divergence 2), the load cases name their state.

## Acceptance criteria

| # | Status | Evidence |
|---|---|---|
| A0 | met | ticket r1.2: decisions answered (VAL-105), D1 (a) text removed, closure section 18; r1.3 the correction pass |
| A1 | met | `feasibility.md`, `fj-confirm.tsv`: G0 all routes reproduced; F1 passes; F2, F3 recorded |
| A2 | met | R1 and R1c: base 280 failures, 0 errors, head 0 failures (`fixed`); T9, T10 pass on both |
| A3 | met after the correction pass | R1 (run 45307575111) showed `test-metab-03d` `new-failure` (A3 stop); R1c (run 45313963471): every other file `same` |
| A4 | met | `r1-check-compare.txt`, `r1c-check-compare.txt`: nothing only on head |
| A5 | met | R2 (run 45319721518): 457 `same`, the new file `only-on-head` |
| A6 | met | `effect.tsv` |
| A7 | met on the branch | NEWS under the unreleased heading; ledger row `VAL-<merge>`; `DESCRIPTION`, `NAMESPACE`, `man/`, `packaging/`, lock files unchanged |
| A8 | met | this file, `runs.tsv` |

## Divergences recorded (for the implementation audit)

1. **7.3, how the remembered start is found** (`omicNormRunStart`, `R/mod_omic_norm.R`). The ticket's pseudocode takes
   the start's save label (`start$identity$parameters`) as its state name, then checks that the name is in the active
   history and that its current reference has the start's entity id. The code looks for the active-history state whose
   current reference has the start's entity id. Why it is needed: the two are the same condition only while the label
   equals the state name, which the ticket left for F1 to confirm ("its save label = its state name") and which a
   recommit or a content-addressed save could break; matching by entity id needs no assumption about labels and cannot
   select a state whose table differs from the one the previous run read. Tested by T3 (a) to (c) and T4.
2. **7.5, the DA fallback's A01 check only when the owner holds a `<prefix>_normalized`**
   (`omicDaFallbackHeadNormalised`, `R/mod_omic_da.R`). The ticket applies the derivation check to every
   `_ruv_corrected`, `_correlation_filtered`, `_norm_complete` head; with no `_normalized` in the owner the check is FALSE
   and the DA would be refused. Why it is needed: the hazard A01 closes (an earlier run's state made current by a Revert
   after Reset and a second run) needs a second run, which always commits `_normalized`, so an owner without one has no
   earlier run to confuse; without the condition the lipid and metab no-contrasts and invalid-preflight cases of
   `test-omic-da-module.R` (owner fixtures holding only `<prefix>_correlation_filtered`, regression scope) would show the
   DA refusal in place of their recorded messages (read before R1; R1 shows the file `same`). Proteomics' name list
   still applies there.
3. **B02, B03, the managed `failed` callbacks** (`omicQcManagedFailed`, `R/mod_omic_qc_tabs.R`; trx's in
   `R/mod_trx_qc_filter.R`). The register's texts show every error passed to `failed`; the code shows a notification
   only for the QC guard's refusal (class `omic_qc_input_normalised`) and otherwise does what it did before (removes the
   working notice; trx had no callback). Why it is needed: a managed step's own failure is already reported by the step
   (`installWorkflowOmicsActions`, its step's `failed`, `R/utils_workflow_omics_actions.R:57-66`), so showing it again
   would duplicate every worker error on managed sessions, a visible change beyond this ticket; the only failure that
   reaches `failed` unreported is a refusal made when a queued action starts (A02), which is what B02 and B03 address.
   Tested by T1 A02 and T2 managed.
- Also: the trx run's refusal is written to the tab's log as well as notified (as lipid/metab's).

## Not checked

- The real worker's busy refusal and start-time QC check under a press that actually waits (stubs only; the worker
  journey presses at idle and stayed `same`).
- R2 is against staging 19997323; the merge with 0b926295 is covered by R1b under 6.3, not by a second full suite.
- Enrichment downstream of an affected DA; the report text of an affected record; F-m, F-p, F-q, F-r (candidates).
- Lipid's recorded log2 argument (D4 (b2)) and the state names older builds wrote.

## Left for the owner's wording pass (implementation audit I02; VAL-41)

Not changed here; the auditor's proposals, verbatim:
- **The DA text for a loaded file** (D4 (b1), the loaded branch of `runOmicDaAnalysis`). Today it shows
  `omicNormGuardTexts()$da` ("…Run the normalisation (and the correlation filter or its Skip) again, then run the
  analysis."), which that adopted session refuses (A07). Proposed separate text: "This session file was saved after
  Quality Control was reverted past the normalisation, so its table ("<state>") was not normalised and cannot be used
  for the differential analysis. Import the data again, run quality control and the normalisation, then run the
  analysis." T8 (b1)'s expectation (`k194DaText`) changes with it.
- **The export text** repeats the correlation sentence. Proposed: "the current table is not the one the last
  normalisation produced (Quality Control was reverted since), so the exported session would not hold the normalised
  table the analysis reports. Press "Run Normalization Pipeline" again, then export."

## For the change report

**What was wrong, and who was affected.** In lipidomics and metabolomics the normalisation tab remembered that a
normalisation had run, but nothing checked that the table in front of it was still that normalisation's result.
- After Quality Control was reverted past the normalisation (or after "Reset to Pre-Normalization"), "Apply Correlation
  Filter", its Skip, "Export Normalized Data" and the differential analysis worked on the quality-control table. The DA
  tested intensities with no internal-standard step, no log2 and no between-sample normalisation, and reported the
  analysis complete. A Revert that stopped inside the normalisation gave a partly normalised table.
- Pressing "Run Normalization Pipeline" a second time without Reset (to try another method) normalised the normalised
  table: log2 of log2 values, then normalisation again; at the default offset 0 every value at or below 1 on the first
  log2 scale became missing. The same happened in a reopened session. In transcriptomics a second run with fewer RUV
  factors stopped with an error and kept the earlier factors.
- Quality-control steps could be applied to a normalised table (thresholds computed on log or normalised values); in
  transcriptomics the expression filter on `trx_ruvs` kept factors estimated on the unfiltered genes.
Affected: any lipidomics, metabolomics or transcriptomics session that took one of these routes in the application.
Not affected: the default path (quality control, one normalisation, correlation or Skip, DA), a run after Reset, a run
after a Revert into Quality Control, every replayed analysis record, the workbooks, and proteomics (KIT-182).

**Code at the base (447ba678).** `R/mod_omic_norm.R:343-378` (`runOmicNormPipeline` reads the head with no state check),
`:457-500` (`applyOmicNormCorrelation` checks only the tab's flags), `:507-537` (`skipOmicNormCorrelation`),
`:558-604` (`exportOmicNormSession` checks only `normalization_complete`); `R/mod_omic_da.R:286-289`, `:299-309`
(`runOmicDaAnalysis`: the head fallback needs only the status "complete"); `R/mod_omic_qc_tabs.R:174-207`, `:356-383`,
`:758-799` (QC steps read the head unchecked); `R/mod_trx_qc_filter.R:52-86`; `R/utils_workflow_trx_norm.R:136-158`;
`R/utils_workflow_omics_actions.R:249-266` (managed QC actions read the head when they start).

**Code at the fix (09f05774).** The guards: `R/mod_omic_norm.R:821` (`assertOmicQcHeadNotNormalised`), `:835`
(`omicNormRunStart`, the run rule: busy refusal, re-run from the previous run's start, refusal of a normalised or
adopted input), `:863` (`omicNormApplyRunStart`), `:883` (`assertOmicNormHeadIsResult`); their calls at
`R/mod_omic_norm.R:356` (run), `:472` (Apply), `:519` (Skip), `:589` (export); `R/utils_workflow_trx_norm.R:140`;
`R/mod_omic_qc_tabs.R:183`, `:380`, `:782` (QC presses) and `:267` (`omicQcManagedFailed`);
`R/utils_workflow_omics_actions.R:264` (managed QC judged when it starts); `R/mod_trx_qc_filter.R:61`, `:89`;
`R/mod_omic_da.R:299` (`omicDaFallbackHeadNormalised`), `:315` (`omicDaLoadedStateRefused`), `:341` (their call).

**What changed in results.** Only for sessions that took the routes above: the step is refused with a message naming
what to press, or (a second run at the tab's own result) the run starts again from the table the first run started
from, so its result equals a single run with the new settings. Every default-path result is unchanged (R2: 457 files
`same`). Analysis records are unchanged and replay as recorded.

**Evidence and its limits.**
- G0 (run 45304509096) executed the defect on the base: four routes on both omics, the DA testing the QC table
  (maximum 1,373,458.9 against 20.38 normalised), 514 rows, status "complete" (`fj-confirm.tsv`).
- The new tests (R1 45307575111, R1c 45313963471): 280 failures at stated values on the base, none at the fix;
  second-run maximum 4.33 on the base against 20.11 at the fix (`effect.tsv`).
- R2 (45319721518): the full suite unchanged apart from the new file.
- Limits: synthetic fixture F172 (520 features, two groups of six), a mechanism check, not an effect size on real
  data; on F172 "ITSD by pattern" finds no standard, so the ITSD step was not exercised (feasibility.md, F3); managed
  sessions were tested with an in-process runner (the real worker only through the unchanged worker journey); how often
  users took these routes is not known. Earlier results can be checked from their analysis record (NEWS, the recipe).

## Merge record (coordinator)

Merged into `v0.9-staging` on 2026-11-24; the ledger row written `VAL-<merge>` above is VAL-139. The merge check R1b (run 45328973335) was on the merge with staging 0b926295; staging has since gained KIT-175 and KIT-199, which share no file with this branch except NEWS and the methods ledger.
