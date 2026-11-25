# KIT-215 outcome: enrichment lists not tested are shown and reported as "Not tested" (VAL-170)

Ticket r1.3 (`tickets/kit/KIT-215.md`, section 18: step 0 on f93dfd85). Plan audit `AUDIT-KIT-215-plan-20261125.md`
(re-checks 1 and 2; step-0 ruling SR-1 to SR-4). Branch `kit-215`, cut from `v0.9-staging` f93dfd85.

## Commits

| Commit | What |
|---|---|
| 01725d76 | tests first: `tests/testthat/test-kit215-not-tested-lists.R` and G0's script `g0.R`; nothing under `R/` |
| 6d002423 | G0 instrument fix (VAL-54): the W1 module fixture's annotation gains its two absent GO aspects as empty columns |
| c9219652 | the fix (`R/`, the four report templates, NEWS, ledger `VAL-<merge>`, one line of `test-workflow-omics-binding-steps.R`, S0-10) |
| f7ee8aa3 | merge of `origin/v0.9-staging` d048898d (ledger: staging's VAL-191, VAL-192 first, this ticket's row after; no `R/`, `inst/` or test file changed on staging since f93dfd85) |
| 66f32932 | correction pass (R2): T7b asserts the refusal's class (the ticket's stated expectation) and its sentence, not the caught error's message |

## Runs (observed)

See `runs.tsv`.

- **G0 + R1, run 45493249015 at 01725d76** (`r1/`, `g0a/`):
  - G0 engine conditions on the runner (R 4.6.1, enrichit 0.2.0, clusterProfiler 4.20.0, gson 0.2.0, GO.db 3.23.1): (i) n = 0 gives NULL, (ii) m' = 0 gives NULL, (iii) a tested term holding a list gene gives an `enrichResult` with 1 row, (iv) m = 0 gives NULL. All four `pass` (`g0a/engine.tsv`), so stop rule 3 did not fire.
  - The script then stopped in the W1 module family (`g0a/console-tail.txt`: the helper's annotation has no MF/CC columns). This is an instrument defect: no family row was produced or judged. Fixed in 6d002423; re-run once as G0b.
  - R1: `test-kit215-not-tested-lists.R` 51 passed, 44 failed, 5 errors (`r1/tests.tsv`, `r1/failures.txt`).
    - Every failing-first row (T1, T2, T2b, T3, T4, T5, T6, T7, T7b, T7c, T10, T11, T12) fails at its stated base value. Examples: "No up-regulated enrichment data"; "Up-regulated: 0"; "No gprofiler2 results available."; "Down-regulated: 0"; outcome NULL; no `not_tested` or `no_list`. Where the head refuses, the base completes with no error.
    - The 5 errors are `conditionMessage(NULL)` on the rows whose base completes (T4, T7, T7b, T7c, T10).
    - T4b (head only) fails on base, as expected.
    - The guards T7d, T8, G1, G3 and G5 pass. Stop rule 1 did not fire.
- **G0b, run 45495664583 at 6d002423** (`g0b/g0.tsv`):
  - t7_complete: both lists tested (m = 97).
  - k175_incomplete: both lists tested (n 11 and 4, m = 47), so KIT-175's incomplete-annotation runs complete as before.
  - k208_lacks_up: the up list is not tested (n = 0, m = 55); down is tested.
  - w1_module (as written): m = 0 for a universe of 6. The real engine would return NULL. The captures call a fake engine, so they do not move (the fixture limit below).
- **R2 + G2, run 45495673217, c9219652 (head) against 6d002423 (base)** (`r2/`, `g2/`):
  - 25 of the 26 files are `same`, including `test-omic-enrich-module.R` (W1/W113: 19 failures, identical on both sides), `test-workflow-omics-binding-steps.R`, `test-gui-entry-proteomics.R`, `test-kit208-gprofiler-identifiers.R` and `test-prot-annotation-uniprot-empty.R`. The failures in the other `same` files are identical on both sides and recorded on the laptop (KIT-156). Stop rule 2 did not fire.
  - R CMD check: findings identical on both sides (3 WARNINGs, 3 NOTEs); nothing only on head.
  - New file on head: 109 passed, 1 failed. T7b's second expectation was my exact message text: dplyr's `mutate()` replaces the raised message by its header ("In argument: ..."), so the refusal ends with that header. The ticket states only the class for T7b, and the class passed. This is a result of an over-specified test expectation, not of the product: one correction pass, 66f32932. R3 runs the full suite, including R2's scope, on it.
  - G1's results digest and G5's 16 report-text digests (four templates, records with and without the outcome, no record; the limits sentence removed) are identical base against head (`r2/digests-*.txt`).
  - G2 (KIT-203's gate g1 on head): environment as KIT-189's; "match on every list" (400 of 400).
- **R3, run 45502798552: full suite on 94d57634 (head: the merge of staging d048898d plus the correction pass) against d048898d (base)** (`r3/`): 478 files `same`. The new file is only on head: 110 passed, 0 failed, so T7b's corrected row passes. No file timed out or crashed. Check findings are identical on both sides (3 WARNINGs, 3 NOTEs). This also re-runs R2's scope after the correction pass.
- **Runner use:** 24 + 13 + 215 + 1,281 = 1,533 job-minutes = 102.2 core-hours (expected 105.67, ceiling 469.33).

## The design of the GO classification (S0-1; step-0 ruling SR-1)

The GO route calls `clusterProfiler::enricher` as before. Only for a list with proteins for which `enricher` returned NULL without an error does it compute the comparison's annotated universe and m from the route's TERM2GENE (`enrichmentGoTestable`: positional second column, terms of 10 to 500 universe genes). It then classifies the list:

- universe 0: refused, VAL-170 (3);
- m = 0: refused, VAL-170 (2);
- no annotated list protein: not tested;
- otherwise m' = 0: tested, no term, VAL-170 (1).

**Grounds.**
- **Engine fidelity.** The classification never overrules a result object the engine returned. A comparison is refused only when the engine itself returned no result for its lists.
- **Tested paths byte-identical.** Every path that reaches `enricher` today reaches it unchanged. That includes the random-number stream: `check_gene_id`'s `sample()` still runs for a C2 list, so P19's note no longer applies. Direct callers of `processEnrichments` are unchanged.
- **The decided outcomes are asserted on the real engine.** On enrichit 0.2.0, NULL arises only for n = 0, m = 0 and m' = 0 (G0 (i) to (iv)). The standing tests run the installed engine through a wrapper that logs whether each call returned NULL:
  - T1: C2's NULL, not tested;
  - T4: C3, refused, every call NULL;
  - T10: C4, refused with the reason, every call NULL; it failed first (base: completed);
  - T8: C5's NULL, no signal;
  - G1: every call returns a result.

  A later engine that changed these conditions would fail these rows; it would not change the output silently.

**Fixture limit, for the end-of-project re-record list.** The W1 `kernel_clusterprofiler` capture (and its W113 case) records a fake engine's result (`fake_enricher`, `tests/testthat/helper-e2e-enrichment-doubles.R`) on a six-protein universe that the real engine could not test (m = 0). That is a fidelity limit of the fixture, not of the product. It is not a reason to re-freeze (ruling SR-1).

## Acceptance criteria

| # | Criterion | Evidence | State |
|---|---|---|---|
| A1 | G0: engine conditions and fixture families | `g0a/engine.tsv`, `g0b/g0.tsv` | met |
| A2 | failing-first rows fail at base values; guards pass; all pass on the fix | R1; R2 (all but T7b's over-specified text, corrected 66f32932); R3 110/0 | met |
| A3 | R2: every file `same` except the new file | `r2/compare.tsv` | met |
| A4 | KIT-203 gate G1 on head | `g2/g1-summary.tsv` | met |
| A5 | R3 full suite on the merge | `r3/compare.tsv`: 478 same, new file only on head | met |
| A6 | NEWS, ledger; no version, `NAMESPACE`, `man/` change | diff of c9219652 against f93dfd85: none of `DESCRIPTION`, `NAMESPACE`, `man/` changed | met |

## For KIT-216 (VAL-175; KIT-215 goes first)

KIT-216, merging second, must:
- carry `not_tested` and `no_list` in the product's record config (KIT-215 commits them in the `enrichment_completed` state's config, beside `outcomes`, each only when it has rows: `saveWorkflowEnrichmentState`);
- return them from the product path of `restoreWorkflowEnrichment` (KIT-215 reads them from the state's config: `getWorkflowStateReference`, `readComponentMetadata(...)$config`);
- keep the tab rule: the fields `enrichment_not_tested` and `enrichment_no_list` are assigned only when present and cleared only when they exist (`holdOmicEnrichOptionalFields`, `R/mod_omic_enrich.R`). A clearing site KIT-216 adds (dropping the product when its table or DA changes) clears both by that rule;
- re-run in its R2: `test-kit215-not-tested-lists.R` (its record accessor `k215EnrichmentRecord` reads the product's record first), `test-kit208-gprofiler-identifiers.R`, `test-prot-annotation-uniprot-empty.R`, `test-omic-enrich-module.R`.

KIT-216 also adds a refusal in `runOmicEnrichAnalysis` (`R/mod_omic_enrich.R`), after KIT-215's edits to `holdOmicEnrichRun` in the same file.

## Limits (observed or by reading)

- The report bullet splits the outcome `enrichment_lists_not_tested` on "; " (KIT-208's format). A comparison display name containing "; " would be split wrongly (ruling SR-4 (b)). No change.
- An error caught after dplyr wraps it (g:Profiler's per-comparison handler) is reported with dplyr's header ("In argument: ..."), not the raised message, and the sentence then ends with two full stops. The run is refused all the same (VAL-170 (4)). Wording is for the owner's pass (VAL-41).
- On g:Profiler, "no significant term" cannot be told from "no recognised name" (VAL-170 (6)). The report says so; plan section 13 row (ticket 12).
- When one GO list raises an error and the other direction has no list, the run is refused as an error. With both lists present, a comparison-level refusal is judged first.
- An absent outcome means "unknown" in a record made before KIT-215, and "every list tested" in a record made after it.

## Not checked

- Real data (C4 and C5 frequency); a later enrichit version.
- KIT-216's final code.

## For the change report

**What was wrong, and who was affected.** In proteomics enrichment, a list that was never tested looked exactly like a list in which no term was enriched, in the tab, the download and the report:
- On g:Profiler, a list whose proteins all lack a gene name. KIT-208 already stopped sending such lists and gave a notice, but the tab, the download and the report still showed them as empty results.
- On the Gene Ontology route (organisms g:Profiler does not serve), a list none of whose proteins has a Gene Ontology annotation the analysis tests.

Three further cases also completed as if nothing were enriched, when nothing could have been tested or the test had failed:
- a GO comparison whose tested proteins have no annotation at all;
- a GO comparison with too few annotated proteins for any term of 10 to 500 of them;
- a list whose test raised an error (both routes; on g:Profiler also an error in a comparison's processing). The error was only a warning.

A direction in which no protein passed the cutoffs read "0". Affected: g:Profiler users with unnamed proteins in one direction; GO-route users with incomplete or mismatched annotation; users of either route whose test failed.

**What changed.**
- Such a list is labelled "Not tested" with its reason in the tab (plot placeholder and summary), the download's summary and the four proteomics report templates (a bullet). It is recorded as the outcome `enrichment_lists_not_tested` on both routes.
- A direction with no list says "no list (no protein passed the cutoffs)".
- The three comparison-level or error cases now refuse the run with a message: classes `workflow_enrichment_no_identifiers` and `workflow_enrichment_test_error`.
- The g:Profiler report paragraph states that a list with no significant term cannot be told apart from one in which the service recognised none of the names.
- No term, p-value or list of a tested list changes. A list that was tested and in which no tested term holds a list protein is still shown as tested, with no term. Records of the refused cases made earlier no longer replay.

**Code references.** At the base f93dfd85:
- `R/func_general_enrichment_orchestration.R:341-417` (g:Profiler lists and handlers), `:564-641` (GO lists and handlers);
- `R/utils_workflow_dia_da.R:678-757` (kernel and KIT-208's judge);
- `R/utils_workflow_enrichment.R:1070-1079`, `:1184-1215`, `:1302-1315`;
- `R/mod_omic_enrich.R:669-680`;
- `R/mod_omic_enrich_panels.R:111-193`, `:288-406`, `:544-640`, `:689-699`;
- `inst/reports/proteomics/DIANN_report.rmd:869-893`.

At the fix c9219652:
- orchestration `:173-199` (`signalEnrichmentListStatus`, `enrichmentGoTestable`), the GO classification `:700-722`;
- `R/utils_workflow_dia_da.R:744-800` (`diaEnrichmentNoIdentifiers`);
- `R/utils_workflow_enrichment.R:1226-1227` (the committed fields);
- `R/mod_omic_enrich.R:690-698`;
- `R/mod_omic_enrich_panels.R:115-143`, `:385`, `:457`, `:616-631`;
- `DIANN_report.rmd:892-911` (the other three templates have the same blocks);
- `R/utils_workflow_analysis_record.R` (map entry `enrichment_lists_not_tested`).

**Evidence.**
- Run 45493249015: the engine's four NULL conditions on the runner's enrichit 0.2.0; every new test failing first at its base value.
- Run 45495664583: the fixture families.
- Run 45495673217: R2, 25 regression files `same`; check findings identical; KIT-203's gate matching on all 400 lists; GO results and report texts of runs without such lists identical base against head.
- Run 45502798552: R3, 478 files `same`.
- Limits: synthetic and E. coli fixtures only; no real data; enrichit 0.2.0 only; how often the GO cases arise in real projects was not measured.
