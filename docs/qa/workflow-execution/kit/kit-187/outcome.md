# KIT-187 outcome: stopped at R2 (a stored reference moves)

State: **stopped before merge at R2**, by section 6.2 ("A `new-failure`, `message-changed` or `compare-error`
anywhere stops the ticket") and 6.2's stored-reference rule. R3 was not run. NEWS and the ledger row are not yet
written (they carry S2's figures and wait for the owner's ruling below).

Commits: tests-only 4a58a3be; fix 07c97b54; merge of `v0.9-staging` 604ec4dd (KIT-186 and the version/NEWS
restructuring; transcriptomics code only, none on S2's path) e6c42d1d. First branch point 98df468b; R2's base 604ec4dd.

## Observed

### R1, run 45272974170, tests-only commit 4a58a3be (`r1-tests.tsv`, `r1-failures.txt`, `s2-lc/`)
- `test-prot-missing-value-counts.R`: 39 passed, 22 failed, 0 errors. Each failure is an assertion at the value the
  ticket states: T1 kept `P_full` only; T2 `P_failOne` removed; T3 the failing set; T4 `P_oneB` removed; T6 the
  55 % row; T7 kept `P_full` only; T8 (i) `exists()` FALSE; T8 (iii) marker absent; T8 (iv) (b) and (c) markers
  absent, and, because the base has no rule option, (b)'s proteins 1 to 5 removed and the recorded rules of (3) and
  (4) "counts" where "rounded_percentage" is expected (8 assertions); T9 `protein_count` 1; T12 five times no error.
  T5, T6's 50/60 rows, T8 (ii) and every replayed state's `identical()` passed on the base.
- `test-prot-qc-peptide-intensity-legacy-counts.R`: T10's percentage call removes `TWOOFTHREE` (1 failure); the
  count call passes.
- `test-prot-summary-report-templates.R`: T11 fails at "exceeded the allowed share" (12 failures, 4 per template).
- The other three files: ok.
- A2: T3's base failing set (`r1-t3-base-set.txt`) is 62 (n, m) pairs and 6 (G, g) pairs, identical to section
  2.2's list for n up to 20 (and to the indicative arithmetic for G up to 7). No correction of 2.2 is needed.
- S2 lifecycle (`s2-lc/`, 111 s): falsifiers 1 and 2 are 0 in every cell and mask; falsifier 3 is 2 in cell A at
  vd5, vd10 and vd25 (the positive control). Fault run (`s2-lc/fault/`, 28 s): stopped with class `k187_falsifier`
  after writing `falsifiers.tsv` with falsifier 2 = 1 in K/vd5; falsifier 1 is also 1 there, since the one dropped
  protein is by definition one "removed by head and kept by base".

### R2, run 45278755270, head e6c42d1d against 604ec4dd (`r2-compare.tsv`, `r2-compare-failures.txt`, `r2-check-compare.txt`, `s2/`)
- Verdicts over the 45 files: 41 `same`; `count-changed` `test-prot-summary-report-templates.R` (T11's 12
  expectations added, none failing); `only-on-head` the two new files (67 and 3 passed, 0 failed);
  **`message-changed` `test-prot-norm-core.R`** (base 272 failed, head 306).
- What moves: the three KIT-111 tests (W1 (a), W2 (W1 (b)), W2 (R7)) compare each case with frozen values written in
  the test file itself (`protNormW1Frozen()`, `tests/testthat/test-prot-norm-core.R:558` onward; the capture's
  `provenance` field is `protNormW1Sha(protNormW1Provenance(snapshot, root))`, `:485`). On head the `provenance`
  field differs for eight cases (`dia_normalization_manual`, `_skip`, `_skip_repeat`, `_refused`, `_retry`,
  `dia_correlation_apply`, `_reapply`, `_skip`), in-process and through the worker: 34 expectations that pass on base
  fail on head. Every other field's failures are identical on both sides (the known laptop-recorded digests).
- Cause (inferred from the diff and the fields, not traced line by line): the normalisation declaration's new
  `fixed = list(missing_value_comparison = "counts")` (D2 (a+), section 7.3) is written into each normalisation
  step's provenance, and the correlation cases hash the head's provenance including the normalisation step. No
  state, product, file or record digest changes on head beyond the base's.
- R CMD check: 0 errors, 3 warnings, 3 notes on both sides; no finding only on head (A4 met).
- S2 (R = 100 per cell and mask, 76 min, workers 2 per side): no falsifier triggered (`s2/falsifiers.tsv`, all 0);
  cells K and E identical base against head (falsifier 2 = 0; head-minus-base added, removed 0). Fixtures: both
  e2e fixtures keep 4 of 5 proteins on the base (P00002 removed: its one value at or below the 1st percentile in
  KO) and 5 of 5 on head. Headline figures, mean with 95 % Monte Carlo interval, 3,000 proteins per draw:
  - Cell A (2 x 3): kept 2,442 -> 2,838 (vd5), 2,060 -> 2,668 (vd10), 1,261 -> 2,086 (vd25); power 0.650 -> 0.735,
    0.548 -> 0.676, 0.333 -> 0.506; FDP 0.050 on both sides at each mask (paired differences within +/-0.001, intervals
    include 0); of the added proteins about 14 % are called, with FDP 0.061, 0.061, 0.057.
  - Cell D (3 x 4, the group-level case): added 48, 116, (vd25 in `s2/cells.tsv`); power +0.020, +0.028.
  - Cells B (4 x 3) and C (5 x 3): added 107 / 252 / 605 and 75 / 196 / 538; power +0.019 / +0.039 / +0.085 and
    +0.012 / +0.029 / +0.069.
  - Cell E (3 and 6): no change; the per-group rule of D1 (b) would change the fate of 43, 105 and 285 proteins.
  - No replicate fell back from Storey's q-value to BH on either side.
  These are simulated data, KIT-185's generator, R = 100 draws per cell and mask; the size on real data will be
  shown by KIT-161 (Cox, 2 x 3); KIT-185's cell P2, also simulated, measures RUV-III-C on the fixed filter.

## Acceptance criteria
| # | Status | Evidence |
|---|---|---|
| A1 | met | R1 log (`r1-failures.txt`); R2 head: both new files pass |
| A2 | met: equal to 2.2 | `r1-t3-base-set.txt` |
| A3 | **not met**: `test-prot-norm-core.R` `message-changed` | `r2-compare.tsv` |
| A4 | met | `r2-check-compare.txt` |
| A5 | met | `s2-lc/`, `s2/falsifiers.tsv`, `s2/cells.tsv` |
| A6 | not run (R3 waits) | |
| A7 | not yet (NEWS, ledger wait); DESCRIPTION, packaging, NAMESPACE untouched by this branch's commits | diff |
| A8 | met | section below |

## Decision needed (the owner's)
The ticket predicted that no stored reference moves (2.5, 6.2); the KIT-111 frozen `provenance` digests in
`tests/testthat/test-prot-norm-core.R` do, because the normalisation step now records the marker. The file already
fails on base (laptop-recorded digests), but 34 expectations that pass on base fail on head. Options, none taken:
(a) re-freeze those eight `provenance` values (an owner-approved exception, by file; section 11 re-record list);
(b) keep the marker on the filter step only (D2 (a) instead of (a+)), so a session saved before this version and
continued after it may not replay exactly (5.4); (c) record the normalisation marker outside the step provenance
the KIT-111 capture hashes (a design change needing an audit).

## A8. Section 2.4 restated, citations re-read at the branch point

Read at 98df468b (the first branch point) and unchanged at 604ec4dd (the merge brought transcriptomics code only:
`R/mod_omic_da.R`, `R/utils_omic_gui_config.R` at the transcriptomics DA rows, `R/utils_workflow_trx_da.R`). Line
numbers are those of 98df468b.

| Route or entry | Where (re-read) | Uses | Affected before this ticket | After |
|---|---|---|---|---|
| Proteomics DIA, label-free, TMT: tab "Protein Intensity Filter" | tab listed for `common <- c("DIA", "TMT", "LFQ")` (`R/utils_omic_gui_config.R:2205`, `:2227`); `runProteinIntensityFilterApplyStep` (`R/mod_prot_qc_protein_intensity.R:121-139`) reaches `computeDIAReferenceProteinIntensityFilter` (`R/utils_workflow_prot_declarations.R:239-298`) | rounded cutoffs (`R/func_general_parameter_helpers.R:254-267`) and `>` (`R/func_prot_qc_filtering_helpers.R:238-241`) | yes | counts (the helper) |
| Worker handler and registry | `computeWorkflowProteinIntensityAction` (`R/utils_workflow_execution.R:2863-2866`); registration `R/utils_workflow_operations.R:778-784`; declaration `diaProteinIntensityFilterOperation` | same function | yes | counts; the declaration fixes the marker |
| Clean-up after RUV-III-C | `R/utils_workflow_dia_products.R:481`, `:492`; the tab's step `R/mod_prot_norm_ruv_helpers.R:241` | the S4 method on `@args` cutoffs | yes | counts; `diaNormalizationOperation()` fixes the marker |
| DIA with limpa | protein tabs `R/utils_omic_gui_config.R:2230-2232` (DPC quantification, normalisation, peptide count) | none | no | unchanged |
| Peptide intensity filter (GUI, rollup) | direct counts (`R/func_prot_qc_peptide_group_filters.R:332`) | integer counts | no | unchanged |
| Peptide filter, deprecated percentage-only call | `R/func_prot_qc_peptide_group_filters.R:385-390` | percentages, tolerance 1e-12 | yes, for a scripted caller | the same integer reading (D6) |
| `removePeptidesWithMissingValuesPercent` | `R/func_prot_qc_peptide_group_filters.R:662`, `:666` | percentages | not at its only caller (50 / 50) | unchanged (F5) |
| `removeRowsWithMissingValues` | `R/func_prot_qc_filtering_helpers.R:106-140` | integer counts | no caller | unchanged |
| Replay of an analysis record | `replayWorkflowAnalysisRecord` loop | the current code | yes | each step with its own rule (marker) |
| Workbooks | `inst/workbooks/proteomics/DIA_workflow_starter.rmd:869-874` (and the two others) call `updateMissingValueParameters(design_matrix, config_list, ...)` | would be the tab's rule | the chunk cannot run as written (F1, inferred, not run) | unchanged (F1) |
| Reports | `inst/reports/proteomics/{DIANN,LFQ,TMT}_report.rmd` | the stored percentage | sentence false at a tie | allowed-share wording (D3) |
| Lipidomics, metabolomics | `R/func_lipid_qc.R:66`, `:73`; `R/func_metab_qc_filtering_helpers.R:34`, `:41`; slider `R/mod_omic_qc_tabs.R:104-105`; declarations `R/utils_workflow_lipid_declarations.R:48-53`, `:97-98`, `R/utils_workflow_metab_declarations.R:52`, `:86`; help `R/utils_omic_gui_config.R:1232-1233` | the user's proportion over all samples, `< proportion`, `na.rm = TRUE` | not this defect | no code change (F3, F8 stand) |
| Peptide imputation | `R/func_pept_imputation_helpers.R:155-160`, `:467` | a user proportion, `<=` | no | unchanged |
| `proteinIntensityFiltering` (S4) | `R/func_prot_s4_qc_filter_methods.R:515` | errors before any rule | no | unchanged |
| Config default | `inst/config/config.ini:95-98` (60 and 60, `ruv_grouping_variable=group`) | explicit percentages | only at float ties | read the same way |
| Report variables | `R/utils_workflow_analysis_record.R:1023-1024` | file-sourced cutoffs | n/a | unchanged |

## Inferred, not observed
- That the `provenance` move comes from the normalisation marker alone (no other provenance writer changed; the
  filter step is not in the KIT-111 cases). Not traced by a run.
- The exactness bound (316) and the earlier-rule equivalence beyond the tested grid (n up to 20, G up to 7).

## Not checked
- R3 (full suite); R2b/R3b on KIT-177's merge (not merged at stop time).
- Other stored references outside R2's 45 files that may hash step provenance (R3 would show them).
- NEWS and ledger wording (not written).
- Real data (none used); Windows (no platform code touched); whether a managed worker can run a different build
  from the controller.
