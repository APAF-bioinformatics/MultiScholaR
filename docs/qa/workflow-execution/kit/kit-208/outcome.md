# KIT-208 outcome: STOPPED at stop rule 2 (a stored capture moved in R2)

Ticket r1.3 (`tickets/kit/KIT-208.md`). Branch `kit-208-gprofiler-gene-names`, cut from `v0.9-staging` 241914cb
(KIT-203 merged). Runs in `runs.tsv`. R3 (the full suite) was not sent: stop rule 2 stops the ticket first.

## Commits

| Commit | What |
|---|---|
| 0a524465 | tests only: `tests/testthat/test-kit208-gprofiler-identifiers.R` (T1 to T9); nothing under `R/` |
| 60e49a57 | the fix (5 files under `R/`), the ledger row `VAL-<merge>`, the stub signature in `test-general-enrichment-shared.R` (D2) |
| 9dc2bd64 | NEWS entry; R1 records |

## What was observed

- **R1 (run 45367973598, 0a524465).** Every failing-first test fails at the ticket's base value, at its first
  stated comparison: T1 and T2 blank on 70 of 70 rows; T5 10 calls, each (0, 0); T6 completed, 10 calls, each (0, 0);
  T7 five (0, 45) then five (10, 45); T9 5 calls with 0 ids. Guards T3, T4, T8 pass. T6 also errors after its stated
  comparisons (`conditionMessage` of a NULL error, because base completed): test mechanics after the comparison, no
  effect on head. Stop rule 1 not met (`r1/`).
- **R2 (run 45370085901, 9dc2bd64 against 0a524465, 44 files, check both sides).** `r2/compare.tsv`:
  - the new file: 26 passed / 20 failed / 1 error on base, 49 / 0 / 0 on head (`fixed`);
  - 41 files `same`, including the five proteomics GUI-entry lanes (`test-gui-entry-proteomics.R`, 4 failures
    identical on both sides), `test-omic-enrich-module.R` (19 identical), `test-prot-annotation-uniprot-empty.R`,
    `test-general-enrichment-shared.R` (the D2 stub), `test-workflow-analysis-record.R`;
  - check: 0 errors, 3 warnings, 3 notes on both sides, identical;
  - **one `new-failure`: `test-omic-da-module.R`, test "(b) in process: prot_load"** (passed 28 to 27). The moved
    leaf: the stored SHA-256 of `da_output/da_proteins_groupB-groupA_long_annot.tsv` in the W1 (b) `prot_load` run
    files, `fbe4b6a9...` (frozen, equal on base) to `180c3306...` on head (`r2/compare-failures.txt`).

## Why the capture moved (by reading; no run beyond R2)

The `prot_load` case loads a session file written by `module_ci_prot_da_write_session`
(`tests/testthat/helper-module-ci-prot-da.R:156-167`): its `uniprot_dat_cln.RDS` is `Entry` plus `gene_names`, and
the object's identifier table also carries `gene_names` (`:70`). That is case 4, the collision part (a) repairs: on
base the TSV holds `gene_names.x`, `gene_names.y` and a blank `gene_name`; on head one `gene_names` and a filled
`gene_name`. The move is therefore the intended effect of (a) on a fixture the plan audit's P03 reading did not list
(it named the module-CI writer test `test-module-ci-prot-da.R:670-697`, not the W1 capture of
`test-omic-da-module.R`). INFERRED; the TSV contents were not dumped.

## Decision needed (stop rule 2, E-4)

Stop rule 2 stops the ticket. Under E-4 the capture `omicDaW1Frozen[["prot_load"]]` (`test-omic-da-module.R:3254`)
is re-frozen only on the owner's word, with a base-against-head control (dump the TSV on base and head and show the
difference is only the gene-name columns). No capture was re-recorded here.

## Acceptance criteria

| # | Met? | Evidence |
|---|---|---|
| A1 | met | R1 fails at the stated values; R2 head passes the new file 49/0 |
| A2 | **not met** | R2: one stored capture moved (`test-omic-da-module.R` prot_load DA TSV digest) |
| A3 | not run | R3 not sent (stop rule 2) |
| A4 | partly | NEWS, ledger row; no version, `NAMESPACE` or `man/` change; this file |

## Deviations (for the implementation audit; accepted as engineering by the coordinator)

- D1: `not_tested` is a data frame (comparison display name, direction, list size), absent when empty, not the
  ticket's character vector: 7's notice needs the list size, the refusal the display name in the kernel.
- D2: the `perform_enrichment` stub in `tests/testthat/test-general-enrichment-shared.R:609` gains the two optional
  arguments; its fixed signature would have errored on the labels. No expected value changed (R2: `same`).
- D3: the refusal's `<n>` is the comparison's DA row count, read in the kernel from `da_results`.

## Wording for the owner's pass (VAL-41; drafts)

Refusal (case 5):
> g:Profiler tests gene names, and none of the proteins of comparison <contrast> has one (0 of <n> proteins named),
> so nothing would be tested. Check that the protein annotation was retrieved for this organism (save the design
> again). This happens when the annotation was not retrieved for these proteins, or when their identifiers are not
> UniProt accessions.

Notice (case 6), one per run, one sentence per list:
> Comparison <contrast>, <up|down>: none of its <n> proteins has a gene name, so this list was not sent to g:Profiler
> and was not tested (its empty result does not mean that no term is enriched).

Signal message (for callers): "g:Profiler was given no gene names to test: <list_size> proteins in the list,
<background_ids> named proteins in the background".

## Not checked

- R3 (full suite) not run; files outside R2's scope not compared.
- The moved TSV's contents were not dumped; the cause is by reading.
- Real g:Profiler, real data, case 5's reach in the app, DIA-NN group keys (VAL-147 (4), hand-off).
- KIT-205's background: T7 and T8 not aligned to it (whichever merges second runs R2b).

## For the change report

- **What was wrong.** (a) When a protein table built in a script carried `gene_names` and the UniProt annotation did
  too, the DA publication's join produced `gene_names.x`/`gene_names.y` and every `gene_name` was blank, so
  g:Profiler was sent nothing (KIT-175 F1: 70 of 70). (b) On the g:Profiler route a list or comparison with no named
  protein was sent empty; `gost` answered NULL, read as "no terms", five times per list, and the run was saved as
  completed. Who: (a) scripts and test fixtures (no application import puts `gene_names` in the identifier table);
  (b) g:Profiler users whose annotation names none of a comparison's proteins (refused now) or none of one direction's
  (reported not tested now).
- **Code.** Base 241914cb: `R/func_prot_da_results_methods.R:66-70` (join), `R/func_general_enrichment_engine_helpers.R:10-127`,
  `R/func_general_enrichment_orchestration.R:234-235`, `:319-329`, `:353-363`, `R/utils_workflow_dia_da.R:681-691`,
  `:770-771`, `R/utils_workflow_enrichment.R:995`. Fix 60e49a57: join `R/func_prot_da_results_methods.R:67-81`;
  `perform_enrichment` arguments `:21-22`, size `:30`, check and signal `:71-86`; `purrr::imap` `:236`, labels
  `:330-331`, `:366-367`; kernel `R/utils_workflow_dia_da.R:689-700`, helpers `:712`, `:729`; run
  `R/utils_workflow_enrichment.R:1024-1037`, helper `:830`.
- **What changed in results.** Script-shaped DA tables regain gene names (the annotation's first, the table's for the
  rest); script tables sharing other columns with the annotation now name them `X.input`/`X` instead of `X.x`/`X.y`.
  A g:Profiler comparison without a named protein is refused; a direction without one is not sent, notified and
  recorded (`enrichment_lists_not_tested`); unnamed proteins per comparison are recorded
  (`enrichment_unnamed_proteins`). A record of a case-5 run made earlier no longer replays. Left unchanged (E-3): the
  legacy writer `.outputDaResultsAllContrastsProteomicsLegacy` and the workbooks' writer `outputDeAnalysisResults`
  keep the default-suffix join.
- **Evidence.** R1 45367973598, R2 45370085901 (`r1/`, `r2/`); recorded `gost` stand-in, no service call. R2 also
  showed the module-CI W1 `prot_load` DA TSV digest move (above): awaiting the owner's re-freeze decision.
