# KIT-208 outcome: the gene-name join, and g:Profiler not asked with no named identifier

Ticket r1.3 (`tickets/kit/KIT-208.md`; decisions VAL-147, re-freeze VAL-161). Branch `kit-208-gprofiler-gene-names`,
cut from `v0.9-staging` 241914cb (KIT-203 merged), merged with staging a1b397b1 (KIT-207, KIT-196). Every run is in
`runs.tsv`. Done and ready for the implementation audit.

## Commits

| Commit | What |
|---|---|
| 0a524465 | tests only: `tests/testthat/test-kit208-gprofiler-identifiers.R` (T1 to T9); nothing under `R/` |
| 60e49a57 | the fix (5 files under `R/`), the ledger row `VAL-<merge>`, the stub signature in `test-general-enrichment-shared.R` (D2) |
| 9dc2bd64 | NEWS entry; R1 records |
| a57f3ce3, 5728b7eb, af23a090 | R2 records (stop rule 2); the control script; C1 records |
| 0e8ffd05 | correction pass: the one leaf of `omicDaW1Frozen[["prot_load"]]` re-frozen, 180c3306... to fbe4b6a9... (VAL-161, C1) |
| 189df0d7 | merge of `origin/v0.9-staging` at a1b397b1 (ledger conflict only: staging's rows first, KIT-208's after) |
| 57aecb31 | R2c records |

## What was observed

- **R1 (run 45367973598, 0a524465).** Every failing-first test fails at the ticket's base value, at its first
  stated comparison: T1 and T2 blank on 70 of 70 rows; T5 10 calls, each (0, 0); T6 completed, 10 calls, each (0, 0);
  T7 five (0, 45) then five (10, 45); T9 5 calls with 0 ids. Guards T3, T4, T8 pass. T6 also errors after its stated
  comparisons (`conditionMessage` of a NULL error, because base completed): test mechanics after the comparison, no
  effect on head. Stop rule 1 not met (`r1/`).
- **R2 (run 45370085901, 9dc2bd64 against 0a524465, 44 files, check both sides).** `r2/compare.tsv`:
  - the new file: 26 passed / 20 failed / 1 error on base, 49 / 0 / 0 on head (`fixed`);
  - 41 files `same`, including the proteomics GUI-entry lanes (`test-gui-entry-proteomics.R`, 4 failures identical
    on both sides) and `test-general-enrichment-shared.R` (the D2 stub);
  - check: 0 errors, 3 warnings, 3 notes on both sides, identical;
  - **one `new-failure`: `test-omic-da-module.R`, "(b) in process: prot_load"** (passed 28 to 27): the stored SHA-256
    of `da_output/da_proteins_groupB-groupA_long_annot.tsv` in its run step, stored `180c3306...`
    (`test-omic-da-module.R:3343`; base reproduces it), head `fbe4b6a9...`, the value stored for `prot_main`
    (`:3464`) (in waldo's diff `-` is the value produced, `+` the stored one; first reported reversed, corrected).
    **Stop rule 2: stopped.**

## What moved: control C1 (run 45377688470, script commit 5728b7eb)

Declared before sending: one scripts job, ceiling 100 job-minutes, charged to the correction row's instrument
allowance; no expected value changed. `control/control-k208.R` runs the `prot_load` case of `test-omic-da-module.R`
in three separate processes (base 0a524465; head 9dc2bd64 twice), loading each tree as the test jobs do, and keeps
the masked text of the DA table and the digest `guiEntryArtifacts` itself returns. Records in `control/`
(`comparison.tsv`, `columns.tsv`, `status.tsv`, the three dumps and their `digests.tsv`). 12.7 job-minutes.

| Check | Base | Head | Holds |
|---|---|---|---|
| precondition: base digest = stored prot_load (`:3343`) | 180c3306... | 180c3306... | yes |
| head digest = stored prot_main (`:3464`) | fbe4b6a9... | fbe4b6a9... | yes |
| head1 and head2 dumps byte-identical (md5) | f3f99496... | f3f99496... | yes |
| row count | 6 | 6 | yes |
| row order (`uniprot_acc` in order) | same | same | yes |
| columns only on base | `gene_names.x`, `gene_names.y` | | |
| columns only on head | | `gene_names` | |
| the other 15 common columns (`uniprot_acc`, `fdr_qvalue`, `raw_pvalue`, `log2FC`, `comparison`, `numerator`, `denominator`, 4 `log2norm.*`, 4 `raw.*`) byte-identical | 15 | 15 | yes |
| `gene_name` blank rows | 6 of 6 | 0 of 6 | |
| head `gene_name` = first word of head `gene_names` | | | yes |

Observed: the move is exactly the collision part (a) repairs and nothing else: base's two name columns become one,
and `gene_name` goes from blank to `GENE<i>`. On head the session-load path's table is byte-identical to the main
path's (`prot_main`'s stored digest). The comparison is of the masked text the test digests (`guiEntryMask` masks
hexadecimal-looking digit runs in the numeric columns), not of unmasked numbers. Cause, by reading: the `prot_load`
session file written by `module_ci_prot_da_write_session` (`tests/testthat/helper-module-ci-prot-da.R:156-167`) holds
an annotation `Entry` plus `gene_names`, and the object's identifier table also carries `gene_names` (`:70`): case 4.
The plan audit's P03 reading did not list this capture.

## Fixtures part (a) can move (identifier table and annotation both carry `gene_names`; by reading)

1. `tests/testthat/helper-module-ci-prot-da.R:70` with the session file of `:156-167`: `test-omic-da-module.R` W1
   `prot_load` **moved in R2** (above; re-frozen under VAL-161); `prot_main` (stored value already the head's); `test-module-ci-prot-da.R:670-697`,
   `test-omic-enrich-module.R`, `test-workflow-analysis-record.R`, `test-omics-enrich-core.R`,
   `test-workflow-integration-input.R`: `same` in R2. `test-omic-summary-module.R` (via
   `helper-module-ci-prot-summary.R`, which writes its proteomics DA table directly, not through the writer): not in
   R2's scope; `same` in R3.
2. `tests/testthat/helper-module-ci-prot-da.R:218` (`module_ci_prot_da_result_table`, a result table that already has
   `gene_names`, written directly by its consumers, not through the join): consumers `same` in R2 and R3.
3. `tests/testthat/helper-go-propagation.R:134` (`goPropRunProject`) and its stored base capture
   `tests/testdata/go-propagation/replay-base/` (record, capture, evidence, product made from the collision-shaped
   table): `test-go-propagation.R` `same` in R2.
4. `tests/testthat/test-enrichment-significance.R:196` and its stored base fixtures
   `tests/testdata/enrichment-significance/` (`s3-*`, `s6-identity.rds`, `t7-*`; clusterProfiler route, no DA-table
   digest): `same` in R2.
5. `tests/testthat/test-prot-annotation-uniprot-empty.R:645` (`k175Project`): `same` in R2.
6. `tests/testthat/helper-workflow-checkpoint-fixtures.R:228-231` (identifiers `gene_names = SYNTHETIC_<i>`): no test
   found that publishes these with a UniProt annotation; the four checkpoint test files `same` in R3.

## Correction pass (VAL-161)

The owner allowed the re-freeze of that one leaf (methods ledger VAL-161, on the evidence of C1); commit 0e8ffd05
changes only `tests/testthat/test-omic-da-module.R:3343` (180c3306... to fbe4b6a9...). No other stored reference was
changed. The branch was then merged with staging a1b397b1 (189df0d7; no shared code file changed on staging since
241914cb).

- **R2c (run 45386281915, 189df0d7 against staging a1b397b1, the same 44 files, check).** 43 files `same`
  (`test-omic-da-module.R` 1,316 passed / 53 failed on both sides: the re-frozen leaf passes on head, the old leaf on
  base); the new file only on head, 49 / 0; no new failure; check identical (0 / 3 / 3). Stop rule 2 not met.
- **R3 (run 45392959369, 57aecb31 = the merge plus records, against staging a1b397b1, every test file, 20 shards per
  side, check).** 465 files `same` (426 `ok` on both sides, 39 failing identically on both, with
  the same messages; their causes were not examined here), the new file only on head
  (49 / 0); no `new-failure`, `message-changed`, `count-changed`, timeout or crash; `test-omic-summary-module.R` and
  the four `test-workflow-checkpoint*` files `same`; check identical (0 / 3 / 3). 1,208.4 job-minutes.

## Acceptance criteria

| # | Met? | Evidence |
|---|---|---|
| A1 | met | R1: every failing-first test fails at its stated base value, guards pass; R2, R2c, R3: the new file 49 / 0 on head |
| A2 | met after the correction pass | R2: one stored capture moved (stop rule 2), re-frozen on the owner's word (VAL-161, control C1); R2c: every file `same` but the new one, which alone asserts the outcomes `enrichment_lists_not_tested` and `enrichment_unnamed_proteins` (no other file's result changed, so no other run's recorded outcomes moved where a test reads them); check identical, nothing only on head |
| A3 | met | R3: 465 `same`, the new file only on head; check identical |
| A4 | met | NEWS (one Fixes entry), ledger row `VAL-<merge>` citing VAL-147 and VAL-161; no version, `NAMESPACE` or `man/` change; this file with "For the change report" |

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

## Budget

R1 0.78, R2 15.43, C1 0.85, R2c 14.14, R3 80.56: 111.76 core-hours used, against 91.1 expected and the 449.33 ceiling
(the correction row, 222.67, holds C1, R2c and R3's re-run share). R3 alone 80.56 against its row ceiling 186.67.

## Not checked

- Real g:Profiler and real data; how often case 5 or case 6 occurs in the application (VAL-147 (4), hand-off:
  DIA-NN multi-accession group keys).
- C1 compared masked text (as the test digests it), not unmasked numbers.
- Replay of a case-5 record made before this change (refused by reading; no test).
- KIT-205's background: T7 and T8 are not aligned to it; whichever of KIT-205 and KIT-208 merges second runs R2b.
- `gost` 0.2.4's behaviour on an empty query is from source (plan audit 1.3), not run.

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
- **Evidence.** R1 45367973598 (failing first), R2 45370085901 (the one moved capture), C1 45377688470 (what moved:
  the name columns only), R2c 45386281915 and R3 45392959369 (no other difference; full suite 465 files `same`)
  (`r1/`, `r2/`, `control/`, `r2c/`, `r3/`); the tests use a recorded `gost` stand-in returning NULL on KIT-173's
  70-protein E. coli fixture, no service call; the effect on real data and real g:Profiler was not measured. One
  test reference was re-frozen on the owner's word (VAL-161): the module-CI session-load case's DA table digest,
  recorded from the collision-shaped table.
