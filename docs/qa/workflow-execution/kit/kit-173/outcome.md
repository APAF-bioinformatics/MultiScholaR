# KIT-173 outcome: Gene Ontology parent terms on the proteomics UniProt route

Ticket revision 1.5. Branch `kit-173-go-propagation` from `v0.9-staging` a60c4e9a. Every run on `run-tests`
(R 4.6.1, reference BLAS, ubuntu-22.04); run ids, commits, minutes and core-hours in `runs.tsv`.

| Commit | Content |
|---|---|
| 36a34cab, 32ffcf04 | F1 script (`f1-feasibility.R`) and its correction (R0b) |
| 536080c5 | tests only: `test-go-propagation.R` (T1 to T11, fixture hashes), `helper-go-propagation.R`, fixtures, T7(a) capture script, evidence scripts; `R/` equal to base |
| 62a41ba1 | T7(a) capture from R1 under `tests/testdata/go-propagation/replay-base/` (no store, see A4); `R/` equal to base. Base of R2 and R2b |
| e910e19f | the fix: `R/func_general_geneset_kernel.R` (new), orchestration, record, operation, workflow, annotation and report changes, workbooks, Rd, `Collate`, NEWS, ledger row |
| a3d6d84f | D9 script |
| 93c3b946 | correction pass after R2: three existing test literals, T2(d)'s fixture, `t7-replay.R`, R2 evidence |

## Observed

- **F1 (R0, R0b)**: every value of section 2.2 reproduced; no stop rule triggered; details and the two departures (no
  real root-only protein, so one is synthetic; GMT committed as ids and genes) in `feasibility.md`.
- **Failing first (R1, 536080c5)**: `test-go-propagation.R` 67 passed, 163 failed, 12 errors, 1 skipped (T7(b), no
  capture yet). Every new test except T11(e) (regression) and the fixture-hash test (which must pass on base) failed on base: by error where the function is
  absent (T1, T2, T3(a), T4, T6, T9 replay), at the stated value elsewhere: T5 median C 0.29 < 0.90 and gain 0 < 0.20;
  T3(b) attribute NULL; T7(c) nine names, no `go_propagation`; T8 no term outside the direct annotations; T9 outcomes
  absent; T10 errored in the unrelated direct-set defect (all-NA aspects) before any of its expectations was evaluated, so its failing-first at the stated value is not shown; it passes on head (12 expectations); T11(a) the overlay has no
  correction source, (a)/(d)/(f) the old text. `r1-...` failures in `r2/` compare files (base side).
- **Fix (R2 a3d6d84f, correction R2b 93c3b946; base 62a41ba1)**: in R2b `test-go-propagation.R` 317 passed, 0 failed,
  1 skipped (T7(b), store not present in the test job). R CMD check base and head: 0 ERROR, 3 WARNING, 3 NOTE, nothing
  only on head (`check-compare.txt`). Of the other 18 files: 13 `same`; `test-workflow-omics-binding-steps.R` and
  `test-workflow-operation-registry.R` `count-changed` with the same failures on both sides (laptop-recorded
  references); `test-omic-enrich-module.R` (6) and `test-omics-enrich-core.R` (3) new failures, all in
  clusterProfiler (taxon 999999) lanes, listed below.
- **T5 (R2 script `t5-agreement.R`, the route's chain, TERM2GENE as passed to `enricher`)**: 2,703 terms in scope.
  Shipped value: median C 0.933 (IQR 0.833 to 1.000), median P 1.000 (0.95 to 1.000); "none": median C 0.073, P 1.000
  (960 terms empty). All three conditions met. Two-relation option (descriptive): median C 0.893 over all terms; 0.935
  over the 1,592 terms with no descendant reached only through a regulates edge; 0.810 over the other 1,111; P 1.000.
  Mapping loss 0 of 6,589 g:Profiler genes; no reference term outside GO.db. `t5-agreement.tsv`, `t5-summary.tsv`,
  `t5-two-relation.tsv`, `t5-lowest.tsv`.
- **T7 (R2b script `t7-replay.R`)**: base 62a41ba1 installed into a private library in the job (its namespace has no
  `propagateGoAnnotations`), `capture-record.R` run on it, then `test-go-propagation.R` on the installed head with the
  capture: all tests pass except T7(b), 9 of 11 expectations passed, **2 failed**. Which two is not known: the test
  console stayed in the job by the script's design and only the counts were uploaded (`r2b/t7-tests.tsv`).
- **T7 instrument run 1 (VAL-54, run 45141146856, commit 7f4f4060)**: the instrument defect of R2b (only counts
  returned) fixed by writing the failing expectations (test, class, first 12 message lines) to `t7-failures.tsv`;
  nothing else changed. T7(b)'s two failures: `evidence$tables` (the `geneID` strings of the up table differ in order
  only, e.g. row 1 has the same 10 accessions) and `evidence$files` (the up and down results TSV SHA-256 differ: b50f9813...
  against dba72b18..., a13410fa... against e1c525b8...). Every other T7(b) expectation passed, TERM2GENE included.
  Reading (inferred, not tested): enrichit writes `geneID` in the order of the `gene` vector, which `processEnrichments`
  takes from the DA table's row order; the TERM2GENE is identical, so the order of the replayed DA table differs from
  the original's. Whether a base-on-base replay shows the same was not run.
- **Effect (R2e, run 45127866711; D8 accepted, VAL-65)**: one proteome, the reviewed E. coli K-12 entries of UniProt
  2026_04 (4,024 annotated proteins, the background for every setting; no root-only protein, `go_root_only_proteins`
  0), synthetic lists, 20 per scenario, the same lists for none / two relations / five relations. Gene-set pairs 21,882
  / 127,808 / 133,420; terms in the testable range (10 to 500 annotated proteins) 314 / 1,242 / 1,249 (947 enter at
  10, 12 leave at 500, five relations). Median terms tested per list: specific 21 / 123.5 / 123.5, mid-level 129 / 638
  / 646.5, broad 206 / 963.5 / 962.5, null 144.5 / 692.5 / 694; median significant (q-value below 0.05): specific 1.5 /
  20 / 19, mid-level 10 / 53.5 / 46.5, broad 12.5 / 129.5 / 129, null 0 / 0 / 0. With parent terms the planted term is
  significant in 20 of 20 lists in every scenario; without them the specific and broad terms are not tested at all
  (fewer than 10 directly annotated proteins) and the mid-level term is significant in 20 of 20. Random lists: 1 of 20
  with any term significant by q-value under each setting (95 % interval 0.1 % to 24.9 %), by Benjamini-Hochberg 0 / 1
  / 1 of 20; without parent terms that list (draw 5) had 120 of its 156 tested terms significant by q-value and none by
  Benjamini-Hochberg, with parent terms (draw 6) 4 by either; per list the two flags differ by up to 120 terms without
  parent terms and by at most 3 with them. Most significant terms under propagation are neither the planted term nor
  its biological-process ancestors or descendants; every molecular-function and cellular-component term counts as
  neither, and each planted list is half random proteins. Descriptive figures for one proteome and these lists, not
  estimates of power or of the false-discovery rate. `effect.tsv`, `effect-draws.tsv`, `effect-sets.tsv`.

## Declared control: attribution of T7(b)'s difference (authorised by the coordinator, 2026-11-23)

Purpose: to learn whether the difference found in T7(b) (the order of proteins inside the `geneID` cells of the two
saved enrichment result tables, and therefore the two results files' hashes, when a record saved on the base is
replayed on the head) is caused by KIT-173's change or exists in replay independently of it. Design, in one scripts
request, changing no test and no product code: (1) capture on the BASE commit, replay on the BASE commit (the base
installed into a private library as the one-job design already does), compare the same eleven expectations; (2)
capture on the HEAD commit, replay on the HEAD commit, compare the same eleven; and, to test the reading of the cause
directly, (3) in each pair record whether the differential-abundance table the enrichment reads has the same ROW ORDER
at capture and at replay (a yes or no and the count of rows whose position differs; no feature-level rows leave the
job), and whether the `gene` vector handed to the enrichment call is the same sequence. Outputs: a small table with one
row per pair and expectation (pass or fail), the two order facts per pair, and the failing expectations' first 12
message lines as before. Ceiling: one job, at most 2 core-hours. Decision it supports, fixed now: if base-on-base shows
the same difference, the cause is in replay as it stands on `v0.9-staging` and not in this ticket: A4 as worded ("a
project saved on the base replays unchanged on the head") is then not met for a reason that predates the ticket, and
what follows (a product fix of replay's row order as its own ticket; or A4 restated as "replays as the base itself
replays it", with the defect recorded) is the owner's decision, not the implementer's or the coordinator's. If
base-on-base and head-on-head are both identical and only base-on-head differs, the change is KIT-173's: report exactly
which code path reorders, and stop; the fix is then a code change after a result, which comes back to the owner under
the one-correction rule. If the control itself fails on mechanics, the instrument-defect rule applies (T7i1 was
instrument run 1 of 3).

Implementation: `t7-control.R` (the job) and `t7-control-pair.R` (one capture and its replay in one R process, on the
library given). The eleven comparisons are T7(b)'s, computed in the script on the same objects (the test file is not
run or changed). The capture and the replay each record, through wrappers that call the real functions, the protein
id column of every DA table `processEnrichments` receives and every `gene` vector handed to `clusterProfiler::enricher`;
only the comparison (same order yes or no, rows at another position) leaves the job.

Result (run C2, run 45145792186, see `runs.tsv`; C1 was an instrument defect, fixed under VAL-54, instrument runs 2 and 3
of 3 with T7i1 the first). Observed, `control/`:

| Pair | Eleven expectations | DA table read by the kernel, capture against replay | Gene vectors handed to enricher |
|---|---|---|---|
| base on base (base in a private library; no `propagateGoAnnotations`) | 11 of 11 pass | 70 rows, same order, 0 at another position | calls 1 (25) and 2 (10): same sequence |
| head on head (installed head) | 11 of 11 pass | 70 rows, same order, 0 at another position | calls 1 (25) and 2 (10): same sequence |

The script exited 1 after writing every output, on its last line (a print indexing a two-column table with `[, 1:4]`);
no output is affected and no further run was made (the instrument runs are used up).

Branch of the declared decision: base-on-base and head-on-head are both identical, and base-on-head differed (T7i1),
which is the branch "the change is KIT-173's". It is not established, because the two designs differ in a way the
control did not hold fixed: in each control pair capture and replay ran in one R process outside testthat, while in
T7i1 the base capture ran in its own Rscript process and the head replay in a second process under
`testthat::test_file` (which sets its own test context, for example the collation locale). Reading of the code
(inferred, not tested): on a record without `go_propagation` the head path differs from the base only by passing
`go_propagation = "none"` to the kernel (`buildProtEnrichProcessEnrichmentsArgs`, `runDIAReferenceEnrichment`'s extra,
the operation kernel), `checkGoPropagation`, and `buildUniprotGoGeneSetsDirect` (the base code verbatim, TERM2GENE hash
identical in T7i1); none of these reads or reorders the DA table or the gene list. `geneID` is written by enrichit's
compiled `ora_cpp`, from the gene list and the sets. No KIT-173 code path that reorders was found; the cross-process
difference is not explained. Which of KIT-173's change or the process arrangement causes it is the owner's question;
the test and the product code are unchanged.

## Declared control 2: the cross pair by collation (owner, VAL-57: "take recommended")

Purpose: to learn whether T7(b)'s difference (protein order inside the `geneID` cells when a record captured on the
base is replayed on the head, run 45141146856) comes from the process arrangement of that run (base capture in its
own Rscript process, head replay in a second process under `testthat::test_file`, which sets its own collation) or
from KIT-173's change. Design: one scripts request, `t7-cross.R` (job), `t7-cross-capture.R` (base capture process)
and `t7-cross-replay.R` (head replay process). The base commit 62a41ba1 is installed into a private library as
before. The cross pair (base capture, head replay; two processes, as in T7i1) is run twice: (a) both processes with
`LC_COLLATE=C`; (b) both processes in the runner's locale. Each variant records the eleven expectations of T7(b), the
DA table's row order at capture and at replay (same order yes or no, rows at another position), the gene vectors
handed to `enricher` (same sequence yes or no, positions that differ), the `geneID` cells that differ (term id,
number of proteins, same set yes or no; no accession leaves the job), and `Sys.getlocale("LC_COLLATE")` in every
process. The job also records the collation of the T7i1 arrangement: an Rscript process in the runner's locale (the
capture) and a `testthat::test_that` context in a second process (the replay). No test and no product code changes.
Commit bound: the commit that adds these scripts (stated in `runs.tsv`). Ceiling: one job, 1.5 core-hours.

Deciding rule, fixed now: identical in both (a) and (b) = the earlier failure was the test's process arrangement;
different in both = the difference is KIT-173's; any other pattern is reported as found. Nothing further is decided
here: no test is adjusted, no product code changes, A4 is not re-judged by the implementer. One run producing a
result; if it fails on its own mechanics with no comparison produced, the mechanics may be fixed and the run sent once
more, both recorded.

Result (run C3, 45149307211; `cross/`). Observed:

| | T7i1 (run 45141146856) | (a) both processes C | (b) both processes runner locale |
|---|---|---|---|
| Collation, capture (base) | C.UTF-8 (Rscript, runner locale) | C | C.UTF-8 |
| Collation, replay (head) | C (inside a testthat test) | C | C.UTF-8 |
| Eleven expectations | 9 pass; result tables and pathway files fail | 9 pass; result tables and pathway files fail | 9 pass; result tables and pathway files fail |
| DA table, capture vs replay | not recorded | 70 rows, same order | 70 rows, same order |
| Gene vectors to `enricher` | not recorded | calls 1 (25) and 2 (10): same sequence | same sequence |
| TERM2GENE, TERM2NAME hashes | identical | identical | identical |
| `geneID` cells that differ | up table rows 1 to 4 shown (message cut at 12 lines) | 11: up GO:0009085, GO:0005829, GO:0042802, GO:0042803, GO:0005737; down GO:0070180, GO:0022625, GO:0002181, GO:0003735, GO:0005737, GO:0005829; each the same set of proteins | the same 11 |

Pattern: different in both (a) and (b). By the rule fixed above, the difference is KIT-173's (withdrawn after
control 3, below: the same difference occurs on the base alone across two processes). Observed beside it: the
inputs this control records (the DA table's order, the gene vectors, the TERM2GENE and TERM2NAME content hashes, which
compare names and values as text in row order) are identical at capture and replay; the reordering happens inside or
after the enrichment call. Which code path does it was not determined here (no further run; nothing decided).

## Declared control 3: base capture, base replay, two processes (owner, VAL-61)

Purpose: to test the reading that the `geneID` order is set by string memory addresses in Rcpp's `unique()`
(enrichit `ora_cpp`), not by KIT-173. Design: `t7-basebase.R` installs the base 62a41ba1 into a private library as
before and runs `t7-cross-capture.R` and then `t7-cross-replay.R` (variant `basebase`), each in its own Rscript process
and both on the base, in the runner's locale; the same columns as control 2 (eleven expectations, DA row order, gene
vectors, differing `geneID` cells, collation and whether the package has `propagateGoAnnotations` in each process).
No test and no product code changes. Commit bound: the commit adding `t7-basebase.R` (in `runs.tsv`). Ceiling: one
job, 1 core-hour. Rule fixed now: if the `geneID` cells differ as in C3, the cause is shown and the record withdraws
the attribution to KIT-173; if identical, stop and report (the diagnosis is then wrong). One mechanics rerun at most.

Result (run C4, 45151486004; `basebase/`). Observed: both processes on the base (no `propagateGoAnnotations`), collation
C.UTF-8 in both; 9 of 11 expectations pass, result tables and pathway files fail; the DA table (70 rows) and both gene
vectors (25, 10) in the same order at capture and replay; the same 11 `geneID` cells as in C3 differ, each the same
set of proteins. By the rule fixed above the cause is shown: the order inside a `geneID` cell is set per R process
(enrichit 0.2.0 `src/enrichit.cpp:24-28, 69-76` iterate an `unordered_set` filled from `Rcpp::unique()`, whose order
follows the hash of each string's memory address, Rcpp 1.1.2 `IndexHash.h:111-117, 232-238`), and occurs without
KIT-173. **The attribution of C3 to KIT-173 is withdrawn.**

## The geneID fix (VAL-61), run R2f

Commits cf4b5fd4 (T12 and scripts, `R/` unchanged) and 19d24b96 (`sortEnrichmentGeneIds()` after each `enricher` call
in `processEnrichments`; T7(b) compares pre-fix records' `geneID` cells as sets and their results TSV files by name;
NEWS; ledger row). Run 45153295682 (base cf4b5fd4, head 19d24b96), observed:

- T12 fails on cf4b5fd4 (112 `geneID` cells not in alphabetical order) and passes on 19d24b96.
- Two Rscript processes, capture and replay both on the head (`t7-headhead.R`, `r2f/headhead-*.tsv`): 11 of 11
  expectations pass, the result tables and the results files' SHA-256 included; no `geneID` cell differs.
- T7(b) as restated, with the base store (`t7-replay.R`): 12 of 12 expectations pass.
- `test-general-enrichment-shared.R` same; `test-omic-enrich-module.R` 6 to 11 and `test-omics-enrich-core.R` 3 to 7
  failures, every one in a clusterProfiler (taxon 999999) lane's capture (the live lane's stand-in table now has its
  `geneID` cells sorted): "(b) in process: clusterprofiler", "... kernel_clusterprofiler", "(KIT-113 W1) ...",
  "(W2, core side) ... clusterprofiler run step", "... kernel_clusterprofiler run step". For the plan's re-record row.

## Full suite (R3, run 45219690402, head 4bc1a361 against base 62a41ba1)

449 files. 444 `same` (36 of them fail identically on both sides: stored references recorded on the laptop and
runner-environment failures, as in KIT-171's R3). Three of the 36 (`test-integration-input-transcriptomics.R`,
`test-kit129-run-path.R`, `test-trx-ruv-designs.R`) fail on both sides because this request did not set the full-suite
template's `sva`, `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` and 120-minute timeouts; they pass with them (KIT-180 R5) and call no
code this ticket changed. `test-go-propagation.R` `fixed` (the new tests: 77 passed and 176
failing or erroring on base, 319 passed and 0 failing on head, 1 skipped: T7(b), which needs the base store and runs in
`t7-replay.R`). Two `new-failure`, every failure a capture of a clusterProfiler (taxon 999999) lane:
`test-omic-enrich-module.R` "(b) in process: clusterprofiler" (1), "(b) in process: kernel_clusterprofiler" (8; one of
them also fails on base), "(KIT-113 W1) each W1 (b) session's committed record ..." (2); `test-omics-enrich-core.R`
"(W2, core side) ... clusterprofiler run step" (1) and "... kernel_clusterprofiler run step" (6). They move with the
propagated sets, the step's outcomes, the declared `go_propagation` and the sorted `geneID` cells (6.2; the plan's
re-record row). Two `count-changed` with the same failure on both sides: `test-workflow-omics-binding-steps.R` (the
KIT-094 characterisation fails on base and head alike; the passes grow by the two KIT-173 assertions) and
`test-workflow-operation-registry.R` (two failures on both sides). No 9606 capture moved. R CMD check: base and head
0 ERROR, 3 WARNING, 3 NOTE; nothing only on head (`r3/r3-check-compare.txt`).

## Acceptance criteria

| # | Status | Evidence |
|---|---|---|
| A0 | met | VAL-28 |
| A1 | met | `feasibility.md`, `f1/` |
| A2 | met except T10 (failing-first by an unrelated error) | R1 (base failing), R2b (head passing); `runs.tsv` |
| A3 | met | `t5-summary.tsv` (three conditions TRUE), the two-relation report, the 40 lowest with causes where identifiable, `results.md` |
| A4 | **met as restated (VAL-61)**: records made after the fix replay byte-identical across processes (R2f head-on-head 11 of 11); the base record replays with every expectation equal and its `geneID` cells equal as sets (T7(b) 12 of 12). History:  the base record replays on head with the same nine declared names, product parameters, service, contrasts, TERM2GENE and TERM2NAME hashes and plot-data digests; the stored result tables differ in one column, `geneID` (each row lists the same proteins in another order), and so the two `*_enrichment_results.tsv` files' SHA-256 differ. Found by instrument run 1 (VAL-54); stopped and reported, the test not adjusted | `t7i1/t7-failures.tsv`, `r2b/` |
| A5 | met: shown to the owner, accepted (VAL-65); written into NEWS and the ledger row | `effect.tsv`, NEWS, ledger |
| A6 | met: full suite base against head, every difference a new test or a clusterProfiler capture (listed); check unchanged | `r3/` |
| A7 | met | diff: NEWS under the unnumbered heading, no version change, Rd by hand, `Collate` |
| A8 | met | T11(a) to (f) pass on head in R2b; (e) passes on base too |

Captures that moved (all clusterProfiler, taxon 999999; for the plan's section 11 re-record row; none edited):
`test-omic-enrich-module.R` "(b) in process: clusterprofiler" and "(b) in process: kernel_clusterprofiler" (3) at
`:3656`, "(KIT-113 W1) each W1 (b) session's committed record ..." (2: the state's `config$outcomes` digest);
`test-omics-enrich-core.R` "(W2, core side) runWorkflowEnrichment gives W1 (b)'s clusterprofiler run step" (product
gains `parameters$go_propagation`) and "... kernel_clusterprofiler run step" (2: TERM2GENE hash). Test literals changed
(not references): `test-omic-enrich-module.R:4261-4262`, `:4295` and `test-workflow-analysis-record.R:1001` compared a
9606 run's declared names with the whole parameter list, which now ends in `go_propagation`; the expected side drops it
(the runs' names are unchanged).

## Inferred, not observed

- The g:Profiler route's computation is unchanged: no 9606 capture moved in R2b; the D9 texts differ only in the
  paragraph.

## Unrelated (recorded, not changed)

- `processEnrichments`' direct-set code (unchanged, now `buildUniprotGoGeneSetsDirect`) fails with "arguments imply
  differing number of rows: 0, 1" when a GO aspect column has no annotation at all (`cbind` of an empty table).

## Not checked

- The per-source correction statement rests on the gprofiler2 vignette; the 2023 and 2020 papers' full texts could not
  be fetched (HTTP 500).
- `R CMD check` at the records commit.
