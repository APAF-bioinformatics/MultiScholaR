# KIT-175 outcome (revision 1.6): done, ready for the implementation audit

Branch `kit-175-uniprot-annotation`. `v0.9-staging` 5f930c33 merged in at af04046f (true merge). Tests-only commits
def420a1 and 0f64e328 (nothing under `R/`), fix 698f37aa, diagnosis script 873d6e8a (docs only). Earlier history:
gate F3 (45118859939) and measurement F4 (45121855739), `feasibility.md`, `f4-measurement.md`. Runs: `runs.tsv`.

## Observed

- **R1** (45281585187, def420a1, tests-only commit): `test-prot-annotation-uniprot-empty.R` 153 passed, 124 failed,
  6 errors; every T1 to T8 block that A2 names fails on base at its stated value (no classed condition, no record
  row, empty cache written or returned, no refusal, no httr2 request made, `uniprotAnnotationNotice` absent); the pins
  T5 (a), (d), (f3) and T6 pass on base (T6's literal column names and values matched base exactly).
  `test-prot-design-core.R` fails on the D6 sentence only; `test-prot-annotation-uniprot-direct-shared.R` passes with
  both doubles. `r1/failures.txt`, `r1/tests.tsv`. F1 ran in the same request (`feasibility.md`, `r1/f1-*`).
- **Instrument defect, T7** (fixed in 0f64e328, condition (ii) of VAL-54): R1 showed `workflow_replay_invalid` on
  both T7 cases because the test passed the record's path where `replayWorkflowAnalysisRecord` takes the record as
  read (a list), as KIT-173's replays do. No expected value changed. Its base value is re-shown in R2 (base side).
- **R2** (45284015220, 698f37aa vs 0f64e328, 29 files, 4 shards per side, check): the new file 319 passed, 0 failed on
  head (base 154/123/6), `test-prot-design-core.R` fixed (37/0), 26 files same, one `new-failure`:
  `test-omic-enrich-module.R` (b) refusals, the frozen digest of the `gprofiler_plot` placeholder in session
  refusals/pre_da (steps run_before_da, da, open; 3 expectations). Check: base and head 0 ERROR, 3 WARNING, 3 NOTE,
  findings only on head 0. `r2/`.
- **R2 diagnosis** (45289730439 head 873d6e8a, 45289724427 base db167236 on branch `kit-175-diag-base`; no result
  judged): that session run alone gives the frozen digest d4fb7d73... on head and on base; the placeholder text is
  "Run enrichment analysis first" with plotly's random element ids masked as `<HEX>` (`r2diag/`).
- **R3** (45291921569, 873d6e8a vs af04046f, all 455 files, 20 shards per side, full-suite template, check): 454 same,
  1 only-on-head (the new file, 319/0); 0 new-failure, 0 message-changed. `test-omic-enrich-module.R` same (11 failing
  identically on base and head, stored references). 35 files fail identically on both sides. Check identical (0
  findings only on head). Longest job 56.6 minutes. `r3/`.

- **Isoform accessions under the field-qualified query** (implementation audit I01, live at release 2026_04, by the
  auditor): `(accession:P05067-4 OR accession:P04637) AND organism_id:9606` returns `P05067-4` and nine isoforms of
  P04637, not the canonical P05067; so an isoform term lost its protein's annotation (silently when the answer fit on
  one page). Corrected in the correction pass: isoform suffixes are removed before the query (see the correction
  section below).

## Inferred (not verified)

- The R2 new-failure is a flake of the test's mask, not a behaviour change: the head code gives the frozen value alone
  (diagnosis runs, both sides) and in the full suite (R3, `same`), and in the same R2 steps the identical
  `clusterprofiler_plot` placeholder kept the frozen value. The masked text's only random content is plotly's element
  id, which `guiEntryMask` masks only at 10 or more hex characters (`tests/testthat/helper-gui-entry.R:60`); a shorter
  id escapes it (inferred). The stored value is correct and is not for re-record. An instrument repair (masking
  plotly's `visdat`, `cur_data` and `attrs` ids by position) is a candidate for the plan's instrument list, not this
  ticket.
- F1's g:Profiler identifiers without annotation came from the fixture's own `gene_names` (feasibility.md).

## Acceptance criteria

| # | State | Evidence |
|---|---|---|
| A1 | met | F1 table and F3/F4 counts in `feasibility.md`, `f4-measurement.md`; 2.3's answers: the dependents grep was re-run at the merged base (29 files in R2); line numbers changed by the merged tickets are cited at the merged base in the ticket's section 7; F1 did not run the summary report, and its g:Profiler rows describe KIT-173's fixture, whose identifier table carries `gene_names` (no import in the app does; I02, I05) |
| A2 | met | R1 (base values) and R2 (base side 0f64e328 against head) |
| A3 | met | T6 passes on base (R1) and head (R2, R3) |
| A4 | met | R2 and R3 above; R2's one difference explained as a flake of the test's mask (not for re-record); R3 0 differences; check identical; KIT-173 merged first: T6's literal values were run on the merged base (R1) |
| A5 | met | NEWS under the unreleased heading; `wording.md` (two load-error sentences added, marked); no version change |
| A6 | met | `docs/methods-decision-ledger.md` row `VAL-<merge>` after VAL-92, citing VAL-30 and VAL-92 |

## Not checked

- The retry loop itself with real responses (httr2's tested behaviour; T8 tests the policy only).
- The bundle's httr2 1.2.3: the tests ran on the runners' httr2 1.3.0 (each run's `library.tsv` artifact); the policy
  fields T8 reads were read in 1.2.3's source by the plan (I07). The lock's `Declared` for httr2 is now stale
  (`indirect`; httr2 is in Imports): to be set at the next lock refresh (I07), not edited here.
- What the real `gprofiler2::gost` 0.2.4 does with an empty query (F1 replaced it).
- The summary report consumer in F1.
- The app in a browser (notices checked through the recording seams only).
- The exact head text of the R2 order-dependent digest (inferred from the standalone run and R3).

## Budget

The R2 diagnosis runs (45289730439, 45289724427; 1.56 core-hours by the API) are charged to the correction row
(ceiling 230, otherwise unused); the correction pass below uses the same row.

## Correction pass (implementation audit I01, I04): stopped after its R2

- Merged `v0.9-staging` 447ba678 (true merge 2c311694; NEWS and ledger conflicts resolved mechanically). Tests-only
  commit 51d8223f: T2b (iv) with the register's text; I04's T3 `not_attempted` row, T5 (d) filter outcome and notice,
  T7 replayed `annotation_sha256`. Fix 622ea91a: isoform suffixes removed before the query; Rd, NEWS, ledger.
- **R2c** (45303243737, 622ea91a vs 51d8223f, the 29 files and check; 15.39 core-hours, correction row):
  - T2b (iv) fails on base (`accession:P10060-2` sent) and passes on head: I01 fixed. T3 not_attempted and T7 pass.
  - T5 (d) fails on both sides with the same values: `organism_filter_applied` is TRUE (expected FALSE) and the "no
    organism mapping" notice is absent. My expected values were written without reading them from a run (the
    register asked to read them from this R2); the observed values are the pin's values.
  - `test-prot-annotation-uniprot-direct-shared.R` new-failure: it asserts that the query holds `Q8N158-2`; with I01's
    fix the query holds `accession:Q8N158`. The assertion encodes the behaviour I01 removes.
  - 26 files same; check identical on both sides.
- Not run: R3. Stopped under the brief's rule (one correction pass; unsure whether these two are instrument defects
  or results, so treated as results). Proposed edits for the coordinator's decision, both test-only: T5 (d)
  `expect_identical(project$outcomes$organism_filter_applied, TRUE)` and drop the notice assertion (or assert its
  absence); in the direct test replace `expect_true(grepl("Q8N158-2", ...))` with
  `expect_true(grepl("accession:Q8N158)", request_log$query, fixed = TRUE))` and
  `expect_false(grepl("Q8N158-2", request_log$query, fixed = TRUE))`. Then R2 and R3 again.
