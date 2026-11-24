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

## Inferred (not verified)

- The R2 new-failure is order-dependent and not a behaviour change: alone the session gives the frozen value on both
  sides, and the full suite (R3) shows the file `same`. The digest covers plotly's random element ids, which the mask
  replaces only when they are ten or more hex characters; the earlier cases of the file change how many random names
  are drawn before this plot. Listed for the plan's re-record list (`test-omic-enrich-module.R`, the refusals case,
  frozen `gprofiler_plot` of pre_da/run_before_da, da, open), not edited here.
- F1's g:Profiler identifiers without annotation came from the fixture's own `gene_names` (feasibility.md).

## Acceptance criteria

| # | State | Evidence |
|---|---|---|
| A1 | met | F1 table and F3/F4 counts in `feasibility.md`, `f4-measurement.md`; 2.3's answers: the dependents grep was re-run at the merged base (29 files in R2); line numbers changed by the merged tickets are cited at the merged base in the ticket's section 7 |
| A2 | met | R1 (base values) and R2 (base side 0f64e328 against head) |
| A3 | met | T6 passes on base (R1) and head (R2, R3) |
| A4 | met | R2 and R3 above; R2's one difference explained (inferred) and listed for re-record; R3 0 differences; check identical; KIT-173 merged first: T6's literal values were run on the merged base (R1) |
| A5 | met | NEWS under the unreleased heading; `wording.md` (two load-error sentences added, marked); no version change |
| A6 | met | `docs/methods-decision-ledger.md` row `VAL-<merge>` after VAL-92, citing VAL-30 and VAL-92 |

## Not checked

- The retry loop itself with real responses (httr2's tested behaviour; T8 tests the policy only).
- Whether UniProt matches `accession:X-2` (isoform accessions in field-qualified terms); the design path strips
  isoforms before the query, the enrichment step's own retrieval does not.
- The summary report consumer in F1.
- The app in a browser (notices checked through the recording seams only).
- The exact head text of the R2 order-dependent digest (inferred from the standalone run and R3).
