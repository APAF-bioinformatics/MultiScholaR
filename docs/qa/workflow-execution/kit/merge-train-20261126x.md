# Merge train 31 of 2026-11-26: KIT-272 alone (report text for panels d to f of Figure 1)

**This train is NOT stacked.** Branched from `origin/v0.9-staging` at `ec27f69e0bbe970ab80fee9865b6b6d86269cb19`
(fetched before building; it holds trains 25 and 27 and the ledger commit VAL-359). It does not contain the other
trains still out (28 to 30 and any later). KIT-272 is the only member.

## Train

Base: `ec27f69e0bbe970ab80fee9865b6b6d86269cb19` (staging).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-272 (the report templates' sentences about the quality-control figure) | `origin/kit-272` 100d7e58 | f211b6675 | `VAL-<train31-1>` (the coordinator numbers it at landing) | AUDIT-KIT-272-impl-20261126 (READY FOR A MERGE TRAIN; I01, I04, I05, I06, I08, I09 applied below) |

Head (bound by the runs): `f211b6675` (`git merge --no-ff origin/kit-272`). The docs commit after it (the register's
record fixes, NEWS wording, this file's predecessors in `kit-272/`) changes no code, template or test; the render
script copy and this file follow in the declaration commit. Neither of those is under test.

### Conflicts met and their resolutions

- **`docs/methods-decision-ledger.md`** only: every staging row kept (through VAL-359); the branch's row follows
  VAL-359 as `VAL-<train31-1>`, and I numbered nothing (the coordinator numbers it at landing). The record-fix commit
  appends the run ids and the two re-sends to that row (I08).
- NEWS, `lipidomics_report.rmd`, `metabolomics_report.rmd` and `test-kit236-report-sentences.R` merged by themselves
  with KIT-255's edits. No conflict under `R/`, `tests/`, `inst/` or `man/`.
- `test-kit236-report-sentences.R` on the merged tree holds both KIT-255's block (the `test_that` "KIT-255 (I01, I02): the
  comment counts each sample set missing once ..." at `:736`) and KIT-272's declared edit (`:494-524`). The merged constants
  block, verbatim:

```r
k272OldBullet <- "• Density boxplots to assess distribution consistency"
k272NewBullet <- c(
    lipid = paste0("• Box plots of the samples' scores on the first two principal components of every lipid with ",
        "at least two values, a separate analysis from the PCA panels, which use the tenth of the lipids with the ",
        "largest ratio of standard deviation to mean"),
    metab = paste0("• Box plots of the samples' scores on the first two principal components of every metabolite ",
        "with at least two values, a separate analysis from the PCA panels, which use the tenth of the metabolites with ",
        "the largest ratio of standard deviation to mean"))
```

### Scope of the change

`git diff --name-only origin/v0.9-staging f211b6675` is exactly KIT-272's change and nothing else: the six templates
(`inst/reports/proteomics/{DIANN,DIANN_limpa,LFQ,TMT}_report.rmd`, `inst/reports/lipidomics/lipidomics_report.rmd`,
`inst/reports/metabolomics/metabolomics_report.rmd`), the new `tests/testthat/test-kit272-prot-box-panel-text.R`, the three
declared-edit files (`test-kit227-report-ruv-sentences.R`, `test-kit235-lipid-metab-report-parser.R`,
`test-kit236-report-sentences.R`), `NEWS.md`, `docs/methods-decision-ledger.md` and `docs/qa/workflow-execution/kit/kit-272/`
(records). Nothing under `R/`, `man/`, `DESCRIPTION`, `NAMESPACE` or `tests/testdata`.

### Path lengths

The longest path new since the base is under 100 bytes with the `MultiScholaR/` prefix (the branch's records, e.g.
`docs/qa/workflow-execution/kit/kit-272/r1b-45728034972/compare-failures.txt`-length names, and this train's own files
are shorter).

## Declaration (written before any run is sent)

**Purpose.** (1) The merge check on the merged tree: KIT-272's 19-file regression scope with R CMD check, base against
head. (2) A render job: the first knit of the edited proteomics templates (the branch's tests call `FormatResults`,
`FormatMethodDetails*` and the placeholder chunk; no knit has run). Two requests; one run each.

**Scope recount (comparison).** The ticket's grep (`git grep -lE 'system\.file\("reports"|"inst", "reports"|inst/reports|
reports", (spec|dir)|reports", "(proteomics|lipidomics|metabolomics|transcriptomics|shared|integration)"' -- tests/testthat
':!tests/testthat/helper-*' ':!data/benchmark'`, run on the merged tree) gives **19 files**, the same 19 as the branch's R2:
test-enrichment-significance, test-go-propagation, test-kit204-prot-da-not-tested, test-kit205-enrichment-background,
test-kit215-not-tested-lists, test-kit220-da-threshold-consumers, test-kit227-report-ruv-sentences,
test-kit234-workbooks-withdrawn, test-kit235-lipid-metab-report-parser, test-kit236-report-sentences,
test-kit245-github-push-safety, test-kit259b-second-correction, test-kit264-trx-export-current,
test-kit272-prot-box-panel-text, test-omic-gui-config, test-omic-summary-module, test-omics-summary-core,
test-prot-summary-report-templates, test-workflow-analysis-record. Between the branch's base 7f17547f9 and this base
the only one of the 19 whose test file changed is test-kit236 (KIT-255 appended a `test_that`: 289 to 293, observed in
train 25's run).

**Run settings, comparison, copied from the branch's R2 request (`kit272-r2-085e57e5.json`):** `env` empty; no extra
packages; `blas: reference`; `shards: 1`; `file_timeout_minutes: 45`; `file_timeouts` empty; `check: true`; no scripts;
`tools_ref` empty; the 19 files in that request's order. Base `ec27f69e0bbe970ab80fee9865b6b6d86269cb19`, head
`f211b6675` (full hash in the request).

**Render request.** One scripts job in the manner of train 23's (`train23-render-460253ba.json`): no tests, no check,
`blas: reference`, `NOT_CRAN=true`, `out` `kit259-render`, `out_pattern` `\.(txt|tsv)$`, script timeout 55 minutes.
Binding: the scripts job reads its script from the `package_commit` tree, and the script copy is in the declaration commit, so the render binds `ed941465d70537595a0da606d285a08d41e622bc` (the declaration commit: f211b6675's code and templates plus docs and NEWS only); the comparison binds `f211b6675` as above.
Script: `docs/qa/workflow-execution/kit/kit-272/render-train31.R`, a copy of `kit-259/render3.R` that differs in one
line: `BASE <- "ec27f69e0"` (was `82a5df02d`), which only changes the commit whose templates the V0 base side renders.
What it covers: the four proteomics lanes (`prot_lfq`, `prot_tmt`, `prot_dia`, `prot_dia_limpa`) run through
`guiEntryRunLane`, each knitted with the head's template and the base's template (V0), with the KIT-259 placeholder
list (V1) and with a list without a combined figure (V2). **Lipidomics and metabolomics are not covered**: the script
has no lane for them and rendering them needs new code, which I was told not to write. Their templates are **not shown
by a knit** (T6, T11 and the declared edits of KIT-235/236 pin their text, as before). The panel-d placeholder label
(`kinds`) is not reached either, because V1 injects the panels a and h only (T7 pins it).

**What the render must show** (declared before the run). Read from the script's own `checks.tsv` (about 27 rows) plus an
offline text check on the returned V0 texts (`<lane>-V0-head.txt`, `<lane>-V0-base.txt`: plain text, compared with
Python on the fetched files; no R):

1. Every template knits without error on the head: 4 lanes captured, V0, V1 and V2 each rendered (no `render error` row
   and no `captured FALSE` row).
2. The head's rendered text carries the new sentences exactly once each, where section 2.2 puts them: in each lane, the
   Results bullet of the variant the lane prints (predicted: the two-stage variant for LFQ, TMT and DIA, whose lanes skip
   RUV; the three-stage variant for DIA with limpa; the variant that the base's V0 text prints the old bullet for),
   the PCA bullet of that variant once, and the Methods phrase "box plots of the principal component scores" once
   (limpa: the paragraph its route prints); and none of the old ones ("Density plots", "distribution of protein
   intensities per sample", "density boxplots") in the head's text. The base's V0 text carries the old ones exactly once
   each and none of the new.
3. The placeholder sentence ("could not be drawn by the application") appears only in the placeholder case: V0 head none
   (`checks.tsv` row "no placeholder sentence" FALSE), V1 one on LFQ, TMT, DIA (letters a and h) and none on DIA with
   limpa, V2 none.
4. The text is otherwise equal to the base's: the head's V0 text with the new sentences put back to the old ones equals
   the base's V0 text (checked offline, whitespace-normalised). **In the script's own table the row "V0 text equal to
   staging's template" is predicted DIFFERS on all four lanes, by design (the text changes)**; every other row as
   predicted: V0 drawings 3, 3, 3, 1; V0 no placeholder; V1 sentence counts 1, 1, 1, 0; V1 the predicted sentence and
   "text otherwise equal to V0" TRUE; V1 limpa "text equal to V0" TRUE; V2 no sentence.

Any other result is recorded as it is. A missed prediction is a result; no re-run.

**Ceilings.** Comparison 340 job-minutes = 22.67 core-hours (R2's: 2 x 80 + 2 x 90; core-hours = job-minutes x 4 / 60),
expected about 112 (R2's 112.25: one shard per side about 31 and 35 test-minutes, check on both sides). Render 60
job-minutes = 4.0 core-hours, expected about 24 (train 23's 23.7). Total ceiling 400 job-minutes = 26.67 core-hours,
expected about 136.

### Predictions, comparison

Base counts verified against the newest committed tables: the branch's R2 `tests-base.tsv` (base = staging 7f17547f9;
none of these files' test files changed since, except test-kit236) and train 25's `merge-train-20261126r-compare.tsv`
for test-kit236 (293 on its head, which is on this base). Where the train's own base run differs, the train's base run is
the reference.

| File | Base p/f/e/s | Head | Verdict |
|---|---|---|---|
| test-kit272-prot-box-panel-text.R | absent | 373/0/0/0 | only-on-head |
| test-kit227-report-ruv-sentences.R | 162/0/0/0 | 186/0/0/0 | count-changed (G1 42 to 66) |
| test-kit235-lipid-metab-report-parser.R | 54/0/0/0 | 58/0/0/0 | count-changed (T6 20 to 24) |
| test-kit236-report-sentences.R | 293/0/0/0 (with KIT-255's test) | 297/0/0/0 | count-changed (G1 24 to 28) |
| test-kit220-da-threshold-consumers.R | 71/2/0/0 (standing) | 71/2/0/0, identical failures | same |
| test-omic-summary-module.R | 1675/28/0/0 (standing) | as base, identical failures | same |
| test-omics-summary-core.R | 676/25/0/2 (standing) | as base, identical failures | same |
| the other 12: test-enrichment-significance (181), test-go-propagation (321/0/0/1), test-kit204 (81), test-kit205 (228), test-kit215 (112), test-kit234 (29), test-kit245 (87), test-kit259b (90), test-kit264 (85), test-omic-gui-config (1911), test-prot-summary-report-templates (123), test-workflow-analysis-record (351) | as listed | as base | same |
| R CMD check | 0 ERRORs, 3 WARNINGs, 3 NOTEs | identical findings, nothing only on head | identical |

Totals expected: 19 files; 15 `same`, 3 `count-changed` (the declared edits' match assertions: +24, +4, +4), 1
`only-on-head`; no `new-failure`, no `message-changed`, no `fixed`.

### Decision and attempts

- Every prediction holds: TRAIN GREEN; landing needs the coordinator's ledger number and status cells.
- Any other verdict is a result: recorded as it is, not re-sent, handed back.
- One comparison run and one render run. A re-send only for an instrument defect with no result judged (VAL-54), and
  only on the coordinator's word. An HTTP 500/502 at dispatch with no run created: started once by hand with the
  README's command, recorded.

## What the coordinator resolves at landing

- The ledger number of KIT-272's row (`VAL-<train31-1>`), and the status cells: VAL-316 "(7) implemented by KIT-272
  (VAL-<n>), merged by train <t>"; VAL-332 "KIT-272 merged by train <t>" (I08).
- Placeholders: `kit-272/outcome.md:65` ("ledger row `VAL-<merge>`", a record's mention of the row).
- Stale lines, listed not edited: `kit-236/wording.md:36` (P18) and `kit-259/wording.md` W8 (add "replaced by KIT-272
  (VAL-332); see kit-272/wording.md").
- NEWS: the two precision points of I04 are applied in the record-fix commit (older-project path names all three
  omics; "new to 0.9" for proteomics only); the owner's wording pass may reword.
- Wording items for the owner are in `kit-272/wording.md` (the lipid/metab caption's "density", "of the lipids" /
  "of the metabolites").
- I07: KIT-266's G3 compares `tests/testdata/kit236/base` `FormatResults` with head byte for byte and fails once this
  lands: KIT-266 needs its declared G3 edit and restated R3 reference (not this train's).
