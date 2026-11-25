# KIT-259: outcome (implementation, ticket r1.3)

**State: implementation audit AUDIT-KIT-259-impl-20261126 done: READY FOR A MERGE TRAIN on its conditions (records
repairs I01, I04, I05, I06, I08 applied here; train scope, section 6; ledger resolved on staging 828dd105).** R2 had
stopped at stop rule 3 on one unpredicted test-code expectation that froze the A1 defect; the coordinator ruled it a
declared edit (ticket r1.3) and the one correction pass ran it. **Part B is a release prerequisite under E-4: both (b)
and (c) left the normalisation saved but unfinished at the default warn (stop rule 2 did not apply).** Open before the
release commit: I02 (code, not applied) and I03 (the owner's question), section 7.

Branch `kit-259` from v0.9-staging c4505c7f. All computation on GitHub-hosted runners (run-tests); no R on the
laptop. Runs and minutes: `runs.tsv`. Returned tables: `g0b-45623233776/`, `r1-45626688355/`, `r2-45629589384/`.
The triage that preceded the ticket (KIT-259 G0, run 45610170354) is recorded on branch `kit-259-g0`.

| Commit | What |
|---|---|
| e96ef0f6 | G0-B script `g0b.R` (predictions written in it before the run); no package code |
| f5001330 | G0-B records |
| e0cb5a99 | R1 tests only: `tests/testthat/test-kit259-pearson-histogram-figure-failure.R`; declared edit `test-metab-01k-s4-plot-pearson-characterization.R:215` |
| 040fc8a2 | R1 records |
| 1d6737c1 | the fix (R/ only) |
| 55278377 | R2 records, NEWS, ledger row, `wording.md`, this file (stopped at stop rule 3) |
| b2703094 | correction, tests only: the declared edit of `test-workflow-saved-exports.R:730-731` |
| 1fc916b2 | merge of origin/v0.9-staging 42846c0e (ledger append collision: both rows kept, staging's first) |
| f3f58167 | correction-pass predictions (`correction-predictions.md`), committed before the run |
| `run-k259-c-base` 545c4557 | the correction pass's base: e0cb5a99 + the edit (3a0f3c04) + staging 42846c0e |
| d09d6fed | R2c records |
| (this commit) | the implementation audit's records repairs (I01, I04 to I06, I08, I09, I14 routing, I02/I03 open, train scope); then the merge of staging 828dd105 |

## 1. Observed

### G0-B (run 45623233776, at e96ef0f6 = c4505c7f plus the script; 16.5 job-minutes)

55 of 57 predictions held (`g0b-base.tsv`: 57 rows, 55 TRUE, 2 FALSE) (`g0b-45623233776/g0b-base.tsv`; every value in `g0b-observed.tsv`).

- (a) 'Color by' and 'Shape by' offer `group, factor1, factor2, batch` on the label-free lane's builder design; not
  `sample_id` (the builder names it `Run`). Design columns: `Run, group, factor1, factor2, factor3, batch, replicates,
  tech_reps, tech_rep_group`. Lane `warn` = 1.
- (b) In process, label-free lane, 'Color by' = `batch` (unset: all NA): Run Normalization committed `normalized`
  (head `normalized`), `ruv_normalized_obj` set, `normalization_complete` FALSE, publication pending
  (`normalization:FALSE`); notice "Error in normalization: 'from' must be a finite number"; the failing figure was the
  post-normalisation correlation (PCA, RLE, density drawn). Then Apply correlation and Skip: "Error in correlation
  filtering: normalization outputs are pending; retry that action"; DA entry: "Error in DA analysis: normalization
  outputs are pending; retry that action"; Reset: head `protein_replicate_filtered`, "Normalisation reset: ...";
  Run again (same 'Color by'): committed again, the same error, unfinished; 'Color by' = `group` then Run again:
  completed with the completion notice. The supplementary in-process run on the DIA owner (unset batch) gave the
  same, and Run again from the stranded state (before Reset) failed the same way.
- (b) Managed (KIT-122's in-test harness), after the lifecycle check (a managed run completed with its notice on the
  worker; a worker failure reached the tab's failed handler with "[error 10] Error in normalization: "): the action
  committed `normalized` on acceptance; the figures action failed ("'from' must be a finite number"); `publish` never
  ran (nothing pending, `ruv_normalized_obj` NULL); Apply and Skip: "normalization must complete before correlation
  filtering"; Run again failed the same way; Reset recovered; 'Color by' = `group` then Run completed.
- (c) 16 groups x 2 samples, default 'Shape by' (`group`): the post-normalisation PCA fails at `ggsave` with
  "Insufficient values in manual scale. 16 needed but only 15 provided." (generator on the S4 object, and the tab in
  process). In the tab: `redraw_pre` FALSE on the first run; head `normalized`, unfinished, publication pending,
  notice "Error in normalization: Insufficient values in manual scale. ...". The failed `ggsave` left a 5,187-byte
  `post_norm_pca.png` behind (generator on the S4 object).
- Not held (2), as stated here: `ip.a.zero_pair_variable` predicted `batch`, observed `factor2,batch`. That value is the gate's own count
  (`g0b.R:248-252`), which drops empty text. The builder leaves an unassigned `factor2` as empty text, not NA
  (observed). The package's pair function treats empty text as a value (`R/func_prot_s4_qc_correlation_methods.R:216`),
  so by reading 'Color by' = `factor2` would put all ten samples in one group (45 pairs across KO and WT), not zero
  pairs; not observed. `batch` (NA) was used, as predicted; `c.c_entry.events` predicted the pre-normalisation
  PCA failing at tab entry, observed nothing drawn at tab entry in the harness (the `selected_tab` reactive given to
  `testServer` did not trigger the tab-entry observer): the tab-entry probe in (c) not observed; not re-sent, no
  decision depends on it.

**Stop rule 2 (did not apply):** both (b) and (c) leave the normalisation committed and unfinished at the default `warn` (0 in the
harness, 1 in the lane; neither promotes a warning). Part B stays; E-4 (as amended, "including (b) alone") makes it a
release prerequisite (both (b) and (c) left it saved but unfinished).

### R1 (run 45626688355, at e0cb5a99 = c4505c7f plus the tests; 13.4 job-minutes)

`r1-45626688355/`. New file: 39 passed, 32 failed; `test-metab-01k`: 16 passed, 1 failed (the declared edit). Every
E row failed by assertion at the base value the ticket states; every G row passed.

| Row | Base value observed |
|---|---|
| A-E1, A-E2 | 0 of 20 pairs drawn (ymax NA); `ggplotGrob` warns "Removed 1 row ..." |
| A-E3, A-E4, A-E5 | error (`simpleError`), native and rebuilt |
| A3-E1/E2 (lipid, metab) | 1 of 2 pairs counted; x lower limit 0 > -1 |
| A3-E4 (lipid, metab) | no subtitle |
| A3-E3 | the declared pattern does not match the base source |
| B-E1 | error; three files, no placeholder, no warning |
| B-E2, B-E6 | not complete, "Error in normalization", Apply refused "pending" |
| B-E3 | not complete, the chain's error |
| B-E4 | the earlier run's image kept (md5 unchanged) |
| B-E5 | not complete, no PCA file in the project |
| A-G1, A3-G1, B-G1, B-G2 | pass |

### R2 (run 45629589384)

Base e0cb5a99 (tests only), head 1d6737c1 (fix); 43 files, 4 shards per side, R CMD check both sides; 358.3
job-minutes. `r2-45629589384/`.

- 40 files `same` (identical failure lists on both sides where base already fails: `test-prot-norm-core.R` 306,
  `test-gui-entry-proteomics.R` 4, `test-workflow-omics-binding-steps.R` 1, `test-workflow-operation-registry.R` 2,
  `test-workflow-gui-diagnostic.R` 3 errors (lobstr, shinytest2 absent)); among them `test-prot-norm-core.R`,
  `test-prot-norm-worker-figures.R`, `test-workflow-saved-qc.R`, the three GUI-entry lanes, as predicted.
- `fixed`, as predicted: `test-kit259-pearson-histogram-figure-failure.R` (39/32 to 85/0: every E row passes on the
  fix) and `test-metab-01k-s4-plot-pearson-characterization.R` (16/1 to 17/0, the declared edit).
- **Not predicted: `test-workflow-saved-exports.R` `new-failure`** (base 271 passed, 0 failed, 1 error; head 261
  passed, 10 failed, the same error). All 10 are one expectation in the block "ten small normalization revert style
  export cycles release action resources", once per cycle: `expect_warning(ggplot2::ggsave(..., styled, ...),
  "Removed 1 row")` (`tests/testthat/test-workflow-saved-exports.R:730-731`). The test draws the correlation product
  (`buildProtNormQcPlotProduct(..., "correlation")`) of the "small" fixture (two groups of six, 30 pairs) and expected
  ggplot2's removed-rows warning, i.e. the censored bar of A1. On the fix the bar is drawn and no warning is raised:
  the test froze the defect A1 repairs. The ticket's section 6 did not list it.
- R CMD check: base and head 0 errors, 3 warnings, 3 notes; findings identical (6 on both, none only on either side).

**Stop rule 3 (an unpredicted R2 verdict): stopped.** R2 found that `test-workflow-saved-exports.R:730-731` required the
A1 defect. It is test code, not a stored reference. The coordinator's ruling (2026-11-26): a declared edit (ticket
r1.3, section 6) asserting the opposite, that saving the correlation product raises no removed-rows warning
(`expect_no_warning(..., message = "Removed [0-9]+ rows?")`; the rest of the block unchanged), then the one
correction pass (R2c below).

### R2c, the correction pass (run 45641244855)

Base `run-k259-c-base` 545c4557 (tests-only e0cb5a99 + the declared edit + staging 42846c0e; R/ without the fix),
head f3f58167 (the fix + the edit + the same staging merge + the predictions file); three files, no check (R2's check
stands; the merge train repeats it); 28.4 job-minutes. `r2c-45641244855/`. Predictions: `correction-predictions.md`,
committed before the run. Every prediction held:

| File | Base | Head | Verdict |
|---|---|---|---|
| `test-workflow-saved-exports.R` | 261 passed, 10 failed (the edited expectation, once per cycle: the warning is still raised), 1 error (lobstr absent) | 271 passed, 0 failed, the same 1 error | `fixed` |
| `test-kit259-pearson-histogram-figure-failure.R` | 39 passed, 32 failed | 85 passed, 0 failed | `fixed` |
| `test-metab-01k-s4-plot-pearson-characterization.R` | 16 passed, 1 failed | 17 passed, 0 failed | `fixed` |

## 2. Acceptance criteria

| # | Criterion | Evidence | Met |
|---|---|---|---|
| A1 | G0-B on base | `g0b-45623233776/g0b-base.tsv`; section 1 | met (2 rows not held, stated) |
| A2 | E rows fail at base by assertion, G rows pass; all pass on the fix | R1; R2 head | met: all 85 expectations of the new file pass on the fix (R2 head) |
| A3 | R2 every file `same` except the new file and the declared `test-metab-01k` | R2 `compare.tsv` | met after the correction pass: R2 40 files `same`, 2 `fixed`, check identical; the one unpredicted file is the declared edit (r1.3) and is `fixed` in R2c |
| A4 | NEWS, ledger, wording; this file with "For the change report" and P16's note | `NEWS.md`, `docs/methods-decision-ledger.md` (VAL-<merge>), `wording.md`, section 4 | met (drafts for the owner's pass) |

## 3. Inferred, and notes

- P16: under `load_all` the render identity's software digest falls back to the version plus
  `body(newDIAStateAdapter)` (`R/utils_workflow_dia_state.R:1747-1756`), which this fix does not change, so a
  development store reused across the fix could serve a pre-fix image; an installed build keys on the binary's
  sha256. Not a user route.
- P17 (for KIT-247's record): after A1 the label-free fixture raises no removed-rows warning, so KIT-247's browser
  lane passes that point at `warn = 2`; any other warning still fails it there (F4). KIT-247 must not count this
  instance as its own fix.
- A2's wording for every-pair-NA (A-E5): the ticket's single A2 text ("No two samples share a value ...") would be
  false when pairs exist but none has a correlation, so that case has its own sentence (`wording.md` W5).
- B2 also covers the DIA-with-limpa pass-through's notice (the pass-through draws through the same generators; its
  completion is wrapped by the same one-notice handler), and the worker's canonical-correlation figure.
- Convention divergence (recorded, ticket section 3): `drawProtNormQcFigure` (`R/utils_workflow_dia_norm.R:385-397`)
  catches any error raised while a QC figure's display data are read or computed (its QC product: state reads, PCA,
  RLE, density or within-group correlations, the product's save and root) or while the figure is drawn and saved. It
  logs the error with its message and writes the placeholder. Re-raised: `workflow_figure_unwritten` and every
  `saved_entity_*` class (including `saved_entity_limit` and `saved_entity_cancelled`). Interrupts are not errors and
  pass. Also caught: R's own allocation error, and the head-ownership conditions `workflow_execution_stale` and
  `workflow_execution_busy` that a product's root can raise in process when another controller holds the project.
  The normalisation's data are never inside: they are computed and committed before the figures, and its tables are
  written outside them. With `qcDir` NULL a failing figure is logged only: no placeholder, no warning (I05).
- The notice is shown after Run Normalization & RUV (in process, managed and the DIA-with-limpa pass-through). Without
  a notice the placeholder is still written: at tab entry and after an aesthetics change (by design), at the DIA-with-
  limpa RUV covariates step (`R/mod_prot_norm_workflow_helpers.R:104-109`, `:290-299`), and on a managed session for a
  pre-normalisation figure redrawn because the input changed (`:628`). There the warning reaches the console only, and
  at warn 2 or above it is an error. The covariates step's canonical-correlation figure
  (`R/utils_workflow_dia_products.R:1784-1786`) is not guarded. (I04; these are limits of the fix, not covered by
  NEWS's "a warning names the figure".)
- I10 (note): the worker logs the reason through logger's default glue formatter (the application sets
  `formatter_paste`, `R/utils_shiny_logging.R:68`; the worker does not), so braces in a reason would be interpolated
  there; a parse failure would fail the worker's figure action as before B2. Not repaired (optional one-line follow-up:
  `logger::skip_formatter`).
- I11 (note): no test shows the partial PNG of a failed `ggsave` replaced by the placeholder's bytes (B-E5 checks the
  header size, which the partial file shares); sound by reading (`writeProtNormFigurePlaceholder` writes to the same
  path). Optional test-only follow-up: compare B-E5's PCA md5 with the writer's output, as B-E1 does.
- I12 (note): the declared edit's pattern `Removed [0-9]+ rows?` also matches stat_bin's non-finite removals: correct
  for its complete-data fixture, and stronger than A1's claim (proteomics NA pairs are still dropped with a console
  warning).

## 4. For the change report

**Legacy (main part).**

- **A1, proteomics within-group Pearson histogram.** What was wrong: the count axis was fixed at 0 to 4. Any 0.001-wide
  bar holding more than 4 sample pairs was removed when the figure was drawn, with a console warning only. Who: every
  proteomics user, every route (DIA iq, DIA with limpa, label-free, TMT), in process and managed, whose within-group
  correlations put 5 or more pairs within 0.001 of each other (likely with tight replicates or many pairs; ticket
  section 5 gives the approximation). It reached the tab, the composite QC figure, its copy in the results package,
  the report and the session archive. Since 0efbb3bf (2024-10-17), which moved the limit from the QC workbook into the package's
  `plotPearson`; in every release tag made after it (v0.80, v0.9, v0.2, v0.3.6.3, v0.5.0-rc1, v0.9.0; the older tags
  v0.73 and v0.75 do not have it). Code: base
  `R/utils_workflow_dia_products.R:84`, `R/func_prot_s4_qc_correlation_methods.R:68`; fix
  `R/func_prot_s4_qc_correlation_methods.R` `protPearsonCountScale` (used at `:73` and
  `R/utils_workflow_dia_products.R:88`). Change: the axis runs to the tallest bar (at least 4) with whole-number
  ticks; panels whose bars hold at most 4 are byte-identical to before (A-G1). Evidence: run 45610170354 (5 of 20
  pairs missing on the label-free fixture, before and after normalisation); R1 (0 of 20 drawn at base), R2 (20 of 20
  at the fix). Unchanged: the bin width 0.001 and the x range from the lowest correlation minus 0.5 (owner item O-1).
- **A3, lipidomics and metabolomics histogram.** What was wrong: the x axis and the bins covered [0, 1], so a negative
  within-group correlation (a failed or near-empty sample, a processing error) was never drawn; pairs with no
  correlation (too few shared values, or a sample with no variation; after KIT-242 every pair of an empty sample)
  were removed without a note. Since c9b950b8 (metabolomics, 2025-05-04) and ce640728 (lipidomics, 2026-01-18). Code:
  base `R/func_lipid_s4_correlation_methods.R:276`, `:285`; `R/func_metab_s4_qc_methods.R:757`, `:766`; fix lipid
  `:278`, `:289`, `:300`; metab `:759`, `:770`, `:781`. Change: the axis and bins extend below 0 to the lowest
  correlation (rounded down to 0.1); a subtitle counts the pairs with no correlation. Correlations at or above 0 with
  none missing draw exactly as before (A3-G1). Bin width 0.01 unchanged.
- **A2, no within-group pair.** A 'Color by' that leaves no two samples in one group (an unset batch/factor, one
  sample per value, only Pool/QC groups) gave an error; now a labelled empty figure says why. Fix
  `buildProtPearsonUnavailablePlot` (`R/func_prot_s4_qc_correlation_methods.R`), called at `:67` and
  `R/utils_workflow_dia_products.R:82`.
- The correlation filters were never affected on any omic: they read the correlation tables, not the figures. No
  number, filter decision or record changes.

**New to 0.9 (appendix).**

- **B, a normalisation figure that cannot be drawn.** The 0.9 reference states commit the normalisation before its
  figures (KIT-111 D3). A figure that failed then left the normalisation saved but unfinished: in process every later
  step (Apply, Skip, DA, session export, Copy) was refused "normalization outputs are pending; retry that action"
  until Reset; managed, Apply/Skip were refused "normalization must complete before correlation filtering"; Run again
  failed the same way; the project kept the earlier run's figure under the same name. Triggers a user reaches at the
  default `warn`: an unset 'Color by' variable (zero pairs; removed by A2) and 16 or more groups with the default
  'Shape by' (the PCA's 15 point shapes; not removed). On main a post-normalisation figure failure stopped the step
  before saving (`main:R/mod_prot_norm.R:1373-1380`) and an RUV figure failure was ignored (`:1795-1799`). Fix: each
  figure is drawn through `drawProtNormQcFigure` (`R/utils_workflow_dia_norm.R`), which replaces a figure that cannot
  be drawn by a placeholder image under its own name and size (`writeProtNormFigurePlaceholder`); the generators
  (`R/mod_prot_norm_qc_generation_helpers.R`), the canonical-correlation figure
  (`R/utils_workflow_dia_products.R:1084`), the worker's stage loop (`:2074`, skipped and logged) and composite loop
  (`:2157`, `:2166`); one warning notice per run (`withProtNormFigureNotice`, `R/mod_prot_norm_workflow_helpers.R:496`,
  `:527`, `:610`). Limits: no notice at the DIA-with-limpa RUV covariates step or for a managed pre-figure redraw, and
  the covariates step's canonical-correlation figure is unguarded (section 3, I04). Report text: the report's Methods
  and Results do not mention a replaced panel (I03, the owner's question). Evidence: G0-B run 45623233776 (base); R1 run 45626688355; R2 run 45629589384.

## 5. Not checked

- The pre-normalisation PCA at tab entry with 16+ groups (G0-B (c): the tab-entry probe in (c) not observed; see
  section 1), and the `redraw_pre` path on the fixed code.
- The RUV covariates route with a failing figure; multi-assay subtitles (A3); byte identity of lipid/metab images for
  non-negative data (A3-G1 checks layer data and scale); I02's labels (by reading).
- The browser (shinytest2) lanes, and any launch at `warn = 2` (F4).
- Whether the worker's and the application's placeholder images are byte-identical (affects only the composite's
  reuse; not tested).
- TMT and label-free managed routes beyond the DIA owner; the DIA-with-limpa pass-through's notice (by reading).
- Windows; real data; the frequency figures of ticket section 5 (an approximation).
- The wording: drafts for the owner (`wording.md`).

## 6. The merge train: scope and predictions (AUDIT-KIT-259-impl-20261126 section 5, I07)

Binding: base staging 828dd105; head = 828dd105 + `kit-259` (this branch after its merge of 828dd105).

Scope, with R CMD check:
- the ticket's two greps at 828dd105 plus its five named files: 44 files (R2's 43 plus
  `test-kit243-summary-current-outputs.R`);
- KIT-244's new file `test-kit244-limpa-pass-through.R` (B2 wraps the pass-through KIT-244 changed);
- I07's dependents, which R2 did not run: `test-prot-design-log2.R`, `test-workflow-gui-memory-lifetimes.R`,
  `test-module-ci-prot-dia-limpa.R`, `test-workflow-analysis-record.R`, `test-workflow-integration-input.R`,
  `test-workflow-saved-da.R`, `test-prot-05c-norm-module-seams.R`, `test-workflow-server-shared.R`.
- Of R2's 43, `test-lipid-12-norm-module-contracts.R` and `test-metab-03b-norm-module-helper-characterization.R` are
  quarantined (0 expectations both sides) and `test-workflow-memory-mechanism-member.R` is all skipped (13).

Predictions (by reading):

| File | Base | Head | Verdict |
|---|---|---|---|
| `test-kit259-pearson-histogram-figure-failure.R` | absent | 85 passed, 0 failed, 0 skipped (R2/R2c head) | new, passes |
| `test-metab-01k-s4-plot-pearson-characterization.R` | 17/0 (old pattern, old code) | 17/0 (new pattern, new code) | same |
| `test-workflow-saved-exports.R` | 271 passed, 0 failed, 1 error (lobstr) | the same | same |
| `test-prot-norm-core.R` | staging's list | identical list; the body-digest test passes | same |
| every other file | | | same |
| R CMD check | | | identical |

Any other verdict is a result: stop rule 3.

Partners merging later (ticket 8.1): KIT-257 (`test-kit257-norm-own-files.R` with this file; one code conflict hunk in
`R/mod_prot_norm_workflow_helpers.R:524-527`, resolved by keeping both: `finish, fail, preQc = ...) else
withProtNormFigureNotice(finish(runDIAReferenceNormalization(...)), showNotificationFn)`); KIT-242 (after its OD-1 an
empty sample's pairs are NA and gain A3's subtitle); KIT-247 (P17).

## 7. Open after the audit (not applied)

- **I02, an open inferred defect (code, not applied).** A3's x-axis labels may print in scientific form, e.g.
  5.551115e-17 for 0, when the lowest correlation is in [-0.3, -0.2), [-0.6, -0.5) or [-0.7, -0.6): `seq(lower, 1, 0.1)`
  gives 5.55e-17 (lower -0.3) or 1.11e-16 (-0.6, -0.7) for the tick meant as 0, and `format()` then turns every label
  scientific (INFERRED; no R). Every pair is still drawn and counted; no number changes. The register's fix, quoted:
  "Replace `breaks = seq(lower, 1, 0.1)` with `breaks = if (lower < 0) round(seq(lower, 1, 0.1), 1) else seq(lower, 1,
  0.1)` (non-negative data keep today's expression, hence A3-G1's bytes)." Its test row, quoted: "New row A3-E5 (lipid
  and metab): an assay whose lowest within-group correlation lies in [-0.3, -0.2); `ggplot2::get_guide_data(p,
  "x")$.label` contains "0" and none matches `"e[-+]"`. Fails on d09d6fed, passes on the fix." Not applied: the
  ticket's one correction pass is used and an extension is the owner's to allow (parked row 45). Before the release
  commit (A3 is a release prerequisite, E-3).
- **I03, the owner's open question.** When B2 replaces a figure, the report's Methods ("the samples were checked ...
  using PCA plots ...", `inst/reports/proteomics/LFQ_report.rmd:636-641`; DIANN `:685`, DIANN_limpa `:581`, TMT `:627`)
  and Results (each panel described, `LFQ_report.rmd:660-668`) stay as written; only the panel says "<Figure>: could not
  be drawn. <reason>". (a) Accept: the labelled panel suffices (records-only: NEWS item 4 gains "The report's text is
  unchanged and does not list the replaced panels; the panel itself says why."). (b) Recommended before release: the
  report names the panels that could not be drawn (code; a later bounded ticket).
- **Routed findings (I06, I14; the owner informed; plan section 13 candidates):**
  - F9: "Proteomics: a 'Color by' column the design builder leaves empty (empty text, e.g. `factor2`) puts every sample
    in one group, so the within-group Pearson histogram draws every pair, across groups, as within-group (by reading,
    `R/func_prot_s4_qc_correlation_methods.R:216`; not observed). Display only."
  - F10 (I14 a): the proteomics histogram groups pairs by 'Color by' (`R/utils_workflow_dia_products.R:302`) while the
    correlation filter groups by 'RUV Grouping Variable' (`R/mod_prot_norm_observer_helpers.R:198`); the same by
    default; after a 'Color by' change the figure used to choose the threshold shows other pairs than the filter tests.
  - F11 (I14 b): proteomics NA pairs are still removed from the histogram with a console warning only, while
    lipid/metab now count them in a subtitle (A3); omics parity.
- Optional follow-ups (coordinator): I04 (wrap the covariates' finish/publish and the managed `preQc` delivery with
  `withProtNormFigureNotice`; guard the covariates' cancor `ggsave`); I05 (re-raise `workflow_execution_stale` and
  `workflow_execution_busy`); I10 (`logger::skip_formatter`); I11 (test-only).
