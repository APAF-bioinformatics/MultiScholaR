# KIT-233 outcome: the lifecycle run L (discovery, test mode)

**State: L instrument proven** by the owner-allowed fourth run (VAL-234), run 45601782272 on driver `76e4b519`, commit
under test v0.9-staging `e0747fc1`. All L-0 to L-4 criteria pass (25 pass, 1 recorded: "no worker code reads test
mode"), and L-5 passes. J (PXD055981 s500) is not started: it waits for the release commit. L-u is not sent.

- **Earlier runs.** The first send and VAL-54's three instrument reruns (on `1809c51f`, `4793c9e8`, `0260ac73`,
  `628545be`, commit under test `67569959`) ended at stop rule 1.
- **The third rerun's defect.** The driver mapped the DA tab's contrast (`KO_vs_WT`) to no saved table, because the
  table is named by the design's contrast (`da_proteins_groupKO-groupWT_long_annot.tsv`). So N2 was not computed and
  L-4 (b) could not run.
- **The fix for the fourth run.** The tab's contrast is mapped to the saved table by its two group names in order. It is
  a mechanics fix: no expected value, threshold, bound or stop rule changed.
- **Which commit L binds, and why.** The fourth L binds staging `e0747fc1`, merged into `kit-233-journey` as `8664431b`.
  Merge trains 3 and 4 changed steps the journey walks: KIT-234 (S1), KIT-220 and KIT-227 (the report sentences behind
  N4), and KIT-212 and KIT-215 (S10 and S9). KIT-237 and KIT-238 run on this tree.
- **Records.** `runs.tsv`; each run's returned tables are in `l<n>-<run id>/` (`l4-45601782272/` is the fourth).
- **KIT-159 conditions (r1.2, section 10).** L-0 and KIT-159's G13 are never cited as evidence for each other. L-0's
  launch and attach results go to KIT-159 as information only. Nothing under `kit-159/` was changed.

### The fourth L (45601782272, `l4-45601782272/`)

- **Lifecycle criteria.** Every criterion of the table below passes, including L-1 "N1 to N4 computed" and L-4 (b).
- **Interface numbers against the saved files:**
  - N1: 314 proteins and 8 samples, equal.
  - N2: `KO_vs_WT` 16, equal to `da_proteins_groupKO-groupWT_long_annot.tsv` (`fdr_qvalue < 0.05`).
  - N3: the label and description, equal.
  - N4: the report's 10 up + 6 down = 16, equal.
  - The planted copy (row 1's q moved across 0.05) gives 15 against 16 shown, so a kind-4 finding is recorded as planted.
- **Step outcomes and timings.** As in the third rerun's table below, with times about 10 to 20 % shorter. S6.1 took
  86.5 s.
- **S9 refused offline**, as the addendum of 2026-11-25 (item 3) now expects.
- **KIT-251 record at S5.5** (one click): outcome `complete`; head `protein_replicate_filtered`; one "already
  running" warning; one working notice. KIT-251's defect (A2) is still present: it is not merged at `e0747fc1`.

## Observed

### Runs (`runs.tsv`)

| Run | Driver | Job min | What it showed |
|---|---|---|---|
| 45523435786 | 1809c51f | 12.1 | No step outcome. `shinytest2::AppDriver$new` skipped with "Reason: On CRAN": sudo's environment reset dropped `NOT_CRAN`. Instrument defect |
| 45525410012 | 4793c9e8 | 24.2 | Fix 1 (`NOT_CRAN=true` in the phase file). L-2 and L-3 passed. Three driver defects: (a) the import's digest signal completed S2 while the organism dialog was still open and unanswered; (b) S6.1's digest signal (`ruv_corrected` in the history) comes about 80 s before the step's own completion notice (L-5 disagreed at S6.1); (c) S8.2 was misreported as "not in the page" |
| 45573306268 | 0260ac73 | 22.6 | Fix 2: expected dialogs are answered before a step can complete; S6.1 is DOM-only; a page that does not answer is recorded as such, and the pane's HTML is kept. L-5 passed. Two driver typos: S8.2's target looked up as an element id instead of a test id (the HTML dump shows the button with `data-testid="prot-da-run"`); the declared output `normalization-norm_log` does not exist |
| 45577106977 | 628545be | 25.5 | Fix 3. Everything passes except N2: the tab's contrast `KO_vs_WT` is mapped to no file. The worker log shows the table is written as `da_proteins_groupKO-groupWT_long_annot.tsv`, under the design's contrast name. L-4 (b) could not run without that table |

Each fix changed driver mechanics only: no expected value, threshold, stop rule, bound or declared comparison. The fix
to S6.1's completion signal and the declared outputs are instrument definitions. L-5 exists to prove these; they were
changed before J, as section 14 allows.

### Lifecycle criteria at the last run (45577106977, `l3-45577106977/lifecycle.tsv`)

| Criterion | Result | Evidence |
|---|---|---|
| L-0 Chrome; shinytest2 0.5.1 | pass | `/usr/bin/google-chrome` |
| L-0 `entry.R` via `Rscript --vanilla -e` passes its trailing arguments, starts, answers | pass | HTTP 200 after about 4 s; no usage error |
| L-0 `MSR-LAUNCHER` names `process` | pass | `{"launcher":"process","reason":"KIT-233 hosted runner: no user D-Bus"}` |
| L-0 legacy RUV route (no `MSR-TEST-FEATURES`; panel absent at S6) | pass | `ruv_design_summary` absent |
| L-0 `AppDriver` attaches by URL | pass | (information for KIT-159 G13, not evidence for it) |
| L-0 digest read; controller in test mode | pass | digest keys listed |
| L-0 no worker code reads test mode | recorded | Not executable by the driver. Corrected wording (KIT-233 addenda of 2026-11-25; KIT-252 plan audit section 6): the worker inherits the application's environment variables under the plain launcher (so `MULTISCHOLAR_TEST_MODE`, if set as a variable), not the R option; no code run on the worker reads test mode. The earlier inference here (option forced off does not switch test mode off in a worker) cited lines that run in the controller, not the worker |
| L-0 network namespace; inner user = outer user (R1-02) | pass | links `lo` only; `runner` = `runner` |
| L-0 20 MB DIA-NN table uploads with no size error | pass | 22.0 MB, "Upload complete" |
| L-0 a tall tab's outputs inside a saved image | pass | Pane images are cut at the viewport height (about 1,065 px); the outputs below are saved as their own element images |
| L-1 every row has an outcome and seconds; reached S10.5 | pass | 27 rows |
| L-1 S3 complete with KIT-175's warning | pass | the persistent annotation warning |
| L-1 organism dialog answered and recorded | pass | rendered values `selected_organism=9606`, `filter_to_organism=FALSE` |
| L-1 screenshots and JSON exist, non-empty; R log; active-tab checks; progress dialogs never terminal | pass | |
| **L-1 N1 to N4 computed** | **fail** | N1 equal (314 proteins, 8 samples); **N2 not computed (driver mapping)**; N3 equal; N4 equal (report: 10 up + 6 down = 16 = tab's 16) |
| L-2 planted refusal | pass | S1 refused "Base directory does not exist: <path>" after 60 s; S2 to S10 `not_reached` |
| L-3 controller killed during S2 | pass | S2 `crashed` (exit -9), not `timeout` |
| L-4 (a) emptied DA table gives kind 3 | pass | planted finding `K233-F01` |
| **L-4 (b) q moved across 0.05 gives kind 4** | **fail** | not run: no saved table was mapped |
| L-5 digest and DOM agree at every step | pass (since rerun 2) | |

### What the application did on the 8-run fixture (last run, `l3-.../steps.tsv`)

| Step | Outcome | s | Shown |
|---|---|---|---|
| S1 project | complete | 4.9 | |
| S2 import | complete | 5.1 | "Detected format: DIA-NN, Confidence 100%"; organism dialog (one organism, 100 %) |
| S3 design | complete | 10.4 | KIT-175's persistent warning: UniProt not reachable (expected offline) |
| S4.1 to S4.7 peptide QC | complete | 13.5 to 16.8 each | rendered order = the application's table (Q-Value Filter, Protein Evidence, Precursor Rollup, Intensity Filter, Sample Quality, Replicate Filter, Imputation); no default refused on 8 runs |
| S5.1 to S5.5 protein QC | complete | 13.3 to 15.0 each | S5.5: warning "A protein replicate filter action is already running" after one click |
| S6.1 normalisation, RUV Automatic, cyclic loess | complete | 111.4 | |
| S6.2 correlation filter at 0.5 | complete | 26.4 | "Proteins remaining: 314, Samples remaining: 8" |
| S7 export | complete | 4.9 | |
| S8.1 load session, S8.2 DA (`~ 0 + group`, q 0.05, lfc 0), S8.3 table | complete | 0.8, 26.9, 1.9 | 16 significant for `KO_vs_WT` |
| S9 enrichment (g:Profiler, taxon 9606) | **refused** | 62.1 | "Error in enrichment analysis: Enrichment needs protein annotation and this project has none. Save the design again to retrieve it, then run the later steps again." |
| S10.1 to S10.5 Save, Copy, Report (docx), Download, Export state | complete | 0.2 to 6.5 | |

### Application findings and observations from L (not fixed; for the coordinator to file)

| # | Kind | What | Evidence | Known? |
|---|---|---|---|---|
| A1 | 1 (now expected) | In **test mode**, offline, enrichment is refused for want of annotation. The ticket expected the test-mode substitute to run (`R/utils_workflow_enrichment.R:967`). The design step stored an empty annotation table (0 rows), not NULL, and the gate refuses an empty table even in test mode | `l3-.../findings.tsv` F02; steps S9 | triaged 2026-11-25: KIT-175's design (an empty table is refused in every mode); S9 is now expected refused offline in L and J (addendum item 3); remedy text = wording row ENR-22 |
| A2 | 2 | One click on Protein Replicate Filter shows "A protein replicate filter action is already running". The harness pressed this button twice on purpose (`helper-e2e-browser.R:1336-1341`) | `l3-.../steps.tsv` S5.5 warnings | KIT-251 |
| A3 | 2 | A worker killed during Precursor Rollup is shown as "Error applying precursor rollup: managed action ended: worker_failed": internal words, no remedy | `l3-.../observations.tsv` (L3w), F07 | new (UI wording) |
| A4 | obs | Import is marked complete before the organism dialog is answered. A later dialog (the design save's progress dialog) replaces an unanswered organism dialog. A user who skips the dialog keeps the default organism without being told | `l1-.../steps.tsv` S2/S3 modals | new |
| A5 | obs | The normalisation tab's RUV state is committed about 80 s before the step finishes. The correlation filter could be pressed meanwhile, and it completed | `l1-.../steps.tsv` S6.1/S6.2 | new (observation only) |
| A6 | obs | `design.tsv` given as the report: detection "Unknown format, 0%". With DIA-NN forced, import refused with "Missing required DIA-NN columns: Protein.Group, Stripped.Sequence, Run" (a clear message) | L2o rows | — |
| A7 | obs | The DA tab shows contrasts as `KO_vs_WT`; the saved tables use `groupKO-groupWT` in their file names. A user looking for the tab's contrast among the files must translate | worker log | new (minor) |

## Inferred

- The N2 failure is a driver mapping defect, not an application defect. The file exists in the worker's output with
  117,636 bytes (worker log). INFERRED that it was published into `results/proteomics/da_proteins`; the published path
  was not listed.
- With a mapping by the two group tokens in order (as `k233MatchContrast` already does for enrichment), N2 and L-4 (b)
  would most likely pass. Not run: the rerun budget is spent.

## For the coordinator

1. Decided by the owner (VAL-234): one fourth L run, made as 45601782272, which proved the instrument. Discovery use so
   far: (84.4 + 21.3) job-minutes x 4 / 60 = 7.05 core-hours of the 64.0 ceiling.
2. KIT-159 (information only): `Rscript --vanilla -e '…; source("packaging/electron/entry.R")' <port> <run folder>
   <pid>` passes its arguments. The plain launcher starts. `shinytest2::AppDriver$new(<url>)` attaches by URL on
   ubuntu-22.04. It needs `NOT_CRAN=true`, or AppDriver skips with "On CRAN".
3. KIT-237 / KIT-238 needs N1 to N3: met in `journey.R` (sourcing runs nothing; `k233LaunchSpec(form, app_dir,
   network)`; `K233_ROUTE_DIA` and `match.fun` step functions). No conflict with r1.2.

## Not checked

- J (PXD055981 s500), its data phase (download, subset, stop rule 2b), the J timings and N2 on real data.
- L-u (user mode, the real `shinyFiles` chooser): not sent. It needs L proven, then one run of about 25 job-minutes. The
  driver already has the chooser click-through (`k233ChooseFile`). No other driver work is needed.
- Whether the DA long table is published into `results/proteomics/da_proteins` in a managed session.
- The pane screenshot's clipping at the viewport height (a chromote property); outputs below are captured as their own
  images.
- No screenshot was withheld by the upload scan in these four runs (no `withheld.txt`).
- The independent screenshot review (7.3), which belongs to the audit; `review-index.tsv` is the hand-off.
