# KIT-101 W3: the entry-point behaviour suite

## Declaration (written 2026-11-07 07:34, before the counted runs)

- **Purpose.** Freeze today's entry-point behaviour of the four public module servers on every committed fixture, so that
  KIT-102 to 109 can show what they change.
- **Build and inputs.** Source tree `bb87a463` (HEAD after KIT-100; `R/`, `inst/` and `DESCRIPTION` unchanged), loaded with
  `pkgload::load_all(export_all = TRUE, helpers = TRUE)` (the series runner, `kit/scripts-078/run1.R`, one process per test
  file); `R_LIBS` = kit_100 (`core/browser-smoke.2jjdqX/library`, installed from `53e66b4a`, whose `R/`, `inst/` and
  `DESCRIPTION` equal `bb87a463`'s); R 4.6.1. Fixtures: the nine committed e2e lanes plus the new `trx_screen`
  (tests/testdata/e2e/manifest.json), the offline UniProt cache (tests/testdata/gui-entry/uniprot_cache/), the e2e
  enrichment doubles.
- **Runs (the experiment budget).** The four test files in parallel, each run once per item:
  1. run 1: record every lane and write provisional oracles to a scratch directory (`MULTISCHOLAR_GUI_ENTRY_FREEZE=1`,
     `MULTISCHOLAR_GUI_ENTRY_ORACLE_DIR=<scratch>`);
  2. run 2: the suite against run 1's oracles, under `MULTISCHOLAR_STRICT_SKIPS=1` (lifecycle (a) and (d));
  3. the mutation check (b): one scratch copy of the source tree with one parameter default changed, the trx lane only,
     against run 1's oracle;
  4. the failure check (c): one scratch copy of the tests with one step's actions replaced by an error, the trx lane only;
  5. the freeze run: the suite writes tests/testdata/gui-entry/<lane>.json (`MULTISCHOLAR_GUI_ENTRY_FREEZE=1`); its
     records are compared with runs 1 and 2.
- **Ceilings.** 900 s per lane (asserted in the test); 4,500 s per test file process (the proteomics file holds five
  lanes).
- **Decision rule.** A lane is frozen only if its run 1, run 2 and freeze-run records are identical. One bounded
  correction pass for the whole ticket; a lane that is still nondeterministic after it is reported with its unstable
  captures and not frozen.
- **Development runs, not counted.** While the lanes were built, each lane was run in process to see its steps (scratch
  `dev/`); the trx lane was run twice in a row as a development probe of the volatile-field policy (the only difference
  was the download fix made between them). These are implementation, not the determinism runs.

## The suite

- **Files.** `tests/testthat/helper-gui-entry.R` (driver, captures, doubles, the lanes' steps);
  `test-gui-entry-{proteomics,lipidomics,metabolomics,transcriptomics}.R` (one `test_that` per lane, no skip: a lane
  that cannot run, or has no oracle, fails); oracles `tests/testdata/gui-entry/<lane>.json`; the UniProt cache
  `tests/testdata/gui-entry/uniprot/` (`build.R` writes it from the e2e FASTA).
- **Driver.** `shiny::testServer` on `mod_proteomics_server`, `mod_lipidomics_server`, `mod_metabolomics_server` and
  `mod_transcriptomics_server`, with `MULTISCHOLAR_TEST_MODE=true` (the import's `fileInput` fallback). Inputs are set
  through the namespaced ids only; each tab's UI defaults are sent first (the KIT-075 `omicsUiInputDefaults`), as the
  browser would. **Browser double:** an `update*Input()` / `updateTabsetPanel()` call is applied to the session's inputs
  after the flush (testServer has no client); choices are recorded so the proteomics design picks the runs the builder
  offers. Notifications and modals are recorded. No MultiScholaR internal is mocked. **No network:** the UniProt cache is
  placed where `getUniprotAnnotations()` reads its cache; gprofiler2/clusterProfiler are the e2e enrichment doubles; the
  HTTP(S) proxy points at a closed port for the whole lane.
- **Lanes and steps** (step ids in each oracle):
  - trx_screen (16): import refusal (a negative count), import, design refusal (save with nothing assigned), design
    (factors, groups, technical-replicate sets, contrast), QC filter, QC finalise, DA refusal (before normalisation),
    normalisation refusal (manual k above the cap), normalisation (TMM/logCPM, grid, RUVs), DA, summary silent refusal,
    RUVIII export, save, copy, report (package template), session export.
  - lipid_canonical, metab_lc, metab_gc, metab_combined (20): import refusal (empty assay name), import, design refusal,
    design, intensity filter, duplicates, ITSD refusal ("No internal standards found"), finalise, DA refusal (load before
    any export), normalisation refusal (export before normalisation), normalisation (ITSD, cyclic loess, RUV skipped),
    correlation filter, session export, DA load, DA, summary silent refusal, save, copy, report, session export.
  - prot_dia (29), prot_dia_limpa (24), prot_tmt, prot_lfq, prot_lfq_fragpipe (20): import silent refusal (no file),
    import (search results, FASTA, format, route), design refusal, design, every QC sub-tab of the route's spec tables in
    order (the e2e lanes' permissive settings on these six-sample fixtures), DA refusal (load before any export),
    normalisation refusal (export before correlation filtering), normalisation (RUV skipped; the limpa route's
    pass-through), correlation filter, RUV covariates (limpa), session export, DA load, DA, enrichment (prot_dia only,
    gprofiler2 through the doubles), summary silent refusal, save, copy, report, session export.
- **Captures after every step:** head state, state history, `tab_status`; a content digest of every committed state
  (`getState` in the owner's scope; S4 slots, tables and numbers rounded to 10 significant digits; environments and
  functions excluded); every output by name (text masked, plot/image alt text, table and UI digests, download outputs by
  name); artifacts under the project (paths masked; sha256 of `.tsv/.tab/.csv/.fasta`; others by path; store internals
  counted); the step's notifications; for a refusal, whether the head and history are unchanged and whether it raised a
  message (a silent refusal must raise none).
- **Volatile-field policy** (masks before hashing): the lane's temporary root and `<PKG>` (the package path, so check and
  load_all compare), times, dates, durations, UUIDs, hexadecimal runs of 10 or more characters (tempfile suffixes,
  htmlwidget/plotly ids), shiny's random tabset ids, process ids.

## Results

| Run | What | Result |
|---|---|---|
| run 1 (07:34, 491 s) | record, provisional oracles | 10 lanes run, 138 expectations pass |
| run 2 (07:43, strict) | against run 1 | trx_screen identical; the other nine **differ** (32 to 107 captures): shiny's random tabset ids inside rendered UI (`qc-intensity-assay_results_tabs`, `quality_control-dynamic_qc_tabs`), plotly/htmlwidget ids in the enrichment plot outputs, and tempfile suffixes in export file names (`export_complete_dia-export-<hex>.json`, `filtered_session_<time>_native-export-<hex>.msdia.tar`) |
| correction (the ticket's one pass) | volatile-field policy | hex mask from 16 to 10 characters; tabset ids masked |
| run 3 (07:52) | record, provisional oracles | 10 lanes, 138 pass |
| run 4 (08:01, strict) | against run 3 | seven lanes identical; prot_tmt, prot_lfq, prot_lfq_fragpipe differ in one capture only: the store file count (`<store files>`, for example 538 against 540) at `session_export`, `da_load`, `da` and `summary_copy` (8 captures) |
| (b) mutation (08:10) | scratch copy, `R/mod_lipid_da_ui.R:38` FDR threshold default 0.05 to 0.1, trx lane against run 3 | **fails, first at step `da`**: 22 pass, 7 fail (`da` and the six summary steps after it); every step before `da` identical |
| (c) failure | scratch copy of the tests, the trx normalisation step's first line `stop(...)` | **error 1, skip 0** ("Error in `step$run(api)`: KIT-101 failure check ..."). The first attempt injected nothing (the step pattern matched two lines; the `sed` address was malformed) and passed; that is an instrumentation failure, recorded, and the check was run with the fault in place |
| (d) strict | run 4 under `MULTISCHOLAR_STRICT_SKIPS=1` | the seven stable lanes pass (35 + 105 + 29 + the prot_dia and prot_dia_limpa expectations) |
| freeze (08:11, 481 s) | writes the oracles | identical to runs 3 and 4 for the seven stable lanes; the three protein-route lanes again differ only in the store file count |

**Frozen (7):** prot_dia, prot_dia_limpa, lipid_canonical, metab_lc, metab_gc, metab_combined, trx_screen
(`tests/testdata/gui-entry/<lane>.json`, `source_sha` `e989e93e` = `bb87a463` plus KIT-101's test-only commits; load
mode `load_all`; R 4.6.1).

**Not frozen (3): prot_tmt, prot_lfq, prot_lfq_fragpipe.** Nondeterministic after the correction pass in one capture:
the count of files in the owner's store (`.multischolar`/duckdb/qs2) after the session export, which varies by one to
three between runs. Every other capture of these lanes was identical in runs 3, 4 and the freeze run. Their tests fail
("no frozen oracle") until a decision: drop the store-file count from the capture (or count only stable entries) and
freeze them. The budget's one correction pass is spent, so that is not done here.

**What the oracles show (today's behaviour, frozen as is; see `findings.md`):**
- lipid, metab and trx reach the summary with every step committed and no error notification.
- proteomics: import, design (UniProt from the cache), the protein-route QC, normalisation, correlation filter and
  export run; the DIA routes' Q-value filter and rollup refuse the e2e DIA-NN fixture (columns it lacks), the DA fails on
  every proteomics fixture ("interactive volcano output is missing"), so enrichment (prot_dia) is refused ("complete the
  differential expression analysis first") and the publication copy refuses ("da outputs are pending"). The proteomics
  DA, enrichment and report contents are therefore not covered by these fixtures.
- Refusals: every tab has a messaged refusal except the summary tab (all omics: silent, `req()`-gated), the proteomics
  import (silent), the trx and protein-route QC tabs (none reachable) and enrichment (whose only step is refused).

**After rebasing onto the KIT-101 addendum** (`871ed080`: the F-01 to F-05 fixes in `R/`), the suite was run once on
the rebased tree (`1911c607`, 09:10, 463 s) against the oracles frozen from `bb87a463`'s code: the seven frozen lanes
pass unchanged (lipid 35, metab 105, trx 29, prot_dia and prot_dia_limpa within the proteomics file's 115); the three
unfrozen lanes fail with "no frozen oracle" (3 fail, 0 error). The addendum's fixes change no capture of these lanes.
The check-mode baseline (`test-baseline.md`) also shows the frozen lanes passing under the installed package.

**Runtime:** lanes 13 to 182 s (ceiling 900 s); the four files in parallel 481 to 500 s.

**Driver limits (not defects):** under `MockShinySession`, outputs assigned a plain `reactive()` for conditional panels
(`design-data_available`, `import-file_loaded`, `summary-report_ready`, ...) read as `<error> unused arguments (self,
name)`; the value is deterministic and frozen as such. Inputs rendered by `renderUI` get no browser default unless the
lane sets them (the proteomics format options and the trx RUVIII `M` are set explicitly).

## Revision 2.2 (2026-11-07): new proteomics fixtures, the store capture (F-13), the F-11 refusal kinds; phase 1, not frozen

**What changed in the suite.**
- **Fixtures (F-07).** The five proteomics lanes run on the new synthetic fixtures (`kit-101/fixtures.txt`): 400
  proteins, 4 (DIA) or 5 (TMT, LFQ) replicates per group, and a planted KO effect.
- **Store capture (F-13, the approved second correction).** `guiEntryArtifacts` no longer counts the owner's store
  files (`<store files>`). It records the store's committed entries (`<store entries>`): each entity directory that has
  a `commit.json`, counted per masked entity path.
- **F-11 refusal kinds (coordinator, after the R fixes on `ba019048`).**
  - The proteomics `import_refusal` and every lane's `summary_refusal` are now messaged refusals.
  - The protein routes (TMT, LFQ) have a new step, `qc_protein_intensity_filter_refusal`: a percentile of 150 is
    refused with "intensity percentile must lie in [0, 100]".
  - The trx lane has no QC refusal (see the findings below).

**Runs** (all on `origin/performance-v2` `ba019048` plus this work, load_all, R 4.6.1, kit_100 library):
- **Interrupted, not counted (pre-reboot).** A first determinism double run on `e0a05f84` (before the rebase) completed.
  A second one on `83160cc9` was killed by a host reboot. The reboot also wiped the scratch records of both. This was an
  instrumentation interruption, not a budget attempt. The completed pre-rebase double run had shown 9 of 10 lanes
  identical and prot_dia unstable in the same capture as below.
- **Verification (each proteomics lane once, in parallel, HEAD `2cfe18a1`).** Every lane reaches import, QC,
  normalisation, DA, enrichment (prot_dia only) and the summary/report. DA tables (`fdr_qvalue`):

  | lane | proteins in DA | samples | FDR non-NA | q < 0.05 | of which planted (P00001–P00040) |
  |---|---|---|---|---|---|
  | prot_dia | 385 | 8 | 385 | 15 | 15 |
  | prot_dia_limpa | 393 | 8 | 393 | 41 | 18 |
  | prot_tmt | 400 | 10 | 400 | 19 | 18 |
  | prot_lfq | 400 | 10 | 400 | 17 | 16 |
  | prot_lfq_fragpipe | 400 | 10 | 400 | 24 | 23 |

  prot_dia enrichment (gprofiler2 through the doubles) gives 10 terms (5 up, 5 down).
- **Determinism, runs A and B (HEAD `2cfe18a1`, the four files in parallel).**
  - 8 of 10 lanes are identical: prot_dia_limpa, prot_lfq, prot_lfq_fragpipe, lipid_canonical, metab_lc, metab_gc,
    metab_combined and trx_screen (at `2cfe18a1` the trx lane still had its pre-import QC refusal).
  - After the trx lane correction (`26f09136`), runs C and D of the trx file are identical too.
  - Lane runtimes: 43 to 369 s (ceiling 900 s). Test-file wall times: 48 to 929 s.
- **Unstable after the F-13 correction (named, not corrected: the correction budget is spent).**
  - **prot_dia**, `states/enrichment_completed`: the state's content digest, 6 captures (the steps from `enrichment`
    to `summary_export`). The `EnrichmentResults` state holds plotly objects (`enrichment_plotly`). Their random
    widget ids (14-hex names under `x$visdat`, `x$cur_data` and the `attrs` keys) differ between runs. The volatile
    mask applies to text, not to list names inside a state digest. A diagnosis probe before the reboot showed this
    as the only difference between two runs. The same capture was unstable in the pre-rebase double run.
  - **prot_tmt**, `artifacts/<store entries>`: 14 captures (from `norm` to `summary_copy`). In run A, one committed
    `dia-qc-product-` entry's request id ended in only 9 hex characters (`e2d23a8b9`). The 10-character hex mask did
    not hide it, so it was listed under its own name instead of being counted with the masked entries (8 against 9).
    The pre-rebase double run had prot_tmt identical, so this is intermittent: it depends on the length of a random
    id.
- **Against the frozen oracles.**
  - lipid_canonical, metab_lc, metab_gc, metab_combined and trx_screen (C/D) differ only in the store capture
    (`<store files>` replaced by `<store entries>`, every step) and in `summary_refusal` (the F-11 message).
  - prot_dia and prot_dia_limpa differ from `import_refusal` on, in every step: new fixtures, the F-11 import message,
    and the run now reaching the Q-value filter, rollup, DA, enrichment and report.
  - prot_tmt, prot_lfq and prot_lfq_fragpipe have no frozen oracle.

**Findings from these runs** (`R/` untouched here):
- trx QC refusal before the import. The expression filter's F-11 message (`mod_trx_qc_filter.R`) is not reachable
  through `mod_transcriptomics_server`:
  - the QC sub-module servers start only after the design (`mod_lipid_qc_server`, `qc_trigger`), so a click before
    the import reaches no observer;
  - when the filter observer is created later (no `ignoreInit`), it runs the earlier click. In run A the filter was
    applied at the design step (head `trx_expression_filtered`), before the lane's own filter step.
  The trx lane therefore has no QC refusal step.
- prot_dia (F-15, being fixed elsewhere). After enrichment, `summary_save` warns "The analysis record
  (analysis_record.json) was not written: component checkpoint is incomplete", and `summary_report` warns that the
  analysis record bundle was not exported. Captured as is.
- prot_dia_limpa, normalisation tab. The tab still warns "Please complete the Quality Control filtering steps before
  accessing normalization." on `norm_refusal`, `norm`, `correlation` and `session_export` after the limpa protein QC
  completes (head `peptide_count_filtered`). The F-10 fix does not cover the DIA_LIMPA route.
- prot_dia_limpa, DA. 41 proteins have q < 0.05, but only 18 of them are among the 40 planted changes. The route ran
  normalisation with RUV-III covariates (k = 2, 8 samples). The other routes found 15 to 23 planted proteins, with at
  most 1 false positive.
- MaxQuant/FragPipe zeros. The import keeps a 0 intensity as a quantity (`log2(0 + 1)`). The LFQ fixtures therefore
  carry no missing values (`kit-101/fixtures.txt`).

**E2E proteomics browser lanes, once each** (HEAD `26f09136`/`e5d04531`, R/ equal to `ba019048`; the KIT-097/W5 way:
`CHROMOTE_CHROME=/usr/bin/brave`, the kit_100 library, ceiling 1,800 s per file, one attempt):

| file | lanes | result | seconds | failed step |
|---|---|---|---|---|
| test-e2e-proteomics-dia.R | prot_dia | pass 13, error 1 | 267 | after the normalisation export: `e2e_assert_prot_filtered_session_artifacts` (helper-e2e-browser.R:675), `filtered_session_data_latest.rds` does not exist |
| test-e2e-proteomics-dia-limpa.R | prot_dia_limpa | pass 11, skip 1 | 92 | skips at the QC stage by design (`e2e_run_prot_dia_limpa_qc`: "no browser path to drive yet") |
| test-e2e-proteomics-enrichment.R | prot_dia (gprofiler2, clusterProfiler) | pass 14, error 2 | 456 | the same session-artifact assertion, before the DA and enrichment |
| test-e2e-proteomics-tmt-lfq.R | prot_tmt, prot_lfq, prot_lfq_fragpipe | pass 48, error 3 | 311 | the same session-artifact assertion, after the normalisation export |

- **The DIA lanes now get past the Q-value filter and the rollup.** Import, design, every peptide and protein QC step,
  normalisation and the export all run in the browser, so F-07's blocker is gone.
- **Every lane then stops at one harness check.** The check expects the legacy filtered-session file. The native export
  writes `filtered_session_<time>_native-export-<hex>.msdia.tar` instead, as the entry lanes record.
- **This is a stale assertion in `tests/testthat/helper-e2e-browser.R`, which is outside this work's scope.** The same
  failure occurred in a complete pre-rebase run on `e0a05f84`.

**F-14.** The ticket's experiment budget says "6 in parallel, the KIT-099 checkpoint's settings", but the KIT-099
checkpoint ran 8 in parallel (`head -1 kit-099/suite.txt`). The number for the checkpoint's settings is therefore 8.
The per-file baseline (`test-baseline.md`) ran 6 in parallel, and stands as recorded.

**Not frozen** (phase 1). The freeze waits for F-15, and for a decision on the two unstable captures above.

## Revision 2.2, phase 2a (2026-11-07): F-19 fixture fix, LFQ missing values, F-22; not frozen

Base: `origin/performance-v2` `7e8c7875` (the F-08 to F-12, F-15, F-17, F-18, F-20, F-21, F-23 and F-26 fixes). Work
on top of it:
- `4325329b`, the fixtures (`fixtures.txt` points 1 and 2 still hold for everything else):
  - **F-19, both DIA-NN lanes.** `Precursor.Normalised` now has each run's loading removed: each run is centred on its
    median log2 quantity, with a residual run error of SD 0.1. `Precursor.Quantity` stays raw and byte-identical to
    before. On prot_dia_limpa the per-run median log2 normalised values are now KO 20.17, 19.94, 20.03, 20.18 and
    WT 20.08, 20.08, 20.13, 20.10; before, they were 19.95 to 20.14 (KO) against 20.14 to 20.48 (WT).
  - **F-18 follow-up.** The MaxQuant and FragPipe lanes carry the KIT-021 recipe's about 8% missing values again,
    written as 0 (336 of 4,100 MaxQuant cells). The generator matches all 51 records.
- `aa629c09`, F-22, the third correction (approved):
  - `guiEntryCanonical` reads an htmlwidget/plotly object's random ids as `<ID>`: hex list names and values, and the
    `id`, `elementId` and `cur_data` fields. A check on two plots with the same data gave equal digests, and changed
    data still changed the digest.
  - The store entries mask a store or request id's whole random suffix (`-<ID>`).

**Verification (each proteomics lane once, HEAD `aa629c09`).** Every lane reaches import, QC, normalisation, DA,
enrichment (prot_dia: 10 terms) and the summary/report. The LFQ lanes keep all 10 samples. The F-15 warnings are gone
from prot_dia. Planted recall against the 40 planted proteins (P00001–P00040):

| lane | proteins in DA | samples | FDR non-NA | q < 0.05 | planted | false positives | BH < 0.05 (false positives) |
|---|---|---|---|---|---|---|---|
| prot_dia | 385 | 8 | 385 | 15 | 15 | 0 | 15 (0) |
| prot_dia_limpa | 393 | 8 | 393 | 19 | 15 | 4 | 16 (1) |
| prot_tmt | 400 | 10 | 400 | 19 | 18 | 1 | 19 (1) |
| prot_lfq | 400 | 10 | 400 | 17 | 16 | 1 | 17 (1) |
| prot_lfq_fragpipe | 400 | 10 | 400 | 25 | 24 | 1 | 23 (1) |

prot_dia_limpa now matches the F-19 diagnosis's centred variants: f1 gave 19 / 15 / 4 and f2 gave 22 / 16 / 6.
Before the fix it gave 41 / 18 / 23.

**Determinism, runs A and B** (HEAD `aa629c09`, the four files in parallel): **all 10 lanes identical.** Lanes ran 49
to 373 s (ceiling 900 s). Test-file wall times were 55 to 1,130 s; the proteomics file in run B shared the host with
other work. The two captures that were unstable in phase 1 are now stable: prot_dia `states/enrichment_completed` and
prot_tmt `<store entries>`.

**Against the frozen oracles:**
- lipid_canonical, metab_lc, metab_gc, metab_combined and trx_screen differ only in the store capture and in
  `summary_refusal` (the F-11 message).
- prot_dia and prot_dia_limpa differ from `import_refusal` on.
- prot_tmt, prot_lfq and prot_lfq_fragpipe have no frozen oracle.

**Tests.** test-e2e-manifest.R: pass 237. test-e2e-fixture-integrity.R: pass 1.

Not frozen: waiting for F-24 (the lipid/metab/trx QC observers).

## Revision 2.2, phase 2b (2026-11-07): the freeze

Tree `ddf42372`: `origin/performance-v2` with the F-08 to F-12, F-15, F-17, F-18, F-20, F-21, F-23, F-24 and F-26 fixes,
plus this work (the revision 2.2 fixtures, the F-13 and F-22 captures and the F-11 refusal kinds). load_all, R 4.6.1,
the kit_100 library. The four files ran in parallel in each run.

| run | what | result |
|---|---|---|
| A, B (21:20, 21:40) | determinism | all 10 lanes identical. Also identical to phase 2a's run A on `aa629c09`, so F-24 changes no capture |
| freeze (21:58) | `MULTISCHOLAR_GUI_ENTRY_FREEZE=1` writes the oracles | every lane's record identical to runs A and B; oracles committed in `05ecdcb4` with `source_sha` `ddf42372` |
| verify (22:13, HEAD `05ecdcb4`) | the four files on the frozen oracles | proteomics 187, lipidomics 35, metabolomics 105, transcriptomics 29 passes; 0 fail, 0 error, 0 skip |
| strict (22:27, HEAD `05ecdcb4`) | the same under `MULTISCHOLAR_STRICT_SKIPS=1` | proteomics 187, lipidomics 35, metabolomics 105, transcriptomics 29 passes; 0 fail, 0 error, 0 skip (proteomics file 891 s) |

| lane | deterministic (A = B = freeze) | frozen | steps | lane seconds (freeze run) |
|---|---|---|---|---|
| prot_dia | yes | yes | 29 | 304 |
| prot_dia_limpa | yes | yes | 24 | 170 |
| prot_tmt | yes | yes (first time) | 21 | 121 |
| prot_lfq | yes | yes (first time) | 21 | 122 |
| prot_lfq_fragpipe | yes | yes (first time) | 21 | 132 |
| lipid_canonical | yes | yes | 20 | 98 |
| metab_lc | yes | yes | 20 | 95 |
| metab_gc | yes | yes | 20 | 60 |
| metab_combined | yes | yes | 20 | 72 |
| trx_screen | yes | yes | 16 | 53 |

**Differences from the previous oracles** (frozen from `e989e93e`), attributed to findings:
- lipid_canonical, metab_lc, metab_gc, metab_combined, trx_screen:
  - every step's store capture changed: `<store files>` was replaced by `<store entries>` (F-13), with the request-id
    suffixes masked (F-22);
  - `summary_refusal` is now a messaged refusal, whose notification reads "Please copy the results to the publication
    directory..." (F-11).
  - Nothing else changed. F-17 (trx filter early click) and F-24 (lipid/metab/trx QC early clicks) were expected to
    change no capture, and none changed: the lanes make no early click.
- prot_dia and prot_dia_limpa: every capture from `import` on. These lanes were not attributable step by step to one
  finding, because the inputs changed:
  - **New fixtures (F-07, F-19).** 400 proteins and 8 runs, with the DIA-NN confidence columns and a normalised
    `Precursor.Normalised`. These change outputs, artifact digests, state digests, heads and histories from the design
    step on. The Q-value filter and rollup now commit `qvalue_filtered` and `precursor_rollup`; before, they refused.
  - **F-11:** `import_refusal` is messaged ("Please choose a search results file and a FASTA file..."), and so is
    `summary_refusal`.
  - **F-08:** the DA publishes (no `saved_entity_outputs_pending`). DA, enrichment (prot_dia) and the publication copy
    then reach their outputs: `tab_status` changes from `da` (prot_dia) and from `correlation` (prot_dia_limpa).
  - **F-15:** `summary_save` and `summary_report` write the analysis record, with no warning.
  - **F-10 and F-20:** the normalisation tab no longer warns that QC is incomplete, on the DIA and the DIA_LIMPA route
    respectively.
  - **F-09:** the DIA_LIMPA protein sub-tabs commit `protein_dpc_quant`, `protein_normalised` and
    `peptide_count_filtered`.
  - **F-23:** expected to change nothing, and nothing changed: the route default stays `none`.
  - **F-13 and F-22:** the store capture, as for the other lanes.
- prot_tmt, prot_lfq, prot_lfq_fragpipe: no previous oracle. They are frozen on the new fixtures, with:
  - the F-11 import and summary messages;
  - the protein-intensity refusal step (percentile 150);
  - F-18: the MaxQuant and FragPipe readers read a 0 as NA, so the LFQ lanes keep all 10 samples with about 8%
    missing values;
  - F-26 changes no capture: the TMT fixture writes a missing abundance as empty, not 0.
- F-12 (check-mode guards) and F-21 (e2e session-artifact check) are outside the entry suite and change no capture.
