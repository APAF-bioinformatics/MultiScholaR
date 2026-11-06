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
