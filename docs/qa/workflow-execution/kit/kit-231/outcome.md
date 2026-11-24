# KIT-231 gate G0: outcome

**Stopped after G0, as instructed.** No fix has been written. The run binds to `origin/v0.9-staging` cf8af245 plus `g0.R`:
- commit 2c6c2498: run 45509141744, an instrument defect;
- commit 6bddfc23: run 45511537888, observed.

The script is `g0.R`. Its tables are in `g0-out/`:
- `g0.tsv`: the P, C and F rows;
- `g0-stale.tsv`: S, the stale calls;
- `g0-calls.tsv`: S, all 2,235 calls examined;
- `g0-dirs.tsv`: S, all 199 `project_dirs` reads;
- `g0-symbols.tsv`: S, the free symbols, recorded and not judged;
- `g0-session.txt`.

Run 1's partial tables are in `g0-run1/`.

**Result against the stop rules.** The chunks do **not** run cleanly, so stop rule 1 does not apply. Every point the ticket predicted from its reading (U1 to U8, `:1540`, F1, F5, and the unused arguments named in section 7) failed as predicted, with R's messages recorded below.

S also lists stale calls that sections 2.2 and 10 do not name (section 4 below). Under stop rule 2 these are recorded for the coordinator's ruling.

## 1. Runs

| Run | Commit | What happened | Job-minutes |
|---|---|---|---|
| 45509141744 | 2c6c2498 | Part P was observed in full (18 rows, kept). The fixture then stopped with `could not find function "ProteinQuantitativeData"`: the generator is not exported, and KIT-200's builder runs under `load_all(export_all = TRUE)`. No C, F or S observation was made. This is an instrument defect (VAL-54). The fix resolves the generator from the namespace inside the copied builder only, so the search path that S scans is unchanged. | 12.3 |
| 45511537888 | 6bddfc23 | Instrument rerun 1 of 3. The script exited 0 after 43 s. P, C, F and S were all recorded. P is identical to run 1. | 8.9 |
| **Total** | | against the ceiling of 60 job-minutes (4.0 core-hours) | **21.2 (1.41 core-hours)** |

## 2. Part P: what R 4.6.1 does with empty paths (observed)

| Expression | Result |
|---|---|
| `list(a = 1)$de_output_dir` | `NULL` |
| `file.path(NULL)`, `file.path(NULL, "x")` | `character(0)` |
| `file.exists(character(0))` | `logical(0)` |
| `logical(0) && TRUE` | `NA` (no error) |
| `logical(0) && FALSE` | `FALSE` |
| `if (file.exists(character(0)) && !FALSE) 1` | error "missing value where TRUE/FALSE needed" |
| `if (file.exists(character(0)) && !TRUE) 1` | `NULL` (no error) |
| `dir.create(character(0))` | error "invalid 'path' argument" |
| `ggplot2::ggsave(character(0), p)` (ggplot2 4.0.3) | error "`filename` must be a single string, not an empty character vector." |
| `ggplot2::ggsave(<temp>.png, plot = NULL)` | **no error; a file is written** |
| `grDevices::png(character(0))` | error "'filename' must be a non-empty character string" |
| `grDevices::pdf(character(0))` | **no error. A pdf device opened and wrote a file named `NA` in the working directory** (closed by G0) |
| `saveRDS(1, character(0))` | error "'file' must be a non-empty character string" |
| `vroom::vroom_write(df, character(0))` | error "argument is of length zero" |
| `writexl::write_xlsx(df, character(0))` | error "is.character(path) && length(path) is not TRUE" |
| `createDAResultsForEnrichment(..., de_output_dir = stop("evaluated"))` | error "unused argument (de_output_dir = stop("evaluated"))". The argument was never evaluated. |
| `fs::dir_copy(<shipped proteomics workbooks>, <existing dst>, overwrite = TRUE)` | **The edited copy was overwritten in place.** Its first line is now the shipped `---`. No nested folder was created, and the user's own extra file was kept. |

Against what the ticket expected:
- `&&` with a zero-length operand gives `NA`, not an error (audit 1.2 said "RUN-SETTLES"). So `force = FALSE` errors and `force = TRUE` passes, as the audit's first branch read it.
- `ggsave(plot = NULL)` does **not** error, which stop rule 5's note asked about.
- `pdf(character(0))` writes a file called `NA` into the working directory. The ticket inferred that "nothing goes to the working directory". That holds for the workbooks' path, because `ggsave` fails before any `pdf()`, but not for `pdf()` itself.
- The `dir_copy` behaviour confirms the NEWS draft's clause "replaces any edits made to those files" (OD-2).

## 3. Parts C and F, per workbook (observed)

The fixture is KIT-200's `ProteinQuantitativeData`: 60 proteins × 12 samples, groups A, B and C, seed 200001, with `P11..P16` +1 in B. `protein_deAnalysisWrapperFunction` ran on `B_vs_A=groupB-groupA` and returned all 16 elements, with 60 rows in `de_proteins_long`.

*Not explained, recorded as observed:* 0 rows in `de_proteins_long` have `fdr_qvalue < 0.05`.

Only `writeInteractiveVolcanoPlotProteomics` was mocked in C, and it was called 0 times. In F only `shiny::runApp` was mocked.

G0 ran only the chunks section 7 names, not whole workbooks. So for each workbook, "first stop" below means the first stop among those chunks, in workbook order. The design chunk comes first.

| Workbook | Design chunk (F) | Chunks run in C, in order: R's message verbatim; files written before the stop | Ticket's prediction |
|---|---|---|---|
| `DIA_workflow_starter.rmd` | `:234` stops: "missing value where TRUE/FALSE needed" (at `R/func_general_applets.R:240`). `runApp` is never reached. | U1 `:1725` and U2 `:1768` give `character(0)`. **34 Output** completes. Its messages are "Processing contrast: B_vs_A", six times "Saving 7 x 7 in image", then "Error processing contrast: B_vs_A" and "Error message: `filename` must be a single string, not an empty character vector.". Written before the error: `publication_graphs/PCA/PCA_plot{,_with_sample_ids}.{pdf,png}` and `RLE/RLE_plot.{pdf,png}`. Nothing went into the DA folder or the working directory. **35 Enrichment** `:1765` stops: "unused argument (de_output_dir = file.path(project_dirs[[paste0(omic_type, "_", experiment_label)]]$de_output_dir))". | as predicted (U1, U2, F1) |
| `DIA_workflow_limpa.rmd` | `:340`: **not observed**. The chunk stopped on "object 'omic_type' not found" because G0's design-chunk environment lacked `omic_type`, which the workbook assigns at `:97`. This is an instrument gap and was not rerun. The direct `RunApplet(..., force = TRUE)` row below covers this workbook's `force = TRUE` call. | U3 `:2021`, `:2064` give `character(0)`. **34** behaves exactly as the starter's: the same messages, PCA and RLE written. **35** `:2061` stops with the same unused-argument message. | as predicted for U3. The ticket's earlier stop at `:177` (F4) was not executed by G0. |
| `DIA_workflow_limpa_starter.rmd` | `:298`: not observed (the same instrument gap). | U5 `:1170` gives `NULL`: the `integration_your_experiment` key is absent, and the set-up keys are `proteomics_your_experiment` only. **17 Audit Trail** stops at `:1167` with "'file' must be a non-empty character string" (`saveRDS`), after 2 of its 9 expressions. Written before the stop: `protein_qc/ruv_normalised_results_cln_with_replicates.{tsv,RDS}`. U4 `:1333`, `:1378` give `character(0)`. **20** behaves as the starter's. **21** `:1375` stops with the unused-argument message. | as predicted (U4, U5). The chunk 1b stop (F6) is confirmed by S, not executed. |
| `DIA_workflow_limpa_experienced.rmd` | no uncommented `RunApplet` call in any chunk | U6 `:1470`, U7 `:1563` and U6 `:1663` give `character(0)`. **20** behaves as the starter's. **21** `:1660` stops with the unused-argument message. **extract-top-de-genes** stops at `:1536`, the `imap_dfr` statement, after 1 of its 5 expressions: "ℹ In index: 1. ℹ With name: B_vs_A. Caused by error in `UseMethod()`: ! no applicable method for 'filter' applied to an object of class "NULL"". Nothing was written. **`:1563-1564` alone** stops at `:1564`, `vroom_write`: "argument is of length zero". | as predicted (`:1540-1541` lies inside the statement that starts at `:1536`; U6, U7) |
| `DIA_workflow_experienced.rmd` | `:128` stops: "missing value where TRUE/FALSE needed". `runApp` is never reached. The 'Raw Data QC' chunk `:223` stops: "unused argument (publication_graphs_dir = publication_graphs_dir)". | U8: `de_output_dir` (`:1147`, `:1175`), `publication_graphs_dir` (`:1148`) and `pathway_dir` (`:1186`) each give "object '<name>' not found". **Output** completes with "Processing contrast: B_vs_A", "Error processing contrast: B_vs_A" and "Error message: object 'de_output_dir' not found". **Nothing was written**, no PCA and no RLE. **Enrichment** `:1172` stops: "unused argument (de_output_dir = de_output_dir)". | as predicted (U8: the stop at `:333`, before PCA; F1, F2) |
| `metabolomics_workflow_starter.rmd` | `:130` (the chunk; the call is at `:133`) stops: 'argument "omic_type" is missing, with no default' | none (section 7 names no metabolomics chunk) | as predicted (F5) |
| `metabolomics_workflow_experienced.rmd` | `:130` stops with the same message | none | as predicted (F5) |
| `integration_workflow_starter.rmd` | no `RunApplet` call | none | none predicted |

The F rows on `RunApplet` itself, with `project_dirs` from `setupDirectories`: the `proteomics_kit231` element `output_dir` is absent.
- `RunApplet("designMatrix", "proteomics", "kit231", force = FALSE)` gives "missing value where TRUE/FALSE needed". `runApp` is reached 0 times.
- With `force = TRUE`: no error, and `runApp` (mocked) is reached once.
- `saveRDS(..., file.path(<output_dir element>, "design_matrix_results.rds"))`, the call at `:387`, gives "'file' must be a non-empty character string". N5's base values are therefore these two messages.

## 4. Part S: stale calls in all eight workbooks (observed; file:line)

**Unused arguments** (`match.call` against the formals; the generic's formals for S4):
- `createDAResultsForEnrichment(de_output_dir =)`, five times:
  - `DIA_workflow_starter.rmd:1765`
  - `_limpa.rmd:2061`
  - `_limpa_starter.rmd:1375`
  - `_limpa_experienced.rmd:1660`
  - `_experienced.rmd:1172`
- `updateProteinFiltering(publication_graphs_dir =)`, twelve times in `DIA_workflow_experienced.rmd`: `:223`, `:247`, `:266`, `:311`, `:349`, `:393`, `:427`, `:554`, `:615`, `:675`, `:892`, `:974`.
- `RenderReport(suffix =)` at `DIA_workflow_experienced.rmd:1239`.
- `peptideMissingValueImputationLimpa(dpc_method =)` at `DIA_workflow_limpa.rmd:876`.

S records unused arguments only. It does not flag missing required arguments, so the metabolomics `RunApplet("designMatrix")` appears in F, not S.

**Directory elements and keys missing** for each workbook's own set-up (10 rows, exactly the ticket's N1 list):
- `$de_output_dir` at:
  - `DIA_workflow_starter.rmd:1725`, `:1768`
  - `_limpa.rmd:2021`, `:2064`
  - `_limpa_starter.rmd:1333`, `:1378`
  - `_limpa_experienced.rmd:1470`, `:1563`, `:1663`
- the key `integration_<label>` at `_limpa_starter.rmd:1170`.

All other directories reads resolve; there are 199 reads in total. The experienced workbook's bare `de_output_dir`, `publication_graphs_dir` and `pathway_dir` appear in the free-symbol scan (`g0-symbols.tsv`, chunks 'Output DE Analysis Results' and 'Enrichment Analysis').

**Functions not found** on the search path after `library(MultiScholaR)` and the dependencies:
- `conflict_prefer` / `conflicts_prefer` at `DIA_workflow_limpa_starter.rmd:94-106` (13 calls; F6).
- **Not in the ticket:** `conflict_prefer` / `conflicts_prefer` also at `DIA_workflow_limpa_experienced.rmd:69-86` (18 calls), in chunk "1 MultiScholaR FIRST INSTALL".
- **Not in the ticket:** `ProteinQuantitativeData(...)` is **not exported**, yet three workbooks call it bare: `DIA_workflow_starter.rmd:789`, `_limpa.rmd:1182` and `_experienced.rmd:518`, each in chunk "18 Protein Data S4 Object Creation" or "Protein Data S4 Object Creation". By this scan, the starter workbook stops there with "could not find function" before it reaches DA. This was not executed. The same fact stopped G0's own run 1.

**Counts per workbook.** These count the stale calls above. They exclude the runner-environment rows and the instrument false positives listed after the table.

| Workbook | Unused argument | Directory | Function not found | Total |
|---|---|---|---|---|
| DIA_workflow_starter | 1 | 2 | 1 | 4 |
| DIA_workflow_limpa | 2 | 2 | 1 | 5 |
| DIA_workflow_limpa_starter | 1 | 3 | 13 | 17 |
| DIA_workflow_limpa_experienced | 1 | 3 | 18 | 22 |
| DIA_workflow_experienced | 14 | 0 (plus 3 bare names, in the symbol scan) | 1 | 15 |
| metabolomics_workflow_starter | 0 | 0 | 0 | 0 |
| metabolomics_workflow_experienced | 0 | 0 | 0 | 0 |
| integration_workflow_starter | 0 | 0 | 0 | 0 |

The metabolomics workbooks score 0 in S because their stop is a missing required argument, which F shows.

Rows in `g0-stale.tsv` that are **not** workbook defects:
- `devtools::install_github` / `devtools::load_all` (one row per workbook; eval=FALSE in the metabolomics ones). devtools is not installed on the runner.
- Seven rows `not_found_in_namespace` for `dplyr::any_of` and `purrr::set_names`. These are an instrument false positive: both are re-exports, which G0's `exists(..., inherits = FALSE)` guard does not see. They are exported, and the export check passed.

## 5. Differences from the ticket's reading (results, not explained away)

1. `ProteinQuantitativeData` is unexported and called bare in three workbooks (starter `:789`, limpa `:1182`, experienced `:518`). Sections 2.2, 2.3 and 10 do not list it. On the starter it would be a stop before DA, after the design step (static; not executed).
2. `conflict_prefer` / `conflicts_prefer` in `DIA_workflow_limpa_experienced.rmd:69-86`, the first-install chunk. F6 names limpa_starter only.
3. `logical(0) && TRUE` is `NA`, not an error; `ggsave(plot = NULL)` writes a file without error; `pdf(character(0))` writes a file named `NA` into the working directory.
4. In the four `project_dirs` workbooks, the DA output step wrote the PCA and RLE plots before the error, as 2.2 says. No file went into the working directory on that path.

Under stop rule 2, items 1 and 2 go to the coordinator for a scope ruling before R1. They are KIT-232 material unless ruled otherwise.

## 6. Coordinator's ruling on G0 (2026-11-25) and hand-offs to KIT-232

G0 is accepted as the gate's result: the workbooks stop as read. **After KIT-231 alone no DIA or metabolomics workbook is shown to run to the end.**

The two stale calls outside sections 2.2 and 10 are ruled NOT KIT-231's; they go to KIT-232, with every other stale call G0 S counted that is outside KIT-231's outcome. Hand-offs, file:line (from `g0-out/g0-stale.tsv`):
- `ProteinQuantitativeData` (unexported) called bare: `DIA_workflow_starter.rmd:789`, `DIA_workflow_limpa.rmd:1182`, `DIA_workflow_experienced.rmd:518`.
- `conflict_prefer` / `conflicts_prefer` (package `conflicted` not attached): `DIA_workflow_limpa_experienced.rmd:69-86` (18 calls) and `DIA_workflow_limpa_starter.rmd:94-106` (13 calls; F6).
- `updateProteinFiltering(publication_graphs_dir =)`, unused argument: `DIA_workflow_experienced.rmd:223`, `:247`, `:266`, `:311`, `:349`, `:393`, `:427`, `:554`, `:615`, `:675`, `:892`, `:974` (F2).
- `RenderReport(suffix =)`, unused argument: `DIA_workflow_experienced.rmd:1239` (F2).
- `peptideMissingValueImputationLimpa(dpc_method =)`, unused argument: `DIA_workflow_limpa.rmd:876` (F4).
- The metabolomics design call `RunApplet("designMatrix")` without the required arguments: `metabolomics_workflow_starter.rmd:133`, `metabolomics_workflow_experienced.rmd:133` (F5; seen in F, not counted by S).
- The free-symbol list `g0-out/g0-symbols.tsv` (58 rows, not judged) as KIT-232's starting list for objects read before assignment.

KIT-231 repairs, of G0's counts, only the five `createDAResultsForEnrichment(de_output_dir =)` arguments, the ten directories reads, the experienced workbook's three bare folder names and `:1540`.

**Part P facts for KIT-232 and KIT-230** (R 4.6.1, run 45511537888, identical in 45509141744):
- `pdf(character(0))` writes a file named `NA` into the working directory, without error.
- `ggsave(plot = NULL)` writes a file, without error (ggplot2 4.0.3).
- `fs::dir_copy(overwrite = TRUE)` overwrites an edited workbook copy in place (the application's call, `R/app_server.R:544`): users' edits to a copied workbook are lost at the next session start (F11).

## Not checked

- The design chunks of `DIA_workflow_limpa.rmd` (`:340`) and `DIA_workflow_limpa_starter.rmd` (`:298`). G0's environment lacked `omic_type`, so those rows show "object 'omic_type' not found", which is not the workbook's behaviour. This was not rerun; the direct `RunApplet(force = TRUE)` row stands for their `force = TRUE` call.
- Whole workbooks were not executed. Chunks between set-up and DA did not run (KIT-232's executed-workbook check).
- shiny's re-raise of an error inside `stopApp()`. With `runApp` mocked the server never runs, so only the equivalent `saveRDS` was observed.
- `processEnrichments` was not run (it is network-dependent and outside G0).
- `loadDependencies()` itself was not called. Its listed packages that were installed were attached in its order; three were not installed and not installed by G0 (`shinyWidgets`, `shinythemes`, `shinycssloaders`; see `g0-session.txt`).
- Why the fixture's B_vs_A comparison gave 0 proteins at `fdr_qvalue < 0.05`.
- The free-symbol list (`g0-symbols.tsv`, 58 rows) was recorded and not judged.
