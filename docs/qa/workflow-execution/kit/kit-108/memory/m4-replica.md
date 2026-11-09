# KIT-108 W5: route A, measured with the replica at the W4 and W5 trees (M4's first half, 2026-11-09)

**Outcome.** Route A (the DA workbook built one sheet at a time, each transient collected before the next) lowers the copy's peak. It leaves every `results_summary` file and the report text identical. It does not come near D8's threshold. On the replica, the copy's resident peak falls by 0.32 GiB (`108p`) and 0.36 GiB (`108q`), and gc's max used falls by 5 to 7%. What stays resident after the copy is unchanged on `108p` (+0.08 GiB) and 0.32 GiB lower on `108q`. Projected onto the GUI members, the report-step peak stays about 0.8 to 0.9 GiB above D8's limit, and teardown stays 0.7 to 1.3 GiB above it. M4's GUI members (`108v`, `108w`) are not expected to meet D8.

## The fix (W5)

The change is in `R/func_general_filemgmt_results_summary_helpers.R` only, inside `copyToResultsSummary`'s DA-workbook block (`:653-677`, and `:726-727`). The signature is unchanged.

- **The loop:** the `purrr::imap` becomes `for (idx in seq_along(da_files))`. The index, the sheet names (`DA_Sheet<idx>`), the index rows and the order are the same.
- **Release per sheet:**
  - `gc(FALSE)` after the read and before `writeData`, so the read's garbage is collected before the write;
  - `rm(data_content); gc(FALSE)` after the write, before the next sheet;
  - `rm(da_wb, enrichment_wb); gc(FALSE)` after both saves.
- **The failure record:** a failed read still warns. Its `failed_copies` record was local to the former closure: the `<-` assignment inside the closure never reached the function's list. So no record is added now either, and the return value is unchanged. This is a pre-existing defect, left as it was (see Findings).
- **Two variants were measured on `108p`:**
  - **W5a** is the lead's literal route A, with release after each `writeData` only. It lowered gc max used by 85 MB and the raw rise by 0.09 GiB. The rise after `gc()` was 0.30 GiB higher.
  - The block's other large transient overlapped too: each sheet's read garbage (about 1 GB at this size) with its own `writeData`. The same mechanism, a `gc(FALSE)` at the read/write boundary, removes that overlap. That is **W5**, the committed fix. It counts as the work item's one correction pass.

The enrichment block was checked and left alone: none of the three projects has an `*_enrichment_results.tsv`, and the enrichment workbook is 7 kB.

## The replica at W4 and W5

**Setup:**
- the replica is `replica.R --load-all <tree> --keep-work` on each member's `candidate-1` project, with `R_LIBS = kit_107` (`browser-smoke.66sPmE`), R 4.6.1 and `TMPDIR` on `/home`;
- W4 is a `git archive` of `1eead052`; W5 is that tree with the fixed file;
- each run went once, in sequence, with no `multischolar*` unit active;
- the rows are in [m4/](m4/) (`replica-rows-<project>-<tree>.tsv`), and the table below is `m4/m4-replica-w4-w5.tsv`.

**The measures,** all over the copy:
- **gc max used:** Ncells + Vcells MB, since the `gc(reset = TRUE)` before the copy;
- **raw rise:** `RssAnon` at `after_copy_raw` minus `before_copy`, the state the application is left in;
- **rise after gc:** the same at `after_copy`, after the replica's `gc()`;
- **peak rise:** `VmHWM` at `after_copy` minus `before_copy`.

| Project | Tree | Copy (s) | gc max used (MB) | Raw rise (GiB) | Rise after gc (GiB) | Peak rise (GiB) |
|---|---|---|---|---|---|---|
| `108p` metab | W4 | 47.1 | 2,320.3 | 1.773 | 1.379 | 1.772 |
| `108p` | W5a (literal) | 49.2 | 2,235.7 | 1.681 | 1.681 | 1.716 |
| `108p` | **W5** | 54.5 | **2,161.6** | **1.457** | **1.457** | **1.457** |
| `108q` lipid | W4 | 59.9 | 2,281.5 | 1.740 | 1.593 | 1.740 |
| `108q` | **W5** | 52.5 | **2,155.1** | **1.274** | **1.274** | **1.376** |
| `099l` trx | W4 | 0.8 | 543.9 | 0.053 | 0.053 | 0.053 |
| `099l` | **W5** | 1.7 | **508.3** | **0.014** | **0.014** | **0.024** |

**The reduction, W5 − W4:**

| Project | gc max used | Raw rise | Rise after gc | Peak rise |
|---|---|---|---|---|
| `108p` | −158.7 MB (−6.8%) | −0.316 GiB (−18%) | +0.078 GiB (+6%) | −0.315 GiB (−18%) |
| `108q` | −126.4 MB (−5.5%) | −0.466 GiB (−27%) | −0.319 GiB (−20%) | −0.364 GiB (−21%) |
| `099l` | −35.6 MB | −0.039 GiB | −0.039 GiB | −0.029 GiB |

- gc used after the copy is the same at both trees (435.7 MB on `108p`/`108q`, 343.3 MB on `099l`), so the fix adds no live residue.
- The render's raw rise is 0.001 GiB on the metab and lipid projects and 0.059 GiB on `099l`, the same at both trees.
- The copy's time moved by −7.4 s to +7.4 s on the large projects and +0.9 s on `099l`. These are single runs, so that is within run-to-run spread, not a measured cost.
- **The floor of this mechanism.** A stage-by-stage run on `108p`'s two DA files, collecting after every read, write and save (`m4/stages.R`, output in `m4/stages-108p.txt`), shows:
  - each `read.xlsx` alone has a transient of about 1.0 GB;
  - the combined workbook's cell store grows by about 0.41 GB per sheet and is held until the save;
  - the lowest peak reachable by releasing per sheet is therefore the second read over the first sheet's store: 1.54 GB of R heap.

  Route A cannot go lower without changing the reader or the workbook's lifetime, and those are P-30's routes.
- **The `108q` replica is not attribution evidence** (D7: +55% against its GUI rise; [report-memory-diagnosis.md](report-memory-diagnosis.md) §2). Its numbers show the direction, not a validated size.

## D8's content checks (W4 vs W5 replica outputs): all identical

Each `results_summary/<omic>` directory the replica wrote was compared file by file. The timestamped `_backup_` directory of the member's earlier output is excluded. The check script is `m4/content-check.R`, and the tables are `m4/content-<project>.tsv`.

- **Workbooks:** sha256 of each sheet's data frame (`openxlsx::read.xlsx` per sheet, then `digest`). This is `kitFingerprintFiles`' `sha256_xlsx_sheets` method, applied per sheet.
- **The report (`.docx`):** sha256 of pandoc's plain text, with the timestamp pattern masked and the replica's own work path masked. The path appears in the template's "directory not found" warnings and differs between the two work copies.
- **Every other file:** sha256 of the bytes.

| Project | Files | Units compared | Identical | Workbook sheets | Report |
|---|---|---|---|---|---|
| `108p` | 39 | 41 | 41 | 4 (`DA_Results_Index`, `DA_Sheet1`, `DA_Sheet2`, `Enrichment_Index`) | identical (33 lines) |
| `108q` | 44 | 46 | 46 | 5 (as `108p`, plus `lipids_num_sig_de_molecules.xlsx`) | identical (33 lines) |
| `099l` | 20 | 22 | 22 | 5 (`DA_Results_Index`, `DA_Sheet1`, `Enrichment_Index`, two correlation-optimisation sheets) | identical (113 lines) |

**Limit:** the replica's metab and lipid reports are short. Each is 33 lines, three of which are the template's warnings that `Publication_figures/{NumSigDeMolecules, Volcano_Plots, Enrichment_Plots}` are not in the replica's `results_summary`. This is the same at both trees. So the report check covers what the replica renders, not a full GUI report. The GUI members' products (M4's second half) cover that.

## D8's thresholds for the GUI members, and the projection

**Thresholds,** from the baselines (I2, the R process's `RssAnon`; [report-memory-diagnosis.md](report-memory-diagnosis.md) §3):

| Member | DA-step peak | Pre-report level | Report-step peak limit (DA peak + 0.25) | Teardown limit (pre-report + 0.25) | Baseline report peak / teardown |
|---|---|---|---|---|---|
| `108v` vs `108p` (metab) | 1.540 | 1.308 | **≤ 1.79 GiB** | **≤ 1.56 GiB** | 2.980 / 2.945 |
| `108w` vs `108q` (lipid) | 1.398 | 1.398 | **≤ 1.65 GiB** | **≤ 1.65 GiB** | 2.815 / 2.669 |

**The projection.** The replica's reductions are applied to the baselines. The replica starts from an empty heap, while the application starts from about 1.3 GiB that it may reuse, so these figures are indicative only.
- **Metab:** report peak about 2.98 − 0.32 ≈ **2.67 GiB** (limit 1.79). Teardown is about unchanged, **≈ 2.9 GiB** (limit 1.56), because the retained rise after gc did not fall.
- **Lipid:** report peak about 2.82 − 0.36 ≈ **2.45 GiB** (limit 1.65). Teardown about 2.67 − 0.32 ≈ **2.35 GiB** (limit 1.65), from an unvalidated replica.

**Plainly: the replica does not suggest the GUI members will meet D8.** Both criteria are expected to miss by roughly 0.7 to 1.3 GiB on both corpora. Route A's floor (§ above) leaves the read transient of about 1 GB and the workbook's store of about 0.8 GB for both sheets in the copy. glibc keeps the freed pages (H2), so teardown follows the peak. Meeting D8 needs one of P-30's routes: B, `malloc_trim` after the copy, for teardown; or C, the copy on the worker, for both. That is the user's decision.

## Test-first and regressions ([m4/tests.txt](m4/tests.txt))

- **The fix's test** fails on the W4 tree: 664 pass, 1 fail. The copy's gc max used is 518.9 MB against 300.0 MB for both sheets read and held. It passes at W5: 665/0, with 2 skips (Q-175 and Q-183, as at W4).
- **W1:** `test-omic-summary-module.R` gives 1703/0, identical.
- **The report helpers:** `test-general-filemgmt-report-helpers-shared.R` gives 85/0/1, the same as the W4 tree.
- **The ten lanes:** proteomics 188, lipidomics 35, metabolomics 105 and transcriptomics 29, all passing (identical to their oracles).

## Findings

- **F108-W5-01 (pre-existing, left as is).** In `copyToResultsSummary`, the `failed_copies[[...]] <- ...` assignments inside the former DA and enrichment `purrr::imap` closures were local to the closure. A failed DA or enrichment read therefore warned but never appeared in the returned failures. The enrichment closure still has this. The two `saveWorkbook` `tryCatch` handlers (`:712`, `:722`) have the same local-assignment defect. W5 keeps the return value unchanged. A fix would change the returned failures, so it is out of W5's content contract.
- **The fix's test** (`test-omics-summary-core.R`, the "(W5, D8)" block) measures in a fresh R (`callr`) whose collection triggers are set high (`R_VSIZE = 3G`, `R_NSIZE = 40M`). Collections then run only where the code asks for one, and gc's max used counts all garbage alive at once. Without that, gc's max used at this size depends on R's adaptive triggers and the process history: in-process trials varied by up to 90 MB between repeats, and the W4/W5 difference was 2 to 60 MB. With it, the results are deterministic: W4 518.9 MB against 300.0 MB for both sheets read and held (fails), and W5 244 MB (passes; this figure comes from the scratch prototype with the same setup, because the passing test prints no number). The residue after a warmed copy is 0 MB against one sheet's 3.5 MB at both trees.
