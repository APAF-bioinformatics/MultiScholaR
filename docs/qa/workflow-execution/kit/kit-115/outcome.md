# KIT-115: defects left by the GUI workover audit (2026-11-11)

**Outcome:** eight of the eleven items are fixed, each with a test that failed on the tree before its fix. One is closed as not a defect. One is stopped. One waits for the user. No oracle, lane capture or stored reference changed.

Every count below is `load_all`, one process per file, with `R_LIBS` at `browser-smoke.kH2HPy/library` (kit_920's).

## Items

| # | Item | Result | Commit | Evidence |
|---|---|---|---|---|
| 1 | F108-04 | **Fixed.** `resolveWorkflowGitInfo` no longer looks in the working directory or `~/Documents/MultiScholaR`. It uses the package's own checkout, or `MULTISCHOLAR_PATH` when the user sets it; otherwise the version block says "Installed Package" with the package version. Reproduced on the installed `kit_920` build first: from an unrelated repository the block named that repository and its branch. | `3f5a1540` | `test-omics-summary-core.R` (new block, 5 expectations); proven on an installed build as well. Dependents: the three `test-general-filemgmt-*` files, 87/0. |
| 2 | A101-09 | **Fixed.** Reproduced: a failed UniProt retrieval (an empty table, F-28) was written to `uniprot_dat_cln.RDS`, which every reader takes as the project's annotations, and the step said "retrieved successfully". Now no file is written and the notice is a warning that says what happened and how to retry. | `b690de16` | `test-prot-design-core.R` 37/0. Dependents: `test-prot-04a-design-module-characterization.R` 109/0/18 skipped, `test-omic-gui-config.R` 1854/0. |
| 3 | A107-07 | **Fixed.** An enrichment replay warns (class `workflow_replay_resource`) when a service it queried is not among the record's `external_resources`, by name and release. A replay against the recorded release says nothing. | `5679f928` | `test-workflow-analysis-record.R` 292/0. |
| 4 | A107-09 | **Fixed.** The status text and the download's summary give the cutoffs of the run whose results the tab shows: the session's Run, or after an Open the completed state's step provenance. They no longer follow the inputs after the Run. The values are held inside the module, so the tab's recorded fields are unchanged. | `d24db76a` | `test-omic-enrich-module.R` 1066/0 (two new blocks: inputs changed after a Run; an Open with other inputs). The first attempt put the values on the tab, which changed seven frozen field captures; the values were moved into the module (the item's one correction pass). |
| 5 | A103-06 | **Fixed.** Reproduced: a contrasts table saved without the `group` prefix ("KO-WT") was reopened with no denominator and saved again as "KO-WT-". It now survives reopening and saving under either formula, and a contrast that cannot be split is written back as it was read. | `5936276c` | `test-omics-design-core.R` 61/0/10 skipped; `test-omic-design-module.R` 160/0. |
| 6 | F112-06, remainder | **Fixed.** After a project swap the step refused every click with "retry pending publication with its original parameters first". A publication left pending by a replaced project is now dropped and the step runs afresh. On the same project the refusal stands. | `29b30a00` | `test-workflow-execution.R` 118/0/10 skipped. |
| 7 | A110-04 | **Fixed.** The validator refuses `norm$run$handler` without `norm$run$settings`. | `b242fe77` | `test-omic-gui-config.R` 1854/0. |
| 8 | Hex text | **Fixed.** "0x1A" was read as 26. It now refuses the import, in every proteomics reader and in the shared numeric-like rule (lipid, metab, trx). | `f177efc4` | `test-prot-import-quantity-text.R` 344/0; `test-workflow-omics-import.R` 24/0/8 skipped. |
| 9a | DIA-NN parquet | **Not a defect.** The parquet path already goes through the check. A test now proves it: clean, a token, a decimal comma, an error cell. | `f177efc4` | as item 8 |
| 9b | `.xlsx` error cells | **Stopped.** `readxl` reads an Excel error cell as NA, so the reader cannot tell it from an empty cell. Telling them apart needs another reader option or a cell-level read, which is a design choice. It is also not testable without a real `.xlsx` file containing error cells. | none | Raised for the user. |
| 10 | Worker error text | **Closed as not a defect at the handoff.** The worker's handed-over error equals the same step's error in process, hint lines included (`test-prot-norm-core.R`, W1 (a)'s frozen value). The difference is in the notice only: in process it shows `$message`, through the worker the full text. The first attempt changed the handoff and broke that equality, so it was reverted. | `6f24978a`, reverted by `8e8a6cdc` | `test-prot-norm-core.R` 1642/0 after the revert. Whether the managed notice should drop the hint lines is a wording choice, raised for the user. |
| 11 | P-02 interim label | **Not started:** waits for the user's word. | none | |

## Checks

- **Dependents:** `test-omic-da-module.R` 1369/0, `test-omics-enrich-core.R` 241/0/1 skipped, `test-prot-norm-core.R` 1642/0.
- **Entry screens:** lipidomics 35/0, metabolomics 105/0, transcriptomics 29/0, proteomics 193/0. The same counts as the full suite at `cc8764a4`: no capture changed.
- **`R CMD check --no-tests`** on a build of `8e8a6cdc`: 3 WARNINGs and 3 NOTEs, the series baseline. Install: done.

## For the user

- **9b:** should an Excel error cell be refused at import? This needs a different read of `.xlsx` files.
- **10:** should the managed failure notice show the first line only, as in process? It would hide the hint lines (for example "choose a smaller k ...") from users of the worker path. Recommendation: keep the hints, and add them to the in-process notice instead. That changes the in-process proteomics lane captures, so it belongs with the next re-freeze.
- **11:** the P-02 interim label.
