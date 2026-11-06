# KIT-096 outcome: the transcriptomics import, design and QC tabs (2026-11-06)

Ticket revision 1.1. Plan: [plan.md](plan.md), written before any code. In-process and module CI only: no app launch, no
kit journeys (KIT-098).

## What was done

- **Reuse by parameterisation (F01).** The lipid import, design-builder and QC modules take one `profile` argument.
  `NULL`, the default, is lipidomics: every lipid call is made as before (same functions, same arguments; the profile is
  passed to a collaborator only when it is not `NULL`). A profile supplies the omic's collaborators at the seams these
  modules already take (`buildModulePanel`, `handleProcessRequest`, `validateImportState`, `writeImportArtifacts`,
  `applyResultToWorkflow`, `runBuilderObserverShell`, the finalise's `get/validate/save/update/report` functions), plus
  labels with lipid fallbacks. Nothing is copied. Metab and proteomics modules are not edited.
- **The profile** (`trxOmicProfile()`, `R/mod_transcriptomics.R`): labels and the `trx` test-id prefix, class
  `TranscriptomicsAssayData`, workflow type `TRANSCRIPTOMICS`, format `trx_counts`, the id and name column labels, the
  preview loader and processing, run order (`identity`: the counts file), the level order, the design check, object,
  provenance and state name, the QC title, sub-tabs, sub-module servers, metadata and state readers, summary and the
  finalise's collaborators.
- **Wiring (F13).** `mod_transcriptomics_ui`/`_server`: the bus, the session owner and the runner install of
  `mod_lipidomics.R`, tabs import, design and QC, the shared stepper (Import, Design, QC), QC enabled when the design
  completes. `app_server.R` dispatches by name and is not edited; the stale "(finding m-2)" comment above the
  transcriptomics tile in `app_ui.R` is removed (the integration tile's comment stays).
- **Import (F05, F06, F14).** One format, the notebook's count table read with vroom (as the notebook and KIT-095's
  journey read it). The user picks the gene id and name columns (`GeneID`, `GeneName` proposed); every other column is a
  sample, in file order. KIT-093's count policy runs at import (`checkTrxImportCounts`): a refusal keeps KIT-093's class
  and names it in the message the tab shows, e.g. `NA or NaN count (...) (trx_counts_na)`; there is no rounding. No
  sanitisation (the control is not offered). No vendor-normalised, internal-standard or second-assay controls. The help
  text states the whole-number policy. The owner is bound to `TRANSCRIPTOMICS` at import; the summary file is
  `transcriptomics_import_summary.tsv`.
- **Design (F04, F06, F07, F10, F15).** The lipid builder and its commit under the profile:
  - a fresh design keeps the counts file's run order (no `mixedsort`) and names;
  - at save the builder's result carries the level order the user entered (`trxDesignGroupLevels`: the factor levels in
    the order added, restricted to assigned groups) and the runs the user removed;
  - the commit runs KIT-093's `validateTrxDesign` against every imported count column except the removed runs (so an
    unassigned run is refused as `trx_design_runs_extra`), by name only, plus at least two runs per group
    (`trx_design_group_singleton`, a child of `trx_design_invalid`); each refusal is shown with its class;
  - the state `trx_raw_data_s4` is `createTranscriptomicsAssayData(list(counts = ...), design, group_levels, roles =
    list(gene_id_column, annotation_id_column, technical_replicate_id = "tech_rep_group"))`, the owner bound to
    `TRANSCRIPTOMICS`, provenance through `recordTrxDesignProvenance`; the completion text is unchanged
    ("Design matrix and contrasts saved successfully!");
  - the "Import Existing Design" button is not offered under the profile (its directory format is the lipid import's);
    the help text says the design is matched by sample name only.
- **QC (F08, F09).** The QC tab under the profile has two sub-tabs:
  - **Expression Filter** (`R/mod_trx_qc_filter.R`, new): KIT-095's `trx_expression_filter` with no library-size
    argument; managed it is `omicsManagedCommit(trxExpressionFilterStep())` then `omicsManagedProduct("trx_qc_figures",
    list(stage = "raw"))`, in process `runTrxReferenceExpressionFilter` and the state read for one action. It shows
    genes before, after filterByExpr and after the zero-variance step, library size and genes detected per run (the
    filter's `filtering_metrics` product), the filtering figure and the raw-logCPM PCA and RLE. Revert reuses the lipid
    intensity sub-tab's revert observer.
  - **Finalize QC**: the lipid finalise under the profile: the store step `trx_qc_finalize` -> `trx_qc_complete`, the
    worker's `trx_filtering_metrics` on a managed session (the profile's runner and operation), a data summary from the
    state's metadata (row count and design; no count read), no per-assay table.
- **Stable ids (F11)**: module `transcriptomics_workflow`; tab input `transcriptomics_tabs` (`import`, `design`, `qc`);
  sub-modules `import`, `design`, `qc` (QC sub-tabs `filter`, `s4_finalize`; buttons `qc-filter-apply_filter`,
  `qc-s4_finalize-finalize_qc`); the builder's shared ids and completion text; test ids `trx-workflow-tabs`,
  `trx-tab-import`, `trx-tab-design`, `trx-tab-qc`, and `trx-import-format`, `trx-import-assay1-file`,
  `trx-import-id-col`, `trx-import-annotation-col`, `trx-import-process`, `trx-qc-apply-filter`.

## Evidence

- **Test first** ([test-first.txt](test-first.txt)): with only the tests and the test data written, every new transcriptomics
  test errored or failed (import 5 errors, design 6, QC 4, app 2 failures and 3 errors), and the trx residency journey and
  both KIT-082 probe tests failed. The manifest test passed then (it validates the data written with it), as did the
  walk's positive control (it tests the walk).
- **Targeted** ([tests.txt](tests.txt), final tree): `test-module-ci-trx-import.R` 96, `-design.R` 71, `-qc.R` 30,
  `-app.R` 22, `test-module-ci-manifest.R` 33, `test-workflow-omics-residency.R` 97, `test-workflow-omics-application-reads.R`
  216 pass; 0 fail, 0 error, 0 skip in each.
  - **Import (MCI-026):** the screen corpus and the whole-number-as-double fixture are accepted with names and file order
    unchanged, the table held by reference, the owner bound to `TRANSCRIPTOMICS`, the source copy byte-identical; NA,
    negative and fractional are refused, the class shown in the one error notification and raised by
    `checkTrxImportCounts`; nothing is held or bound and the tab does not complete.
  - **Design (MCI-027):** the fresh design is `truth$column_order`; the committed object's inputs keep those names and that
    order; the level order entered is recorded (both orders tested); technical replicates partition the runs as the
    corpus's biological samples; the record is valid. **Identity with KIT-095's journey setup:** on the same corpus and
    design, the GUI owner and `trxDeclarationsOwner(..., technical_replicate_id = "tech_rep_group")` have identical history,
    head, workflow id, state object (its store attributes aside), state metadata (class, slots, assay names and rows,
    design), `design_builder` provenance (operation, parameter names, values) and inputs (`trx_counts_1`, `design_file`:
    role, bytes, sha256). Refusals by class: `trx_design_runs_missing` (a design run with no count column),
    `trx_design_runs_duplicated`, `trx_design_level_unknown`, `trx_design_group_empty`, `trx_design_runs_extra` (an
    unassigned run), `trx_design_group_singleton`; the head stays `initial`. Removed runs leave design and counts; three
    groups are accepted in the order entered.
  - **QC (MCI-028):** in process the filter commits `trx_expression_filtered` and the finalise `trx_qc_complete`; the
    counts at each step equal `trxExpressionFilter`'s masks; metrics equal `trxFilteringMetrics`; the PCA and RLE equal
    `trxPcaStats`/`trxRleStats` of the raw logCPM; the record is valid. **Mask:** on the KIT-040 `trx_edge_12` corpus
    (hash-checked), imported and designed through the tabs, the filtered count matrix (runs put in the reference's design
    order by name) has the reference's `count_matrix_filtered` sha256, and the filterByExpr count equals the reference's.
    **Managed:** the worker ran `trx_expression_filter`, `trx_qc_figures` and `trx_filtering_metrics`, every action was
    delivered, and what the tab shows equals the in-process values; the raw PCA and RLE PNGs were staged to `gene_qc_dir`.
  - **App (MCI-029):** `app_server`'s environment finds `mod_transcriptomics_ui`/`_server`; the UI renders the module, not
    the fallback; the stable ids above are asserted; the server registers the bus and the owner and opens QC after the
    design.
  - **Residency (KIT-075 walk):** the trx journey under `testServer` on `mod_transcriptomics_server` (import, builder
    inputs, filter, finalise) reaches `trx_raw_data_s4`, `trx_expression_filtered`, `trx_qc_complete` with no
    `TranscriptomicsAssayData` and no table with the corpus's sample columns reachable, and the bus holding references
    only, after every step. Positive control: a held object and a held count table are found.
  - **KIT-082 probe:** on a managed session the QC tab (filter and finalise, every output read) materialises none of the
    owner's count tables and leaves nothing reachable; in process the same tabs read them (the control). Positive control:
    a variant that reads and holds the state is caught by both the probe and the walk.
- **Dependents** ([dependents.txt](dependents.txt)): 43 files. First pass two failures, both fixed in the one correction
  pass (the observer passed `profile = NULL` to a seam mock; the generator looked for the corpus tool under its output
  root). Afterwards every file passes with 0 fail and 0 error: lipid import/QC/design seams and characterisations, the
  lipid and metab module CI (6 + 6 files), the lipid and metab declaration tests (259 each), the frozen dump
  (`test-workflow-omics-declarations.R` 68), the KIT-094 characterisation (51), KIT-095's binding test (308),
  `test-module-ci-fixture-integrity.R` (2559). **No lipid or metab test file was modified.**
- **Generator:** `generate-fixtures.R --root <tmp>` reproduces `fx/transcriptomics` byte for byte.
- **R CMD check** ([check.log](check.log), `--no-tests --no-manual --no-build-vignettes`, HEAD `f874cd09`, kit_095 library):
  3 WARNINGs and 3 NOTEs, the sections of KIT-095's baseline; the dependencies WARNING and the code-problems NOTE bodies
  are identical to KIT-095's; Rd usage and code/documentation checks OK.
- **Install** ([installed-test.txt](installed-test.txt)): `kit_096` = `$P/core/browser-smoke.1VgwvX` (source `f874cd09`),
  installed through the kit runner from kit_095's library (RUVSeq, EDASeq, ShortRead, Rsamtools, aroma.light carried);
  `test-module-ci-trx-design.R` against the installed build: 71 pass, 0 fail, 0 error, 0 skip. Path in
  `$P/kit/scripts-096/install.txt`.
- **Full suite:** not run (the per-ticket rule: targeted, dependents, check, install; the next checkpoint runs it).

## Divergences and notes

- **Outside the file list:** `DESCRIPTION` Collate gains the two new files (needed to load them; as KIT-095); `man/`
  usage and arguments of the five Rd files whose exported functions gained `profile` (otherwise a codoc mismatch); the
  roxygen `@param profile` lines match.
- **Conventions (F12):** the shinydashboard/`wellPanel` shell, the managed runner in place of `ExtendedTask` and the
  `fileInput` fallback are kept, to match the other omics (the user's instruction) under the scope rule of
  `conventions/R.md`. New functions are about 50 lines or fewer; the existing lipid functions keep their length (scope rule).
- **Not changed:** `fixture_packs.json` and `helper-module-ci-fixtures.R` (the MCI-002 catalog is the vendor-format class
  matrix; the transcriptomics fixtures are KIT-040's, pinned by their manifest and `fixture.json`). Regenerating the
  catalog with `--root` already differs from the committed one in its paths (`fixtures/` vs `fx/`); that predates this
  ticket and the integrity test passes.
- **Run order (hand-off to KIT-098):** a GUI design keeps the counts file's order (the ticket's rule); KIT-040's reference
  is in the design file's order, so a GUI journey's matrices are the reference's with the runs permuted. The mask test
  puts runs in the reference's order by name; KIT-098's identity check needs the same alignment, or a decision.
- **Config:** the design state's config holds the builder's formula field (`deAnalysisParameters$formula_string`), as the
  lipid commit's does; KIT-095's journey saves an empty config. Not part of the identity compared above; the DA formula
  is realised from the design (KIT-095).
- **Contrasts (hand-off to KIT-097):** the builder's contrasts table is the lipid shape (`contrasts` without the name,
  `friendly_names`, `full_format`); KIT-095's DA takes `trxContrastsTable` entries (`name=expression`), i.e. `full_format`.
- **Source copy name:** the import copies the file as `import_sources/assay1_<file>`; KIT-095's journey names it
  `counts.tsv`. Role, size and sha256 of `trx_counts_1` are identical; only the recorded basename differs.
- **Level order** is the order the factor levels were added, restricted to assigned groups; an empty declared group can
  therefore only reach the commit from a result built outside the builder (tested at the commit).
- **Logs:** the reused lipid functions keep their lipid log lines (for example "Lipidomics QC: Initializing sub-module
  servers", "Reverted lipid intensity filter"); user-facing text is the profile's.

## Hand-offs

- **KIT-097:** the module and the profile; the design state `trx_raw_data_s4` with `technical_replicate_id =
  "tech_rep_group"` and `args$group_levels`; `trx_qc_complete` as the QC head; the contrasts note above.
- **KIT-098:** the stable ids above; the run-order note.

Rollback: `git revert` of this ticket's commits; the evidence is kept.
