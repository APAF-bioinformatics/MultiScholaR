# Harness revision v1.11: inventory of deferred items (2026-11-19, DRAFT for the user)

Branch `kit-harness-v1.11` from `origin/kit-129-integration` 039ce71d. Every item below was deferred by a ticket to
"the harness revision", "protocol v1.11" or "the single re-record". Kinds:

- **H**: harness-bound, in this branch's change set (committed, not executed).
- **P**: package-side, not harness. It lands in its own commit on the laptop day, before the re-record.
- **R**: a stored reference to re-record.
- **G**: a gate option to switch.
- **D**: open for the user (see `re-record-plan.md`, section 7).

Line numbers are on this branch's tip unless a ticket's own anchor is quoted.

## 1. Harness-bound (H): the change set

| # | Item | Source | What changes (file:line) | Commit |
|---|---|---|---|---|
| H1 | The step maps' RUV and DA rows follow the RUV design presets | KIT-129 s4, v1.11 item 1 | `docs/qa/workflow-execution/kit/kit-001/step-map{,-protein,-limpa,-lipid,-metab,-trx}.json`, the RUV and DA rows only (new kernels with file:line, writes, holders, notes; unit `whole_matrix`, no flag) | e65aee22 (b1) |
| H2 | `kitJourneySteps` and the driver follow the presets | KIT-129 s4 (`workflow_kit.R` kitJourneySteps; driver lines 260-290, 620-627, 775), v1.11 item 1 | `tools/profiling/workflow_kit.R`: `kitRuvDesign` (a binding without `ruv_design` reads as before) and `kitJourneySteps(route, ruv_mode, ruv_design)`; the single-journey binding check passes the settings. `tools/profiling/workflow_kit_driver.R`: `kitSelectRuvDesign`, plus the DIA_LIMPA RUV, classic normalisation, lipid/metab and transcriptomics normalisation branches. `tests/testthat/test-workflow-kit-harness.R`: one new block. | 60d3b200 (b2) |
| H3 | `checkCheckpointSummaryFiles` expects KIT-061's integration input | KIT-061 merge gate; KIT-129 s4; U5 (pin lift); v1.11 item 3 | `tools/profiling/run_workflow_checkpoint_proof.R:1905` (the integration rows) and `:1945-1952` (the check): the four files under `integration/proteomics/`, checked with `readIntegrationInput`. The RUV table names at `:1911-1914` are unchanged, per KIT-126's W0 table (the corrected matrix stays the RUV state's table). `docs/qa/workflow-execution/kit/bindings-generator.R:1588-1590`: the CORE check skips the lifted file and binds its hash (`core_pin_lifted`). | a91cfeee (b3) |
| H4 | The nine transcriptomics digest expectations off the laptop; the corpus manifest covers counts and design, not `truth.json` | KIT-129 "Added 2026-11-17 (B-01 closure)"; KIT-134 B-01; v1.11 item 4 | `tests/testthat/helper-workflow-trx-declarations.R`: `trxDeclCorpusFiles` (the rule on every machine; replaces the run-tests switch `MULTISCHOLAR_TRX_MANIFEST_NOTE`), plus `trxDeclDigestOffLaptop` and `expectTrxDeclDigest`. `tests/testthat/test-module-ci-trx-summary.R:65`. `.github/workflows/module-ci.yml`: env `MULTISCHOLAR_TRX_DIGEST_RULE=values`, and a full-gate step that runs `compare-results.R` and fails outside the rule. | c2fef76f (b4) |
| H5 | The harness readers take both folder layouts | KIT-052 (the harness-bound names table; option A, 2026-11-19); v1.11 item 6 (new) | `workflow_kit.R`: `KIT_*_RE`, `kitStoreDirectory`, `kitRequestDirectory`, plus the readers formerly at :1101, 1275, 1322, 1534, 1744-1748, 1804, 1937-1938, 1971-1972. `workflow_kit_coordinator.R:240-250`. `workflow_checkpoint_dia_driver.R:881, 1051, 5079` (`:1046` `.dia-native-import-` unchanged). `workflow_core.R:481-485`. `workflow_perf106.R:171-268, 486, 1006`. `tests/testthat/helper-gui-entry.R:168-169`. `tests/testthat/helper-workflow-bulk-rollup.R:10`. | b3dc89aa (a1) |
| H6 | The package writes `.n-<8>`, `.we`, `.wx`, `w-<16>`; earlier projects still open | KIT-052 decision 2026-11-19 (option A) | See section 2 (the package side of H5/H6). | 4a7c211f (a2) |
| H7 | The protocol text | KIT-129 s4 (five items, approved U11) + KIT-052 option A | `protocol-v1.11-draft.md` (this folder). It goes into `protocol.md` as section 13 on the laptop day; that changes the protocol hash. | not in `protocol.md` yet |
| H8 | The bindings regenerated once: protocol hash, harness file hashes, the six step-map hashes, `ruv_design` per journey, `core_pin_lifted` | KIT-129 s4, W9 | `bindings-generator.R` (the structure is ready, H3). New ticket keys for the B0-B5 members and their builds are written on the laptop day; the installs do not exist yet. | laptop |
| H9 | The lifecycle pair through the new input (complete, forced failure) before B4's journeys | KIT-129 s4, W9 | No code. A run. | laptop |

## 2. Package side of the folder names (H6, commit 4a7c211f)

| Name | Before | Now | Old projects | Sites |
|---|---|---|---|---|
| Session folder | `tempfile(".dia-native-")` (27 characters with `/`) | `.n-<8>` (`newWorkflowSessionFolder`) | open as before (the path is recorded); `removeDIAReplacedSessionFolder` accepts both | `R/utils_workflow_dia_sessions.R` (three creators, the cleanup rule) |
| Store | `.workflow-entities` | `.we` (`workflowEntityDirectory`) | an existing `.workflow-entities` opens and is never renamed | `R/utils_workflow_saved_entities.R` (newWorkflowEntityStore), `utils_workflow_analysis_record.R:687`, `utils_workflow_execution.R:326, 2347` |
| Request/staging root | `.workflow-export` | `.wx` (`workflowExportStagingRoot`) | an existing root is used, so pending journals are found | `utils_workflow_saved_exports.R` (workflowExportPaths, workflowStagedRequest, the Q-value journal), `utils_workflow_dia_sessions.R` (findDIANativePendingExport), `utils_workflow_execution.R` (the completion check) |
| Worker request id | `wfx-<32>` | `w-<16>` | `isWorkflowRequestId` accepts both; `workflowRequestUnitKey` keeps the 12-character unit names | `utils_workflow_execution.R` (id, three validators, two unit names) |
| Native archive member prefix | `.workflow-entities/` | the store folder's own name (`.we/` for new stores) | archives with `.workflow-entities/` open; manifest schema unchanged; an older build refuses a `.we` archive loudly | `utils_workflow_dia_sessions.R` (member check, required entries, staging, open) |

Not renamed: `.dia-native-import-`, `.dia-native-legacy-`, `.dia-gui-adoption-`, the GUI snapshot id `wfx-<32>` (a message
name, not a folder), the scratch `tempfile("wfx-")` under `tempdir()`. No ticket names them.

Tests: 15 test files adapted (test code that assumed the old names on a freshly written store). Hand-built old layouts
are kept as compatibility tests. Three figure-path fixtures were lengthened so their premise (a path of 260 or more
characters) still holds. A new file, `test-workflow-folder-names.R`, has five tests.

## 3. Package-side (P): not harness, land before the re-record

| # | Item | Source | What changes | When |
|---|---|---|---|---|
| P1 | Stop writing `D66_debug.log` in the lipid/metab reports | KIT-146 decision 3 | `inst/reports/lipidomics/lipidomics_report.rmd:84`, `inst/reports/metabolomics/metabolomics_report.rmd:81` (and the lines that write to `d66_log_file`) | own commit before step 6 of the plan; moves R1 (four entry oracles list the file) |
| P2 | The KIT-145 fixture swap and manifest | KIT-145 hand-over 1-3 | `tests/testdata/e2e/pending-kit-145/<lane>/` over `tests/testdata/e2e/<lane>/` for `metab_lc`, `metab_gc`, `metab_combined`, `lipid_canonical`; `tests/testdata/e2e/manifest.json` (`files_sha256` as listed in KIT-145, `fixture_recipe`); remove `metab_lc/lcms_features.tsv`, `metab_combined/combined_features.tsv`; generator text | with G2, one commit |
| P3 | The transcriptomics report template's static RUV paragraph follows the design | KIT-129 sixth pass, open item ("transcriptomics report template unchanged (static paragraph)") | `inst/reports/transcriptomics/transcriptomics_report.rmd:59-63` | before B5; else trx report text moves twice |
| P4 | KIT-135/136 integration-input writers (lipid, metab, trx) | KIT-129 s5 "Order with KIT-135 and KIT-136: after B5, or their oracles join this re-record"; KIT-135/136 report | `origin/kit-135-integration-inputs` 6501e008 (3 commits on 0cb92e0e) | D2: merge before B5 (recommended) or after |

## 4. Gates to switch (G)

| # | Gate | Where | Switch | Source |
|---|---|---|---|---|
| G1 | `MultiScholaR.ruv_design_panel` | `R/utils_omic_gui_config.R:1322-1324` (`omicGuiRuvDesignEnabled`) | default `FALSE` to `TRUE` in its own commit (B5, "the default switch"; Recommended pre-selected, technical replicates on designs that declare them, M2/M4, U2) | KIT-129 s4 B5 |
| G2 | `MultiScholaR.cyclic_loess_minimum` | `R/func_omics_norm_shared.R:19-21` (`omicsCyclicLoessMinimumActive`) | default `FALSE` to `TRUE`, plus `tests/testthat/test-omics-cyclic-loess-minimum.R` (first block, two "rule off" blocks set the option `FALSE` explicitly) | KIT-145 hand-over 4 |

## 5. Stored references to re-record (R), by cause

| # | Reference | Causes (the mechanism each difference is attributed to) |
|---|---|---|
| R1 | GUI-entry oracles, all ten (`tests/testdata/gui-entry/*.json`) | H6 (`<store entries>` keys: `.dia-native-<ID>/.workflow-entities` to `.n-<ID>/.we`; every lane); KIT-128 (trx_screen from design on); KIT-126/127/G1 (RUV and DA rows); KIT-129 UI (the norm tab's UI digest, every lane); KIT-061 already re-frozen for proteomics; P1 (`D66_debug.log` row: lipid_canonical, metab_lc, metab_gc, metab_combined); P2+G2 (the four feature-table lanes from import on); P4 (new `integration/<omic>/` rows: lipid_canonical, metab_lc, metab_gc, metab_combined, trx_screen) |
| R2 | Frozen keys in `test-omic-da-module.R` (335 `<store entries>` keys, plus 14 DA-field hashes from KIT-128) | H6; KIT-128 |
| R3 | `test-omic-enrich-module.R` (159 keys; `state_entries` literal at :3688) | H6 |
| R4 | `test-omic-summary-module.R` (49 keys; 12 `createWorkflowArgsFromConfig` hashes; the inline `omicSummaryW1Frozen` captures, which gain `integration/<omic>/` files; a stale S4 digest that predates KIT-135) | H6; KIT-128; P4 |
| R5 | `test-workflow-figure-paths.R:747-751` (Summary archive limit `c(before = 64, after = 125)`, becomes about 77/138) | H6 |
| R6 | KIT-128's list: MCI-027.2 (`test-module-ci-trx-design.R`), the KIT-073 dump (`test-workflow-omics-declarations.R`), MCI-028.2 (`test-module-ci-trx-qc.R`), `test-omics-summary-core.R` (12) | KIT-128 |
| R7 | Frozen digests per KIT-129 s4: `test-omic-norm-module.R`, `test-prot-norm-core.R`, `test-workflow-saved-da.R`, `test-trx-norm-module.R`, `test-omics-da-core.R`, `test-trx-kernels.R`, `test-trx-publication.R`; fixtures `kit073-omics-declarations-frozen.rds`, `kit094-characterisation.rds`; the integration-input tests | KIT-126/127/128, G1 |
| R8 | Module CI: `tests/testdata/mci/refs`, `helper-module-ci-oracles.R`, the lipid/metab/proteomics/trx module-CI files; the nine trx digests on the laptop (value-bearing ones re-frozen by the causes above) | KIT-128, G1, H4 |
| R9 | e2e lanes (`tests/testdata/e2e`, seven lanes) and the report-template doubles (`tests/testdata/e2e/report_templates`) | P1, P2, P3, G1, G2 |
| R10 | Kit journey records and in-process oracles (protocol 7.1): DIA_LIMPA (m, r), lipid and metab (a, b, s, v), trx (c), LFQ with RUV (f); LFQ and TMT (l, t) verified, not re-recorded; DIA full (q) store against store | KIT-126/127/128, G1, H6 (store paths only) |
| R11 | KIT-040 notebook references (`docs/qa/workflow-execution/kit/kit-040/`): no longer describe trx_screen after KIT-128; kept as the Legacy check on a corpus without technical replicates | KIT-128; D3 (U6) |

## 6. Read and not added (no ticket defers them to this revision)

| Item | Source | Why not added |
|---|---|---|
| The GUI-entry recorder's mask (`helper-gui-entry.R:60`, `[0-9a-f]{10,}`) hides a run of 10 or more digits, so table digests do not see digits past the ninth | KIT-134 outcome :699, :705 ("a weakness of the recorder, not changed") | Not deferred. Changing it would move every entry oracle; the user decides (D5). |
| About 30 stale module names in the lipid/metab/limpa step maps | KIT-117 outcome ("for a follow-up") | Not deferred to this revision; documentation that no code reads (D6) |
| The dead `reclaim_memory` launch option | KIT-117 outcome | Not deferred |
| An independent check of the default transcriptomics route (a direct `RUVSeq::RUVs` + edgeR call on the journey's inputs) | v1.11 item 5 | The protocol text is in the draft. The script is not written: it waits for U6 (still open) and the corpus variant without technical replicates (D3). |
| `.dia-native-import-`, `.dia-native-legacy-`, `.dia-gui-adoption-` stay long | KIT-052 lists only the four names | A session adopted from a legacy or imported one gets less of the saving (D4) |

## 7. Counts

| Kind | Count |
|---|---|
| H, harness-bound (in the change set: H1 to H6; on the laptop: H7 to H9) | 9 (6 committed) |
| P, package-side before the re-record | 4 |
| G, gates | 2 |
| R, stored-reference families | 11 |
| D, read and not added | 5 |
