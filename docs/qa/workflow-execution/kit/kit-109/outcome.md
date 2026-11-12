# KIT-109: removals and closure of the GUI workover series (2026-11-12)

**Outcome:** done. Every work item is done, with no numeric result changed, and every exit run passes. The one regression found at the exit (W6, check mode) is corrected. All counts below are `load_all`, one process per file, unless a row says otherwise.

## Work items

| W | Item | Result | Commits | Evidence |
|---|---|---|---|---|
| W0 | Sites | Every site confirmed by reading at `35769a54` ([sites.md](sites.md)). The grep missed multi-line calls and stubs, which the tests then found; `sites.md` was corrected. | `1f1641c3` | |
| (entry) | QC digests | Two DIA_LIMPA QC UI digests changed by KIT-116's P-02 label (`08e9b25e`), a file KIT-116 did not run. They were re-frozen with that attribution (the file at `08e9b25e~1` is 158/0, at `08e9b25e` 155/3). | `347e5dd5` | `test-omic-qc-module.R` 157/0 |
| W1 | `profile` removed | Transcriptomics' registry entries mount the shared `mod_omic_*` modules under `omicGuiConfig("transcriptomics")` and its `norm`, `da` and `summary` groups. The eight lipid wrappers lose `profile` and `trxOmicProfile()` is deleted. `prepareLipidSummarySessionState` takes the state reader, workflow type and template as arguments. Eight Rd pages were edited by hand (no `document()`). The tests' comparisons with the profile's values are now those values, frozen at `35769a54`. | `f02e5183` | trx entry lane 29/0; module-CI trx files, `test-omic-*-module.R` and `test-omic-gui-config.R` at their entry counts, with no capture change. In `R/`, the grep finds only the DIA rollup's own `profile` argument (a different one, kept) and one dated provenance comment (`utils_omic_gui_config.R:1959`, "at 5be34ee5", `@noRd`). In `tests/`, it finds only absence assertions and history comments. |
| W2 | Dead functions | `trxSummaryIntegrationReady` and `buildLipidSummarySessionState` are deleted, and their seam tests are already Q-170/Q-171. The four test-only validators are kept: each is a one-line per-omic entry point whose siblings are live ([sites.md](sites.md)). | `a204b763` | `test-omics-summary-core.R` 701/0/2 |
| W3 | A106-05 | Proteomics' DA group drops `status_keys` and `status_text`. Five fields stay because the validator requires them: `columns`, `summary`, `table_filter`, `lfc$tested`, `views$heatmap$kernel`. | `d498b1b8` | `test-omic-gui-config.R` 1852/0 |
| W4 | P-25 | Proteomics' empty `heatmap_manual_save_warning` slot is deleted. The proteomics DA HTML digest changed by that node only. The proteomics entry oracle did not change: a freeze changed only its metadata, so it was reverted. Lipid, metab and trx render the note and are untouched. **The note's wording is wrong for lipid (audit, P-25) and stays as it is.** | `eae95778` | `test-omic-da-module.R` 1369/0 |
| W5 | A101-08 | One method-to-normaliser lookup (`switch`) is now used by both the protein table and the DPC EList in `func_prot_s4_norm_methods.R`. | `23fda3fc` | `test-prot-norm-core.R` 1690/0, identical; eleven dependents at their entry counts |
| W6 | Guards and CI | `impact-map.json` gains a `shared-gui-layer` rule (`mod_omic_*`, `func_omics_*`, `utils_omic_gui_config.R`) and a `transcriptomics-source` rule. The four broad rules and `detect-impact.R`'s fallback now include transcriptomics. `run-module-ci.R` maps the trx scenarios to their files (none were mapped before). The globals guard now scans the shared files: a per-family check replaces the stale `> 100` count, and the `<<-` allowlist is updated. **Exit correction:** in check mode the per-family check failed six times where the old count failed once, so it was made one expectation. | `4c56d2f0`, `a15d64c0` | `test-ci-impact-*` pass (new cases failed first); module CI `--omic transcriptomics` 10 scenarios pass; globals guard 14/0 |
| W7 | Records | Two stale comments are fixed. Stale statuses in `kit-101/findings.md` are corrected, append-only: F-06/F-07 `88661446`, F-13 `958273df`, F-14 `60cc2605`, F-22 `6314065a`/`724c10be`. A101-07 is answered and the convention divergences are recorded (both below). | `63ce6d6c`, `1f1641c3` | |
| | Bindings | kit_922 = `browser-smoke.iZxvdt` at `63ce6d6c`; members 922a, b, c, l, t from 112a/b/c/l/t. | `e7559e28`, `2033ed7e` | |

**Moved out at W0:** the diagnostic exit boundary is in a `tools/profiling` harness file, so it goes to KIT-117.

## A101-07, answered by reading

- **Q1, DIA_LIMPA QC plots:** no defect. The plots read `protein_quant_table`, which equals E except in RUV covariates mode, where the difference is by design. Side note: the pre- and post-normalisation plots on DIA_LIMPA both show the already-normalised state. This is a labelling issue; it is recorded here and not changed.
- **Q2, a literal 0 in `Precursor.Normalised` on DIA (IQ):** no defect, because iq drops zero cells.
- **For the user:** IQ quantifies from `Precursor.Quantity` (`Peptide.RawQuantity` → `Peptide.Imputed`) with `normalization = "none"` (KIT-025's fixed parameters, `utils_workflow_prot_declarations.R:68`, `:84-86`, `:419-420`). A literal 0 in `Precursor.Quantity` therefore enters the imputation and pulls its mean down. Whether 0 should be treated as missing there is a methodology question. It belongs with KIT-114 (e), zeros as missing.

## Convention divergences, recorded and not refactored

- **A104-12:** `mod_omic_qc_tabs.R:601`, `:644`, `:647` find sample columns with `is.numeric` rather than from the design.
- **A111-07:** the norm dispatchers repeat the in-process gate and its resolution, with 9 to 12 parameters each.
- **A112-08:** the DA module detects the session class by name (including `MockShinySession`), and one observer is over 50 lines.
- **A113-10:** `runOmicEnrichAnalysis` is about 105 lines.

Each one works and is tested, and refactoring it is outside a removal ticket.

## Closure of the GUI workover series

27 of the 28 omic × stage cells are on the shared frame. **The one exception is proteomics' normalisation tab (P-32; user, 2026-11-11: keep the exception).**

## Exit runs

- **Full suite, both modes, at `63ce6d6c`** ([exit/suite.md](exit/suite.md)): `load_all` 36,582 passed, 0 failed, 2 error blocks (the two diagnostic ones), 0 with no result. Check mode: 31,465 passed and 154 failed. The 154 is the entry's 152, plus five from the W6 regression (corrected), minus the three re-frozen QC failures. There is no new file error, new file skip or file with fewer blocks run.
- **`R CMD check --no-tests`:** 3 WARNINGs, 3 NOTEs (baseline). **Install:** kit_922.
- **e2e lanes** ([exit/e2e.txt](exit/e2e.txt)): seven of seven pass on kit_922, from a worktree at `63ce6d6c`, one attempt each. Lipid 88/0 (130 s), metab LC-GC 190/0 (283 s), trx screen 16/0 (89 s), proteomics DIA 42/0 (290 s), DIA_LIMPA 86/0 (216 s), enrichment 98/0 (583 s), TMT/LFQ 153/0 (292 s). The counts equal KIT-116's.
- **Module CI `--omic all`** ([exit/module-ci-all.log](exit/module-ci-all.log)): passed, 110 scenarios over 34 test files, the seven `module-ci-trx-*` files among them. The skips are the quarantined and guarded ones already in the register.
- **Journeys** ([exit/journeys/](exit/journeys/)): five members on kit_922, one attempt each, all accepted under `exactness_and_completion`: 922a lipidomics, 922b metabolomics, 922c transcriptomics, 922l protein LFQ, 922t protein TMT. Each is exact against its oracle.

## Cleanup

All eight `audit-r*` branches had their commits on `performance-v2` (`git cherry`). For `audit-r5`'s `f66db461` the patch id differs, but its changed lines are identical to `a5704be3`. The branches and their six worktrees were removed. To restore, the tips were: r1 `d85e1461`, r1b `0b316584`, r2 `edadb433`, r2b `29aeb4c1`, r3 `ccf11ddf`, r4 `79e00d7d`, r5 `b0bbc0b4`, r9 `8f2bef6c`.

## Acceptance (revision 1.3)

- [x] Every work item is done (W2's four validators were kept on W0's reading).
- [x] The W1 and W2 static greps are empty, apart from the hits listed under W1.
- [x] The entry-point suite is unchanged: in the exit suite (`load_all`), lipid is 35/0, metab 105/0, proteomics 193/0 and trx 29/0, which covers all ten lanes against their oracles. No oracle changed: W4's freeze changed only metadata and was reverted.
- [x] Full suite in both modes: no file has fewer tests run, a new failure or a new file-level skip. The check-mode globals regression is corrected (`a15d64c0`).
- [x] Journeys exact; e2e lanes pass; module CI passes; `R CMD check` at baseline; install.
- [x] The closure names P-32 as the one exception.
- [x] The `audit-r*` branches and worktrees are removed.

## Next

KIT-117 (harness revision), with KIT-114 alongside, then KIT-118 and KIT-119. Each needs the user's go.

## Correction, 2026-11-12

"P-32; user, 2026-11-11: keep the exception" overstates the record. The user did not choose P-32 specifically: "keep the exception" was the audit table's recommendation, adopted with the rest of the table ("i want everything fixed according to your recommendations"). On 2026-11-12 the user asked for the deep join to be planned: KIT-120. The closure above is the state at this ticket's date.
