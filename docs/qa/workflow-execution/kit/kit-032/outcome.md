# KIT-032 outcome: metabolomics declarations and the stage-3 closure

Status 2026-10-31: **completed; stage 3 closed**. Code 4819f030 (item 1) + d0b99a45 (bindings, the record on the install) + 8e5897ee (the exact member and its replay) + d3801500 (the lifecycle members) on `performance-v2`; base e4572f01 (KIT-031 closed). Install `kit_032` = `core/browser-smoke.Tmn2Ur` (4819f030) over the kit_031 library. No counted member. Plan: [plan.md](plan.md) (the planning gate, revision 2.1).

## Stage-3 closure statement

Lipidomics and metabolomics are on the kit. Both modules hold the on-disk session owner from session start (KIT-070), write their session files as the owner's publications (KIT-071), and run every step as a registered in-process declaration with complete provenance (KIT-031, KIT-032). For each omic, one screen member through the GUI to the report is **exact** against an in-process record made on the same build, the member's analysis record validates with no legacy step, and the record replays exact from the design state. The proteomics routes are unchanged throughout the stage: the KIT-028 replay is exact and the LFQ, TMT and DIA_LIMPA records are identical to the kit_065 oracles at the closure build. The kit's rule held: the two new routes needed declarations, registry entries, record and replay branches and module code only.

**Claim limits.** No memory claim is made for these omics: the stage ran lifecycle members only, on synthetic screen corpora (24 samples, two assays of about 200 features); the peaks below are those of the application on these corpora, not a prediction for real data. Enrichment is not on either route (no GUI exists; KIT-066). The report templates' shared logic is KIT-067's. Complete-journey qualification and the summary rule stay with KIT-050. The `ceiling()` removal is in both intensity filters and changed no product on any lane or corpus recorded here.

## What changed (item 1)

The metabolomics twin of KIT-031's item 1 ([KIT-031 outcome](../kit-031/outcome.md) describes the pattern):

| Where | Change |
|---|---|
| `R/utils_workflow_metab_declarations.R` (new) | twelve declarations on `MetaboliteAssayData` (`metabolite_data`, `design_matrix`): `metab_intensity_filter` (row; the per-assay quantile pre-pass; called by the method's own name, as the module calls it, so the `@args` key is `metaboliteIntensityFiltering`), `metab_duplicate_resolution` (feature group, through `resolveMetabDuplicateAssayData`), `metab_qc_finalize` / `metab_norm_complete` / `metab_reset` (store), `metab_itsd_normalisation` (complete sample; the mean-factor pre-pass), `metab_log2` (row), `metab_between_sample_normalisation`, `metab_ruv` (the module's optimisation and correction steps in one), `metab_pearson_pairs`, `metab_pearson_selection`, `metab_da_output` (product only; the manifest `metab.da_results.v1`); the wrappers `runMetabReference*`, `recordMetabDesignProvenance`, `runMetabReferenceDaAnalysis`, `replayMetabolomicsStep` |
| `R/utils_workflow_lipid_declarations.R` | the omic-neutral wrapper core named `runWorkflowOmicsReferenceOperation`, `commitWorkflowOmicsStoreStep`, `markWorkflowOmicsStoreStep`; the lipid names kept as aliases (one body each) |
| registry, outputs, record, overlay, harness | the metab loop (`omic_type = "metabolomics"`), defaults from the shipped `[metaboliteIntensityFiltering]` keys, labels, `workflowSectionOperationsByType$METABOLOMICS`; `metab.da_results.v1`; the metab RUV/ITSD outcomes, the replay branch and DA re-run; the overlay's metab sources; `kitOperationIdentityProducts` covers `metab_*_science` |
| metab module sites | the intensity observer, the duplicate workflow (the apply step publishes without a second commit on the owner), the QC finalise, the ITSD, log2 and between-sample steps, the RUV shell, the correlation dispatcher (the outcome told the state is committed), skip and reset, the DA shell, the design builder and design import; a test double keeps the direct path |
| summary and template | `analysis_record.json` beside `study_parameters.txt` at Save Workflow Args; the record bundle at the report; `metabolomics_report.rmd` overlays the record |
| `metaboliteIntensityFiltering` | `quantile(...)[1]`; the `ceiling()` removed |

Divergence recorded: the metab normalisation module splits RUV into an optimisation step and a correction step; on the owner both run as the one declaration `metab_ruv` (one commit, as before), and the QC step runs on its result. Tests: new `test-workflow-metab-declarations.R` (223, the KIT-031 test's shape on a metab fixture), the registry lists, `test-metab-02t` (the pinned threshold source) and `test-metab-03b` (the outcome mock) updated to the contract; [affected-tests.txt](affected-tests.txt).

## Results

| Gate | Evidence | Result |
|---|---|---|
| Item 1: metab module-CI | [module-ci-comparison.txt](module-ci-comparison.txt) | design, QC, DA lanes byte-identical; norm and summary differ only in owned-export journals and date/version lines; no data product changed |
| Item 1: R CMD check | [check-4819f030.log](check-4819f030.log) | baseline (3 WARNINGs, 3 NOTEs; the code note identical to KIT-031's) |
| The recorder against the GUI | [recorder-vs-kit031-member.txt](recorder-vs-kit031-member.txt) | equal to the kit_031 metab member state by state (9 journey and 3 intermediate states) at the first run; the member ran with the `ceiling()`, so its removal changes no row on the metab screen corpus |
| **Acceptance: the metab screen member exact** | [kit032-journey-032m-stage.json](kit032-journey-032m-stage.json) | **accepted, exact**: 0 oracle mismatches, 0 canonical-form differences, 10/10 steps; 241 s, 2.60 GiB peak |
| **Acceptance: its record and replay** | [member-record-replay-032m.json](member-record-replay-032m.json) | record valid, 10 steps, 0 legacy; replay from the design state **exact** (10 states, the DA re-run); a tampered input refused |
| Lifecycle members on kit_032 | [kit032-lifecycle-lipid-stage.json](kit032-lifecycle-lipid-stage.json), [kit032-lifecycle-metab-stage.json](kit032-lifecycle-metab-stage.json); [verify-member-kit032-lipid.txt](verify-member-kit032-lipid.txt), [verify-member-kit032-metab.txt](verify-member-kit032-metab.txt) | both observed at the first attempt (lipid 244 s / 2.53 GiB, metab 243 s / 2.61 GiB); every commit through the GUI names its declaration for both omics |
| Projects against the kit_031 members' | [member-project-comparison-kit031-vs-kit032.txt](member-project-comparison-kit031-vs-kit032.txt) | both omics: every DA table, publication workbook and `study_parameters.txt` identical; the metab project gains the record bundle; the lipid store's entities name the new build |
| **Closure: both stage-3 replays** | [member-record-replay-031l-on-kit032.json](member-record-replay-031l-on-kit032.json), [member-record-replay-032m.json](member-record-replay-032m.json) | the KIT-031 lipid member's record replayed under kit_032 **exact** (10 states); the KIT-032 metab member's **exact** (10 states) |
| **Closure: proteomics unchanged** | [replay-check-kit032-closure.json](replay-check-kit032-closure.json), [closure-gates.txt](closure-gates.txt) | KIT-028 replay **exact, 16 states, 265 s**; LFQ (11), TMT (11) and DIA_LIMPA (14) records on kit_032 **identical** to the kit_065 oracles |
| **Closure: full suite** | [suite-summary.txt](suite-summary.txt) | **30,025 pass, 0 fail, 2 errors (the two known diagnostic files), 307 skips, 352 files**, one process per file on a clean worktree of d3801500: green by the rule recorded since KIT-017. Two single-process attempts (cross-file leakage, a runner mismatch) are kept and not counted |

## Memory and disk on the screen corpora (reported, not claimed)

| Member | Route | Build | Corpus input | Peak (GiB) | Resident at end (GiB) | Store on disk |
|---|---|---|---|---|---|---|
| kit-journey-031l | lipidomics, declared | kit_031 | 112 KB | 2.54 | 1.87 | 22 MB, 79 entities |
| kit-journey-032m | metabolomics, declared | kit_032 | 148 KB | 2.60 | 1.91 | 21 MB, 77 entities |
| lifecycle (kit_068) | lipidomics, in memory | kit_065 | 112 KB | 2.29 | 1.60 | none |
| lifecycle (kit_068) | metabolomics, in memory | kit_065 | 148 KB | 2.25 | 1.68 | none |
| kit-journey-065m | DIA_LIMPA screen | kit_065 | 28 MB | 2.12 | 1.28 | 187 MB |
| kit-journey-065 | DIA full (counted) | kit_065 | 300 samples | 7.36 | 2.85 | 4.8 GB |

The lipid and metab peaks sit in the normalisation action (one tab action: ITSD, log2, between-sample normalisation, the automatic RUV sweep and the QC figures), about 1.6 GiB above the 0.85 GiB the application holds through import and QC. The in-memory modules peaked about 0.25 GiB lower on the same corpora: the owner's modules hydrate whole objects per read (KIT-070's recorded cost), which a component-streaming read path would remove; that is not in the stage-3 scope.

## Backward compatibility and rollback

Public API, S4, checkpoint and saved formats unchanged; the proteomics routes identical; the lipid names of the wrapper core kept. Rollback: `git revert` of d3801500, 8e5897ee, d0b99a45 and 4819f030 (and KIT-031's commits for the lipid route).

## Handoff

KIT-066 (real lipid and metab corpora, enrichment decision), KIT-067 (the templates' shared logic), KIT-050 (complete-journey qualification and the summary rule), and stage 5 (integration: the declared DA manifests `lipid.da_results.v1` / `metab.da_results.v1` as the inputs' source).
