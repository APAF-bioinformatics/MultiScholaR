# Review of the KIT series after the DIA stage, and the stage-2 re-cut (2026-10-19)

User request: pause after KIT-033 and evaluate the planned remainder of the series against what the DIA refactor actually produced, name the generalised functions that worked, and say where the codebase falls short of best practice. Read-only on the code; the outcome is a re-cut of stage 2 (tickets in `tickets/kit`, gitignored; this file and the plan rows are the tracked record).

## What stage 1 produced

Proven on the 300-sample DIA journey and the kit-journey-028 project, with no omic name in their bodies:

- the class declaration registry (`workflowClassDeclaration`; `PeptideQuantitativeData`, `ProteinQuantitativeData`, `LipidomicsAssayData`, `MetaboliteAssayData`, the result classes);
- the operation runner and registry (`runSavedOperation`, `runProjected`, `registerWorkflowOperation`; 16 DIA data steps registered, one place names an operation);
- the three partitioners (`partitionByRows`, `partitionBySamples`, `partitionByFeatureGroup`; cap-driven batch counts);
- the managed worker (request, poll, publish, acknowledge; 13 module call sites of `installWorkflowManagedStep`);
- the harness (journey driver, per-step populations and oracles, protocol v1.4, the ledger discipline, the per-ticket test rule);
- the analysis record (head provenance filled by every registered operation; derive, validate, replay 16/16, render).

Not generalised:

| Finding | Evidence (2026-10-19) |
|---|---|
| The storage and owner layer is DIA-named | 204 functions named `*DIA*` in 18 files; the generic kit files total about 780 lines against about 7,400 in `utils_workflow_dia_{state,products,da,qc_rollup,sessions}.R`; stage 3's first call would be `newDIAStateAdapter` and `getDIAActiveManager` |
| The GUI runs two code paths | 55 `usesDIAReferenceStates(` sites in 29 `R/mod_prot_*.R` files, each with an in-memory `else`; `WorkflowState$new()` constructed in `mod_proteomics.R`, `mod_lipidomics.R`, `mod_metabolomics.R` and four design helpers; 11 direct `$saveState(`, 11 `$getState(`, 19 `$getHistory(` calls in proteomics modules; the design helper gates the kit path on `workflowType == "DIA"` |
| The owner attaches only for DIA-NN with iq on Linux | `initializeDIANativeImport`: `sysname == "Linux"`, `format == "diann"`, `method != "limpa"`; `adapter$setWorkflowType` refuses anything but `"DIA"`; `ingestDIANativeLegacySession` refuses non-DIA and Limpa sessions |
| Platform identity, not capability | shell pipelines for hashing (`bash`/`dd`/`sha256sum`), durability (`sync`) and archive (GNU `tar`); the worker through `systemctl --user` and cgroups; the P14 native rollup x86_64-Linux (with a `legacy` fallback) |
| Globals | 70 `.GlobalEnv` references in `R/`, 32 of them in the lipid and metabolomics kernels |
| `partitionByPlex` has no consumer | zero `plex` references in proteomics modules; a TMT plex is a column set of one wide table |
| The record misstates three steps | post-hoc audit `analysis-record-series/AUDIT-20261019-implementation.md`: IM-01 (fixed constants dropped, `@args` dumped for IQ, P14 and the design), IM-02 (realised RUV absent), IM-03 (open vocabulary), IM-06 (report prints `NA` versions), IM-07 (SI in kit vocabulary), IM-09 (no build or numerical environment) |

## Manuscript Figure 2 and Methods (held until the series is done; user 2026-10-19)

1. Panel B slot 2 and the "Operation runner" paragraph say the identity hashes the build. It does not (`diaComponentScientificIdentity`: operation, parameters, scientific slots, column roles); the build is recorded in provenance so a product survives a rebuild.
2. "Failure isolation" says each operation runs in a separate worker. Normalisation and RUV, correlation, the DA contrasts and enrichment run in process under `runProjected`.
3. Panel A, proteomics column 3, omits the protein QC steps, one of which (the missingness filter, dataset-wide quantile) is whole-matrix.
4. The metabolomics, lipidomics, transcriptomics and MOFA2 lanes are design, not verified units; the step map confirms DIA only. A caption sentence should say so. Limpa has no declaration until KIT-046.

## Decisions (user, 2026-10-19)

- Submission deadline 31 December 2026; Figure 2 and Methods fixes after the series.
- Stage 2 re-cut (this document); order: record remediation (KIT-038, 039), LFQ/TMT (KIT-020 rev 2.0, 042, 043, 021 rev 2.0), Limpa (KIT-046, 047), portability then the single path (KIT-044, 048), naming (KIT-045), then stage 3.
- `partitionByPlex` dropped.
- The adapter is ported to every platform (KIT-044) before the in-memory path is retired (KIT-048); macOS and Windows qualification stays with KIT-051/052.
- Limpa declared now (KIT-046); `APAF-bioinformatics/fast_limpa_dpc` (`fastlimpadpc`, Rust, aligned to about 1e-8 in `E` and 1e-4 in standard errors, not bit-identical) ingested as a declared, recorded backend (KIT-047) with a user decision on the acceptance rule for a tolerance-aligned alternative kernel before it closes.
- The naming ticket approved (KIT-045) under AGENTS.md's renaming rule as the recorded debt of the KIT-017 checklist.
- The protocol's replay member kind stays undecided (AR-19); nothing counted depends on it.

## Naming inventory for KIT-045 (defined in `R/`, none exported; counts are matching lines in R/ excluding the definition, tests/, tools/)

| Function | R | tests | tools |
|---|---|---|---|
| `getDIAActiveManager` | 102 | 18 | 0 |
| `withDIAStateScope` | 55 | 74 | 9 |
| `newDIAStateAdapter` / `closeDIAStateAdapter` | 10 / 16 | 107 / 230 | 16 / 37 |
| `usesDIAReferenceStates` | 88 | 11 | 0 |
| `getDIAStateReference` / `selectDIAStateReference` | 54 / 14 | 138 / 12 | 15 / 1 |
| `readDIAStateMetadata` / `resolveDIAStateValue` / `materializeDIATable` | 55 / 28 / 13 | 61 / 6 / 46 | 22 / 4 / 10 |
| `commitDIAReferenceState` / `publishDIAStateHead` / `setDIAStepProvenance` | 14 / 29 / 6 | 19 / 6 / 8 | 4 / 1 / 0 |
| `readDIAExecutionComponents` / `checkpointDIAExecutionResult` / `inheritDIAExecutionSource` | 23 / 6 / 8 | 22 / 2 / 0 | 6 / 0 / 0 |
| `saveDIAOwnedValue` / `saveDIATableComponent` | 14 / 16 | 6 / 13 | 2 / 6 |
| `saveDIADisplayTable` / `readDIADisplayTablePage` / `diaDisplayTableFilter` | 8 / 7 / 4 | 7 / 10 / 2 | 0 |
| `rootDIAProduct` / `diaComponentScientificIdentity` | 10 / 6 | 0 / 7 | 0 |

`DIAStateAdapter` is an S3 class on an environment (54 `inherits` sites in R/; five test files set the class string directly). `saveWorkflowTableComponent` already exists (`utils_workflow_saved_table_codec.R`), so `saveDIATableComponent` becomes `saveWorkflowOwnedTable`. The CORE-pinned files and `workflow_kit.R` call `newDIAStateAdapter`, `closeDIAStateAdapter`, `readDIAStateMetadata`, `publishDIAStateHead`, `materializeDIATable`, `getDIAStateReference`, `isDIAComponentCheckpoint`, `assertDIAExportReady`, `openDIANativeSession`, `readDIANativeLatestPath`, `resolveDIANativeSessionView`; every DIA name stays callable through an alias.

## Best-practice debt outside the re-cut (recorded, not scheduled)

- The `cli` import warning and the `.ruv_replicate_block` note in `R CMD check` are the package's own; the Makevars and hidden-file items belong to the vendored extension.
- Fixed interpretive prose in the report templates (IM-08), database identity and upstream search provenance (IM-04/05), method citations and `inst/CITATION` (IM-11): paper and release decisions.
- The open test quarantine entries Q-02, 03, 06, 07, 08, 10, 12, 13, 17, 18, 19 (`quarantined-tests.md`).
- `checkParamsObjectFunctionSimplify` (call-stack parameter resolution) in 62 lipid/metab kernel call sites: adapted to, never rewritten (README "Kernel contracts").

## Limpa as its own workflow (later on 2026-10-19)

User redirection after the re-cut: the DIA workflow must not use limpa; limpa is its own workflow, set up like DIA on the on-disk architecture, following the checkpoints of `APAF-bioinformatics/e_35036_MorpheusTex_GarryNunn_20260329_b/scripts/proteomics_isoseq/DIA_workflow_limpa_isoseq_tissue.Rmd` (its ev, juice and metagenomics twins differ only in project constants) and running on `fastlimpadpc`. The chunk-by-chunk map is `limpa-workflow/script-map.md`.

What limpa is in the package today: an import radio value the owner refuses (`initializeDIANativeImport` returns without attaching), a rollup-tab branch (`runProteinLimpaRollupApplyStep`) that aborts on an owner and otherwise runs in memory, a template choice from `use_limpa`; the GUI's DA path never reaches `dpcDE` (the branch lives in `daAnalysisWrapperFunction`, which has no caller); `generateLimpaQCPlots` and `peptideMissingValueImputationLimpa` have no callers in `R/`.

Decisions (user): `fastlimpadpc` default with limpa selectable and recorded; the DIA route keeps its replicate-mean imputation as an optional step and is advertised as "DIA, iq MaxLFQ (optional replicate-mean imputation)", the new route as "DIA with imputation leveraging limpa (Li et al. 2025)". Design decisions taken here: one proteomics module with a new type `DIA_LIMPA`; RUV as covariates (W appended to the design, the matrix not corrected, the adjusted matrix only for QC plots) as the script does; the script's capped optimiser enters the package as `findBestNegCtrlPercentageCapped` (no "Custom" in names); the between-sample normalisation step exists but defaults to `none` as the script runs; the package's `DIANN_limpa_report.rmd` is the route's template although the script rendered the plain one (deviation recorded); old in-memory limpa sessions are refused by name. The exactness question of the morning's KIT-047 dissolves: a new workflow has no prior oracle, so `fastlimpadpc` is the declared kernel, exactness is on-disk versus in-process on the same kernels, and determinism (two runs, two thread settings) and alignment to limpa are recorded. Tickets: KIT-049, KIT-046 rev 2.0, KIT-047 rev 2.0, KIT-041, KIT-053, between KIT-021 and KIT-044.

Steps the script does not run, which the new route therefore omits: peptide intensity filter, minimum peptides per sample, peptide replicate filter, peptide imputation, IQ, protein-level accession cleanup (done at the peptide level instead), protein intensity, duplicate and replicate filters, RUV-III-C correction of the matrix.

