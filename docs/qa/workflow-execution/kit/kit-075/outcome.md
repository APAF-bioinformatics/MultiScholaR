# KIT-075 outcome: the lipid and metab tabs hold references

Status 2026-11-02: **completed**, with one recorded divergence (the import's file-selection preview) and a memory result that met the stop rule but not the expectation. Stage 3b, ticket 4 of 10.

- **Commits** on `performance-v2`: code f4791128 and 83b11a5a, help pages 8a445a22, bindings 9c138c58, and this evidence commit.
- **Base:** a0001f4a (KIT-074 closed).
- **Install:** `kit_075` = `core/browser-smoke.VuMtXT` (83b11a5a), over the kit_074 library. `browser-smoke.4lGlmC` (f4791128) was an intermediate install, used only for the R 4.6.2 control. The help-page commit changed no code, so the install stands.
- **Ticket:** `tickets/kit/KIT-075.md`, revision 1.1. Plan: [plan.md](plan.md).

## What changed

1. **Owner lookup (item 1).** `getWorkflowOwnedInput`, `getWorkflowOwnedPayload`, `setWorkflowOwnedPayload` and `savedReactiveVal` take an explicit `manager`; NULL keeps the session lookup, so proteomics callers are unchanged. The receipt probe also counts a list of tables.
2. **Readers in the core.** `readWorkflowOmicsState` reads a state's components for one action. `readWorkflowOmicsMetadata` returns the slots, assay names and row counts (from the stored tables' shapes) and the design, reading no assay table.
3. **Wrappers (item 2).** `hydrate` now defaults to FALSE; the argument stays because the KIT-031/032 recorders pass it. The DA wrapper no longer attaches `theObject`.
4. **Holders (item 3).**
   - `normData`'s slots hold state references.
   - The DA session carries the head's reference and slots (`state_metadata`), as proteomics does.
   - `lipidDaCurrentS4` / `metabDaCurrentS4` return metadata.
   - The ITSD tables' renders keep only their candidate tables.
   - Frames whose outputs keep them drop the object first (QC finalise, ITSD analysis, summary).
5. **Metadata and per-action reads (items 4 and 5).** The QC tab, auto-initialisation, DA status and preflight class checks, assay names, and the summary's `@args`, counts and samples use metadata. QC plots, duplicate detection, ITSD analysis, the preflight's alignment check, and the Glimma volcano and heatmaps read per action.
6. **Revert (item 6).** The GUI callers use `selectWorkflowStateReference`; `revertToState` is unchanged. The reset recommits the pre-normalisation state by reference.
7. **Import and design tables (item 3, user decision 2026-11-02).**
   - **The owner binds at the import** (`bindWorkflowOmicsOwnerAtImport`), the proteomics pattern. Binding replaces the pending store, so anything saved before it would be lost. The design step's bind is now a no-op.
   - `data_tbl` and `data_cln` on the bus are references.
   - The design builders read the tables per use (a function, not a caching reactive). They hold `data_cln`, `initial_state`'s tables and the returned result by reference (`omicsBuilderTableVal`, `state_owner`).
   - `design_matrix` stays a plain sample table, as in proteomics.

## Results

| Gate | Evidence | Result |
|---|---|---|
| Residency walk, fail first (S4) | [residency-before.txt](residency-before.txt) | On a0001f4a: 19 failures and 2 errors. The errors are the owner-lookup tests. |
| Residency walk, fail first (tables) | [residency-tables-before.txt](residency-tables-before.txt) | On f4791128: 44 failures. Tables found on the bus and in the builder's `data_cln`, `initial_state` and `result_rv`. |
| Residency walk, final | `tests/testthat/test-workflow-omics-residency.R` | **82 pass, 0 fail**: both omics, every step, with no assay object, no plain assay table, and the bus holding references. A unit test is the walk's positive control. |
| Exactness set (1)(3)(6) | [exactness-gates.txt](exactness-gates.txt) | On R 4.6.1: `lipid75` and `metab75` valid and **identical** to `lipid31` and `metab32`; KIT-028 replay **exact** (16 states); LFQ, TMT and DIA_LIMPA **identical** to the kit_065 oracles. |
| Exactness set (2) | [075l](kit075-kit-journey-075l-stage.json), [075m](kit075-kit-journey-075m-stage.json) | Both accepted and exact, 0 oracle mismatches. `075l`: 241 s, 2.53 GiB peak. `075m`: 241 s, 2.56 GiB peak. No sampler gaps. |
| Exactness set (4) | [075l](member-record-replay-075l.json), [075m](member-record-replay-075m.json) | Both replays exact (10 states); a tampered input is refused. |
| Exactness set (5) | [lipid](kit075-kit-lifecycle-lipid-stage.json), [metab](kit075-kit-lifecycle-metab-stage.json) | Both observed, launched one at a time. |
| Full suite | [full-suite.txt](full-suite.txt) | 355 files, one process per file: **30,127 pass, 0 fail, 0 error**. Two files needed attention; see Deviations. |
| R CMD check | [check-item1.log](check-item1.log) | Baseline: 3 WARNINGs, 3 NOTEs, code note identical to KIT-074's. |

## The declared memory observable

Anonymous peak in GiB, standard (lifecycle) runs:

| Omic | Build | Whole slice | Normalisation | Correlation | DA | Report |
|---|---|---|---|---|---|---|
| Lipid | kit_032 | 2.46 | 2.46 | 2.35 | 1.85 | 1.89 |
| Lipid | kit_074 | 2.48 | 2.48 | 2.34 | 1.86 | 1.90 |
| Lipid | **kit_075** | **2.44** | **2.44** | 2.33 | 1.80 | 1.83 |
| Metab | kit_032 | 2.54 | 2.54 | 2.33 | 1.91 | 1.94 |
| Metab | kit_074 | 2.54 | 2.54 | 2.33 | 1.91 | 1.95 |
| Metab | **kit_075** | **2.50** | **2.50** | 2.33 | 1.93 | 1.97 |

- **Stated outcome, as declared:** total peak is at or below `kit_032` on both omics (lipid 2.50 vs 2.53 GiB, metab 2.57 vs 2.61), so the stop rule did not trigger. This is reported and handed to KIT-080.
- **The expectation was not met.** "Near or below the in-memory build" (2.29 / 2.25 GiB) did not happen; the peaks are about 0.2–0.3 GiB above it.
- **The drop is not distinguishable from noise.** At 0.02–0.04 GiB it is smaller than the run-to-run variation already seen (the lipid peak has ranged from 2.46 to 2.68 GiB). No memory benefit is claimed from these screen runs.
- **Reading, not proven:** the peak sits in normalisation, while the kernels run in the app process with the object they need. References can't remove that; moving the kernels out (KIT-078) is the mechanism that can. Removing the resident holders is what this ticket verifiably achieved, by the walk.

## Deviations and discoveries

- **R changed under the series.** A system upgrade on 2026-11-02 at 11:59 moved R 4.6.1 → 4.6.2 (and onetbb, among 114 packages). Under 4.6.2 every stored reference stopped matching, on unchanged kit_074 code too: saved-product hashes everywhere, and DIA_LIMPA's state summaries from the design step on. kit_075 was identical to kit_074 on all five records there ([environment-control.txt](environment-control.txt); the failing KIT-028 replay is kept as `replay-check-kit075-r462.json`). The user returned R to 4.6.1 and pinned it with `IgnorePkg`; every result above is on 4.6.1, and the DIA_LIMPA match returned with it. **Whether the 4.6.2 DIA_LIMPA difference is in values or only in storage was not established; examine it before R is next upgraded.**
- **Scoped divergence: the import's file-selection preview stays plain.** `localData$assay1_data` and `assay1_import_result` hold the first assay until the next file selection. About 20 render and selector sites read them, and they exist before the owner can be bound. The walk does not exercise them, because it drives the import through its processing helper; so "no plain table reachable" does not cover this holder. Worth revisiting with KIT-078.
- **Bind at import, not at design.** This changes KIT-070's binding point (user decision 2026-11-02, "do it in KIT-075"). I first described this option as having the larger blast radius; on reading the code it was the smaller one. The store directory now appears at import.
- **Gate plan agreed with the user.** Journeys and lifecycle ran once, on the final build. After f4791128 the diff touches omics files only, plus the help pages. Proteomics was still run once on the final build, since its only earlier runs of this code were on R 4.6.2.
- **The ticket's "state names"** in `normData` became state references, the proteomics `correlation_filtered_obj` pattern.
- **`design_matrix` is not in the receipt probe.** Proteomics keeps it plain on purpose, and the proof tool fails on any resident assay holder, so adding it would break proteomics.
- **Kept on purpose:** the session export RDS and the integration RDS still use the owner's full read. They are per action, and the component reader's attributes would otherwise enter the files.
- **Test stand-ins.** 54 seam tests with a stand-in manager list use `omicsStandInReaders()`; wrapper mocks use `omicsMockCommit()`. The readers themselves are exercised on real owners by the declaration tests, the lanes and the walk. A non-S4 state in the finalise lanes is now refused by the component reader before the class check; the assertion accepts either message.
- **Displays.** The walk does not render the lazily drawn Glimma and heatmap outputs. The declaration tests check that their source read equals the state's data and leaves the results otherwise unchanged, and the exports' heatmaps ran in the journeys. Rendered images were not compared.
- **Reset not in the journeys.** The recommit by reference is covered by the module tests only.
- **Note for KIT-078.** `revertToState` cancels a managed worker action and `selectWorkflowStateReference` does not. The omics have no worker yet; KIT-078 must add the cancel.
- **Full-suite first run.** `test-omics-audit-carry-overs.R` had 4 failures (a stand-in manager, fixed). `test-general-filemgmt-report-helpers-shared.R` gave no result line under the parallel run and passes alone in 9 s; the cause was not established.
- **My error, fixed:** f4791128 committed 123 testthat failure snapshots by accident. 83b11a5a removes them and ignores the folder.
- **Pre-existing, left alone:** the `copyToResultsSummary` copy warning (KIT-077).

## Backward compatibility and rollback

- Saved sessions open (`test-workflow-omics-sessions.R`, the journeys' DA load); legacy session files still adopt.
- `revertToState` is unchanged; proteomics holders are unchanged.
- No saved format changed. Rollback: `git revert` the ticket's commits.

## Handoff

KIT-076 gets a DA tab that holds no object. KIT-077's figure sites read per action. KIT-078 inherits `readWorkflowOmicsState` for its worker reads, the cancel note above, and the memory question.
