# KIT-077 outcome: lipid and metab figures as store entities, and the heatmap repair

Status 2026-11-02: **completed**. Stage 3b, ticket 6 of 10. Every figure file on both journeys is a stored, rooted asset. The per-contrast heatmaps now exist, and two copy defects are repaired. No memory benefit is claimed; the costs are below.

- **Commits** on `performance-v2`: heatmap repair 224ff3b9, mechanism e9807ab7, bindings 5c46896f, evidence after them.
- **Base:** c9f2cac9 (KIT-076 closed).
- **Install:** `kit_077` = `core/browser-smoke.FXMKdB` (e9807ab7), over the kit_076 library. Later commits change tests and documents only, so the install stands.
- **Ticket:** `tickets/kit/KIT-077.md`, revision 1.2.

## What changed

1. **Heatmap repair (user decision 2026-11-02).** `generate<Omic>DAHeatmap` returns `list(plot, ...)` and the exports drew the list. The error was caught and logged, so a run looked fine while no per-contrast PNG or PDF was written, `all_heatmaps_combined.pdf` was 0 bytes and was copied to `Publication_figures` as such, and every failure left a graphics device open. Two test mocks returned a bare `Heatmap`, which hid it.
2. **The omics figure product** (`R/utils_workflow_omics_figures.R`, new). A hit writes the stored bytes; a miss renders, saves and roots. PCA, RLE, density and correlation keep their plot data; every other figure keeps bytes only. With no omics owner in reach, it renders to the file as before. The proteomics figure functions are untouched.
3. **Sites:** per-assay QC, filtering QC and grids, the composite, and the DA volcanoes, heatmaps, NumSig and combined PDFs.
4. **Copy defects.** The metab copy now looks for the name the module writes and the report reads. The lipid fallback, which could only queue a missing file, is deleted. The fixtures that seeded the wrong name are corrected.

## Results

| Gate | Evidence | Result |
|---|---|---|
| Heatmap, fail first | [heatmap-fail-first.txt](heatmap-fail-first.txt) | With the real generator: 6 of 7 expectations failed per omic before the repair (no files, a 0-byte PDF, a leaked device); 0 after. |
| Heatmap, real run | `077l`, `077m` | 4 PNG and 4 PDF new per member (counting the `Publication_figures` copies), the combined PDF non-empty. One PNG inspected: KO and WT separate, replicates mostly pair. |
| Copy defects, fail first | [copy-defects-fail-first.txt](copy-defects-fail-first.txt) | Metab 2 failures, lipid 1, before the repair; 0 after. |
| Metab composite, real run | `077m` vs `076m` | In `QC_figures` now, absent before. The metab report grew from 11,463 to 334,514 bytes; the lipid report's size is unchanged. |
| Every file from the store | [077l](figure-compare-077l.txt), [077m](figure-compare-077m.txt) | **69 of 69 and 67 of 67 figure files are byte-identical to a stored asset; 37 assets each, all rooted.** Every baseline path present; nothing empty. |
| Plot data | same, [pca-difference.txt](pca-difference.txt) | RLE, density and correlation (36 products) equal the baseline. All 12 PCA products differ from my baseline, which passed `shape_variable = NULL`; with the app's parameters, all 12 are identical to `plotPca` on the old build. |
| Reopen | `test-workflow-omics-figures.R` | A reopened owner writes every QC file byte-identical with the plotting functions refusing to run. A second DA export writes the stored bytes. |
| Exactness | [exactness-gates.txt](exactness-gates.txt), [077l](kit077-kit-journey-077l-stage.json), [077m](kit077-kit-journey-077m-stage.json), [replay 077l](member-record-replay-077l.json), [replay 077m](member-record-replay-077m.json) | 077l and 077m accepted exact against `lipid77` / `metab77`; records valid and replay exact; records identical to KIT-076's. KIT-028 replay exact (16 states); LFQ, TMT, DIA_LIMPA identical to kit_065. Both lifecycle members observed, launched one at a time. |
| Tests and check | [affected-tests.txt](affected-tests.txt), [check-item1.log](check-item1.log) | 42 dependent files: 6,235 pass, 0 fail, run before the last test was added (the figure file: 74 pass). R CMD check at baseline, code note unchanged. No full suite: the next checkpoint is KIT-078. |

**Determinism** ([figure-determinism.txt](figure-determinism.txt)): PNG is byte-identical across runs when the plot has no randomness; PDF is not (timestamp). KIT-076 showed volcano PNGs vary too. So equality is judged on plot data and file presence, as the ticket allows.

## Costs (reported)

| Lifecycle, GiB anonymous | Whole slice | Normalisation | DA | DA step, s | Total peak |
|---|---|---|---|---|---|
| Lipid kit_076 | 2.47 | 2.47 | 1.88 | 11.5 | 2.53 |
| Lipid **kit_077** | **2.53** | 2.53 | 1.86 | 18.7 | **2.61** |
| Metab kit_076 | 2.51 | 2.51 | 2.00 | 12.5 | 2.58 |
| Metab **kit_077** | **2.57** | 2.57 | 1.95 | 19.2 | **2.65** |

- **Memory:** the lifecycle peak rose 0.06 GiB on both omics, in normalisation, where the QC figures are written. The journeys disagree (lipid 2.47 to 2.45, metab 2.50 to 2.56), and the lipid peak has ranged over 0.2 GiB, so this is not separated from noise. It is now above stage 3's total peaks (2.53, 2.61). KIT-080 judges.
- **Time:** the DA step takes about 7 s more on both omics. That covers drawing the heatmaps, which previously failed at once, and saving the figures. The two were not separated.
- **Disk:** a member's store grows from about 23 MiB to about 65 MiB, since every figure is now kept twice: as a file and as an asset. This scales with assays and contrasts, not samples, and was not budgeted in the ticket.

## Deviations and discoveries

- **The owner is found by omic on the session** (`session$userData$workflow_omics_owners`, set at attach), not passed in. KIT-075 rejected a single lookup slot because one session can run several omics; a key by omic does not have that problem, and the figure sites are exported kernels with no manager argument. **For KIT-078:** a bare reopen has no `omic_type` until a session attaches, so the worker must set it or figures silently render without saving.
- **Identity uses what the figure draws, not the source state's id** as the ticket words it, because the DA tab's source need not be the head. It is stricter: different data can never hit.
- **DA plot objects are still built before the store is asked;** only the render is skipped on a hit. The QC sites plot nothing on a hit. DA figures do not recompute any statistics either way.
- **No rebuild function.** Plot data is saved for the four types, but figures are served from bytes; nothing rebuilds a plot from the data yet.
- **Filtering figures are not written inside the app,** as KIT-030 recorded: the modules pass no directory. So the journeys never exercise those sites. Tests cover both, with a directory given. Whether the app should write them is a separate question, not taken up here.
- **Lipid copy repair rests on its test.** Members do not capture the copy step's output, so the baseline's copy-failure list could not be taken as the ticket asks; the `QC_figures` files themselves are in the baseline.
- **Tests changed on purpose:** the fixtures that seeded `composite_QC_figure.*`, and the two heatmap mocks. Each had been hiding a real defect.
- **Out of scope, as the ticket states:** the Glimma interactive volcanoes.

## Backward compatibility and rollback

File names and paths are unchanged; the new files are the heatmaps that should always have been there. Proteomics figure identities and copies are unchanged (`test-workflow-dia-products.R`, the proteomics gates). Rollback: `git revert`; the stored figures become unreferenced and are ignored by the reverted build, and the heatmap defect returns.

## Handoff

KIT-078 gets figure products a worker can produce, with the `omic_type` note above. KIT-067 inherits omics figures that are entities, with their files still written.
