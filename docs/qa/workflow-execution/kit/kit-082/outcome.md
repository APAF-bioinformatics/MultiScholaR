# KIT-082 outcome: the omics application reads no assay outside the managed steps

Status 2026-11-04: **completed**. Stage 3b. On a managed session the lipidomics and metabolomics tabs no longer read a whole assay table in the application process. **On the screen corpora this did not lower the application's memory.** Peak controller memory is unchanged within noise, so this ticket supports no memory claim. KIT-084 tests whether it matters at larger sizes.

- **Commits** on `performance-v2`: code caece7fd; records abb63470; driver correction 963c4333; bindings d46f1699 (closure) and 5f5b5fe6 (rerun).
- **Install:** `kit_082` = `core/browser-smoke.UhTReX` (caece7fd). No package code changed after it.
- **Ticket:** `tickets/kit/KIT-082.md`, revision 1.1.

## What changed

| Read | Now | Sites |
|---|---|---|
| No table | Zero-row prototypes, stored shapes, design, slots | QC finalise validation; QC data summary |
| Id column only | One column | Unique feature counts; duplicate detection; study parameters |
| Page by page | Nothing whole kept | QC stats missingness; ITSD candidate tables |
| Rows the kernel uses | Matching rows only | QC internal-standard analysis; DA heatmap's top-N |
| Worker, commits nothing | Registered `<prefix>_` operation | `filtering_metrics`, `final_pca`, `volcano_page`, `session_export`, `integration_export` |

The runner's queue is now first in, first out, and a queued action dropped by a failure hears the failure.

## Results

| Gate | Evidence | Result |
|---|---|---|
| Probe fails on the base | [probe-on-base.txt](probe-on-base.txt) | Fails on 44cf3abf: 8 assay tables read on each omic, and no iframe. |
| Probe passes, with its control | `test-workflow-omics-application-reads.R` | 0 assay tables read on either omic on the changed tree; the in-process tabs read above 0. |
| Each site equals the whole-object read | same file | 97 pass, including more than one page. |
| Worker sites in the browser | stage JSONs | All five operations dispatched and complete on all four omics members, the final PCA twice. |
| Omics records | [exactness-gates.txt](exactness-gates.txt) | Identical to the `kit_081` oracles on both omics. |
| Proteomics | same; `gui-082.txt`; replay | LFQ, TMT and DIA_LIMPA identical to the `kit_065` oracles. KIT-028 replay exact. `kit-gui-082` observed: 0 fingerprint mismatches against `kit-gui-079`, 13 of 13 pages. |
| Journeys | `kit-journey-082l`, `kit-journey-082m` | `accepted`, exact. |
| Member records | [082l](member-record-replay-082l.json), [082m](member-record-replay-082m.json) | Replay exact on 10 states each; tampered input refused. |
| Lifecycle | `kit-lifecycle-lipid`, `kit-lifecycle-metab` | Both `lifecycle_complete_observed`. Largest sampler gaps 0.06 s. |
| Figures | [082l](figure-compare-082l.txt), [082m](figure-compare-082m.txt) | 99 of 99 and 91 of 91 files are stored, rooted assets. The tool's `RESULT differs` is the 12 PCA products KIT-077 explained, and only those. |
| Dependents | `dep082` | 175 files: 15,690 pass, 0 fail, 0 error, 86 skip, including the module-CI lanes as tests. |
| Package | [check.log](check.log) | Baseline: 3 WARNINGs, 3 NOTEs. |

## The application's memory

Peak controller anonymous memory by step in GiB, lifecycle members. The `kit_081` metabolomics row is KIT-081's baseline. **Those baseline members never drew the volcano** (see the failed attempt below); the `kit_082` rows do.

| Member | import | design | intensity | duplicates | QC finalise | normalisation | correlation | export | DA | report | peak |
|---|---|---|---|---|---|---|---|---|---|---|---|
| metab, `kit_081` | .51 | .53 | .63 | .64 | .65 | .57 | .57 | .61 | .74 | .88 | 0.89 |
| metab, `kit_082` | .51 | .53 | .63 | .65 | .65 | .58 | .58 | .62 | .80 | .93 | 0.93 |
| lipid, `kit_082` | .51 | .53 | .65 | .66 | .67 | .60 | .59 | .76 | .76 | .87 | 0.87 |

The lipid `kit_081` peak was 0.91. The journeys peak at 0.93 (lipid) and 0.92 (metab). The whole-slice peaks are 1.82 to 1.86 GiB, against 1.85 to 1.86 on `kit_081`.

**Reading this:** the screen assays are about 300 rows by 24 samples, so a whole read was never large enough to see in this trace. The peak here is not the matrix. What makes it up (packages, plot objects, the display) was not measured by this ticket.

## What got slower

Both builds on AC in low-power mode, journey members, seconds. Each row is `kit_079` then `kit_082`, lipid / metab.

| Step | `kit_079` | `kit_082` |
|---|---|---|
| QC finalise | 6 / 5 | 9 / 8 |
| Normalisation | 56 / 57 | 59 / 61 |
| Correlation | 9 / 9 | 12 / 12 |
| Session export | 1 / 1 | 9 / 9 |
| DA | 15 / 15 | 23 / 24 |
| Report | 3 / 3 | 7 / 7 |
| **Whole journey** | **148 / 146** | **179 / 179** |

So the journey takes about 30 s longer, roughly a fifth. Part of the DA and report increase is the driver doing more: it now draws the volcano and saves the study parameters.

Worker action by site, request to last journal write ([site-durations.txt](site-durations.txt)): filtering metrics 2.7 s; final PCA 3.6 to 3.8 s; session export 3.7 s; integration export 3.0 to 3.3 s; volcano page 6.6 to 7.4 s. Each of these tabs was synchronous and now shows a wait while the worker runs. The non-worker sites were not timed separately.

## Decision for the user: item 8, the DA tables on the tab

`readWorkflowOmicsDaResults` still hydrates the per-contrast tables in the application, which carry per-sample values. They are not an assay state, so the probe does not count them. Measured on the screen journeys ([da-tab-read-size.txt](da-tab-read-size.txt)): **0.54 MiB lipid, 0.80 MiB metab**, for two products. At this size paging them gains nothing measurable. Recommendation: decide after KIT-084 shows the size at scale.

## Failed attempt, kept

- **`attempt1-volcano-combined`:** the first closure's members all completed, but none drew the volcano. The assay selector opens on Combined, whose banner satisfied the driver's wait. I found it because `<prefix>_volcano_page` was absent from the dispatches, after telling the user the driver exercised it. This was an instrumentation failure of the driver. It is corrected in 963c4333: the driver selects the first assay, waits for the plot itself, and fails the step on a plot error. The four omics members were rerun once; this was the ticket's one correction. The attempt's evidence is in [attempt1-volcano-combined/](attempt1-volcano-combined/) and the runtime directories carry the same suffix.

## Deviations

1. **Baseline build.** The ticket asks for the trace against `kit_079`, which has no controller trace. It is reported against `kit_081`, which changed no package code and whose records are identical. Step times are against `kit_079`, the last members on AC.
2. **Kernel access.** The heatmap kernels gained a one-line data access (`read_source`) so the top-N rows are read once chosen. No arithmetic changed.
3. **`kit-gui-082` was not rerun** after the driver correction: the corrected step is the omics DA step, which the proteomics member does not run. `R CMD check` ran on d46f1699; since then only the driver and the bindings changed.
4. **Module-CI lanes** ran as tests. The content-hash comparison (`kit-069/mci-hash.R`) was not run.

## Limits

- **Tests only, not in the browser:** the QC internal-standard analysis, the heatmap, the ITSD observer's page build. The driver does not open them.
- **Fixtures only:** the probe's zero count is on test fixtures, under `testServer`. The members show the worker operations ran; they do not count table reads.
- **Pixels:** the tab's final PCA is grobs built on a png device, which drew byte-identically to `print(plot)` in the test. If Shiny's device differs, pixels could differ.
- **Time** is one run of each member. An earlier run of the same build took up to 262 s, so step differences of a few seconds are not reliable.

## Rollback

`git revert`. The worker operations commit nothing, so stored states are unaffected.

## Hand-off

KIT-084 runs on this build. KIT-066 and KIT-080 take it for the real corpora.
