# KIT-082 outcome, revision 1.3: the omics application reads no assay outside the managed steps

Status 2026-11-05: **completed**. Stage 3b. This replaces [the revision 1.1 outcome](rev-1.1/outcome.md), whose "no fall in memory" came from comparing members run on different drivers. The audit's findings and their closure are in [the audit](../KIT082-AUDIT-20261104.md) and [the remediation](../KIT082-REMEDIATION-20261104.md).

**On the screen corpora, with one driver on both builds, the application's peak falls by about 0.2 GiB: 1.22 to 1.02 GiB on lipidomics, 1.16 to 0.94 GiB on metabolomics. The total for the whole slice does not fall, and each journey took 13 to 46 s longer.** These are uncounted members on a synthetic corpus, one control run per omic, so this supports no claim about real data. KIT-080 decides that.

- **Install:** `kit_082` = `core/browser-smoke.UzND1n` (c7557ac9). No package code has changed since: `git diff c7557ac9 HEAD` on `R/`, `DESCRIPTION`, `NAMESPACE`, `src/`, `man/` and `inst/` is empty. After it, only the corpus tool, the bindings generator and the records changed.
- **Control:** `kit_081` = `core/browser-smoke.6nbYY3` (dbd5b88b), on the final driver.
- **Superseded installs, evidence kept:** UhTReX ([rev-1.1/](rev-1.1/)), cSzZKS ([rev-1.2a-monitor-race/](rev-1.2a-monitor-race/)), pnUSaD ([rev-1.2b-monitor-fix-regression/](rev-1.2b-monitor-fix-regression/)).

## What changed

| Read | Now | Sites |
|---|---|---|
| No table | Zero-row prototypes, stored shapes, design, slots | QC finalise validation; QC data summary |
| Id column only | One column | Unique feature counts; duplicate detection; study parameters |
| Page by page | Nothing whole kept | QC stats missingness; ITSD candidate tables |
| Rows the kernel uses | Matching rows only | QC internal-standard analysis; DA heatmap's top-N |
| Worker, commits nothing | Registered `<prefix>_` operation | `filtering_metrics`, `final_pca`, `volcano_page`, `session_export`, `integration_export`, and from revision 1.2 `da_tab_results` |

Revision 1.2 also moved the DA tab's rows to the worker. The tab now holds its rows without the per-sample columns, per-assay tables or fits, and the download is written from the store one contrast and assay at a time.

Revision 1.3 changed three things:
- **The monitor** decides from one read of a worker's memory counter. This shared code is in scope by the user's decision. It gives up one thing: a live worker whose counter can never be read has no `nonworker_reservation`. The slice counters, early stop, host floor and memory limit remain.
- **Failure paths:** the fixes for A-04 to A-08 (see the remediation).
- **The driver** now observes each display by what it shows.

## Results

| Gate | Evidence | Result |
|---|---|---|
| Probe fails on the base | [probe-on-base.txt](probe-on-base.txt) | Fails on 44cf3abf: 8 assay tables read on each omic. |
| Probe passes, with its control | `test-workflow-omics-application-reads.R` | 0 assay tables read on either omic; the in-process tabs read above 0. The cells kept by each store read stay below the whole table's. |
| Displays in the browser | [displays.txt](displays.txt) | Lipid on `kit_082`: "Found 12 internal standards"; static volcano ink 0.0278; heatmap ink 0.6571; download 365 rows, 24 sample columns. Metab: found 12; ink 0.0392; ink 0.6615; 564 rows, 24 sample columns. The ITSD tables and results table each showed a data row, or the step would have failed. |
| The checks can fail | [made-to-fail.txt](made-to-fail.txt) | `082c` and `082d` (no internal standard in the corpus) fail at `qc_finalize` on the app's error. `082e` and `082f` (nothing significant) fail at `da`: heatmap "drew nothing", ink 0. |
| Worker sites in the browser | stage JSONs | All six operations dispatched and complete on all four `kit_082` omics members, the final PCA twice. |
| Omics records | [exactness-gates.txt](exactness-gates.txt) | Identical to the `kit_081` oracles on both omics. |
| Proteomics | same; [gui-082.txt](gui-082.txt); [replay](replay-check-kit082.json) | LFQ, TMT and DIA_LIMPA identical to the `kit_065` oracles. KIT-028 replay within budget. `kit-gui-082` observed on the final build: 0 fingerprint mismatches against `kit-gui-079`, 13 of 13 pages. |
| Journeys | stage JSONs | `082l`, `082m`, and on `kit_081` `082a`, `082b`: all `accepted`, exact, no oracle mismatches. Control and candidate have the same state names and row counts at every step, and their two oracles are identical (the omics records row). Entity ids are not comparable: they differ between any two members, even on one build. |
| Member records | [082l](member-record-replay-082l.json), [082m](member-record-replay-082m.json) | Replay exact on 10 states each; tampered input refused. |
| Lifecycle | stage JSONs | Both `lifecycle_complete_observed`. Largest sampler gap 0.07 s, no swap. |
| Figures | [082l](figure-compare-082l.txt), [082m](figure-compare-082m.txt) | 99 of 99 and 91 of 91 files are stored, rooted assets. The tool's `RESULT differs` is the 12 PCA products KIT-077 explained, and only those. |
| Monitor | [monitor-window-stress.txt](monitor-window-stress.txt) | On 30 real starting and stopping units, aborts: second-look 30, the original check-then-read 11, one-read 0. No member's session died on this build. |
| Dependents | `dep082`, on c7557ac9 | 181 files: 16,737 pass, 0 fail, 0 error, 97 skip, including the module-CI lanes as tests. |
| Rerun after the corpus tool changed | [final-targeted.txt](final-targeted.txt) | Harness 525, execution 106 (10 skip), application reads 204, residency 82: 0 fail, 0 error. |
| Package | [check.log](check.log) | Baseline: 3 WARNINGs, 3 NOTEs. |

## The application's memory, like for like

Peak anonymous memory by step in GiB, final driver on both builds. Members `082a`, `082b` are the control. The candidates are the journey and lifecycle members.

**Application (controller):**

| Member | import | design | intensity | duplicates | QC finalise | normalisation | correlation | export | DA | report | peak |
|---|---|---|---|---|---|---|---|---|---|---|---|
| lipid, `kit_081` | .51 | .53 | .64 | .66 | .66 | .59 | .59 | .64 | **1.09** | **1.22** | **1.22** |
| lipid, `kit_082` journey | .51 | .53 | .64 | .66 | .67 | .59 | .59 | .76 | **.87** | **1.02** | **1.02** |
| lipid, `kit_082` lifecycle | .51 | .53 | .63 | .65 | .66 | .59 | .70 | .76 | **.88** | **1.01** | **1.01** |
| metab, `kit_081` | .51 | .53 | .64 | .65 | .66 | .59 | .72 | .76 | **1.03** | **1.16** | **1.16** |
| metab, `kit_082` journey | .51 | .53 | .65 | .65 | .65 | .58 | .58 | .61 | **.83** | **.94** | **.94** |
| metab, `kit_082` lifecycle | .51 | .53 | .63 | .66 | .66 | .59 | .59 | .65 | **.82** | **.94** | **.94** |

**Whole slice (application plus worker):**

| Member | import | design | intensity | duplicates | QC finalise | normalisation | correlation | export | DA | report | peak |
|---|---|---|---|---|---|---|---|---|---|---|---|
| lipid, `kit_081` | .85 | .87 | 1.24 | 1.28 | 1.01 | 1.77 | 1.30 | .99 | 1.83 | 1.57 | 1.83 |
| lipid, `kit_082` journey | .85 | .87 | 1.24 | 1.28 | 1.29 | 1.77 | 1.38 | 1.57 | 1.83 | 1.59 | 1.83 |
| lipid, `kit_082` lifecycle | .85 | .87 | 1.24 | 1.27 | 1.28 | 1.77 | 1.37 | 1.57 | 1.83 | 1.60 | 1.83 |
| metab, `kit_081` | .85 | .88 | 1.25 | 1.27 | 1.01 | 1.77 | 1.42 | 1.11 | 1.81 | 1.51 | 1.81 |
| metab, `kit_082` journey | .85 | .87 | 1.24 | 1.28 | 1.28 | 1.76 | 1.36 | 1.38 | 1.86 | 1.54 | 1.86 |
| metab, `kit_082` lifecycle | .85 | .87 | 1.24 | 1.27 | 1.29 | 1.77 | 1.37 | 1.38 | 1.87 | 1.53 | 1.87 |

**Reading this:**
- **The fall is at DA and report.** The application falls about 0.2 GiB at both steps, and the two candidates of each omic agree within 0.01 GiB. The steps up to normalisation do not move.
- **Noise.** Runs of one build differ by up to 0.11 GiB at correlation (lipid .59 and .70), and lipid's export is .12 above its control while metab's is below. Differences of that size are not a result. There is one control run per omic.
- **The cost moved.** The slice's peak is unchanged on lipid and 0.05 to 0.06 GiB higher on metab. QC finalise and the session export are higher in the slice (up to 0.58 GiB at lipid's export) because their worker now runs alongside the application. Both stay below the journey's peak, at normalisation and DA.
- **Not the matrix.** The screen assays are about 300 rows by 24 samples, so the fall is not their size. What makes it up was not measured. KIT-084 shows how each number scales.

## Time

Seconds, AC in low-power mode, all on the night of 2026-11-04, with other work on the host throughout (load about 3). One run per member.

| Step | lipid `kit_081` | lipid `kit_082` | metab `kit_081` | metab `kit_082` |
|---|---|---|---|---|
| QC finalise | 15 | 20, 21 | 10 | 15, 15 |
| Normalisation | 120 | 109, 131 | 108 | 110, 125 |
| Correlation | 21 | 22, 22 | 21 | 27, 22 |
| Session export | 5 | 15, 17 | 3 | 14, 16 |
| DA | 44 | 52, 65 | 50 | 59, 52 |
| Report | 6 | 11, 11 | 8 | 15, 13 |
| **Whole journey** | **310** | **323, 356** | **292** | **327, 325** |

So the session export takes 10 to 13 s longer, and QC finalise and the report about 5 s longer each, consistently on both omics. DA is slower by 2 to 21 s. Normalisation runs 109 to 131 s on one build and code that did not change, so the whole-journey difference is rough. These times are not comparable with `kit_079` or earlier closures: the same member has run about twice as slowly on this host with the cause not established. They are reported, not a decision.

## Failed attempts, kept

- **Rev 1.1, volcano:** the driver's wait passed on the Combined banner, so no member drew the volcano, after I had told the user the driver exercised it ([rev-1.1/attempt1-volcano-combined/](rev-1.1/attempt1-volcano-combined/)).
- **Rev 1.2a, monitor race:** a session died when a worker exited between the monitor's check and its read.
- **Rev 1.2b, my fix for it (d52303f4):** worse than what it replaced; sessions died at worker start.
- **Rev 1.3, first `082c`:** this made-to-fail member used a bad pattern and was accepted, because the analysis falls back to its own patterns. The gate stopped the closure (`GATE kit-journey-082c did not fail at qc_finalize`). The design was wrong, so the member moved to corpora with no standard. Its runtime directory has the suffix `.rev13a`, and its line stays in `displays.txt` above the empty one.

## Overstated commit messages

- **7deb06ff** says the members use corpora with no standard. A script chain committed it after its edit had stopped partway.
- **d344c1df** says "complete", but the lipid corpus still failed to generate.

**f6a2f001** is the first commit that completes this, and it was gated on its checks. The screen corpora stayed byte-identical. The members on the no-standards corpora ran after it (it was committed at 23:34, they started at 23:47).

## Deviations

1. **Baseline build.** The ticket asks for the trace against `kit_079`, which has no controller trace. It is reported against `kit_081`, whose records are identical to `kit_079`'s.
2. **Kernel access.** The heatmap kernels gained a one-line data access (`read_source`), so the top-N rows are read once chosen. No arithmetic changed.
3. **`082l` was not rerun.** It ran in the first closure of revision 1.3 on the same install and driver, and the later commits changed only the no-standards corpora.
4. **Module-CI lanes** ran as tests. The content-hash comparison (`kit-069/mci-hash.R`) was not run.
5. **Corrections.** Revision 1.3 was a further correction pass past the ticket's one, on the user's instruction ("fix literally every issue"). Valid results were not retried. The reruns were of members whose driver, monitor or made-to-fail design had changed.
6. **Reference GUI.** `kit-gui-082` compares with `kit-gui-079`, its bound reference. It took 780 s here, against 319 and 604 s earlier: not comparable.

## Limits

- **Made to fail, in the browser:** the empty internal-standard analysis and the empty heatmap. The other failure branches (volcano warning, error and output-error ends; a table with no data row; a rejected download; study parameters saved with warnings) were read in review and not forced.
- **Fixtures only:** the probe's zero count is on test fixtures, under `testServer`. The members show the displays and worker operations ran; they do not count table reads.
- **Still in the application by design:** the DA tab's rows without samples (0.08 and 0.12 MiB here, [da-tab-read-size.txt](da-tab-read-size.txt)); the download reads contrast products at click; the heatmap, QC internal-standard analysis and ITSD pages read rows there. How these scale is KIT-084's question.
- **Pixels:** the tab's final PCA is grobs built on a png device, which drew byte-identically to `print(plot)` in the test. If Shiny's device differs, pixels could differ.
- **Proteomics composite:** still uses the old reader (`R/mod_prot_norm_qc_support_helpers.R:225`). Recorded, out of scope.
- **Full suite:** not run here. The standing rule takes it at the next checkpoint, KIT-080.

## Rollback

`git revert`. The worker operations commit nothing, so stored states are unaffected.

## Hand-off

KIT-084 runs on this build, under the declaration in its ticket. KIT-066 and KIT-080 take it for the real corpora.
