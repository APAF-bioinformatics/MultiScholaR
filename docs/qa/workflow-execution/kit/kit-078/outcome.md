# KIT-078 outcome: the lipid and metab steps on the worker

Status 2026-11-03: **completed**. Stage 3b, ticket 7 of 10. Every lipid and metab compute declaration runs in a worker under a managed application and in process otherwise, with equal products. During a managed step the application materialises no assay table. No memory claim is made here; that decision is KIT-080's under protocol v1.6.

- **Commits** on `performance-v2`: items 1-2 d76c7604, item 3 bcf6a512, items 5-7 98d72ea8, the composite on the worker b5307e91 (user decision), verifier repair 3a7bf5ea, bindings and evidence around them.
- **Base:** 4aed0d98 (KIT-077 revision 1.3).
- **Install:** `kit_078` = `core/browser-smoke.mdSh5c` (b5307e91). Gate installs, kept: `IFYQXu` (items 1-2), `h4c3hq` (item 3). `Gw3jSa` (0d028452) was superseded when the composite moved.
- **Ticket:** `tickets/kit/KIT-078.md`.

## What changed

1. **The `omics` role.** The nine compute declarations per omic run on it. The store steps (QC finalise, norm complete, reset) and the v1 DA stay in process. A worker step with no `@args` section declares `section = NULL`: the recorded change to the registry test's rule.
2. **The worker leaves the head alone.** `computeWorkflowOmicsAction` reopens the store; what the handler saves is collected (`deferred_roots`, `deferred_owned`) and `completeWorkflowOmicsAction` records it and commits in one head.
3. **Three publication operations per omic**, committing nothing: `<prefix>_qc_figures`, `<prefix>_composite_figure` and `<prefix>_da_publication`. Their files are staged and published through the action-output journal, so a failed publication is retried from the same request, never recomputed.
4. **Filtering figures** follow the proteomics pattern: the worker computes the metrics from the matrix, the tab draws the figures from them, and the filtering progress stays on the bus.
5. **The tabs.** One runner per omic (`R/utils_workflow_omics_actions.R`) drives a tab's actions over `installWorkflowManagedStep`, which gained a per-apply operation and a mode with no browser acknowledgement. The normalisation is sequential actions under one progress display, the correlation two, the DA one per contrast and then its export. With no managed application the same tabs run the same declarations in process.
6. **Pearson selection** reads its declared `pairs` entity through the store, in both processes; it was the undeclared live `correlations`.
7. **An emptied omics assay** is now saved as the kernel's result. `inheritWorkflowExecutionSource` takes its fields from the class declaration (the protein class's are the same three; the peptide class keeps the old list).

## Results

| Gate | Evidence | Result |
|---|---|---|
| Worker, fail first | [worker-fail-first.txt](worker-fail-first.txt) | Before: `operation adapter unavailable`. After: complete, both omics. |
| Worker equals in process | `test-workflow-omics-worker.R` | Every state, value, correlation, DA table and provenance entry identical on both omics' fixtures. |
| Hydration probe (item 7) | same | 0 assay tables and 0 assay reads in the controller during managed steps, both omics, with a manual ITSD selection. Positive control: the same tabs in process do read the matrix. |
| Journeys on the worker | [078l](kit078-kit-journey-078l-stage.json), [078m](kit078-kit-journey-078m-stage.json) | Accepted exact against the in-process records. Dispatches named in each summary: 15 requests, all complete. |
| Member records | [078l](member-record-replay-078l.json), [078m](member-record-replay-078m.json) | Replay exact, 10 states; tampered input refused. |
| Oracles | [exactness-gates.txt](exactness-gates.txt) | `lipid78`, `metab78` identical to KIT-077's. LFQ, TMT, DIA_LIMPA identical to kit_065. |
| KIT-028 replay | [replay-check-kit078.json](replay-check-kit078.json) | Exact, 16 states. |
| Proteomics on the chain (`078p`) | [baseline](kit-gui-077-decision.txt), [items 1-2](gate-items12.txt), [item 3](gate-item3.txt), [closure](gate-closure.txt) | Observed each time: 13 of 13 dispatches complete, 0 mismatches. |
| Figures | [078l](figure-compare-078l.txt), [078m](figure-compare-078m.txt) | 99 of 99 and 91 of 91 files are stored, rooted assets. The tool's `RESULT differs` is the 12 PCA products KIT-077 explained, and only those. |
| Repeat export | [078l](repeat-export-078l.txt), [078m](repeat-export-078m.txt) | 0 plots built, 11 of 11 figures byte-identical to stored assets, no asset added. |
| Lifecycle members | [lipid](kit078-kit-lifecycle-lipid-stage.json), [metab](kit078-kit-lifecycle-metab-stage.json) | Observed. |
| Full suite | [full-suite.txt](full-suite.txt) | 355 of 358 files in the parallel run, 30,773 pass. Each other result classified below. |
| R CMD check | [check.log](check.log) | Baseline: 3 WARNINGs, 3 NOTEs, code note unchanged. |

**The full suite's other results:** the 3 `ci-impact` files have no result line under the parallel run, as at KIT-075, and pass alone. The heatmap test's 2 failures came from my runner's log threshold (below), and it passes 12 of 12 without it. The harness test's error was my verifier repair, now fixed, and it passes 483 of 483. `perf101` and `ruv-pearson` are the two known diagnostic errors.

## Memory and time (reported; no claim)

| GiB anonymous | kit_077 lipid | kit_078 lipid | kit_077 metab | kit_078 metab |
|---|---|---|---|---|
| Whole slice, lifecycle | 2.50 | **2.90** | 2.52 | **2.91** |
| Total peak, lifecycle | 2.62 | **2.99** | 2.63 | **3.00** |
| After normalisation, before correlation | 2.35 | **0.94** | 2.34 | **0.93** |
| Session export | 2.34 | 1.11 | 2.34 | 0.95 |
| Report | 1.92 | 1.24 | 1.93 | 1.23 |

- **The application's memory fell; the slice's total rose.** After normalisation the slice holds only the application, so its 2.35 → 0.93 GiB is the application's own drop. The slice total rose about 0.4 GiB because the application and a worker now coexist, and a worker carries the package floor.
- **The whole-slice peak is the composite's worker.** It takes about 2.0 GiB beside the application's 0.93, and every other step stays under 1.5 GiB. The `whole_slice_anon` rule KIT-080 proposes would refuse this on the screen corpus. Whether real sizes change that is what KIT-080 measures. Reading the PNGs more cheaply is the next lever; it is not done here.
- **Moving the composite (user decision):** before the move, on `Gw3jSa`, the application rose from 1.1 to 2.55 GiB assembling it and stayed at 2.35.
- **Time is not comparable.** Both journeys and most of the lipid lifecycle ran on battery in `power-saver`, as [power-state.txt](power-state.txt) records. Those members took much longer than kit_077's on mains, though how much of that is the power state is not established: 078l 418 s and 078m 402 s, against 277 s and 265 s. The metab lifecycle, on mains but still in the `low-power` profile, took 220 s against 268 s. Each worker adds about 7 s of start-up. KIT-080 should time this on mains.
- **Disk.** The lifecycle project grew from 114 to 196 MiB (lipid) and from 111 to 188 MiB (metab). 46 MiB of each is request directories under the store, 41 of it the composite's journal copy. They are not cleaned up.

## Host and run record

- **Sampler gaps: cause not established.** Members before the 15:16 reboot had no gaps. Between that reboot and the one at 18:52 they had bursts, up to 0.88 s against the frozen journey limit of 0.5 s, on battery in `power-saver`. I recorded battery as the cause at the time, and that was wrong or at best incomplete: the final members, after the 18:52 reboot, had gaps under 0.1 s, also on battery in the same profile for the first two. The limit was not touched. The first worker journey (on `Gw3jSa`) is kept as `kit-journey-078l.battery-Gw3jSa`: complete, 0 mismatches, its decision invalid on that rule ([journey-078l-battery.txt](journey-078l-battery.txt)). Its gaps fell while the composite ran in the application, which has since moved; the proteomics members' gaps fell in the first seconds of their runs, so that does not explain them all.
- **Proteomics gate on item 3.** Three launches: one killed by the reboot, one invalid on the sampler rule, one observed. All are kept.
- **Baseline's reference.** The bare `kit-lifecycle-complete` stage no longer exists, so `kit-gui-077` was run first and verified against kit_065's member.
- **My runner.** The suite ran with `logger::log_threshold(ERROR)`, which I added when rebuilding the runner after a reboot cleared `/tmp`. `save_heatmap_products` returns its last log call's value, which produced the 2 false failures. A threshold can only add failures of this kind, so the run was not repeated.
- **Scripts** now live in `<proof>/kit/scripts-078/`.

## Deviations from the ticket

1. **More worker operations than the ticket names.** The ticket lists the compute declarations; the probe needs the QC figures and DA export out of the controller too, and the user decided the composite. Hence the three publication operations.
2. **No per-action browser acknowledgement.** A tab's several actions share one display, so each closes on delivery (status `delivered`). Journeys wait on the tab's own completion text, as before.
3. **`acceptWorkflowActionResult` is unchanged.** It never checked the assay's shape; a test proves the list-of-references checkpoint passes it.
4. **`assertWorkflowExportReady`'s state list is unchanged.** Its action-output guard already refuses dependent work while an omics publication is pending (test); the list drives a proteomics-only resume.
5. **The driver is unchanged.** It already waits on each action's text.
6. **The managed module test** is in `test-workflow-omics-worker.R`, driving the real tabs through the real worker entry and acceptance, rather than on the peptide fixture in `test-perf107-audit.R`.
7. **`078p` is the stage `kit-gui-078`:** the harness's stage pattern takes a ticket number only. Its verifier reports the interactive volcano page rather than judging it: the page's bytes differ between members of one build (kit_065, `bdd543d3cd`).
8. **A managed DA runs on the head only.** A worker's input is the head, so a tab left on an earlier state is told to reload the session; in process it could run on either.

## Limits

- **The application still reads the matrix outside managed steps.** These are per-action reads, not holders. They are the next target if the objective is meant to cover them, and KIT-081's controller trace can size them:
  - duplicate detection;
  - the ITSD candidate tables;
  - QC finalise and its summary;
  - the final-QC PCA on the normalisation tab;
  - the session export;
  - the DA tab's interactive volcano and heatmap;
  - the summary's integration file.
- **The DA tables** are hydrated in the application for the tab's display (KIT-076's list); they carry per-sample values for the contrast's features. Not an assay state, and the probe does not count them.
- **The probe is a test-level observable**, on the fixtures and with the service layer mocked. The real journeys prove the dispatches and the exactness, not the absence of a read.
- **A failed QC-figure action fails the pipeline**, as in process; the composite and the DA export only warn.
- **Not inspected by eye:** the composite and the grids (files, stored assets and sizes only).

## Rollback

`git revert`. Requests written by this build are not resumed by the reverted one. Stored states are unchanged in form, so a reverted build reopens them.

## Handoff

KIT-079: a worker that passes the variant to the omics declarations. KIT-081: the controller trace should separate the application from the workers, since the slice's total now moves the other way. KIT-080: the whole-slice rule against the composite's 2.0 GiB worker, and timing on mains.
