# KIT-036: Self-contained pages without pandoc

Status 2026-10-16: **accepted and complete.** Stage `kit-journey-036` decision `accepted`; `kit_036` (`core/browser-smoke.h6QeDJ`, source `6c9c4901` = `kit_017`'s source plus KIT-037 and this ticket) is the control from here (protocol section 2), and stage 1 ends on an accepted build. Revision 1.1 (the page build, not the worker).

## The cause (KIT-017 post-hoc audit PH-01)

`htmlwidgets::saveWidget(selfcontained = TRUE)` writes the page with its dependencies on disk and then runs pandoc over the whole document to embed them. On the DA volcano page (14.7 MB for 300 samples, almost all of it the widget's data payload) pandoc peaks at about 3 GiB resident for several seconds, in whichever process saves the page. That transient was the DA step's peak on both KIT-017 builds; the worker's baseline on top of it is why the pair read 1.129 × on the DA step and 1.077 × on the whole slice while resident end read 0.927 ×.

## The mechanism (one)

`saveSelfContainedWidget(widget, file)` (`R/func_general_plotting_qc_utility_helpers.R`): the widget is rendered with `htmltools` (its head, body and resolved dependencies); every dependency's stylesheet becomes a `<style>` element and every script a `<script>` element, read once from the package's files; the page is written as one UTF-8 file. The Glimma volcano widget's eight dependencies are 0.94 MB of plain script and style with no asset references and no element terminators in their contents. A dependency with attachments, a file whose contents would close its own element, a non-widget, or the option `multischolar.widget_pandoc = TRUE` take the pandoc writer as before. The DA volcano save (`writeDIAInteractiveVolcanoOutputs`) and the two enrichment page saves use it; nothing else changes.

## Evidence before the pair

| Check | Result |
|---|---|
| [page-writer-probe.R](page-writer-probe.R) → [page-writer-probe-kit_017-counted-store.json](page-writer-probe-kit_017-counted-store.json): the DA volcano page of the KIT-017 counted candidate's store (2,476 proteins × 297 samples), each writer in a fresh R process, child processes sampled at 50 ms | inline: 9.4 s, R peak 1.13 GiB (0.65 before the widget), **no child process**, page 14,695,098 bytes. pandoc: 15.8 s, R peak 1.13 GiB, **pandoc child peak 2.84 GiB**, page 14,690,907 bytes. |
| `test-workflow-dia-products.R` (new: the inline page carries the same JSON payload as the pandoc page and every dependency verbatim, references nothing by path, and pandoc is not called; the pandoc writer behind the option) | 149 pass |
| `test-general-enrichment-shared.R` (the enrichment saves, with their stand-in writer) | 103 pass, 2 skips |
| Rule 7 of the README's pre-launch section | the mechanism's cost at the step measured on the real page; the owner of the old peak (pandoc) and of the new one (the widget build in R) named |

Expected from the KIT-017 traces: the candidate's DA step about 3.5 GiB against the control's 5.41 (the worker baseline and the widget build remain; the 2.8 GiB child is gone), the whole-slice peak back at the rollup step, resident end unchanged near 0.93 ×.

## Screen lifecycle run (2026-10-16 10:02–10:25, not counted; stage archived as `kit-journey-036-screen-lifecycle`, decision kept as [screen-lifecycle-decision.rds](screen-lifecycle-decision.rds))

Both members complete and verified: control `kit_018` 566 s, 2.14 GiB anonymous, 1.82 GiB resident end, DA step 1.93 GiB; candidate `kit_036` 557 s, 2.14 GiB anonymous, 1.37 GiB resident end, DA step 2.14 GiB. Whole slice 1.000 (KIT-017's screen read 1.128), resident end 0.755, DA step 1.109 (KIT-017's screen read 1.32): the pandoc transient is gone from the step and what remains above the control is the worker's baseline on an 8-sample page (+0.21 GiB). The candidate's `da_publication` completion record reads `status complete`, `browser_acknowledged TRUE`, `measurement_valid TRUE` (KIT-037's live check; every KIT-017 run had read `browser_ack_timeout`). The screen decision (`rejected`, the DA step on 8 samples) is the lifecycle observation, not the claim.

## Pair (counted; 2026-10-16 10:36–13:51; tracked summary [kit-journey-036-stage.json](kit-journey-036-stage.json); 300-sample corpus, no launch record, 14,400 s ceiling, one session, control first)

| Member | Build | Outcome | Elapsed | Slice peak | Anonymous peak | Resident end | Max gap | Verified |
|---|---|---|---|---|---|---|---|---|
| control-1 | `kit_018` | complete, 18 steps | 4,913 s | 7.47 GiB | 5.66 GiB | 3.04 GiB | 170 ms | every state equal to the KIT-005 oracles |
| candidate-1 | `kit_036` | complete, 18 steps | 4,741 s | 7.50 GiB | 5.69 GiB | 2.80 GiB | 130 ms | every state equal to the KIT-005 oracles; `da_publication` acknowledged (KIT-037) |

Decision (`kitJourneyTicketDecision`, rule `resident_end`, protocol v1.4): **`accepted`.** Resident end 0.921 (≤ 0.95; the KIT-019 mechanism's claim, kept). Whole-slice anonymous peak 1.005 (≤ 1.02; the journey's peak is the P14 rollup on both builds, 5.66 and 5.69 GiB). Per step, every guarded step within 1.02: the DA step **0.675** (5.41 → 3.65 GiB), the step the KIT-017 pair rejected at 1.129; q-value 0.951, IQ 0.945, accession 0.950, protein intensity 0.951, duplicates 0.950, protein replicate 0.949, normalisation 0.992, correlation 0.967, session export 0.995; the highest guarded ratio is the rollup at 1.005. Import and design reported (0.998, 1.007). Time, reported: the candidate ran 172 s faster, the DA step included.

What changed at the DA step and why it is enough: the pandoc child process (2.84 GiB on this page, measured in the probe) is gone, so the step's peak is now the worker's baseline plus the widget build in R on top of the application's plateau, 3.65 GiB, 1.76 GiB below the control's in-process publication. The resident-end reduction is unchanged from KIT-017's rerun (0.927 there, 0.921 here), as expected: this ticket moved nothing that stays resident.

Against the production-default baselines (KIT-005 `default-1` and KIT-011's `kit_006` control, both at the memory guard in imputation at 15.40 GiB after 8 steps): the control build completes the journey at 5.69 GiB anonymous, 0.37 × the guard reading, and ends at 2.80 GiB resident.

## Regression (standing rule)

Full suite at `6c9c4901` from a clean worktree ([logs](check-logs/)): 343 files, 29,666 pass, 0 fail, 140 skip; the one error is the PERF-101 snapshot test, the known runner artefact (passes alone). `R CMD check` 3 warnings, 3 notes, nothing new against KIT-002's record.

## Claim limits

The page's bytes differ from pandoc's (pandoc minifies whitespace slightly); the widget data and dependency contents are identical, and the page is a reported product, not an oracle. Pages with dependency attachments (none in this package's widgets today) still go through pandoc.
