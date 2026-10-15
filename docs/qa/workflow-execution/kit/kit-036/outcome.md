# KIT-036: Self-contained pages without pandoc

Status 2026-10-16: **implemented; install and pair pending.** Revision 1.1 (the page build, not the worker). Control `kit_018`; candidate `kit_036` = `kit_017`'s source plus KIT-037 and this ticket.

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

## Pair

_Filled from stage `kit-journey-036`._

## Claim limits

The page's bytes differ from pandoc's (pandoc minifies whitespace slightly); the widget data and dependency contents are identical, and the page is a reported product, not an oracle. Pages with dependency attachments (none in this package's widgets today) still go through pandoc.
