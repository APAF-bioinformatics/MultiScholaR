# KIT-065: proteomics closure runner

Status 2026-10-28: **in progress**. Code 9497cd26 (the package code f55bef1b: `RenderReport`'s environment; the harness 6b8d2581 + 6f65c571: the report step). Install `kit_065` = `core/browser-smoke.InuLC5` (earlier candidate q23r8g, before the render fix). No counted member yet: `kit:journey-065:candidate:1` is declared below and awaits the user's go.

## Mechanism (one): the journey ends at the report

- `kitJourneySteps()` gains `report` (kind `report`, the key the driver dispatches on) on the DIA, DIA_LIMPA and protein routes: the summary tab's publication copy (`copy_to_publication`, which bootstraps the workflow arguments when they were not saved) and the report render (`generate_report`) with the route's shipped template, waited on the "Report generated successfully!" notification. The step is reported-only in the pair guard (no memory claim) and has no state: the record keeps the head the DA left and the rendered file's path and size.
- The KIT-021 and KIT-053 recorders run the same two actions in process after the DA (`bootstrapProtSummaryCopyFallbackStudyParams`, `ensureProtSummaryReportTemplate`, `copyToResultsSummary`, `RenderReport`) and write `report = list(template, path, bytes)` into the record. The step maps gain the row. The bindings bind the records as `limpa65_screen`, `lfq65_fixture`, `tmt65_fixture`, the four screen members and the counted member.
- Enrichment is not in the journeys (user decision 2026-10-28, revision 1.1): the kit corpora are synthetic, the engines need real accessions and the network, and g:Profiler's results move with its releases. KIT-066 proves the enrichment tab on real data after stage 3.

## What the report step found, and the repair

The first GUI run of the lifecycle member completed the step in 0 s without a document: the step carried kind `store` and fell through the driver's `store` case (fixed: kind `report`). The second run failed at the render: every one of the six report templates' setup chunk looks up `getProjectPaths()` by name (`exists(..., mode = "function")`), and the kit launches the application without attaching the package (`getExportedValue("MultiScholaR", "run_app")`), so `rmarkdown::render` in `new.env(parent = globalenv())` could not see it ("getProjectPaths function not found. Ensure MultiScholaR package is loaded"). The in-process recorders attach the package, which is why they rendered. Repair, once, in `RenderReport`: the render environment's parent is the package namespace (global objects such as `project_dirs` remain reachable through the namespace's parent chain); the two `RenderReport` tests pin it. The six copied setup chunks are the subject of KIT-067 (the user's request that the templates' build logic be synced; the LFQ template is the reference).

## Evidence

### The four screen journeys through the GUI on `kit_065` (InuLC5), each through the report

| Member | Stage archive | Steps | Verification | Elapsed | Peak (anon) | Report |
|---|---|---|---|---|---|---|
| DIA lifecycle (`kit-lifecycle-complete`, screen) | `kit-lifecycle-complete.kit065r-dia-102150` | 19 complete | 0 mismatches against KIT-005 (18 steps); `lifecycle_complete_observed` | 306 s | 2.11 GiB | `DIANN_report_*.docx`, 3.1 s |
| DIA_LIMPA (`kit-journey-065m`, screen) | `kit-journey-065m.kit065r-limpa-102630` | 14 complete | accepted, 0 mismatches (`limpa65_screen`) | 221 s | 2.12 GiB | `DIANN_limpa_report_*.docx` |
| LFQ (`kit-journey-065l`, MaxQuant fixture) | `kit-journey-065l.kit065r-lfq-102940` | 11 complete | accepted, 0 mismatches (`lfq65_fixture`) | 127 s | 1.79 GiB | `LFQ_report_*.docx` |
| TMT (`kit-journey-065t`, PD fixture) | `kit-journey-065t.kit065r-tmt-103249` | 11 complete | accepted, 0 mismatches (`tmt65_fixture`) | 121 s | 1.77 GiB | `TMT_report_*.docx` |

Summaries: [lifecycle-complete.json](lifecycle-complete.json), [kit-journey-065m-stage.json](kit-journey-065m-stage.json), [kit-journey-065l-stage.json](kit-journey-065l-stage.json), [kit-journey-065t-stage.json](kit-journey-065t-stage.json). Each member's project holds the rendered document (two files: the report and its copy under the publication directory). The records were made on InuLC5 by the recorders (`oracles-*.rds` here): LFQ and TMT 11 states, DIA_LIMPA 14, the in-process report render in about a second each.

### Replays, adoption, check, suite, fallbacks

- KIT-028 DIA replay: the package code of 9497cd26 is 1ab2318b's plus `RenderReport`'s environment; the replay on 1ab2318b was exact (16 states, 174 s, KIT-064); the no-tools replay below runs on this code.
- DIA_LIMPA replay against the KIT-065 record (`KIT_LIMPA_TICKET_DIR=kit-065`, the q23r8g library, the record of tdDG9e): **exact, 11 states, DA replayed, 145 s** ([replay-check-screen-kit065.json](replay-check-screen-kit065.json), [analysis_record-limpa-screen-kit065.json](analysis_record-limpa-screen-kit065.json)).
- Legacy RDS adoption and recovery (`test-workflow-native-session.R` 725/2, `test-workflow-saved-recovery.R` 57/2) at 6b8d2581; `test-workflow-kit-harness.R` 435/0 with the report step; `test-general-filemgmt-report-helpers-shared.R` 85/1 with the render-environment assertions.
- `R CMD check --no-tests --no-manual` at 9497cd26: CHECK.
- The closure full suite on a clean worktree of 9497cd26 (the InuLC5 library): SUITE.
- Capability fallbacks (`MULTISCHOLAR_DISABLE_TOOLS=bash,dd,sha256sum,sync,tar,cp`, KIT-044's recipe, in process on this build): the three recorders without the tools (LFQ 11 states, TMT 11, DIA_LIMPA 14, each through the report) are **identical state for state** to the tool-enabled records (class, rows, components, support, summary, products); the KIT-028 DIA replay without the tools **exact, 16 states, 144 s** on this code ([replay-check-kit065-notools.json](replay-check-kit065-notools.json), [analysis_record-origin-kit065-notools.json](analysis_record-origin-kit065-notools.json)).

### The counted member

Declared: `kit:journey-065:candidate:1`, stage `kit-journey-065`: the DIA route's whole journey on the 300-sample DIA-NN corpus on `kit_065` (InuLC5), the 18 recorded steps judged against the KIT-005 oracles and the report step completion-checked; ceiling 14,400 s; one attempt; memory decides, time is reported; the harness decision reported verbatim (the state-summary rule carried to KIT-050 applies, as it did to KIT-048's member). COUNTED.

## Closure statement

CLOSURE.
