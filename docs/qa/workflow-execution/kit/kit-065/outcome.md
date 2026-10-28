# KIT-065: proteomics closure runner

Status 2026-10-28: **completed** (proteomics closed on the kit). Code 9497cd26 (package code f55bef1b; harness 6b8d2581 + 6f65c571). Install `kit_065` = `core/browser-smoke.InuLC5`. Counted member run on the user's go of 2026-10-28.

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
- `R CMD check --no-tests --no-manual` at 9497cd26 (git-archive export, `_R_CHECK_FORCE_SUGGESTS_=false`, the InuLC5 library): 3 WARNINGs, 3 NOTEs, the baseline.
- The closure full suite on a clean worktree of 9497cd26 (the InuLC5 library): **29,084 pass, 0 fail, 2 errors (the two known diagnostic files), 307 skips, 346 files**: green by the rule recorded since KIT-017.
- Capability fallbacks (`MULTISCHOLAR_DISABLE_TOOLS=bash,dd,sha256sum,sync,tar,cp`, KIT-044's recipe, in process on this build): the three recorders without the tools (LFQ 11 states, TMT 11, DIA_LIMPA 14, each through the report) are **identical state for state** to the tool-enabled records (class, rows, components, support, summary, products); the KIT-028 DIA replay without the tools **exact, 16 states, 144 s** on this code ([replay-check-kit065-notools.json](replay-check-kit065-notools.json), [analysis_record-origin-kit065-notools.json](analysis_record-origin-kit065-notools.json)).

### The counted member

Declared: `kit:journey-065:candidate:1`, stage `kit-journey-065`: the DIA route's whole journey on the 300-sample DIA-NN corpus on `kit_065` (InuLC5), the 18 recorded steps judged against the KIT-005 oracles and the report step completion-checked; ceiling 14,400 s; one attempt; memory decides, time is reported; the harness decision reported verbatim (the state-summary rule carried to KIT-050 applies, as it did to KIT-048's member).

Run 2026-10-28 11:22, stage archive `kit-journey-065.kit065-counted`, summary [kit-journey-065-stage.json](kit-journey-065-stage.json):

| Measure | KIT-065 (InuLC5) | KIT-048 (2GA8g6) | KIT-036 candidate |
|---|---|---|---|
| Journey | **complete, 19 of 19 steps** (the report rendered in 27 s) | complete, 18 of 18 | complete |
| Elapsed | 3,334 s | 2,927 s | 4,740 s |
| Peak anonymous | 5.94 GiB (cgroup 7.36 GiB) | 5.96 GiB | 5.69 GiB |
| Resident at end | 2.85 GiB | 2.91 GiB | |
| Harness decision | **`invalid`**, rule `verifier_accepted` | `invalid`, the same rule | accepted |

The verifier's 18 mismatches are every state's `summary` and the normalisation step's identity product, the same set as KIT-048's member; the member's store against the bound KIT-005 oracles, field by field: 17 states differ in `summary_sha256` only, the `initial` placeholder in its recorded class as well, **no state differs in class, rows, components or support**. The cause is the one recorded in KIT-048 and carried to KIT-050: the summary hashes the state's whole `config`/`args`, which have gained shipped default keys since the KIT-005 control's build. The elapsed time is longer than KIT-048's by about 400 s, inside the step-level variance seen between members (imputation 550 s against 384 s, normalisation 1,155 s against 1,048 s; the GUI process's own variance, recorded in KIT-057); time is reported, not judged. The ticket's criterion, completion with the decision, peaks and times recorded verbatim, is met.

## Closure statement

Scope: the proteomics work of the KIT series on the kit host (x86_64 Linux, 32 GiB, the P14 engine), the synthetic DIA-NN corpora (8-sample screen, 300-sample full), the MaxQuant and Proteome Discoverer fixtures, and the KC real LFQ corpus through KIT-054's lane only.

- **Four routes, end to end through the GUI, on one build (`kit_065` = InuLC5, source f55bef1b):** DIA (19 steps), DIA_LIMPA (14), LFQ (11) and TMT (11), each from import through design, QC, normalisation, correlation, session export, differential abundance and the rendered report. Every state the journeys commit verifies against an in-process record or the KIT-005 oracles; the DIA route's 300-sample member completes with every component identical to KIT-005's.
- **The architecture the manuscript describes holds:** the state owner attached from session start (KIT-048), every step a declared operation on saved entities, the managed workers for the memory-heavy steps, the records and replays (KIT-028 exact on this code, the DIA_LIMPA replay exact against its record), the capability fallbacks without the GNU tools identical state for state (KIT-044, re-proven here), and the omic-agnostic names stage 3 builds on (KIT-045).
- **Memory, the series' first claim:** the 300-sample DIA journey runs in 5.94 GiB anonymous peak with 2.85 GiB resident at the end (KIT-001's pre-kit lifecycle exceeded 15 GiB on the same corpus; KIT-011's accepted candidate 5.51 GiB; KIT-036 5.69 GiB; KIT-048 5.96 GiB). DIA_LIMPA on the same corpus 5.35 GiB (KIT-057). Time is reported: 3,334 s here.
- **What is not claimed:** enrichment (synthetic corpora, network engines; KIT-066 on real data after stage 3); the kit's state-summary rule as a cross-build exactness test (KIT-050 carries it; every counted member since the shipped defaults grew reads `invalid` on summaries with exact components); Windows (KIT-051/052); the lipidomics and metabolomics routes (stage 3, with the audit findings carried there); the report templates' shared logic (KIT-067).
- **Open for the user:** the PCA input threshold (BUG-18, left as recorded), the `qvalue()` fallback and `pivot_wider` kernels (scientific decisions), the step-map unit confirmations still pending from KIT-047/041/053, the replicate rule.
