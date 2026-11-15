# GUI journey: audit and readiness (2026-11-15)

Requested by the user on return: fix the cleanup, check the work against best practice and the intended architecture ("the GUI just renders the work done by R operations, on disk"), give executive calls on the parked decisions (`parked-20261113.md`, last section), and judge whether the GUI journey is complete for all omics and ready for the showcase datasets.

Two sources: an independent review of the code this journey changed (`git diff 63ce6d6c..e8cbd388 -- R/`, a separate reviewer, read-only), each finding then verified by the implementer; and the latest journeys on the current build (kit_933), for where the work actually runs.

## 1. Where the work runs (journeys on kit_933)

Every science step on every route is dispatched to the worker (the journeys' dispatch records: 17 operations on DIA, 12 on DIA_LIMPA, 9 on LFQ/TMT, 18 on lipid, the same on metab, 14 on trx). What the app process itself does is its own memory rise per step:

| Route (corpus) | Import | Design | QC steps | Normalisation | DA | Report |
|---|---|---|---|---|---|---|
| Lipid, metab (screen) | +0.17 | +0.01 | +0.08 to +0.11 | **+0.00** | +0.11 to +0.16 | +0.01 to +0.05 |
| Transcriptomics (screen) | +0.18 | +0.02 | +0.00 to +0.03 | **+0.00** | +0.00 | +0.10 |
| DIA_LIMPA (screen) | +0.37 | +0.19 | +0.04 to +0.18 | +0.02 | +0.04 | +0.01 |
| LFQ (fixture) | +0.18 | +0.03 | +0.06 to +0.17 | **+0.58** | +0.00 | +0.01 |
| DIA, IQ (full, 5 million precursor rows) | **+1.95** | **+1.33** | +0.00 to +0.20 | **+1.00** | +0.37 | +0.27 |

GiB, application process only (`controller_anon`), peak above the step's start.

**Reading.** Lipid, metab and transcriptomics meet the rule: the app renders and holds references. Proteomics meets it for QC and for the DA computation, with these exceptions:

| # | Exception | Size | Status |
|---|---|---|---|
| E1 | **Proteomics normalisation draws its QC figures in the app** (the S4 is hydrated and PCA, RLE, density and correlation computed there: on tab entry, in every publication half, on an aesthetics change, and the skip-correlation path wholly). Lipid and metab dispatch the same figures to the worker. | +1.0 GiB on the full DIA corpus; the app's highest point of the journey (2.64 GiB) | pre-existing (KIT-111 D3); KIT-120 mounted it on the shared frame unchanged. **KIT-122.** |
| E2 | **Import and design run in the app** (the report is parsed and stored; the design step builds the first S4). | +1.95 and +1.33 GiB on the full DIA corpus | by design since KIT-048; released afterwards (KIT-121). Not ticketed. |
| E3 | **Enrichment runs in the app and keeps plot objects** in tab state and in the stored state (about 180 MiB per contrast). | per contrast | KIT-119 (decision 9). |
| E4 | Proteomics' report and session export, and the S4 compatibility download, run in the app. | +0.27 GiB (report) | recorded; small. |
| E5 | The proteomics tab keeps small display objects in state between actions (a RUV result with a grob, three plot objects). No assay is held. | small | recorded; goes with KIT-122. |

## 2. The independent review's findings

| # | Sev | Finding | Verified | Disposition |
|---|---|---|---|---|
| F1 | medium | KIT-121 did nothing off glibc Linux: `reclaimWorkflowMemory` returned "unsupported" before collecting. | yes (`src/workflow_memory.c`, the pinned test) | **fixed**: every platform collects; only the page release is skipped (status `collected`). Test failed first. |
| F2 | medium | The RUV summary's mode and grouping are lazy arguments read once, so after a second run in another mode the summary still names the first. Pre-existing, carried onto the frame. | yes | KIT-122 (display defect; the right source is the result, not the inputs) |
| F3 | medium | = E1 above. | yes, and measured | KIT-122 |
| F4 | low | The frame looked up `omicGuiConfig(omic_type)` even for an omic with its own server, so the exported `mod_prot_norm_server` aborted for an `omic_type` label other than "proteomics". | yes | **fixed**: no lookup in that branch; such a group must name its QC directory and assay (validator). Test failed first. |
| F5 | low | The RUV tab's three static outputs are created through a `renderUI`. | yes | KIT-122 (a library entry) |
| F6 | low | One state object carries two vocabularies; two fields are written only by dead code. | yes | KIT-122 |
| F7 | low | Every session logs "=== METABOLOMICS NORMALIZATION MODULE STARTED ===", now proteomics too. | yes (a lipid test pins the line) | KIT-122 |
| F8 | low | The no-op passed to the observer registrar leaves an inert observer. | yes | KIT-122 |
| F9 | low | The reclaim reports "busy" while an owner has worker clusters. Not shown to occur in the app. | not checked: the journeys show the release after import and design; whether a later request (after normalisation) is carried out once worker groups exist was not measured | to measure in KIT-122's runs |
| F10 | low | Reset replaces two summary outputs with constant text, so they stay at the placeholder after a later run. Pre-existing, outside the diff. | yes (by reading) | KIT-122 |
| conv. | low | Two short anonymous functions written `function(...)`; one roxygen `context` line incomplete. | yes | **fixed** |
| records | | `kit-121/outcome.md` said the reclaim follows "import, design, QC or normalisation" for every owner; `sites.md` listed a live function as retired; `kit-120/outcome.md` omitted a display change. | yes | **corrected** (below and in those files) |

Reviewer's verdicts: the diff **conforms to the architecture with named exceptions** (no new data holding, no new in-app compute); KIT-121 is **safe to keep enabled** (guarded, after delivery, failures caught, inert where no request site exists). No double registration of observers or outputs was found.

## 3. Record corrections

- **KIT-121's reach is narrower than its outcome said.** The reclaim runs after proteomics' import, the design builder's save, the q-value display and a normalisation publication. The lipid, metab and trx owners have no request site (their apps already sit flat, section 1), and the DIA_LIMPA normalisation pass-through does not request it.
- **KIT-120:** `registerProtNormRenderOutputs` is live, not retired. The stage plots changed from 400 px to the frame's 300 px image slots. Proteomics sessions now log the frame's start line (F7).

## 4. Is the GUI journey complete, and is it ready for the showcase datasets?

**Complete: yes, for all four omics, as a GUI.** All 28 omic x stage cells run on the shared frame. The full suite is clean in both modes (36,775 passed, 0 failed), the seven browser lanes pass, one journey per route is exact under the current harness, and the architecture rule holds for every science step, with the exceptions of section 1, none of which gives a wrong result.

**Ready for the showcase datasets: yes for a shakedown now; not yet for the numbers that go in the paper.**

- **Start the shakedown now.** Real data will find what synthetic data cannot, and three things have never run on real data through this GUI: enrichment in a whole session, memory at real scale, and RUV on real missingness. Finding those early is worth more than waiting.
- **Do not take publication numbers until these land**, because each changes results on real data:
  1. **The DIA-NN zero rule** (decision 1): real reports contain zeros.
  2. **Lipid and metab zeros as missing** (decision 7): the largest effect measured (14 to 49 points of power).
  3. **RUV control selection with missing values** (P-03): on real LFQ data RUV may refuse every setting today. Until it lands, "Skip RUV" works.
  4. The enrichment fixes (decision 8) for any dataset that uses Reactome by organism name or ids with a colon.
- **Expect during the shakedown:** the app at about 2.6 GiB on a DIA dataset of the full corpus's size (E1); about 180 MiB kept per enrichment contrast (E3); the "Use Precursor.Normalised" box doing nothing (decision 4).

**Order recommended:** KIT-118 (zeros, RUV controls, the label and record items) while the shakedown runs; then KIT-119 (enrichment) and KIT-122 (proteomics normalisation figures on the worker); then the publication runs (KIT-088) and the real-data gate.
