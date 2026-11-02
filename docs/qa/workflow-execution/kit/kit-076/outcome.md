# KIT-076 outcome: lipid and metab DA per contrast on a shared per-assay fit

Status 2026-11-02: **completed**. Stage 3b, ticket 5 of 10. The numbers are unchanged under the frozen rule; no memory benefit is claimed.

- **Commits** on `performance-v2`: code 812a5943, bindings 57f96a9e, evidence b7a29cca and the commits that follow it.
- **Base:** eec0de61 (KIT-075 closed).
- **Install:** `kit_076` = `core/browser-smoke.FGYOUX` (812a5943), over the kit_075 library. The commits after it change tests and documents only, so the install stands.
- **Ticket:** `tickets/kit/KIT-076.md`, revision 1.1. Plan: [plan.md](plan.md).

## What changed

1. **Kernel halves, moved not rewritten.** `fit<Omic>DaAssayModel` and `apply<Omic>DaContrasts`. `runTestsContrasts*DA` calls both, so `runLipidsDA` / `runMetabolitesDA` return exactly what they did, including the all-contrast `fit.eb`. The assay input, the per-contrast table and the counts also moved out of the orchestrators, so both paths build the same tables.
2. **One declaration per contrast.** `lipid_da_contrast` / `metab_da_contrast`, on one saved `da_fit` per assay (`newWorkflowOmicsDaFitIdentity`, new). `<prefix>_da_output` stays registered for v1 records; the wrapper no longer runs it.
3. **Manifests v2**, keyed by contrast, each contrast holding per-assay tables and fit references. v1 stays declared and readable.
4. **Consumers fed from the manifest.** `readWorkflowOmicsDaResults` rebuilds the kernel's list, so the exports, Glimma, volcanoes, heatmaps, selectors and tables needed no edit.
5. **Record and replay.** One entry per contrast and assay, labelled by contrast, with `failures`. The v1 reading is unchanged. Replay accepts either operation and runs per contrast.
6. **Failure isolation.** A failing contrast fails alone, recorded in the manifest and the record, with a warning in the tabs.

## Results

| Gate | Evidence | Result |
|---|---|---|
| Old vs new, frozen rule | [da-comparison.txt](da-comparison.txt) | **0 violations on all 10 captures.** Largest differences below. |
| Public kernels | same | `runLipidsDA` / `runMetabolitesDA` lists `identical()` old vs new on every capture. |
| Exports | same | TSV bytes identical on `~ 0 + group`; workbook content within the rule everywhere; same file set. |
| Figures | [figure-control.txt](figure-control.txt) | 94 of 110 differ old vs new, and the same 94 differ between two runs of the old build: 0 candidates for a real change. |
| Journeys | [076l](kit076-kit-journey-076l-stage.json), [076m](kit076-kit-journey-076m-stage.json) | Both accepted and exact against `lipid76` / `metab76`. 240 s each; 2.53 and 2.56 GiB peak. |
| Member records | [labels](member-records-labelled.txt), [076l](member-record-replay-076l.json), [076m](member-record-replay-076m.json) | Valid, labelled by contrast and assay, replay exact (10 states); a tampered input is refused. |
| v1 records | [031l](member-record-replay-v1-031l.json), [032m](member-record-replay-v1-032m.json) | The kit_031 / kit_032 members' records replay exact on `kit_076`. |
| Older session | [older-session-open.txt](older-session-open.txt) | kit_075 members with a v1 head: the record is valid, the Open succeeds, and a re-run gives tables identical to what the old build saved. |
| One failing contrast | `tests/testthat/test-workflow-omics-da-contrasts.R` | The others complete with exactly what a run without it gives; the failure is in the manifest and the record; the tabs warn. |
| States and proteomics | [exactness-gates.txt](exactness-gates.txt) | All 10 states identical to stage 3 on both omics. KIT-028 replay exact (16 states); LFQ, TMT, DIA_LIMPA identical to the kit_065 oracles. |
| Lifecycle | [lipid](kit076-kit-lifecycle-lipid-stage.json), [metab](kit076-kit-lifecycle-metab-stage.json) | Both observed, launched one at a time. |
| Tests | [affected-tests.txt](affected-tests.txt) | 39 dependent files: 5,437 pass, 0 fail. Two tests added after that run pass (the new file: 124 pass). No full suite: the next checkpoint is KIT-078. |
| R CMD check | [check-item1.log](check-item1.log) | Baseline: 3 WARNINGs, 3 NOTEs, code note identical to KIT-075's. |

**Largest difference, by design:**

| Design | Contrasts | Lipid | Metab |
|---|---|---|---|
| `~ 0 + group` (`screen`, `multi`, `multi_treat`, `e2e`) | 1 to 3 | 0 | 0 |
| `~ 0 + group + batch` (`nonorth`) | 2 | 7.11e-15 (exports), 6.22e-15 (tables, fits) | 5.11e-15 (exports), 4.44e-15 (tables, fits) |

The ticket expected about 4e-16 on non-orthogonal designs; the observed value is about 15 times larger and still about 140 times inside the 1e-12 limit. The non-zero result also shows the comparison can detect a difference.

## Memory (reported, not claimed)

Anonymous peak in GiB, lifecycle members:

| Omic | Build | Whole slice | Normalisation | Correlation | DA | Report |
|---|---|---|---|---|---|---|
| Lipid | kit_075 | 2.44 | 2.44 | 2.33 | 1.80 | 1.83 |
| Lipid | **kit_076** | **2.47** | 2.47 | 2.35 | **1.88** | 1.91 |
| Metab | kit_075 | 2.50 | 2.50 | 2.33 | 1.93 | 1.97 |
| Metab | **kit_076** | **2.51** | 2.51 | 2.33 | **2.00** | 2.03 |

The peak is still in normalisation, where this ticket changes nothing, so the whole-slice moves are noise. The DA step rose 0.07 to 0.08 GiB on both omics, one contrast each. That may be the rebuild reading its products or may be noise; it was not separated. It stays below the whole-slice peak, and KIT-078 moves this work out of the app process.

## Deviations and discoveries

- **The named inputs could not falsify the claim.** The screen corpora and e2e lanes have one contrast each, so one contrast alone is the same computation as v1. Multi-contrast scenarios were added on the same data before any new code ran (group x batch, with and without treat, and two contrasts on the non-orthogonal design). The rule itself is unchanged.
- **Contrast-term validation moved to the contrast half**, not the fit half as the ticket lists. Item 6 needs an invalid contrast to fail alone, and a saved fit must not depend on the contrasts.
- **Added behaviour: if every contrast fails, the run stops** with their messages and the head keeps its previous result. This keeps the lipid preflight's run-level stop now that it runs per contrast. It is also a change for metab, where such a run used to finish "complete" with empty tables.
- **`fit.eb` in the rebuilt list is a coefficient view** (`list(coefficients = )`), since Glimma reads only its names and row order. The full per-contrast fit stays in the saved product and is what the comparison judged.
- **Recorded, not fixed:** a fit-half failure still drops the assay silently, as v1 did.
- **Ids the fit envelope refuses** (duplicates, NA, empty) keep the fit unsaved and recomputed per contrast; results are unchanged (tested).
- **Older-session route.** The check's copied project takes the adoption path, as an older session file does. Reopening an owner at a v1 head is not separately exercised on a real store; the module tests cover the reopen, and nothing on that path reads the manifest.
- **Tests changed on purpose.** The KIT-073 frozen dump now excludes the new declaration and compares tables, counts and warnings, not the v1 manifest; those still match the dump frozen before KIT-073. The module harnesses mock the wrapper, not the kernel.
- **Pre-existing defect, left alone: per-contrast heatmap PNGs are never written, on both omics.** `generate<Omic>DAHeatmap` returns a list and the exports call `ComplexHeatmap::draw` on it (`R/func_lipid_da_export.R:265`, `R/func_metab_da_results_output.R:267`). The error is caught and logged, so the run looks fine. Old and new builds both log it 34 times in the capture, and kit_075's journey members have only `all_heatmaps_combined.pdf` in `Heatmaps/`. Noted in KIT-077, which edits those lines.

## Backward compatibility and rollback

- `runLipidsDA` / `runMetabolitesDA` unchanged for callers; exported tables unchanged; v1 manifests and records readable.
- Rollback: `git revert` the ticket's commits. Stores written with v2 manifests are then readable only through their files: a reverted build has no v2 declaration, so the record falls back to treating the manifest as proteomics'. Re-run the DA after a revert.

## Handoff

From here the exactness set's oracles are `lipid76` / `metab76`. KIT-077 gets per-contrast DA results and the heatmap defect. KIT-078 moves the per-contrast operations to the worker, and the failure path above must survive that move. KIT-061 reads the v2 manifests.
