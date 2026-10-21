# KIT-021: LFQ and TMT journeys: exactness and completion

Status 2026-10-21: **completed** (ticket revision 2.2; protocol v1.5 section 7). Both counted single journey members ran on the bound build `kit_021` (`core/browser-smoke.l7OM2e`, source `88df9314`) under the user's "GO" of 2026-10-20 and were **accepted**: complete, every state exact against the in-process oracles recorded on the same build and corpus, peaks and resident ends recorded, no ratio claimed.

## Counted members

| Member | Corpus | Decision | Elapsed | Peak anonymous | Resident end | Mismatches |
|---|---|---|---|---|---|---|
| `kit:journey-021:candidate:1` ([stage](kit-journey-021-stage.json)) | `lfq_300` (MaxQuant proteinGroups, 300 samples × 8,000 groups + 50 contaminant + 50 decoy rows) | accepted, exact | 2,994 s | 3.59 GiB | 2.61 GiB | none |
| `kit:journey-021t:candidate:1` ([stage](kit-journey-021t-stage.json)) | `tmt_300` (Proteome Discoverer protein export, 300 samples in 30 TMT10 plexes × 8,000 proteins) | accepted, exact | 8,832 s | 4.04 GiB | 3.02 GiB | none |

Per-step elapsed (s) and peak anonymous memory (GiB), LFQ / TMT: import 14 / 15 (1.13 / 1.16), design 552 / 578 (1.39 / 1.56), accession 229 / 446 (3.25 / 3.32), protein intensity 19 / 34 (2.20 / 2.32), duplicates 32 / 61 (2.50 / 2.62), protein replicate 12 / 24 (2.62 / 2.62), normalisation and RUV 1,821 / 7,392 (2.51 / 3.71), correlation 34 / 41 (1.88 / 2.09), session export 7 / 8 (1.69 / 1.89), DA 175 / 144 (3.59 / 4.04). The automatic RUV sweep dominates both (the TMT sweep at 300 samples took two hours; the same kernels run in process took 5,563 s when the oracle was recorded). Elapsed times are reported, not judged. The counted TMT member ran while `R CMD check` and the KIT-054 install (about 15 minutes of compile) shared the host; its peaks are its own slice's.

## Acceptance

1. **Both members complete with every state verified against its oracle** (ten steps each, `populations_complete`, `oracle_mismatches` empty; the member's rendered volcano page and QC plot images are reported as additional products, as on the DIA route). Peaks, per-step peaks and `resident_end_bytes` recorded above; no ratio.
2. **The LFQ contrast tables equal the control build's in-memory output** ([compare-check-lfq_300.json](compare-check-lfq_300.json)): [control-check.R](control-check.R) ran the control build's (`kit_036`, `6c9c4901`) in-memory kernels on `lfq_300` from import through DA (1,091 s; the control build's reader keeps the flagged rows its sampled type guess misses, dropped by prefix in the check and recorded); [da-tables.R](da-tables.R) read the owner path's saved DA tables back from the oracle run's store; **3 of 3 table hashes equal** (`da_proteins_long`, `da_proteins_wide`, `contrasts_results_table`, 8,000 rows, contrast KO vs WT). The counted member reproduces the oracle run hash-identically, so the chain owner path → oracle → control is closed. RUV on the synthetic corpus resolved k = 60 at 1 % (no structured unwanted variation to find; identical on both paths).
3. Elapsed times reported, never judged (above).

## What changed (harness and corpus; one reader repair)

- **Corpora** (`tools/profiling/kit_protein_corpus.R`, seeded): `lfq-300` in MaxQuant proteinGroups format with the dotted headers the fixtures carry, flagged rows interleaved with one contaminant and one decoy leading the table; `tmt-300` in Proteome Discoverer protein-export format (`Abundance: F<plex>: <channel>, <label>`); the two e2e fixtures as `lfq-fixture` / `tmt-fixture`; manifests with sha256s; bound through `bind_corpus` with a file list. Design samples are the import's sanitised run names.
- **Protein route on the kit**: `kitJourneySteps("protein")` (import → design → accession → protein intensity → duplicates → replicate → normalisation and RUV → correlation → session export → DA), the protein step map (`kit-001/step-map-protein.json`, bound by hash), suffixed ticket keys (`021`, `021t`) for whole-journey single members with in-process oracles, the coordinator taking the corpus format and a per-corpus normalisation override (the five-protein fixtures skip RUV: no valid negative-control percentage), the driver opening the LFQ/TMT QC tabset (`qc_tabs_tmt`) and waiting for the normalisation text of the bound RUV mode, the parsed bindings keeping the harness block.
- **Oracle recorder** ([oracle-record.R](oracle-record.R)): an in-process run of the bound build through the GUI's own entry points, the driver's population shape (every step snapshotted, per-step products by first appearance, rendered assets left out), the DA contrasts from the design builder's contrast table.
- **Reader repair** (`R/func_prot_import_readers.R`, `88df9314`): `importMaxQuantData` types `Potential.contaminant` and `Reverse` explicitly. Found by the corpus check: vroom guesses a column's type from a sample of rows across the file, so a sparse "+" column whose sample carried no mark was read as logical with every mark lost, and 50 decoys among 8,100 groups survived the import, collapsing to one unnamed protein at the accession cleanup and failing the normalisation on both builds. Regression test with sparse marks (the control build keeps both flagged rows; the repaired reader drops them). KIT-054 found the companion defect (the reader knows only the dotted header names) and repairs it.

## Lifecycle and corrections (all archived under `<proof>/kit/runtime/`)

Fixture lifecycle members ([LFQ](kit-journey-021-lfq-fixture-lifecycle-stage.json), [TMT](kit-journey-021t-tmt-fixture-lifecycle-stage.json)): accepted, exact, 128 s and 176 s. Six earlier attempts were instrumentation failures, each repaired in the harness or the recorder and archived (`attempt1-harness` route helper outside the helpers environment; `attempt2-harness` the DIA QC tabset id on an LFQ owner; `attempt3-harness` the automatic-mode completion text with RUV skipped; `attempt4-oracle-type` the GUI's integer percentile against the recorder's double; `attempt5-oracle-shape` the RUV percentage type and the rendered volcano page in the oracle; `attempt6-oracle-contrast` the DA contrast from a pairwise combination instead of the builder's table). No counted member was retried.

## Checks, install and suite

- Full suite at the checkpoint before the counted launches, on the tree at `88df9314`: 344 files, **30,261 pass, 0 fail, 0 error, 140 skip** ([summary](full-suite-88df9314.txt)); the earlier checkpoint run (2026-10-20, before the reader repair) 30,258 pass, 0 fail.
- `R CMD check --no-manual --no-tests` at `03ca18b6` and at `88df9314` (clean worktrees): 3 WARNINGs, 3 NOTEs, the baseline ([log](check-88df9314/r-cmd-check.log)).
- Installs: `kit_021` = `core/browser-smoke.0WrMZI` (`03ca18b6`), replaced by `core/browser-smoke.l7OM2e` (`88df9314`) after the reader repair; bindings regenerated; the control stays `kit_036`.
- Harness tests: `test-workflow-kit-harness.R` 364 then 378 pass (with KIT-054's keys); reader tests `test-prot-import-maxquant-direct-shared.R` 18 pass, `test-module-ci-prot-import.R` 90, `test-prot-01-import.R` 15, `test-prot-01c-import-module-contracts.R` 410.

## Deviations and claim limits

- Protocol v1.5 governs: exactness is against in-process oracles on the same build (the control build cannot produce the owner path's states); the control-build cross-check closes the chain for the LFQ DA tables. The TMT DA tables were not cross-checked against the control build in this ticket (not an acceptance criterion); KIT-054 runs that check on `tmt_300`.
- The synthetic corpora carry no unwanted-variation structure, so RUV's automatic sweep picks the extreme of its range; exactness holds regardless. Real-data behaviour is KIT-054's.
- The reader repair is an R change in a harness ticket, made because the corpus check reproduced a defect that silently keeps decoys on real MaxQuant output; it is one line of column typing with a test, declared here.

## Hand-off

KIT-054 (next): the real MaxQuant export beside the in-memory path on `kit_036`, with the reader's spaced-header repair; it reuses this ticket's recorder, cross-check scripts, corpus binding and the counted TMT member. KIT-049 onward: the route-parameterised driver and the single-journey binding form (`whole_journey`, suffixed keys) are the form KIT-053 takes for the `DIA_LIMPA` journey.
