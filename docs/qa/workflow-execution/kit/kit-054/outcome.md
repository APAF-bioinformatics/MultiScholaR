# KIT-054: Real MaxQuant LFQ dataset (and the synthetic TMT corpus): the owner path beside the in-memory path

Status 2026-10-21: **completed** (ticket revision 1.2; user request and decisions of the same day). The user's real MaxQuant export ran the whole LFQ journey on the on-disk owner path (the KIT-021 build with the two reader repairs) and on the legacy in-memory path (the control build `kit_036`); the synthetic TMT corpus ran the legacy lane beside KIT-021's counted member. All three counted members were **accepted**; the lanes' tables are equal; no science difference was found. The data handling rule held throughout: no row, identifier or sample name appears in any committed file, log excerpt or message (every evidence file scanned before commit; the fingerprints of the sensitive corpus carry column lists as a count and a hash).

## Members (counted, ceiling 14,400 s each)

| Member | Lane | Build | Corpus | Decision | Elapsed | Peak anonymous | Resident end |
|---|---|---|---|---|---|---|---|
| `kit:journey-054:candidate:1` ([stage](kit-journey-054-stage.json)) | new (owner path) | `kit_054` = `core/browser-smoke.elSszS` (`cc661b14`) | `kc_lfq` (real MaxQuant export, 6 runs, 2 groups) | accepted, **exact** against the in-process oracles (10 states) | 361 s | 2.18 GiB | 1.51 GiB |
| `kit:journey-054l:candidate:1` ([stage](kit-journey-054l-stage.json)) | legacy (in memory) | `kit_036` = `core/browser-smoke.h6QeDJ` (`6c9c4901`) | `kc_lfq_dotted` (the same export, header renamed) | accepted (completion and measures; no store) | 319 s | 3.60 GiB | 2.64 GiB |
| `kit:journey-054m:candidate:1` ([stage](kit-journey-054m-stage.json)) | legacy (in memory) | `kit_036` | `tmt_300` (synthetic PD export, 300 samples) | accepted (completion and measures) | 7,697 s | 6.57 GiB | 4.66 GiB |
| `kit:journey-021t:candidate:1` (KIT-021, reused) | new (owner path) | `kit_021` = `core/browser-smoke.l7OM2e` | `tmt_300` | accepted, exact | 8,832 s | 4.04 GiB | 3.02 GiB |

Elapsed times reported, never judged (the TMT in-memory RUV sweep 6,802 s against the owner path's 7,392 s; the legacy LFQ lane ran a little faster than the owner path on six samples). No ratio is claimed.

## Real-data results (aggregate only)

The export holds 1,455 protein groups over 6 runs (3 per group, no technical replicates) with 51.8% of intensities missing; the import drops the flagged rows the reader is meant to drop (none survive). Through the journey: 1,455 after accession cleanup, 390 after the intensity filter at the bound defaults (33.333% missing per group, 1st percentile), 390 after duplicate removal and the replicate filter, RUV automatic chose k = 3 at 3% controls, 384 after RUV, 384 after correlation filtering (none removed at 0.5), one contrast F vs E with 384 DA rows. The same counts on both lanes.

## Lane equality

- **DA tables** (acceptance 2): the new lane's three DA tables read back from its oracle store equal the control build's in-memory kernels' tables on the dotted copy, **3 of 3 hashes** ([compare-check-kc_lfq.json](compare-check-kc_lfq.json), [control-check-kc_lfq_dotted.json](control-check-kc_lfq_dotted.json), [da-tables-kc_lfq.json](da-tables-kc_lfq.json)). For TMT, KIT-021's owner path on `tmt_300` against the control build in memory: **3 of 3 hashes**, 7,922 rows, RUV k = 2 at 4% ([compare-check-tmt_300](../kit-021/compare-check-tmt_300.json)).
- **Published files** ([lane-compare-kc_lfq.json](lane-compare-kc_lfq.json), [lane-compare-tmt_300.json](lane-compare-tmt_300.json)): on both corpora 49 files are written by both lanes; **34 are byte-identical**, and that set holds every table (the DA long table and workbook, the normalised matrix before RUV, the RUV-corrected matrix, the replicate-filter output, the significance table, the volcano PDF and PNG, the design, the cleaned long table, the contrasts, the annotation caches). The 15 differing files are the same on both corpora and are classified from their paths and sizes: twelve QC plot renderings a few hundred bytes apart (correlation, PCA and RLE at three stages, the cancor plot, the composite figure and the significance-count PNG: rendered from matrices that are equal), the project `config.ini` (the new build writes the worker-cluster specification portably where the old wrote environment placeholders), `qc_params.RDS` and `ruv_optimization_results.RDS` (serialised R objects embedding plot environments). The new lane writes 155 files the legacy lane cannot (the session export archive and its bookkeeping, per-step filtering QC, the interactive volcano page, display products); the legacy lane 3 files of its own session format. **No science difference.**
- The user's outputs: the two LFQ lanes' project outputs were copied (store excluded) to `~/Projects/APAF/Projects/KC/kit-054-outputs/<lane>/`, the new lane's `results_summary` filled and its report rendered as the Summary tab does ([fill-summary.R](fill-summary.R)); the user read the report and the template defects found there are KIT-055's (closed the same day; the report re-rendered on that build).

## What the ticket changed

- **Reader** (`R/func_prot_import_readers.R`, `R/func_prot_import_detection.R`, `f87e5978`): `importMaxQuantData` reads genuine MaxQuant output (spaced column names "LFQ intensity <experiment>", "Potential contaminant", "Majority protein IDs" mapped to the dotted names the fixtures carry; the detector matches its markers against spaced headers); regression test with a spaced-header twin of a dotted table. Found by the real export: the shipped reader stopped with "No MaxQuant intensity columns found" on both builds (a pre-existing defect, not a refactor regression). The legacy lane, whose build keeps the old reader, ran the header-renamed copy (`kit-054/dotted-copy.R`: the header line renamed, the rows the original bytes, 41 of 85 columns renamed).
- **Binding** (`kit-054/bind-real-corpus.R`, `bindings-generator.R`): a real corpus bound by path and sha256 from its own directory with a manifest of counts only (`synthetic = false`, `sensitive = true`); the fingerprints of a sensitive corpus redact column lists (the recorder and the generator).
- **Harness** (`workflow_kit.R`, `workflow_kit_coordinator.R`): a corpus names its FASTA and design files in the binding; a protein-route single member bound without oracles is the in-memory lane (no store, no populations, exactness not judged). Three keys `054`, `054l`, `054m`.
- **Scripts**: `lane-compare.R` (publication files by hash, counts and paths only), `fill-summary.R`; the KIT-021 recorder, control check, DA table reader and join take the FASTA and design names from the binding, an evidence directory and a control corpus key.

## Checks and installs

- `R CMD check --no-manual --no-tests` at `f87e5978`: 3 WARNINGs, 3 NOTEs, the baseline ([log](check-f87e5978/r-cmd-check.log)). Install `kit_054` = `core/browser-smoke.elSszS` (`cc661b14`). Reader tests: `test-prot-import-maxquant-direct-shared.R` 24, `test-module-ci-prot-import.R` 90, `test-prot-01-import.R` 15, `test-prot-01c-import-module-contracts.R` 410, the detection characterisations 8 and 79; harness tests 382.
- Attempts archived under `<proof>/kit/runtime/kit-journey-054*.attempt1-harness`: the first launches failed on two harness slips (the coordinator looked for the FASTA under the synthetic corpora's fixed name; the decision demanded store populations of an in-memory member), repaired and relaunched; the stopped legacy TMT run is the second of those.

## Deviations and claim limits

- No real TMT data exists (user): the TMT side-by-side is on the synthetic `tmt_300` with KIT-021's counted member as the new lane.
- The legacy lane is the in-memory path on `kit_036`, which carries the series' two intentional scientific changes (declared technical-replicate blocks for RUV-III and DA; the quantile ceiling removed on the DIA route), so the comparison isolates the architecture, not the whole history since the pre-series package.
- The intensity filter removed three quarters of the groups at the bound defaults on six samples with half the intensities missing; a science observation for the user, not a kit finding.

## Hand-off

KIT-048 (the in-memory path retired) has this side-by-side as its regression evidence on real data. KIT-050 takes the real-corpus binding form. The reader's two repairs (sparse flags, spaced headers) ship with `kit_054` and every later build.
