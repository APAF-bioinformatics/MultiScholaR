# Audit: the DIA_LIMPA route's full path on the kit (2026-10-23)

Mode: `audit_and_repair` (user, 2026-10-23: "audit the full path and get this to the production standard I expect"). Scope: every step of the `DIA_LIMPA` journey as the KIT-053 counted member ran it (`kit:journey-053:candidate:1`, `kit_053` = `BTcMbt`, source `81a69944`, the 300-sample DIA-NN calibration corpus), against the DIA route's counted members on the same corpus (`kit-journey-036`: control `kit_018` 4,912 s, peak anonymous 5.66 GiB; candidate `kit_036` 4,740 s, 5.69 GiB). The question: does each step use the kit's on-disk mechanisms, and where does the route's peak come from. Findings are published here before any repair; repairs are KIT-057's.

Evidence: the stage summaries (`kit-053/kit-journey-053-stage.json`, `kit-036/kit-journey-036-stage.json`), the managed actions' completion records in each member's store (`.workflow-export/wfx-*/completion.rds`: the worker unit's `MemoryPeak`), the components hydrated on the counted store (measured in process on `hpvhkD`), and standalone kernel runs on the counted store's data (`VmHWM` of a fresh R process; synthetic corpus).

## 1. The path, step by step

| Step | Operation | Unit (step map) | Variant | Worker | Hydrates | Worker peak | Time | Member peak (anon) | DIA route's counterpart (worker peak, time) |
|---|---|---|---|---|---|---|---|---|---|
| import | import | store | — | app | DIA-NN report (5,846,578 rows) | — | 42 s | 2.88 GiB | same (shared) |
| design | design_builder | store | — | app | — | — | 102 s | 3.30 GiB | same |
| qvalue | peptide_qvalue_filter | feature group | `streamed` | worker | batches | 1.82 GiB | 169 s | 4.88 GiB | 2.07 GiB, 295 s (shared) |
| rollup | precursor_rollup | complete sample | `p14` | worker | batches | 3.58 GiB | 43 s | 5.89 GiB | 3.78 GiB, 58 s (shared) |
| **peptide_log2** | peptide_log2_matrix | row | **`whole` only** | worker | `peptide_data` 5,521,570 × 9 (320 MiB) + design | **4.09 GiB** | 71 s | 6.62 GiB | no counterpart (the DIA route pivots inside IQ) |
| **peptide_accession** | peptide_accession_cleanup | feature group | **`whole` only** | worker | table (320 MiB) + matrix 60,000 × 300 (142 MiB) + design | **8.95 GiB** | 409 s | **11.42 GiB** | protein_accession_cleanup on the 2,500-row protein table: 0.42 GiB, 21 s |
| **dpc_quant** | protein_dpc_quant | whole matrix | **`whole` only** | worker | table + matrix + design | **8.38 GiB** | 5,564 s | 10.89 GiB | protein_iq_rollup (`whole`): 3.44 GiB, 233 s; peptide_imputation (`partitioned`): 2.70 GiB, 666 s |
| protein_normalisation_limpa | protein_between_sample_normalisation | whole matrix | `whole` | worker | protein table 2,500 × 301 (6 MiB) | 1.33 GiB | 21 s | 3.67 GiB | normalization-ruv in process |
| peptide_count | protein_peptide_count_filter | row | `whole` | worker | protein table + owned EList | 1.22 GiB | 18 s | 3.61 GiB | protein filters 0.5–0.7 GiB |
| correlation | pearson-pairs, pearson-selection | whole matrix | in process | app | protein table | — | 37 s | 3.09 GiB | same |
| ruv_covariates | protein_ruv_covariates | whole matrix | in process | app | protein table | — | 250 s | 2.99 GiB | normalization-ruv in process |
| session_export | export | store | — | app | — | — | 11 s | 2.40 GiB | same |
| da | DA-contrast-output, da_publication | whole matrix | in process + worker | app + worker | protein table + owned EList | 1.28 GiB (publication) | 141 s | 3.94 GiB | same shape (limma) |

The member's peak is the application process plus the active worker (`whole_service_peak_bytes` 11.59 GiB at the accession step; the application's own resident set about 2.5 GiB throughout). The route's shared steps and its protein stage are on the kit's mechanisms or on 6 MiB tables; the three KIT-047 peptide steps are not.

## 2. Where the memory goes (standalone kernel runs on the counted store's data, `hpvhkD`)

| Kernel | Hydrated input | Process after hydration | Process peak (VmHWM) | Time | Reading |
|---|---|---|---|---|---|
| `calcPeptideMatrix` + `log2TransformPeptideMatrix` on the rollup state | table 320 MiB | 1.27 GiB | **4.28 GiB** | 53 s | the whole-table `pivot_wider` costs about 3 GiB of working set over a 320 MiB input: the row unit is paid for as one matrix |
| `chooseBestProteinAccession` (peptide method) on the log2 state | table + matrix | 1.23 GiB | **10.41 GiB** | 369 s | the FASTA join to every row, the whole-table `group_by`/`summarise` and the re-pivot each copy the 5.5-million-row table |
| `limpa::dpc` (the fit) on the matrix alone | matrix 142 MiB | 0.68 GiB | **0.74 GiB** | 19 s | the global statistic is cheap: the fit needs the matrix only |
| `limpa::dpcQuant` on the first 500 of 2,500 proteins (their rows of the matrix) | matrix rows of 500 proteins | 1.43 GiB | **6.87 GiB** | 1,924 s | 18,525 matrix rows: the cost is inside limpa's per-protein loop (`peptides2Proteins`, run with the kernel's `chunk = 1000`), not in the fit or the hydration; the chunk size is the first lever (measured next: 200 proteins at chunk 1,000 and 100, memory, time and a digest of the result) |

## 3. Findings

| ID | Severity | Observed | Claim affected | Smallest sufficient repair | Blocks |
|---|---|---|---|---|---|
| F-01 | high | `peptide_log2_matrix` (unit row) runs only a `whole` variant: the worker hydrates the 5.5-million-row table and pivots it in one `pivot_wider`; worker peak 4.09 GiB (standalone 4.28 GiB) for a 142 MiB result. The DIA route's row-unit steps stream batches through the bounded table writer (`persistDIAFilteredSamples`, `computeDIAReferenceQvalue`). | The route's memory bound; the step map's "row" unit is declared but not realised. | A `streamed` production variant: batches of whole peptide rows (the feature key's rows) pivoted per batch and appended through the table writer; `whole` kept as the control arm; hash-exact. | execution of a bounded journey |
| F-02 | high | `peptide_accession_cleanup` (unit feature group) runs only `whole`: the FASTA join, the duplicate aggregation and the re-pivot over the whole table; worker peak 8.95 GiB, the member's journey peak (11.42 GiB). The DIA route's feature-group steps run on `partitionByFeatureGroup` (batches of whole groups, admission recorded and refused above the cap, merge rules). | As F-01; the route cannot run on a 16 GiB host beside the application with margin. | A `partitioned` production variant in two passes on the partitioner: the per-row accession choice (needs only the row's candidates and the FASTA scores), then the aggregation by the chosen accession in batches of whole chosen-accession groups, the pivot per batch; `whole` kept as the control arm; hash-exact. The kernel's helper and aggregation called as they are, per batch. | execution |
| F-03 | high | `protein_dpc_quant` hydrates the 320 MiB table (the kernel reads it only to rebuild the feature-key map and to count peptides per protein) and runs `dpcQuant` whole: worker peak 8.38 GiB, 5,564 s (80 % of the journey). The fit (`limpa::dpc`) alone peaks at 0.74 GiB in 19 s; `dpcQuant` estimates four hyperparameters on the whole matrix once, then quantifies protein by protein (`peptides2Proteins`, internal `chunk` of 1,000 proteins): 500 proteins alone peak at 6.87 GiB, so the loop's working set follows the chunk, which the kernel fixes at 1,000 (`[limpaDpcQuant] chunk`, a declared parameter). | The route's memory bound and its time; the step is the reason for the 10,800 s worker deadline. | Measured 2026-10-24 (200 proteins, 13,580 matrix rows, the same proteins): chunk 1,000 → 1,707 s, VmHWM 6.48 GiB; chunk 100 → 1,707 s, 5.94 GiB; result digests identical. **The chunk is not the lever**: the working set is the per-protein fit itself (`peptides2Proteins` with standard errors on 300 samples), which the kit cannot change without changing limpa. What remains to the kit is the hydration: the feature-key map and the peptide counts per protein from a streamed two-column prelude instead of the 320 MiB table, the matrix alone hydrated; the global statistics once as `dpcQuant` computes them; `whole` kept as the control arm; hash-exact. | execution; the deadline decision revisited |
| F-04 | medium | The worker deadline of 10,800 s (`81a69944`) was set from the `whole` variant's time. With F-03 repaired the time is unknown. | The per-step deadline's justification. | Measure the bounded variant's time; keep or drop the entry on the measurement (KIT-057 item 5). | decision |
| F-05 | low | The step-map rows of the three steps record the unit without a variant; `in_journey` is `false` on every route row although the rows are in the journey (noted by the format audit of 2026-10-22). | The step map's reading. | The rows gain `variant` and the production mechanism; `in_journey` corrected. | documentation |
| F-06 | low | The application process's resident set (about 2.5 GiB) is the same on both routes and is outside this audit's question. | — | none here (KIT-044/048 territory). | — |

Not findings: the shared steps (import, design, q-value, rollup, correlation, export, DA) and the route's protein stage run on the kit's mechanisms or on 6 MiB tables with worker peaks of 1.2–3.6 GiB; the oracle, replay, record and verification machinery held (the counted member was exact).

## 4. Decision design for the repair (KIT-057)

- Exactness is the gate: each production variant hash-identical to `whole` on the screen corpus in process and on the full corpus through one counted control/candidate pair (control `kit_053`, the whole variants; candidate `kit_057`), the kit's pair rule (KIT-005/017).
- Memory targets, per worker on the full corpus: log2 ≤ the DIA q-value worker + 1 GiB (≈ 3 GiB); accession ≤ the DIA route's worst peptide-group worker (3.78 GiB); DPC ≤ 4 GiB; the journey peak ≤ the control's and ≤ 5.7 GiB.
- Time is reported, never judged; no ratio.
- Terminal outcomes: accepted (exact and within the targets); rejected (exact but above a target: the step is named, the next mechanism is a new ticket); invalid (not exact: the variant is withdrawn, `whole` stays the default).

## 5. States

`implementation_verified`: no (nothing repaired yet). `execution_complete`: the audit's measurements, yes. `decision_supported`: the findings and targets above. `promotion_authorized`: no. `claim_scope`: the 300-sample synthetic DIA-NN corpus on the kit host; the per-worker peaks are systemd `MemoryPeak` of the worker unit; the standalone peaks are `VmHWM` of a fresh R process on the same data.
