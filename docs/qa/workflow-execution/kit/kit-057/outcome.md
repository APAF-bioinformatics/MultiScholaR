# KIT-057: DIA_LIMPA peptide stage on the kit's memory mechanisms

Status 2026-10-24: **in progress**. Audit: [AUDIT-20261023-limpa-path.md](AUDIT-20261023-limpa-path.md) (findings F-01..F-06, published before repair). Control build `kit_053` = `core/browser-smoke.BTcMbt` (the KIT-053 counted member: journey peak 11.42 GiB anonymous on the 300-sample corpus; the DIA route's 5.7 GiB).

## Repairs

| Finding | Step | Mechanism | Commit | Exact against `whole` | 300-sample worker peak (own process, VmHWM) | Time |
|---|---|---|---|---|---|---|
| F-01 | `peptide_log2_matrix` | `streamed` variant on the feature-group partitioner: a prelude over the identity columns (the feature-key map and the sample order as the kernel builds them on the whole table), the unchanged kernels per batch of whole protein groups, the matrix rows placed by name, the kernel's updated quantity column streamed as a change map over the source rows | `10594b77` | screen and full corpus: matrix, table, design, config, args, slots identical ([screen](variant-check-peptide_log2-tree.json), [full](variant-check-peptide_log2-full-tree.json)) | **4.84 → 2.10 GiB** | 71 → 219 s |
| F-02 | `peptide_accession_cleanup` | `partitioned` variant: the method keeps the DIA-NN group key and replaces only `Protein.Ids`, and its aggregation groups by (sample, group, sequence), so one pass by protein group suffices; the helper's mapping computed once on the distinct group/candidate pairs (identical to the whole-table mapping, checked); the unchanged method per batch; the aggregated rows saved in the method's sorted order through the order index; the matrix rows placed in the whole pivot's order | `10594b77` | screen and full corpus identical ([screen](variant-check-peptide_accession-tree.json), full: see below) | **10.54 → 2.38 GiB** | 384 → 680 s (the first measurement, 3.59 GiB and 2,381 s, had a quadratic prelude that accumulated every page's feature triples; the prelude now keeps one row per feature per page, re-measured exact: [full](variant-check-peptide_accession-full-tree.json)) |
| F-03 | `protein_dpc_quant` | **No kit mechanism bounds this step without changing the science or the kernel** (measurements below). The `whole` variant stays; the step is named as the route's floor and the decision is the user's | — | n/a | 8.38 GiB (unchanged) | 5,564 s (unchanged) |

The registry defaults are the production variants; `whole` stays registered as the control arm and the in-process entry points follow the defaults. The fixture test (`test-workflow-prot-declarations.R`, 280) compares each variant with `whole` on the rolled fixture with a FASTA table so the cleanup applies.

## F-03: the DPC step's floor (measured 2026-10-24, the 300-sample store, VmHWM of a fresh process)

| Measurement | Time | Peak | Reading |
|---|---|---|---|
| limpa `dpcQuant`, 200 proteins (13,580 rows), chunk 1,000 vs 100 | 1,707 s both | 6.48 vs 5.94 GiB | identical results; the chunk is not the lever |
| hydration: the table + matrix + design (today) vs the matrix + design with the key map from the identity columns streamed from the saved table | — | 2.88 vs 1.05 GiB after hydration | the trim saves about 1.8 GiB of base, but see below |
| limpa `dpcQuant`, the first **2** proteins (904 peptide rows), trimmed hydration | 378 s | **5.62 GiB** | the working set is the per-protein fit of a large protein group (about 450 peptides × 300 samples) inside limpa (`peptides2Proteins` with standard errors); it is not a global statistic and no partition of the kit's lowers it |
| limpa `dpcQuant`, 200 proteins (9,677 rows), trimmed hydration | 1,477 s | 6.58 GiB | the step's floor with limpa on this corpus: about 6.6 GiB per worker, plus the application's 2.5 GiB |
| **fastlimpadpc, the whole kernel as the worker runs it** (table, matrix and design hydrated; all 2,500 proteins; one Rayon thread; the tree after the switch) | **168 s** | **3.69 GiB** | the step on the shipped engine: within the ticket's 4 GiB target and under the policy's 1,800 s with margin, so the `worker_deadline` entry is dropped (item 5); limpa on the same store 5,564 s and 8.38 GiB |
| fastlimpadpc `dpcQuant`, the same 200 proteins, `RAYON_NUM_THREADS=1` | **24 s** (fit 117 s) | **1.22 GiB** | bounded and fast, but a different engine: its DPC fit differs from limpa's (intercept −0.994 / slope 0.0079 against −1.081 / 0.0117 here; KIT-047 measured max abs ΔE 0.28 log2, median 0.0004, on the screen corpus), so its output is not exact against the `whole` variant |

Why the hydration trim is not shipped: the kernel (`proteinMissingValueImputationLimpa`) counts peptidoforms per protein as the rows of the long table (`dplyr::n()` over the 5.5-million-row table) and stores the summary in the state (`peptide_counts_per_protein`, read by the peptide count filter). Any hydration short of the whole table changes that fingerprint; keeping it exact needs the count computed outside the kernel and handed in, a kernel contract change outside this ticket (kernels are never rewritten here). The saving would be about 1.8 GiB on a step whose floor is limpa's 6.6 GiB regardless.

Consequence for the ticket's acceptance: with limpa the candidate's journey peak is the DPC step's (the application plus the whole-matrix worker, 10.89 GiB on the KIT-053 member) and cannot reach the DIA route's 5.7 GiB. The decision design's terminal outcome for F-03 is **rejected** (exact, above the target; the step named; the next mechanism a new decision). The three options are the user's: (1) the route's engine to fastlimpadpc (a scientific change: KIT-041 kept limpa the default pending the user's engine decision; the fit and the quantities differ as measured); (2) a kernel contract change so the DPC worker hydrates the matrix only (about 1.8 GiB saved, the floor unchanged); (3) the limpa floor accepted and recorded as the route's property on 300-sample corpora.

## Decision 2026-10-24: the route's engine is fastlimpadpc (the user: "confirm the switch")

The user's direction of 2026-10-19 named fastlimpadpc for the limpa route; KIT-047's determinism rule had kept limpa the shipped default while you were away (fastlimpadpc's fit was not bit-identical across runs under the host's default Rayon pool; it is on one thread), logged as an afk decision awaiting confirmation. Confirmed today. What changes (one bounded item inside KIT-057): the resolver's default is fastlimpadpc where installed (limpa otherwise, and selectable by name); the shipped `config.ini` key; the DPC tab lists the default first; the kernel pins `RAYON_NUM_THREADS=1` before the engine's first call unless the environment sets it, and records `engine_threads` in the results; the route's oracles are re-recorded on the new engine (the exactness gate is each production variant against its own `whole` run on fastlimpadpc; the limpa-versus-fastlimpadpc difference stays a characterised number: KIT-047 max abs ΔE 0.28 log2, median 0.0004, on the screen corpus); the `worker_deadline` entry is reconsidered on the measured DPC time; the counted pair runs on the new engine (declared below, the user's go).

## Counted pair: declaration (runs only on the user's go)

- Purpose: the ticket's acceptance on the full corpus through the GUI path: the candidate exact against its oracles, its journey peak and every worker's peak no higher than the control's, and its journey peak ≤ 5.7 GiB (the DIA route's on this corpus).
- Members: `kit:journey-057:control:1` = `kit_053` (`BTcMbt`, 81a69944: limpa, the whole variants; its own counted member was 6,922 s and 11.42 GiB) and `kit:journey-057:candidate:1` = `kit_057` (the pass-2 install from bdc3e815: fastlimpadpc, the production variants). Corpus: the 300-sample DIA-NN calibration corpus; normalisation binding as KIT-053's (cyclic loess, RUV as covariates, correlation 0.5).
- Oracles: the control against KIT-053's full-corpus fingerprints (`limpa_full_<step>`), the candidate against kit-057's, recorded in process on the candidate build (`limpa57_full_<step>`).
- Rule: the kit's pair rule (`whole_slice_anon`: candidate ≤ 0.95 × control on the whole-slice anonymous peak, every guarded step ≤ 1.02 × the control's), plus the ticket's absolute 5.7 GiB on the candidate's journey peak, judged in the outcome; no ratio promoted.
- Ceiling 14,400 s per member; the driver's browser budget 7,200 s per member (CORE-pinned; the control used 6,922 s under KIT-053). One attempt; one correction rerun for a demonstrated instrumentation failure only. Expected host time about 2 h (control) + 0.5 h (candidate).
- Terminal outcomes: accepted (exact, within the rule and the 5.7 GiB); rejected (exact, above a bound: the step named); invalid (not exact, or a member incomplete: the cause named, no retry without a diagnosis).

## Open

- Journey time against the driver's 7,200 s browser budget (CORE-pinned): the KIT-053 journey used 6,922 s; the two variants add about 420 s (log2 +148 s, accession +296 s) on the standalone measurements, so the candidate's journey is expected near the budget unless the DPC step (5,565 s) comes down; reported, not judged.
- F-03: the user's decision (above).
- The counted control/candidate pair: with limpa its outcome on the journey peak is predictable (the DPC step); whether to spend the pair (about 4 h of host time) is part of the decision.

## Incident

2026-10-23 19:46: systemd-oomd killed the VS Code scope (30.4 of 33.4 GB, 20 GB swap) while three heavy jobs ran together (the accession whole-variant measurement, 10.5 GiB; limpa's chunk measurement; a test file that runs limpa). My scheduling error; nothing lost (the work was on disk). Rule since: one memory-heavy job at a time on this host.
