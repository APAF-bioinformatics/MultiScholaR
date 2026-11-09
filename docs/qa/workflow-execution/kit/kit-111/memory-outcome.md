# KIT-111 D8: the memory and time measurement (2026-11-10)

**Decision: "no worse, not better".** The four normalisation-tab operations now run on the worker, and every result is identical. The application's (controller's) `normalisation_ruv` peak on the 300-sample DIA corpus did not fall (2.620 → 2.626 GiB; "better" needed ≤ 2.370). Every "no worse" rule holds on the large pair. On the screen members it holds for 111l and 111t, and for 111m after the correction pass's rerun (111n). As declared, the move is kept and contingent **P-33** is raised. It is not tried here.

All members are recorded and uncounted, with one attempt each under the host rule (`power-state.txt`). Every member is valid: largest sample interval ≤ 0.31 s, swap 0, slice stopped. GiB = 2^30 bytes, controller anonymous memory per step unless named.

## The large pair (111p in process on `kit_108`; 111q on the worker on `kit_111`)

- **Exactness:**
  - Both verifier decisions are `rejected` on `oracle_mismatch`, with the same 18 mismatches, the state summaries plus `normalisation_ruv/products/normalization_science`, which are KIT-065's set. The ticket expected the wording "`invalid` on `verifier_accepted`"; the decision is recorded verbatim.
  - Store against store (`memory/compare-stores.R`, `memory/compare-111p-111q.{tsv,txt}`): 16 states, 98 comparisons, **0 different**.
- **Dispatches:** 111q = 111p's 13 + `normalization-ruv`, `pearson-pairs`, `pearson-selection` (each 1/1).

| Rule | 111p | 111q | Limit | Result |
|---|---|---|---|---|
| **Better:** `normalisation_ruv` controller peak | 2.620 | 2.626 | ≤ 2.370 | **not met** |
| No worse, application: every step's peak | — | largest rise +0.093 (`da`); `normalisation_ruv` +0.005 | +0.10 | met |
| No worse, slice: whole-member anonymous peak | 5.726 | 5.813 (+0.087) | +0.25 | met |
| No worse, time: `normalisation_ruv` | 366.4 s | 390.4 s (1.066×) | ≤ 1.10× | met |
| No worse, time: `correlation` | 17.8 s | 19.0 s (1.071×) | ≤ 1.10× | met |
| Reported: `correlation` peak | 1.774 | 1.806 | — | — |
| Reported: `resident_end_bytes` | 2.654 | 2.643 | — | — |
| Reported: whole journey | 1,902 s | 2,132 s (1.12×) | not judged | — |

**Why "better" was not met (the traces, `candidate-1/memory-trace.tsv`, sampled every 50 ms):**
- In 111p the controller sat at 1.95 to 2.00 GiB while the RUV kernel ran (0 to ~330 s), dipped to 1.77, then rose to 2.62 in the step's last ~30 s.
- In 111q the controller is flat at 1.82 GiB while the worker runs the kernel (0 to ~350 s; the slice is 2.80 to 2.87 GiB, so the worker holds about 1.0 GiB). It then rises to 2.63 in the step's last ~30 s.

So the kernel left the controller (−0.13 to −0.18 GiB while it runs). But the step's peak is the in-process publication half (`publishDIAReferenceNormalization`: the files, QC figures and display of the 300-sample result), which is the same on both paths. This is exactly contingent P-33: "the controller's step peak is held by the work left in process". KIT-112's W0 checks reachability before any code change for this reason.

**Time:** the whole journey is 1.12× slower, but most of that is in steps KIT-111 did not change (`imputation` 371 → 468 s, `qvalue` 160 → 178 s). The ticket times only the two moved steps. During 111q, two read-only agents were authoring and reviewing the KIT-112 ticket (file reads and `git`, no R). That is recorded here as host activity in the measured window. The host's load average stayed at 0.84 to 1.10.

## The screen members (against `108l/m/t`, reused)

| Member | Exact | Moved operations | Largest step rise | Slice | `resident_end` | Elapsed |
|---|---|---|---|---|---|---|
| 111l | yes | 3 × 1/1 | +0.002 | 1.702 → 1.749 | 1.098 → 1.101 | 129.5 → 139.7 s |
| 111t | yes | 3 × 1/1 | +0.001 | 1.702 → 1.750 | 1.098 → 1.099 | 125.7 → 136.3 s |
| 111m | yes | 3 × 1/1 | **+0.133** (`peptide_log2`, `peptide_accession`; +0.10 to +0.13 from `rollup` on) | 1.973 → 2.096 | 1.239 → 1.366 | 230.1 → 250.5 s |
| 111n (111m again, the correction pass) | yes | 3 × 1/1 | +0.027 (`peptide_log2`); the rest ≤ +0.007 | 1.973 → 1.966 | 1.239 → 1.238 | — |

**111m and the correction pass:**
- 111m's rise begins at `rollup`, before any operation KIT-111 moved runs in that journey.
- Its trace shows the same transient spikes as 108m's (about +0.09 per step), on a floor that was never released: 108m's falls from 0.857 to 0.730 GiB early in `rollup`, while 111m's stays near 0.875.
- Across `101m` to `110m` (nine builds) the `rollup` peak was 0.823 to 0.874 GiB.
- No D3 code path explains a retention there: the runner installed at the tab's server start is idle until used, and the new input guard reads a reference.
- So the correction pass reran the member once, as 111n (bound in `05b938ef`, 111m's binding). 111n is within the spread on every step and equals 108m at the end.

**Read:** 111m's +0.13 GiB was a one-off where freed memory was not returned to the system, not a retention KIT-111 causes. The application rule holds on the rerun. Both members are kept and reported.

**Screen times** are reported, not judged. `correlation` is 6 → 11 s on 111l and 111t: a worker start-up on a five-protein fixture.

## Lifecycle (M2) and the failure half

111f failed at `normalisation_ruv`: `normalization-ruv` 1/0, "production action failed: … managed action ended: worker_failed", store head `protein_replicate_filtered`. See `proof/chain-b.md`.

## Contingent rows raised

- **P-33 (raised):** move the proteomics normalisation tab's publication half (files, QC figures, display) to the worker, as KIT-082 did for lipid and metab, as its own measured ticket; or accept. The evidence is above. The kernel move alone cannot lower this step's controller peak, because the peak is in the publication.
- **P-34 (not raised):** no "no worse" rule fails after the correction pass.
