# KIT-111 W0 (M1): the large in-process baseline `111p`

**Run.** Started 2026-11-10 00:41, on `kit_108` (`X6BzNV`, `7b52f535`), at `7c7c599e`. The member is the DIA route's whole journey on the synthetic 300-sample corpus (`full`), bound as KIT-065's counted member with the build replaced. It ran through `member.sh` with no other `multischolar*` unit live; the power state is in [../power-state.txt](../power-state.txt). The stage summary is [kit111-kit-journey-111p-stage.json](../kit111-kit-journey-111p-stage.json).

**Validity: valid** (`kitMeasurementValid`, journey rule):
- the largest sample interval is 0.31 s, against the 0.5 s limit (two gaps above 250 ms);
- swap is 0;
- the slice stopped;
- all 19 steps are complete, populations complete, and verification complete.

**The decision, verbatim:** `rejected`, rule `exactness_and_completion`, failing rule `oracle_mismatch`. The verifier's 18 mismatches are every state's `summary` and `normalisation_ruv/products/normalization_science`. This is **the same set KIT-065 and KIT-048 recorded** (`kit-065/outcome.md`, the shipped-default keys in the summary hash), with no canonical-form differences.
- **The wording differs from the ticket's expectation** (`invalid`, `verifier_accepted`, KIT-065 under protocol v1.4). Today's protocol (v1.9) reports the same verifier result as `rejected`/`oracle_mismatch`.
- D8's rule needs the set to be the same on `111p` and `111q`, with the stores identical by `compare-stores.R`. The wording does not change that.

**Figures** (GiB, anonymous):

| Step | Controller (app) peak | Slice peak | Seconds |
|---|---|---|---|
| imputation | 2.122 | 4.509 | 371.0 |
| **normalisation_ruv** | **2.620** | 2.982 | **366.4** |
| **correlation** | 1.774 | 2.135 | 17.8 |
| da | 2.115 | 3.391 | 100.0 |
| report | 2.350 | 2.711 | 26.1 |

Whole member:
- elapsed 1,902 s (KIT-065: 3,334 s);
- slice anonymous peak 5.73 GiB (`peak_anon_bytes`, in the rollup step);
- `resident_end_bytes` 2.65 GiB.

Every step's figures are in the stage summary.

The host's load at the end was 4.92, against 1.55 at the start. The W1 agent was told not to run R while the member was live; its report says whether it kept to that.
