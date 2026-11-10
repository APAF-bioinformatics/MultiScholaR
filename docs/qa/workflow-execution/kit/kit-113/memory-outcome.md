# KIT-113 D8: the memory measurement (2026-11-10)

**Decision: "better" is unreachable at W0, and the probe projects the move to be worse for the application.** As D8 declares, the move is made and judged on the "no worse" rules only. Those are met on every journey. Contingent **P-43** is raised beside **P-40**.

**No managed session ran an enrichment.** No journey reaches enrichment (D1), so the move's effect on a real application was not measured. The only figures for it are the W0 probe's, at statement level in one process.

## W0: the reachability probe (M0p; [w0/enrichment-heap.md](w0/enrichment-heap.md))

**Method:** one counted attempt, 29 min, inside the 30-min ceiling, on a synthetic DIA owner. It ran an in-process DA with its publication, then each statement group of D1. Kernels were offline: the live `processEnrichments` over the e2e doubles, and the deterministic double. No network call was made.

**The decision case** (live kernel, 16,000 proteins × 300 samples × 3 contrasts; the warm pass, declared before the run):

| Figure | GiB |
|---|---|
| In process, every D1 group in today's order | 1.075 |
| Managed controller, the groups that stay | 1.732 |
| In process − managed (rule: ≥ +0.25) | **−0.657** |

The cold pass gives −0.947.

**Why the move raises the controller's figure:**
- **The kernel's results carry `ggplotly` widgets,** about 31 MiB each.
- **In process,** the live object is 16 MiB, because it shares its environments with the loaded packages.
- **On the managed path,** the controller must hydrate the worker's portable copy, which serializes to 187 MiB and is 535 MiB in memory. The state commit then adds about +1.2 GiB on both paths, on top of that.
- **This is not data-scaled.** The results are the same size at 400, 4,000 and 16,000 proteins.
- **What does move is data-scaled, but smaller than the commit:** the DA tables (155 MiB at 16,000 × 3), the annotation match (+0.6 to +1.0 GiB transient) and the kernel (+0.67 GiB).
- **With the deterministic double,** which stores no plots, the move would lower the figure by 1.03 to 1.08 GiB at 16,000 proteins.

**The probe's limits, as D8 states them:**
- it measures R heap in one process at fixture scale, not a managed controller's anonymous memory;
- g:Profiler's real result size (terms, not samples) is not probed;
- the organism filter was off, taxon 9606 only;
- the worker was emulated in the same process.

## The journeys (M1, against `112l/m/t`)

All three are accepted and exact, with dispatches identical to `112l/m/t`'s: added none, removed none, changed none. That is as expected, because no journey reaches enrichment. Every member is valid: largest sample interval 0.06 s, swap 0, slice stopped. Per-step figures are in [memory/steps.tsv](memory/steps.tsv).

| Member | Slice peak (rule: ≤ baseline + 0.25) | Largest step change (reported) | `resident_end` | Elapsed |
|---|---|---|---|---|
| 113l | 1.751 → 1.750 (−0.001) | +0.032 (`report`) | 1.079 → 1.075 | 165.6 → 137.2 s |
| 113m | 2.062 → 1.957 (−0.105) | +0.008 (`import`); `peptide_count` −0.143 | 1.329 → 1.213 | 300.5 → 236.2 s |
| 113t | 1.749 → 1.751 (+0.002) | +0.032 (`report`) | 1.078 → 1.078 | 167.0 → 139.6 s |

"No worse, journeys" is met. The step changes are inside the screen spread KIT-112 recorded; KIT-113 changes no step the journeys run.

## The managed path's other costs (from the code, reported)

- **One dispatch per Run,** about 2.5 s on the screen (KIT-111's `correlation`).
- **The worker's deadline:** 1,800 s (P-41).

## Contingent rows

- **P-40 (raised, the user's):** measuring an enrichment on a real managed session needs a harness revision, a journey step with an offline, data-scaled kernel input. W0's figures are its evidence.
- **P-43 (raised):** "better" is unreachable at W0, and the probe projects the managed controller about 0.66 GiB higher at the enrichment step than in process. The options are the ticket's:
  - accept (the placement is the user's decision);
  - revert the placement (KIT-113's commits touch only its four `R/` files);
  - measure it anyway with P-40's harness ticket.
- **P-44 (not raised):** no journey failed "no worse" or exactness.
