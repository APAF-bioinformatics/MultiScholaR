# KIT-113 findings and lead decisions

## Lead decision at W0 (2026-11-10)

**The line map** (`w0/lines-at-start.txt`) does not stop W0. Every named function exists. Every named body is byte-identical to `525cbd78` apart from moving. Only KIT-112's anticipated changes touch the referenced lines.

**D8's reachability probe** (`w0/enrichment-heap.md`; one counted attempt, 29 min, inside the 30-min ceiling) gives **"better" UNREACHABLE at W0.**
- **The figures:** on the decision case (the live kernel over the offline doubles, 16,000 × 300 × 3), the in-process figure is 1.075 GiB and the managed controller's is 1.732 GiB. The difference is −0.657 GiB against the +0.25 GiB rule. The cold pass gives −0.947 GiB.
- **The move is projected to make the application worse, not just "not better".** The kernel's results carry `ggplotly` widgets (about 31 MiB each). Serialized portable results are 187 MiB for 3 contrasts. Hydrated in the controller, they are 535 MiB, against 16 MiB for the live object, which shares its environments with the loaded packages. The state commit then adds about +1.2 GiB on top, on both paths. The figure is not data-scaled: it is the same at 400, 4,000 and 16,000 proteins.
- **The deterministic double stores no plots.** With it, the move would lower the peak by 1.03 to 1.08 GiB at 16,000 proteins.

**Decision: the frozen protocol is followed.** D8's declared outcome for "unreachable at W0" is that the move is made, as the user's placement decision ("Join + move backend to worker", KIT-111 to KIT-113). It is then judged on the "no worse" rules only, and contingent **P-43** is raised beside P-40. P-43's wording anticipated this cause ("the results carry the plots' captured frames").
- The outcome will say plainly that the probe projects a higher application peak on the managed path, about +0.66 GiB at the enrichment step.
- **The journeys cannot see it,** because none reaches enrichment.
- P-43 recommends the user's choice between keeping the move and reverting it, with this evidence.
- **Reverting is cheap:** KIT-113's commits touch only the files in its `files_to_modify`.

**Pre-run declaration by the W0 agent** (D8 is silent on it): the warm pass decides. Both passes decide the same way.

## Findings

(none yet)
