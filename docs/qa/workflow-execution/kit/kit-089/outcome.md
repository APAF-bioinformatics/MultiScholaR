# KIT-089 outcome: the intensity filter on `partitioned` by default

**The lipid and metabolite intensity filter now runs `partitioned` by default, with or without the worker, and the record says so. Both small journeys produce exactly the same content as the earlier accepted ones: every state and every judged product. The verifier itself marked both journeys `invalid`. Its product fingerprints hash nested tables by their file bytes, and a partitioned write stores the same content in different bytes. Fixing that rule is a protocol decision for the user, and KIT-088 needs it.**

User direction, 2026-11-06: "Switch the intensity filter to partitioned on both omics now." Build `kit_089` = xD1Tin, c8d79db2.

## The change

- **Registry** (`R/utils_workflow_operations.R`): `lipid_intensity_filter` and `metab_intensity_filter` default to `partitioned`. Every other step stays `whole`.
- **In-process path** (`runWorkflowOmicsReferenceStep`, `R/utils_workflow_omics_declarations.R`): passes the registry default when it isn't `whole`. Before this, a run without the worker would have run `whole` and recorded `whole` regardless of the default. Product identity is unchanged, because `variant` is not a declared parameter.
- **Not changed:** `whole` stays selectable. The journey's peak sits in normalisation, so this lowers the step, not the journey.

## Evidence

| Member | Decision | Time | Whole slice | Peak | Intensity step | Recorded variant |
|---|---|---|---|---|---|---|
| `089u` lipid screen, missing values | invalid (`verifier_accepted`) | 181 s | 1.94 GiB | 1.83 GiB | 11 s, 1.24 GiB | `partitioned` |
| `086u` the same on `kit_087`, `whole` | accepted, exact | 186 s | 1.94 GiB | 1.80 GiB | 11 s, 1.24 GiB | `whole` |
| `089v` metab screen, missing values | invalid (`verifier_accepted`) | 194 s | 1.95 GiB | 1.85 GiB | 9 s, 1.24 GiB | `partitioned` |
| `086v` the same on `kit_087`, `whole` | accepted, exact | 180 s | 1.95 GiB | 1.85 GiB | 10 s, 1.24 GiB | `whole` |

At 24 samples this step is too small for the variant to change its memory. KIT-084's `084p` measured it at large: .81 and .77 of `whole`.

**The verifier's findings.** Against `lipid86s` and `metab86s` (recorded on `kit_086`), no state, component, support entry, summary or row count differed. The only mismatches were `products/names` at intensity, duplicates, normalisation, correlation and DA.

**Why the product names differ.** A step's result product holds a reference to the state it saved. The harness's canonical form (`kitCanonical`, `tools/profiling/workflow_kit.R`) hashes a table nested in a product by its file bytes, while a state's own tables are hashed by content, page by page. The partitioned writer stores the same content as different parquet bytes, with different table metadata and dependency lists. The later steps' products refer back to those tables, so the difference carries forward. On the intensity product, the counts are identical (206, 206 → 189, 188); only the referenced bytes differ.

**Content comparison (`content-compare.R`, `content-compare.txt`).** This compares each new journey with the earlier accepted one, `086u`/`086v`, which was exact against the same `kit_086` oracles. The harness file is unchanged. In the script's own copy of the canonical form, a nested table is hashed by its pages, as a state's own components are.

| Pair | States | Judged products | Rendered figure pages |
|---|---|---|---|
| `089u` vs `086u` | 10 of 10 identical (class, rows, components, support, summary, products) | 70 of 70 identical | 65 and 65, not judged |
| `089v` vs `086v` | 10 of 10 identical | 64 of 64 identical | 59 and 59, not judged |

Rendered figure pages differ between two members of the same build (KIT-078), so the harness never judges them across members either.

**Recorded variants (`variants.txt`):** `partitioned` on the intensity state in both new journeys, `whole` on duplicates and log2; `whole` everywhere in `086u`.

## Acceptance, as written and as met

- **Met:** the worker and the in-process path both run and record `partitioned`, with the same result as `whole`.
- **Not met as written:** "Both screen journeys accepted and exact". The protocol's decision is `invalid` and stands as recorded.
- **What was shown instead:** content-identical to journeys that were accepted and exact against the old build's records. This doesn't loosen the comparison: it judges every state and every non-figure product, by content.

## Decision for the user (needed before KIT-088)

KIT-088 compares the final build, which runs `partitioned` on this step, against `kit_086`. Every such journey will hit the same verifier result. The options:
1. Correct the harness's canonical form so a nested table is hashed by content, and re-fingerprint the oracle stores the comparisons use. This is a change to the frozen protocol.
2. Keep the protocol and judge those members with the content comparison, as here.

## Checks

- **Test first (`test-first.txt`):** the new test failed before the change (4 of 12) and passes after.
- **Dependents (`dependents.txt`):** 63 files, 10,096 pass, 0 fail, 0 error, 120 skip.
- **`R CMD check` (`check.log`):** at baseline, 3 WARNINGs and 3 NOTEs.
- **Attempts (`power-state.txt`):** the first run's dependents had 9 failures in 3 module tests. Those tests fake the intensity kernel, which under `partitioned` sees only an empty template, so they now pin `whole`, the path they characterise. That was this ticket's one correction pass; the rerun passed.
- **`fastduckdb`:** in Suggests. Partitioned paths need it, as proteomics' defaults already do. Without it, the omics journey test already errored before this change.
