# KIT-128: what KIT-129's re-record should expect to move

Written 2026-11-17 on branch `kit-128` (held until KIT-129, decision D6). Read-only note; nothing was re-recorded.

## On designs with technical replicates (all three kit corpora: lipid, metab, trx)

- **Lipid, metab:** the DA step and everything after it on `lipid_screen` / `metab_screen`: t, P, q, BH, B, the
  significance calls and counts, the saved `da_fit` entities (new identity, the block is in it), the DA manifests'
  model record, DA figures, exports, report text. logFC and AveExpr are expected unchanged (balanced, complete).
- **Trx:** everything from the design state on, on `trx_screen`: the design state (one sample per biological unit,
  named by `tech_rep_group`; `args$technical_replicates_summed`), filter counts, TMM/logCPM, grid, RUVs, DA, figures,
  exports. The KIT-040 notebook references describe library counts and no longer the route on these corpora.

## On every lipid and metab DA, technical replicates or not (audit K128-4)

The tables are identical on designs without technical replicates (shown on the branch by digest), but the objects
around them gain fields:

- every lipid/metab fit (`MArrayLM`, so `contrasts_results[[assay]]$fit.eb` and the saved `da_fit` products) carries
  `block_reason`, `uses_duplicate_correlation` and `consensus_correlation`;
- every lipid/metab DA manifest's `model` carries `block_reason` and `uses_duplicate_correlation`.

So product hashes of those fits and manifests move on every lipid and metab journey, not only on the corpora with
technical replicates. Trx fits (`DGEGLM`) and proteomics records do not gain them.

## Refusals that are new (K128-1, K128-2)

- A design whose every group is one biological unit measured more than once is refused in every omic (the user's rule
  of 2026-11-17). No kit corpus is such a design; no stored reference is expected to move for it.
- The shared design builder refuses replicate numbers shared inside a group unless the runs were marked by
  'Assign Technical Replicates'. The kit driver marks them that way, so its journeys are not refused.

## Tests left failing for the re-record

- `test-module-ci-trx-design.R`, MCI-027.2 (the journey setup on library counts).
- `test-workflow-omics-declarations.R`, the KIT-073 frozen dump (its lipid/metab DA tables, 6 units by 2 injections).
