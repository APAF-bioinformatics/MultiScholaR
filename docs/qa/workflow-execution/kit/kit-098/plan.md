# KIT-098 plan: blocked at the run-order planning item (2026-11-06)

**Status: incomplete. This file authorises no harness, bindings or protocol change and no run.** The ticket requires the
run-order item to be settled before the notebook comparator is fixed. It cannot be settled without a user decision, so
the rest of the plan (step list, comparator mapping, settings, fault) is not written yet.

## Planning item from KIT-096: run order

A GUI-built design keeps the counts file's run order (KIT-096, `run_order = identity` at `R/mod_transcriptomics.R:44`).
The notebook reorders runs to design-file order (`count_matrix[, design_matrix$Run]`). KIT-040's corpus shuffles the count
columns on purpose.

**(a) Aligning runs by name does not give `identical` results.** Evidence: [run-order.R](run-order.R),
[run-order.txt](run-order.txt). KIT-093's kernels (the functions KIT-095's declarations call) ran on the bound `trx_screen`
corpus on `kit_097` (R 4.6.1, reference BLAS/LAPACK 3.12.0), once per design order:

- **Design-file order:** 21 of 21 products are `identical` to KIT-040's full reference (the control).
- **Counts-file order (the GUI's), each product aligned to design-file order by run name:** 13 of 21 are `identical`.
  - Equal: the filter masks, filtered counts, effective library sizes, TMM samples table, raw logCPM, preliminary order,
    k cap, grid summary, best k, controls, RUVs counts and the control mask.
  - Not equal:
    - TMM logCPM (2.7e-15);
    - W: the first RUVs factor's sign flips under the permutation, and W_2 moves by 2.4e-14;
    - the RUVs-corrected logCPM (1.2e-14);
    - the DA design (it carries W);
    - the edgeR table and topTags (1.1e-08);
    - RUVIII newY (4.0e-14);
    - M (equal values; the alignment's row subset drops `model.matrix` attributes, an artefact of the projection).

  The cause is the sample order, not the path: edgeR's sums over samples and the SVD inside RUVs depend on it.

**(b) Design-file order in the design state needs code in `R/`.** The lipid and metab builders order a fresh design by
`gtools::mixedsort` of the sample names (`R/mod_lipid_design_builder_state_helpers.R:395`,
`R/mod_metab_design_builder_server.R:237`). On `trx_screen`, mixedsort order is not design-file order either, so parity
with lipid/metab would not reach the notebook's order. Only a design read from a file in the notebook's order would. KIT-096
did not offer "Import Existing Design" under the profile, so that also means changing `R/`. Both are outside this ticket's
`files_to_modify`. This is a scope question.

**(c) A declared tolerance** needs the user's go. Not proposed here.

## Decision needed (user)

Choose how the journey store is compared with the notebook:

1. **(b) Change the design state's run order in `R/`.** This could be a new transcriptomics ticket, or KIT-096 reopened.
   Options are to read an imported design file in its own order, or to order by an explicit rule the notebook comparison
   can reproduce. It needs a new build before KIT-098's oracles and journeys.
2. **(c) Declare a tolerance** for the eight products above. This weakens the independent check, so it needs your go.
3. **(d) A second notebook reference in the GUI's order.** Render the unchanged notebook with KIT-040's driver on the same
   counts, with `design_matrix.tab` rows permuted to the counts file's order. The design has the same content, only its
   row order changes. This needs no `R/` change and no tolerance. It is a new KIT-040-style render (about 30 s) and a bound
   input variant, so it needs your go, and possibly the reference's owner, KIT-040.

The same-build oracle comparison (path equivalence) does not depend on this decision. The in-process recorder uses the
GUI's run order. But if (b) is chosen, the build changes, so no oracle or member runs before the decision.
