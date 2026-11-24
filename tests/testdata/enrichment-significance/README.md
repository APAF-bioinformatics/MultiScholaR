# KIT-203 fixtures: the base commit's outputs of the clusterProfiler route

Written on a GitHub-hosted runner by `docs/qa/workflow-execution/kit/kit-203/gate-k203.R` (G0, run 45326277295) and
`capture-k203.R` (run 45329062402, package commit 63b95f90, whose `R/` equals `v0.9-staging` 9453e79b), each capture
inside `testthat::test_that()`; the S3 and T7 captures were byte-identical in a second process (run 45326277295's
capture). Inputs: KIT-173's committed UniProt fixtures (`tests/testdata/go-propagation/`, UniProt 2026_04, CC BY 4.0).
Read by `tests/testthat/test-enrichment-significance.R` (S3, S4, S6, S8). Not to be re-recorded: they are the base.

| File | SHA-256 | What |
|---|---|---|
| `e0-accessions.tsv` | `d5ae6b2572f77d31fee4e2ae48e018dc25520f8f4a84313d125b8a0a0417ecbe` | KIT-173's null draw 5 (126 E. coli accessions), as KIT-189's check E0 drew it; written by gate-k203.R g0 |
| `s3-none-up.rds` | `4eff27bf4a9eb38efbe34f9dbe3242a8efa68f7b22eae8d0eef9fc760cc9a6ca` | the same run: `list(result = @result, plot_inputs)` of the stored up result |
| `s3-none-up.tsv` | `89a382359bacb9942bfe2092a9ada1faeebb0802f0885da44285f76b8f8008cc` | S3's list through `processEnrichments()` on the base commit, `go_propagation = "none"`: the written up TSV |
| `s6-identity.rds` | `1805791905a8179cee1a3a4172574992e15ed546ef6e65af54550e45ba3f5db6` | KIT-173's T7 project on the base commit (real kernel, the tab's default propagation): identity, declared names and values, outcome names |
| `t7-iapor-down.rds` | `703f305739fe6e5a2bedd2b27f426fed77150066bb8f3d1dc5b220f685839502` | KIT-173's T7 fixture with T12's table, `go_propagation = "is_a_part_of_regulates"`, down list: `list(result = @result, plot_inputs)` |
| `t7-iapor-down.tsv` | `71d0f08c76fdfa3691838843ca2d8a45322667f1259c03884545fdd2c9eb6cda` | KIT-173's T7 fixture with T12's table, `go_propagation = "is_a_part_of_regulates"`, down list: the written TSV |
| `t7-iapor-up.rds` | `f6683cb13f5ac409557809f1ad1a56dd3505404ad9fcf197d9b412119631b14e` | KIT-173's T7 fixture with T12's table, `go_propagation = "is_a_part_of_regulates"`, up list: `list(result = @result, plot_inputs)` |
| `t7-iapor-up.tsv` | `c1e2269eb9a3efca201eb82b63baaeb3d290b16186b87606759a09b81ca94a2a` | KIT-173's T7 fixture with T12's table, `go_propagation = "is_a_part_of_regulates"`, up list: the written TSV |
| `t7-none-down.rds` | `91f9d030fe3f1905dd98f169e832988db63412e119c3a1c7ee44a7f7ce06d49a` | KIT-173's T7 fixture with T12's table, `go_propagation = "none"`, down list: `list(result = @result, plot_inputs)` |
| `t7-none-down.tsv` | `57ea14bfaf446f2bfcfd2f8ae1a8eaa43cbe0b3c771a9e5e9a2f8a84c5b7aa7f` | KIT-173's T7 fixture with T12's table, `go_propagation = "none"`, down list: the written TSV |
| `t7-none-up.rds` | `00d84ec56cace7a85e091934c2775f59c6ca13818aebb06d16fef3d0a7064bbd` | KIT-173's T7 fixture with T12's table, `go_propagation = "none"`, up list: `list(result = @result, plot_inputs)` |
| `t7-none-up.tsv` | `e3fd1a05684d486cdda737a7a7fab1a76448be5acba286d6526b319074be784c` | KIT-173's T7 fixture with T12's table, `go_propagation = "none"`, up list: the written TSV |
