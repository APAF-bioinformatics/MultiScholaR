# KIT-193 models: each Lean definition against the R lines it transcribes

For the statements gate (KIT-193 r1.2 section 5.3). Branch `kit-193-lean-proofs`; the Spec files are
`proofs/lean/MultiScholaRProofs/Spec/*.lean`. Code read on `v0.9-staging` at 51959495 (the branch point; no
change under `R/` since 147adde8), on `kit-187-filter-by-count` at e6c42d1d and on `kit-191-enrichment-test-options`
at 2f54469e (both with `git show`). KIT-190's marker row is read as ticket text (KIT-190 r1.2 7.4): its branch does
not carry KIT-187's table yet. Nothing here has been run against R: every pair below is a plan, written and run
after its product ticket merges (section 8).

Kinds of pair (5.5): `package` (a MultiScholaR function or constant, bound by digest, 8.2) or `library` (a locked
dependency, bound by version). Planted defect: the stated wrong version of one shipped function, through
`testthat::local_mocked_bindings()`, that the pair's planted case must catch (`expect_failure()`). The three fixed by
the ticket are marked (ticket); the others are proposals for the gate.

## 1. Definitions against code

| Spec file | Lean definition | Transcribes | Abstracts or omits |
|---|---|---|---|
| FilterCount | `failsRule P n k` | `exceedsRoundedPercent` (`kit-187:R/func_prot_qc_filtering_helpers.R:304-306`), `roundedPercentNumerator` (`:319`) | doubles (exact below 2^53, stated) |
| FilterCount | `allowedShare nstar m` | `max_missing / n_reps` (`R/func_general_parameter_helpers.R:254-258`); group level `(G - g)/G` (`:265-266`) | `round(…, 3)` and `round(1000 * c)`: replaced by `RoundedOK` |
| FilterCount | `RoundedOK P s` | the hypothesis the shipped cutoff must meet | R's decimal rounding (the pair checks it, (i)) |
| FilterCount | peptide floor in `peptide_floor_stmt` | `allowedByRoundedPercent` (`:313-315`) | double division and `floor` (stated) |
| ReplayRule | `Marker`, `TableOK` | `workflowReplayRuleMarkers` (`kit-187:R/utils_workflow_analysis_record.R:650-652`) plus KIT-190's row | the table is a constant (digest) |
| ReplayRule | `ruleLoop`; `rule` (existential) | `workflowReplayRule` (`:657-664`) | value decoding (`decodeWorkflowRecordValue`): values compared after decoding |
| ReplayRule | `inForce`, `compute` | `workflowFixedInForce` (`:686-693`) through `setWorkflowStepProvenance` (`R/utils_workflow_dia_state.R:729`), every writer | R's `NULL` for an empty list |
| ReplayRule | `computeManagedNewDesign` | KIT-190 7.4: the managed writer removes `ruv_control_grouping` (`R/utils_workflow_execution.R:3327-3332` as cited) | read as ticket text |
| ReplayRule | `build` | the builder's fixed entries and registry fallback (`:421-432`) | other entries (`others`) abstract, any place, source not `fixed`/`registry` (`rule` does not depend on order) |
| ReplayRule | `markersFixed`, `replayOptions`, `withOpts` | `workflowReplayRuleOptions` (`:675-680`), `withr::with_options` at the wrapper (`:785`) | per-step option scope (R's `on.exit`) |
| ReplayRule | `helperRule` | the helpers' `getOption(option, current)` with only the earlier token selecting the earlier rule (KIT-187 7.3) | the kernels themselves |
| GoClosure | `parentRel`, `keptEdges` | `goAncestorPairs` edge filter (`R/func_general_geneset_kernel.R:79`) | GO.db's parent map as a given list |
| GoClosure | `bfsStep`, `bfsLoop`, `bfsStart`, `bfs` | the `while` loop (`:80-91`) | a pass budget `fuel` in place of `while`; termination is `bfs_terminates_stmt` |
| GoClosure | `kept`, `propagate`, `propagateFin` | `propagateGoAnnotations` (`:128-138`) | data frames as finite sets; order |
| GoClosure | `build` | root removal (`:204-208`) | pooling over aspects (a union) |
| QvalueBh | `bhSorted` | `stats::p.adjust(p, "BH")` | tie order of R's `order` (does not matter) |
| QvalueBh | `bh` (closed form) | the same, `bh_sorted_eq_stmt` | none |
| QvalueBh | `qv`, `flag` | `qvalue.R:111-119` with a given `pi0`; the flag (`R/func_general_enrichment_orchestration.R:544`, in `processEnrichments`) | `pi0`'s estimation |
| QvalueBh | `gridVal`, `gridWt`, `fdr2`, `mixVal`, `mixWt`, `fdrKeptMix` | KIT-189 2.1's laws | the hypergeometric p-values of a real list |
| Permutation | `pval`, `unchangedDraws` | `integrationEnrichmentPermutation` (`kit-191:R/utils_integration_enrichment.R:1631-1656`) | `T` abstract; exact arithmetic; RNG uniform |
| Permutation | `reflect`, `reflectedDraws` | fixture C's reflection count | none |
| Permutation | `uniformDraws` | `vapply(seq_len(B), function(i) sample(n))` | Mersenne-Twister as uniform (idealisation) |
| RuvAlgebra | `RuvIII`, `alpha`, `Yc`, `ac`, `W`, `Ystar` | `ruv::RUVIII` (ruv 0.9.7.1, the locked version, `R/RUVIII.R:35-42`, read at the statements gate) | `U` abstract (svd's output), constrained by named hypotheses |
| RuvAlgebra | `centre`, `Wc`, `adjusted` | `computeRuvCovariates` (`R/func_prot_norm_optimization_helpers.R:1179-1183`); `computeOmicsRuvAssayFactors` (`R/func_omics_ruv_factors.R:170-173`) | none |
| RuvAlgebra | `restore` | `restoreRuvFeatureLevels` (`R/func_prot_norm_optimization_helpers.R:1133-1140`) on complete data | missing values |
| RuvAlgebra | `residop` | `ruv::residop` | none |
| RuvDownstream | `betaHat`, `resid`, `rss`, `tStat` | `lm.fit`/`lmFit` normal equations | limma's empirical Bayes |
| RuvDownstream | `betaHatW`, `residW` | `dpcDE` weights (`R/func_prot_da_model_methods.R:316-380`); `duplicateCorrelation` (`R/func_prot_da_model_stats.R:259-261`) with a fixed correlation | the correlation's estimation |
| RuvDownstream | `H1_stmt`'s `X * T⁻¹` | `ruvExactContrastsFit` (`R/func_general_ruv_factor_count.R:383-397`, `limma::contrastAsCoef`) | limma's `T` (any invertible `T` with first row `c`) |
| RuvIIIC | `RuvIIIC`, `Y0`, `rows`, `Uon`, `ctlOn`, `aOn`, `WOn`, `ystar` | RUV-III-C per target (KIT-185 2.1; RUVIIIC at c534151f, `src/RUVIIIC_Varying.cpp:151-165`, `:229-303`) | `P_R` and the eigen-solver abstract (no statement uses their form); closed-form and all-missing branches |
| RuvIIIC | `ystarT`, `foldErr` | KIT-188 2.4 item 1's error | none |

## 2. Items: statements, pairs, planted defects

| Item | Spec statements (witness) | Pair file | Kind; calls | Planted defect |
|---|---|---|---|---|
| T1 | `never_stricter`, `exact_lcm`, `exact_316`, `balanced`, `share_mono`, `peptide_floor`, `counterexample`, `counterexample_any_rounding` (witness each) | `test-kit193-pairs-filter.R` | package: `exceedsRoundedPercent`, `roundedPercentNumerator`, `allowedByRoundedPercent`, `updateMissingValueParameters` | (ticket) `roundedPercentNumerator <- function(cutoff_percent) 2 * round(1000 * cutoff_percent)`: the table differs at n* = 3, m = 2, n = 3, k = 1 |
| T2 | `rule_first_match`, `totality`, `faithful`, `without_marker`, `managed_new_design`, `round_trip`, `locality`, `unfixed_default` (witnesses) | `test-kit193-pairs-replay.R` | package: `workflowReplayRule`, `workflowFixedInForce`, `workflowReplayRuleOptions`, `buildWorkflowAnalysisRecord`; constant `workflowReplayRuleMarkers` | (ticket) `workflowReplayRule` returning the current token for an entry of source `registry` |
| T3 | `bfs_terminates`, `bfs_correct`, `propagate_fin`, `idempotent`, `monotone`, `extensive_kept`, `not_extensive`, `labels_mono`, `root_removed`, `all_absent`, `bfs_all`, `aspect_closed` (witnesses) | `test-kit193-pairs-go.R` | package: `goAncestorPairs`, `propagateGoAnnotations`, `buildUniprotGoGeneSets`; constants `goPropagationRelations`, `goAspectRoots` | proposed: `goAncestorPairs <- function(terms, parents, relations) data.frame(term = terms, node = terms)` (no closure): the chain graph's rows differ |
| T4 a-c | `bh_sorted_eq`, `qflag_contains_bh`, `flag_none`, `qflag_extra`, `bh_ge_p`, `bh_mono`, `bh_ties`, `bh_stepup` (witnesses) | `test-kit193-pairs-qvalue.R` | library: `qvalue::qvalue`, `stats::p.adjust`; package: `processEnrichments` (the flag, `:544`) | proposed: the route's stored q-values replaced by `pmin(1, 2 * p.adjust)` (pi0 = 2) through a mock of the package function that builds the `@result` table; the nesting assertion (ii) must fail |
| T4 d-i | `grid_mass`, `grid_superuniform`, `fdr_n20`, `fdr_n20_strict`, `fdr_n25` | same file | library: `stats::p.adjust` (table `bh-counterexample.tsv` from the model) | proposed: the R enumeration applies BH to the whole family where the kept family is asked (a planted helper inside the test file, since base `p.adjust` cannot be mocked); the 0.075 check must fail |
| T4 d-ii | `kept_family_formula` (0 < a < 1), `kept_family_unit_a_counterexample`, `kept_family_instance` (witness: the paper's instance) | none (no R object; the figure is KIT-189's) | none | none |
| T5 | `perm_finite_core`, `identity_clause`, `perm_tolerance_mono`, `perm_measure`, `perm_floor`, `perm_floor_reflect`, `bh_floor`, `bh_floor_strict`, `bh_floor_sentence_counterexample` (witnesses) | `test-kit193-pairs-permutation.R` | package: `integrationEnrichmentPermutation` (through its exported caller), `integrationEnrichmentSetT`; library: `stats::p.adjust` (5c') | proposed: `integrationEnrichmentSetT` wrapped so that permuted columns' t is scaled by 1e-6 (only unchanged draws counted): fixture A's `#{p ≤ 1/3} ≤ 4608` must fail |
| R-A1, R-A2 | `A1`, `A2` (witness `A2_witness`) | `test-kit193-pairs-ruv.R` | library: `ruv::RUVIII` | proposed: `ruv::RUVIII` mocked (`.package = "ruv"`) to return `newY = Y`; `Y - newY = Wα` must fail. If testthat cannot mock it there, the gate decides the replacement |
| R-A3 | `A3` | same | package: `computeProtRuvModelFactors`, `computeOmicsRuvAssayFactors` | proposed: `computeOmicsRuvAssayFactors` mocked to return status `ok` with a zero W at `k = 0` |
| R-B1, R-B2, R-E1 | `B1`, `B2`, `E1` (witness `B1_witness`), `E2` (instance) | same | package: `computeRuvCovariates`; library: `RUVIIIC::RUVIII_C_Varying` | proposed: `restoreRuvFeatureLevels <- function(corrected, original) corrected` (B2's exact-c shift must fail) |
| R-C1 | `C1` | same | package: `restoreRuvFeatureLevels`, `computeRuvCovariates` | (ticket) `restoreRuvFeatureLevels <- function(corrected, original) corrected` |
| R-D1, R-D2 | `D1`, `D2` (witness `D_witness`) | same | package: `computeRuvCovariates` | proposed: `computeRuvCovariates` wrapped to return `2 * W` for the scaled input; "W unchanged" must fail |
| R-F1, R-F2, R-F3 | `F1`, `F2`, `F3_span`, `F3_fwl` (witnesses `F_witness`, `F3_fwl_witness`) | same | package: `computeRuvCovariates()` (`$W`, `$adjusted`), `computeProtRuvModelFactors`; library: `limma::lmFit` | proposed: `$adjusted` replaced by `Y - 2 * W α` through a wrapper of `computeRuvCovariates`; equal group coefficients must fail |
| R-G0 | `G0` (witness) | same | package: `ruvIII_C_Varying` (P1), `computeRuvCovariates()$adjusted` (P6), the lipid method (L1) | proposed: the corrected matrix replaced by the uncorrected one; rank `n - g - k` must fail (it is `n - g`) |
| R-G1, R-G2 | `G1` (witness), `G2`, `G2_sign` | same | library: `limma::lmFit` on shipped `W`, `$adjusted` | proposed: as R-F1 |
| R-H1 | `H1`, `H1_weighted` (witness) | same | package: `ruvExactContrastsFit`; library: `limma::contrastAsCoef` | proposed: `ruvExactContrastsFit` mocked to return `contrasts.fit`'s approximation on features with missing rows |
| R-J1, R-J2 | `J1` (witness), `J2` (witness), `J2_contribution` | same | library: `RUVIIIC::RUVIII_C_Varying`, `ruv::RUVIII` (upstream's own equality) | proposed: the zero-filled matrix replaced by mean-filled in the replicated upstream check; equality must fail |
| R-K1 | `K1`, `K1_gls` (witnesses) | same | library: `RUVIIIC::RUVIII_C_Varying`, `ruv::RUVIII` (k = 9 and NULL) | proposed: `ruv::RUVIII(k = NULL)` result replaced by `newY` at `k = 8` (wrapper in the file); equality must fail |
| R-L1, R-L2 | `L1`, `L1_invariance` (witness), `L2` (witness) | same | library: `RUVIIIC::RUVIII_C_Varying`; package: the shipped corrected values | proposed: `e_j` compared with the change of a contrast whose weights do not sum to 0 |
| R-M1 | `M1` (witness) | same | library: `edgeR::glmQLFit` | proposed: `W Q` replaced by a non-invertible `Q` (rank 1); equal fitted values must fail |

## 3. Statements not yet frozen (departure from ticket 5.3 item 1, recorded in r1.3)

Each gets its own new Spec file and its own statements review before its proof; the files frozen now keep their
hashes. None is listed as proved anywhere until then.

| Item | Why not in this pass |
|---|---|
| R-I1 (priority 3, last of the ticket) | KIT-184's removed factors (5.4) are not written; the model (`D_c`, the rank rule at `1e-8 · s_1`, the averaging map) is stated from KIT-184's merged code, before its proof starts |
| R-J3 (priority 3) | not written in this pass; the gate read the branch conditions at c534151f: all missing when `min(residual dimensions, controls complete on the target's runs) < k` (`:168-171`), when fewer than k eigenvalues converge (`:252-263`), and on any exception in the eigen step (`:304` onward); otherwise finite exactly on the target's observed runs (`:186-198`, `:286-301`) |
| R-N1 (priority 4) | needs R's `rank(…, ties = "first")` and `order()` semantics on the three first-pass functions read and transcribed; not cheap in this pass |

## 4. The gate's own reads (ticket 2.3 and 14): done at the gate

- RUVIIIC at the locked commit c534151f and `ruv` 0.9.7.1 were read by the statements gate
  (`AUDIT-KIT-193-statements-20261124.md` section 5); the docstrings cite those lines (S07).
- RUVSeq: no statement depends on it; R-M1 is linear algebra about any W.
- No digest in `bindings.tsv` is filled: rows are written when each pair first runs on merged code (8.2).
