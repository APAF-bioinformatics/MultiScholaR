# KIT-193 models: each Lean definition against the R lines it transcribes

For the statements gate (KIT-193 r1.2 section 5.3). Branch `kit-193-lean-proofs`; the Spec files are
`proofs/lean/MultiScholaRProofs/Spec/*.lean`. Code read on `v0.9-staging` at 51959495 (the branch point; no
change under `R/` since 147adde8), on `kit-187-filter-by-count` at e6c42d1d and on `kit-191-enrichment-test-options`
at 2f54469e (both with `git show`). KIT-190's marker row is read as ticket text (KIT-190 r1.2 7.4): its branch does
not carry KIT-187's table yet. Each pair is written and run after its product ticket merges (section 8). Run so far:
target 3's pair (P0a, run 45302637982) and target 1's two pair files (P1a, run 45323167019); the others are plans.

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
| QvalueBh | `bhSorted` | `stats::p.adjust(p, "BH")` | tie order of R's `order` (does not matter for p ≥ 0) |
| QvalueBh | `bh` (closed form) | the same for p ≥ 0, `bh_sorted_eq_stmt` | none |
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
| T1 | `never_stricter`, `exact_lcm`, `exact_316`, `balanced`, `share_mono`, `peptide_floor`, `counterexample`, `counterexample_any_rounding` (witness each) | `test-kit193-pairs-filter.R` and `test-kit193-pairs-filter-table.R` (helpers `helper-kit193-filter.R`) | package: `exceedsRoundedPercent`, `roundedPercentNumerator`, `allowedByRoundedPercent`, `updateMissingValueParameters`, `removeRowsWithMissingValuesPercentHelper` (driven by the pair since P1a; digest due from the next R run, second proofs audit Q10) | (ticket) `roundedPercentNumerator <- function(cutoff_percent) 2 * round(1000 * cutoff_percent)`: the table differs at n* = 3, m = 2, n = 3, k = 1 |
| T2 | `rule_first_match`, `totality`, `faithful`, `without_marker`, `managed_new_design`, `round_trip`, `locality`, `unfixed_default` (witnesses) | `test-kit193-pairs-replay.R` | package: `workflowReplayRule`, `workflowFixedInForce`, `workflowReplayRuleOptions`, `buildWorkflowAnalysisRecord`; constant `workflowReplayRuleMarkers` | (ticket) `workflowReplayRule` returning the current token for an entry of source `registry` |
| T3 | `bfs_terminates`, `bfs_correct`, `propagate_fin`, `idempotent`, `monotone`, `extensive_kept`, `not_extensive`, `labels_mono`, `root_removed`, `all_absent`, `bfs_all`, `aspect_closed` (witnesses) | `test-kit193-pairs-go.R` | package: `goAncestorPairs`, `propagateGoAnnotations`, `buildUniprotGoGeneSets`; constants `goPropagationRelations`, `goAspectRoots` | proposed: `goAncestorPairs <- function(terms, parents, relations) data.frame(term = terms, node = terms)` (no closure): the chain graph's rows differ |
| T4 a-c | `bh_sorted_eq`, `bh_sorted_neg_ties_counterexample`, `qflag_contains_bh`, `flag_none`, `qflag_extra`, `bh_ge_p`, `bh_mono`, `bh_ties`, `bh_stepup` (witnesses) | `test-kit193-pairs-qvalue.R` | library: `qvalue::qvalue`, `stats::p.adjust`; package: `processEnrichments` (the flag, `:544`) | proposed: the route's stored q-values replaced by `pmin(1, 2 * p.adjust)` (pi0 = 2) through a mock of the package function that builds the `@result` table; the nesting assertion (ii) must fail |
| T4 d-i | `grid_mass`, `grid_superuniform`, `fdr_n20`, `fdr_n20_strict`, `fdr_n25` | same file | library: `stats::p.adjust` (table `bh-counterexample.tsv` from the model) | proposed: the R enumeration applies BH to the whole family where the kept family is asked (a planted helper inside the test file, since base `p.adjust` cannot be mocked); the 0.075 check must fail |
| T4 d-ii | `kept_family_formula` (0 < a < 1), `kept_family_unit_a_counterexample`, `kept_family_instance` (witness: the paper's instance) | none (no R object; the figure is KIT-189's) | none | none. `kept_family_unit_a_counterexample` is about the formula's range (a = 1), not about BH. The statements that show BH over the kept terms has no general guarantee are `fdr_n20`, `fdr_n25` and `kept_family_instance`, each under a constructed law. KIT-189's measurement on the real route is a different result: BH over the kept terms failed in none of its 18 settings (VAL-123) |
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

Notes on hypotheses and witnesses (proofs audit P04):
- R-C1: the conclusion also holds without `Hac` (with Lean's junk inverse both `W` and `W_c` are 0); `Hac` is the
  code's refusal (`solve()` fails on a singular matrix), kept by ticket 5.1.
- R-G2: `IsUnit (Xᵀ * X).det` is implied by the full-rank `[X, C]` and is not used; `G2_sign` is the instance with a
  non-zero γ (γ = 1; ratios 4/3 and 2/3).

Further notes (second proofs audit Q02):
- R-B1, R-B2, R-E1: `Hset` and `Hac` are not used by the proofs (the shift argument needs `HU` and `Hone`); each
  conclusion also holds without them. `Hset` is true of every replicate matrix the code builds; `Hac` is the code's
  refusal. Kept by ticket 5.1.
- R-D1: `Hac` is not used (without it both sides are `Y`). R-K1 (i): `Hset` is not used.
- `kept_family_formula`: `0 < m` is not used (at `m = 0` both sides are 0).
- R-K1's witnesses have one replicate set and `Y_c Y_c′ = I`; R-J1's has `k = 0` (second proofs audit, Q01).
- R-I1, `I1_ctl_rank`: witness `A2_witness` with F = `univ` (`ctl ⊆ univ` for every computation); with F = {0},
  `witnessE2` (as `I1_witness` uses it: its rank clause is `I1_ctl_rank` at F = ctl = {0}).


## 3. R-I1, R-J3, R-N1: their own statements review (departure from ticket 5.3 item 1, recorded in r1.3)

Each has its own Spec file; the files frozen at the gate keep their hashes. Drafted on 2026-11-24 (b32c94ca), reviewed
in `tickets/kit/audits/AUDIT-KIT-193-proofs2-20261124.md` section B, edited as its Q06 to Q08 write (the optional
statements included) and committed for the freeze at 2ab4c189 with their three `SPEC.sha256` lines replaced (Q09);
frozen when the auditor confirms that commit by reading. None is proved yet.

| Item | Spec file, statements | What it transcribes | Pair: kind | Planted defect (proposed) |
|---|---|---|---|---|
| R-I1 | `Spec/RemovedFactors.lean`: `I1_UtWc`, `I1_centred_diff`, `I1_centred_diff_rows`, `I1_span`, `I1_rank`, `I1_ctl_rank`, `I1_average` (witnesses `I1_witness`, `I1_average_witness`) | KIT-184 r1.4 5.2-5.4's algebra over the RUV-III model: `D` (uncorrected minus corrected, each feature centred over the runs) is `W_c α`, also on any selection of the runs (centred rows of `W` times `α`); its span lies in `col(W_c)`, equal when `rank α_F = k`, and then of rank exactly `k`; `rank α_F = k` whenever `F` holds the controls; averaging runs (`integrationAverageRuns()`, `R/utils_workflow_integration_model.R:647`) keeps each feature's removed part in `span[1, A W_c]`. Not stated: a feature missing only in a run removed after RUV (the writer's status `approximate`); the writer's numerical rank rule (`1e-8 · s_1`, status): it waits for KIT-184's merged code | package: `integrationInputRemovedFactors()`, `integrationAverageRuns()` (after KIT-184 merges) | proposed: the difference not centred per feature; the span check must fail |
| R-J3 | `Spec/RuvIIICPattern.lean`: `RuvIIICBr` (residual dimension, eigen-step flag, closed-form values abstract), `J3_pattern`, `J3_counts` (witness `J3_witness`) | RUVIIIC at c534151f: all missing when `min(d_j, #C_j) < k` (`:168-171`), closed form at `d_j = k` (`:173-184`), general branch otherwise, all missing when fewer than `k` eigenvalues converge (`:252-263`) or the eigen step throws (`:304`); values on the observed runs only (`:186-198`, `:286-301`); callers drop all-missing features (`R/func_prot_s4_norm_methods.R:264-265`, `R/func_lipid_s4_ruv_methods.R:713-715`). `J3_counts` is the missing-value part of the clean-up's count only; the intensity threshold it recomputes (`R/func_prot_qc_filtering_helpers.R:243-247`) is not covered | library: `RUVIIIC::RUVIII_C_Varying` on the missing-value fixture; package: the clean-up's counts | proposed: a corrected value written at a missing run; the pattern check must fail |
| R-N1 | `Spec/RuvControlRules.lean`: `sweepControls`, `rankFirst`, `firstPassControls`, `orderPos`, `trxFirstPassControls`; `N1_stable`, `N1_ties`, `N1_first_pass`, `N1_trx`, `N1_rules_differ` (two parts) (witness `N1_stable_witness`) | `getNegCtrlProtAnovaHelper` (`R/func_prot_norm_optimization_helpers.R:565-582`; the rule `q_table_order`, which stays for proteomics and for the replay of lipid and metabolite records made before KIT-190; KIT-190's `p_id_two_groups` is not modelled), `protRuvFirstPassControls` (`:1380`), `omicsRuvFirstPassControls` (`R/func_omics_ruv_factors.R:54`), `trxRuvDesignControls` (`R/func_trx_ruv.R:296-297`), read on `v0.9-staging` at 19997323; R's `rank(ties = "first", na.last = "keep")`, `order()` and `x[1:n]` as stated assumptions | package: `getNegCtrlProtAnovaHelper` (equal q-values, tied blocks, a count of 0, no candidate; after KIT-190 merges, the equal case under the option `dropped`), `protRuvFirstPassControls`, `omicsRuvFirstPassControls`, `trxRuvDesignControls` on an oracle table `first-pass.tsv` | proposed: `ties.method = "last"`; the table must differ |

## 3a. Pairs run in P2/P0b (run 45335802998): departures from the proposals above

Every planted defect above was built and caught, except as follows (mechanics; no expected value changed;
`outcome.md` P2/P0b): R-A3 reaches k = 0 through `ruv_k_rule = "num_sv_be"` (the user rule refuses k = 0); R-F1's
"`Y - 2 * W alpha`" uses the uncentred `RuvIII.W` (with `W_c` the defect stays in col(W_c) and cannot fail); R-M1's
rank-1 map is one column `W q` (a two-column rank-1 design is refused by edgeR); 4 (i) adds a planted
`qvalue::qvalue` ignoring pi0 (4 (ii), with `processEnrichments`, waits for KIT-203); R-N1's first-pass planted
defect is `protRuvFirstPassControls` with `ties.method = "last"`. The target 5 oracle (`Oracle/Permutation.lean`)
counts over lists with decidable equality; its bridge lemma to the Spec's classical `unchangedDraws` and
`reflectedDraws` is owed to the Proof files (ticket 5.5).

## 4. The gate's own reads (ticket 2.3 and 14): done at the gate

- RUVIIIC at the locked commit c534151f and `ruv` 0.9.7.1 were read by the statements gate
  (`AUDIT-KIT-193-statements-20261124.md` section 5); the docstrings cite those lines (S07).
- RUVSeq: no statement depends on it; R-M1 is linear algebra about any W.
- `bindings.tsv` rows are written when each pair first runs on merged code (8.2). Filled: T1's four (P1a's script);
  T1's helper, target 3's five, target 5's, target 4's libraries and group R's (run 45335802998, `bindings-t3.R` and
  `bindings-p.R`). Pending: target 2's five (KIT-190) and T4 `processEnrichments` (KIT-203).
