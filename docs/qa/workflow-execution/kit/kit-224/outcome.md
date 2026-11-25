# KIT-224 outcome (r1.2; stopped at R2 by stop rule 2, before the full suite)

Branch `kit-224` from `origin/v0.9-staging` cf8af245. Commits: cea72da6 (tests first, nothing under `R/`), 12205bd4 (the
fix, NEWS, ledger), 5f261d89 (stored references accepted by file, OD-4), the provenance re-dump of
`test-omic-norm-module.R` and these records. Decisions: owner, 2026-11-25 ("yes to all, and extra arm"; OD-0 to OD-5).

## Observed

**G0 and R1** (run 45508928818, tests-only commit cea72da6; `g0/`, `r1/`). Every E row failed on base at its stated
value; E4 (ii), G1, G2 passed (stop rule 1 not triggered):

| Row | Base value |
|---|---|
| E1 lipid, complete | largest level change 10.23 (pos), 11.59 (neg) log2 |
| E1m lipid, 5 % missing | 7.56, 5.98 |
| E2 metab | 10.23, 11.59 (identical to lipid) |
| E3 lipid_ruv Automatic (KIT-172 seed 172101) | 22.37; record lists no marker |
| E4 (i), (iii), (iv) | marker absent; the two rules' corrected tables identical |
| E5 transcriptomics | file and integration input: largest gene-mean change 2.431; `correction$levels` absent |
| E6 (T1b seed 172102, x 2^10) | 22 calls of 516 change, max abs(delta t) 0.4205, max abs(delta logFC) 6.4e-13 (equal to KIT-172 r2) |

G0's site E5 stopped in the KIT-040 corpus manifest check (laptop-recorded digest; the runner switch
`MULTISCHOLAR_TRX_MANIFEST_NOTE=1` was not set). R1, on the same commit, gives E5's base values (above), so G0 was not
re-sent.

**R2** (run 45510911807; 12205bd4, tests unedited, against cea72da6; `r2/`). R CMD check: identical findings on both
sides (0 only on head). The new file: 69 pass, 0 fail (`fixed`). `test-kit129-run-path.R`: `same` (711 pass).
`test-prot-ruv-levels.R` (G3), `test-trx-kernels.R` (G4): `same`. Moves:

| File | Verdict | What moved | Approved (OD-4)? |
|---|---|---|---|
| test-lipid-12c-ruv-s4-characterization.R | new-failure 10 | W1 (8), forwarding (2) | yes |
| test-metab-03a-ruviii-c-varying-characterization.R | new-failure 10 | W1 (8), forwarding (2) | yes |
| test-kit127-omics-ruv-factors.R | new-failure 4 | :241-242 | yes |
| test-omic-norm-module.R | new-failure 12 | the 12 `provenance` fields only (the dumps differ in exactly those 12 lines, `norm/base-head-diff.txt`) | yes |
| test-integration-input-transcriptomics.R | new-failure 2 | :80, :120 | yes |
| test-workflow-trx-declarations.R | message-changed | the export file's SHA-256 at two paths | yes |
| test-module-ci-trx-summary.R | message-changed | the export file's SHA-256 | yes |
| test-workflow-omics-declarations.R (KIT-073) | 22 -> 23 failing, message-changed | registry `fixed`; states lipid/metab `_ruv_corrected` and `_correlation_filtered`; correlations; DA results. No earlier state | deferred; R1-01 holds |
| test-workflow-omics-binding-steps.R (KIT-094) | message-changed | `records`: the RUV step's parameters gain the marker (9 -> 10) | deferred; R1-01 holds |
| **test-prot-norm-core.R** | message-changed (306 failing on both) | the pinned digests of the `lipid_ruv` and `metab_ruv` declaration results | **no** |
| **test-omics-summary-core.R** | message-changed (25 failing on both) | trx_main and trx_managed: the digest of `integration/transcriptomics_normalized_data.tsv` | **no** |
| **test-omic-summary-module.R** | message-changed (28 failing on both) | the same file digest, trx_main and trx_managed | **no** |

Every other file: `same`.

**S1 units, criterion (vi)** (`s1u/`). Passed. Base: 26 of 84 rows change calls (up to 22; max abs(delta t) 1.247), the
same as KIT-172 T1b. Head: 0 of 84 (max abs(delta t) 2.9e-11, max abs(delta logFC) 1.1e-10).

**S1 synthetic** (`s1/`). Instrument defect, no verdict. The declared-path part (v) raised
`omics_ruv_correction_failed` on every `k185` object on both sides. It shared one tryCatch with the journey, so all 50
`k185` results were lost, and with them (iv) and (v). On law `level` (49 objects): (i), (ii) and (iii) passed in every
row. At cutoff 0, 36 of 49 objects changed at least one call; at cutoff 0.5, none did. The correlation filter kept the
same samples in all 49. Object 16 (`level`) failed on the head only, downstream of RUV:
`Assay 'LCMS_Neg' has sample column(s) missing from design: A_6`. The correlation filter removed a sample that was
missing in one assay only (see the unrelated finding below).

## Inferred

- The three unapproved moves are the decided change seen through stored references the ticket did not list. The
  lipid/metab RUV result digests and the transcriptomics export file digest move because the values were restored.
  All three files already fail on runners on base, so only the messages change.
- The fix is unit-invariant on complete data: (vi), E6.

## Acceptance criteria

- A1: met (G0 plus R1 for E5).
- A2: met. R1 failed at the stated values; R2 head passes.
- A3: partly met. (vi) is met; (i) to (iii) are met on law `level` only; (iv), (v), (a) and (b) await the S1
  instrument rerun.
- A4: not met. Three unapproved files moved (stop rule 2). The edited tests (5f261d89, re-dump) have not yet been run.
- A5: not run.
- A6: NEWS, ledger row and M8's status are done; `man/`, `NAMESPACE` and `DESCRIPTION` are unchanged.

## Hand-off to KIT-188 (OD-5 (c))

The un-repaired legacy behaviour is reached by running the legacy route (Automatic or Manual) under
`options(MultiScholaR.ruv_feature_levels = "shifted")` (the replay marker's earlier token, `workflowReplayRuleMarkers$ruv_feature_levels$earlier`). Under it, `ruvIII_C_Varying` for lipid and metabolite data returns RUVIIIC's values exactly (E4 (ii); S1 (iv) once rerun). The option must be set in the process that computes the RUV step. Unset, or `"restored"`, gives the repaired path.

## Unrelated finding (not fixed)

The lipid correlation filter (`filterSamplesByLipidCorrelationThreshold`, through `binding$filter_correlation`) can remove
a sample from the design while one assay keeps its column. The object then fails validity at the DA
(`R/func_omics_design_shared.R:142`). This was seen once in S1, on the head; it does not come from this change.

## Not checked

The full suite (R3, shared with KIT-225); the edited tests on the fix; S1 (iv), (v), (a), (b); Windows; real data.

## For the change report

- **What was wrong.** Legacy lipidomics and metabolomics RUV-III-C (Automatic, the shipped default, and Manual) shifted
  every feature by a constant, s_j = colMeans(W)' alpha_j. That constant depends on the unit of measurement. limma-trend's
  prior therefore read shifted mean abundances, so moderated t, p, q and calls depended on the unit. Fold changes did not.
- **Who was affected.** Every lipid or metabolite analysis that applied RUV, on `main` (metabolomics since 2025-05-06,
  lipidomics since 2026-01-18) and on staging. The transcriptomics RUVIII export file and the RUV-corrected integration
  input carried the same kind of shift; the transcriptomics DA did not.
- **Code.** Base cf8af245: `R/func_lipid_s4_ruv_methods.R:711`, `R/func_metab_s4_ruv_methods.R:711` (no restoration);
  `R/utils_workflow_trx_declarations.R:949-950`, `:1014`. Fix 12205bd4: the same lines plus `restoreRuvFeatureLevels`
  under the marker `ruv_feature_levels`; `R/func_omics_ruv_factors.R:361-364` (the declared path restores once);
  `trxRuviiiExportRestoredValues`.
- **Evidence.** KIT-172 T1b and S1 (vi): a change of unit changed calls in 26 of 84 complete-data rows before the fix
  and in 0 of 84 after (run 45510911807). E6 (seed 172102): 22 calls before, 0 after.
- **Limits.** Missing values with zero-filling, and the Automatic k chosen from the uncentred cancor plot, still depend
  on the unit (KIT-185, KIT-188, KIT-210).
