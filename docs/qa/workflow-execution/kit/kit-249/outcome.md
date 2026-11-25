# KIT-249 outcome (G0, R1, the fix, R2: stopped at stop rules 3 and 4)

Ticket: `tickets/kit/KIT-249.md` r1.3. Branch `kit-249` from `v0.9-staging` at **16112917329a63ab062655f04db5bf1dd7df9729**
(merge train 3). G0 branch `kit-249-g0` (not merged). OD-1 was decided by the owner on 2026-11-25: (a), older records
recompute on replay (VAL-237). Staging e0747fc1 (merge train 4, KIT-218 merged) was merged into `kit-249` (47ad12d0);
E7 (b) for `daAnalysisWrapperFunction` was added as a tests-only commit (db30cf75); the fix is 82470db1.
**R2 stopped the ticket at stop rules 3 and 4: a stored reference moves (section "R2" below); the owner rules by file (OD-2).**

## What was run

| Run | Commit | What | Job minutes |
|---|---|---|---|
| G0 45590869970 | dcbce79f (16112917 + `tools/kit249/g0.R`) | scripts job, installed package, no package change | 10.2 |
| R1 45592793437 | 5f072f85 (tests-only: `tests/testthat/test-kit249-prot-da-qvalue-fallback.R`) | the new file alone | 12.1 |
| R2 45601895895 | head 82470db1 (fix) vs base db30cf75 (tests-only + E7 (b) on staging e0747fc1) | 77 files, 6 shards per side, R CMD check | 394.0 |
| Correction 45614191837 | head 144d2c5d (fix + G1b, staging 63143aa3) vs base 968804de (same without the fix) | the same 77 files, 6 shards per side, R CMD check | 417.4 |

Runner: ubuntu-22.04, R 4.6.1, qvalue 2.44.0, limma 3.68.4, reference BLAS (`g0-45590869970/provenance.txt`).

## Observed

**G0** (`g0-45590869970/g0-routes.tsv`, per protein `g0-proteins.tsv`, the matrix `g0-matrix.tsv`; c shift 0, so
the fixture rule of section 8 was not needed):

| Route | B_vs_A max p | qvalue() / bootstrap | non-NA fdr_qvalue | planted sig. (of 6) | flagged | C_vs_A |
|---|---|---|---|---|---|---|
| G0-a runTestsContrasts topTable | 0.756 | both error | 0 of 30 | 0 | yes | Storey, 30 of 30, 6 of 6, 0 nulls |
| G0-b runTestsContrasts topTreat (lfc 0) | 0.754 | both error | 0 of 30 | 0 | yes | the same |
| G0-c S4 method, LFQ object | 0.754 | both error | 0 of 30 | 0 | yes; count table "Not significant: 30" | the same; "Significant and Up: 6" |
| G0-d application kernel in process, DIA owner (publish on) | 0.754 | both error | 0 of 30 (da_proteins_long 0 of 30) | 0 | manifest names B_vs_A; record `fdr_method` "BH"; notification (the handler's text, printed by G0 when the manifest flags a comparison; the tab itself was not driven) "Some contrasts used the existing BH fallback after q-value estimation failed." | the same; record "qvalue" |
| G0-e protein_deAnalysisWrapperFunction (completed) | 0.754 | both error | 0 of 30 | 0 | yes | the same |
| G0-f lipidomics / metabolomics exported S4 | 0.754 | both error | 0 of 30 | 0 | yes (kernel) | the same |
| G0-c1 control: adjustProtDaResultTable | 0.754 | both error | 30, identical to BH | 6 | `qvalue_failed` TRUE | n/a |
| G0-c2 control: applyLipidDaContrasts | 0.756 | both error | 30, identical to BH | 6 | yes | Storey |

- qvalue()'s message on B_vs_A: "missing or infinite values in inputs are not allowed" (the message of section 2.1's logs); the bootstrap: "missing values and NaN's not allowed if 'na.rm' is FALSE".
- qvalue(p, lambda = 0.5) gives pi0 0.467 on B_vs_A.
- BH on the same raw p-values calls 6 of 6 planted and 0 of 24 nulls on every route (smallest null p 0.085 topTreat, 0.118 topTable; largest planted p 4.9e-20 and 1.1e-17).
- `daAnalysisWrapperFunction` was not run: KIT-218 has not merged at 16112917.

**R1** (`r1-45592793437/tests.tsv`, `failures.txt`): 76 passed, 49 failed, 0 errors. Every failure is an assertion of R1's
failing-first list at G0's value: E1 (ii) 30 missing, (iii), (iv), (v) 0 of 6; E1b (ii) to (iv) for B_vs_A and C_vs_A
(qvalue replaced by a double that errors); E2 the same as E1; E3 (ii), (iii), (v) and the count table 0 "Significant and
Up"; E4 on DIA, LFQ and TMT owners (results (ii) to (v), da_proteins_long (ii) and (v)); E5 managed the same; E6 second
half; E7 (ii), (v); E8 lipidomics and metabolomics (ii), (iii). Every guard passed: E1 (i) and (vi), E2 (i) and (vi), E6
first half (B_vs_A "BH", C_vs_A "qvalue"), G1 (C_vs_A equal to `qvalue::qvalue(raw_pvalue)$qvalues`, unflagged, 6 of 6
significant) in E1, E2, E4 (three owners) and E5.

**R2** (`r2-45601895895/compare.tsv`, `check-compare.txt`, `tests*.tsv`, `saved-da-new.txt`): 77 files.
- 75 `same` (including `test-module-ci-prot-da.R` and `test-prot-07d-da-model-stats-shared.R`; neither falsifies a fix that
  assigns BH on the success path: 07d's double sits in a test quarantined as Q-11 and module-CI's asserts no
  success-path value. The falsifier is G1: on C_vs_A Storey's q is not bit-identical to BH on any route (G0:
  `identical_q_bh` FALSE, `q_equals_storey` TRUE), so G1's exact identity fails under BH, at rounding level; `test-kit218-da-wrapper-tail.R`; `test-kit220-da-threshold-consumers.R`; `test-workflow-analysis-record.R`).
- The new file `fixed`: base 82 passed, 54 failed (including E7 (b), which therefore fails first on the merged tree (committed 42 seconds before the fix and first run in R2, beside it; its base values are G0-e's and KIT-218 A1's, and the fix was not changed for it));
  head 136 passed, 0 failed. G1 holds on head (E1, E2, E4 on three owners, E5).
- **`test-workflow-saved-da.R` `message-changed`: base 780 passed / 171 failed, head 759 / 192.** Every difference is in
  the KIT-112 W1 (a) and W2 (W1 (b)) frozen captures (`savedDaW1Frozen()`, `:764` onward), frozen on the laptop at
  64e7c2ae (2026-11-10), a commit with this defect.
  - 21 expectations fail only on head: `entities` in nine W1 (a) cases and their nine W2 twins (one or two
    `rendered_asset` more, nothing else), and `calls` in three W1 (a) cases (`saveDIAStateValue` 10 -> 11 in
    `dia_one_contrast`, 7 -> 8 in `dia_two_contrasts`, 8 -> 9 in `dia_retry`; W2 does not compare that count).
  - 87 expectations that already fail on runners on both sides change their message: `products`, `dependencies` and
    `files` of those nine cases, and `manifest`, `head_manifest` and `display` of the six fallback cases.
  - The six cases whose comparisons fell back (frozen `fdr` "BH": `:786`, `:818`, `:850`, `:885`, `:920`, `:1091`)
    now carry Benjamini-Hochberg q-values: each fallback comparison publishes an interactive volcano page (a saved
    `rendered_asset`, one more state save) where the defect wrote a 'not available' note, and its contrast products,
    manifests and display change.
  - The three later cases at TREAT 0.5 (`:957`, `:996`, `:1035`; `fdr` "qvalue") move only in fields that accumulate
    over the fixture's store and folder (`entities`, `products`, `dependencies`, `files`). Every changed entry belongs
    to an earlier fallback click. Their own products, dependencies, manifest, head manifest, display and record are
    identical on both sides.
  - `limpa_one_contrast` (`:1135`, "BH" through `adjustProtDaResultTable`, which always kept the fallback) does not
    move: a control.
  - `states`, `fits`, `record` (without the fdr fields) and `fdr` are identical on both sides in every case. In W2 the
    worker equals today's managed path on head in every compared field.
  - **The stored captures encode the defect's result for the six fallback cases** (no q-value, no page).
- R CMD check: identical on both sides (6 findings; none only on head).
- Longest job 45.3 minutes (projected 45.7).

**Correction pass 45614191837** (declared in `r3-declared.md` before it was sent; `r3-45614191837/`): head 144d2c5d (the
fix, G1b, staging 63143aa3 merged), base 968804de (the same tree without the fix), the same 77 files (the grep at
144d2c5d gives the same 70; KIT-239's, KIT-214's, KIT-226's and KIT-241's new files do not match it), 6 shards per side,
R CMD check; 417.4 job-minutes, longest job 47.0.
- The new file `fixed`: base 90 passed / 54 failed, head 144 / 0. G1b passes on both sides. Declared 88/54 and
  142/0: my count of G1b's expectations was wrong (8, not 6: each branch's completion guard was not counted); the
  failing set (54) and the head's zero failures are as declared.
- `test-workflow-saved-da.R` `message-changed`, 171 -> 192: the same 21 newly failing and the same 87 changed leaves as
  R2 (I01), every failure message identical to R2's on each side; `limpa_one_contrast` unchanged.
- 75 files `same`; R CMD check identical (6 findings, none only on head).
- So R2's verdict stands on staging 63143aa3: the fix moves nothing but the fallback comparisons; the ticket waits for
  the owner's word on parked row 36. The ticket's one correction pass is now used.

## Inferred

- The defect is as section 2.1 reads it on every claimed route, including the scripted wrapper and the lipid/metab
  exported S4 method (stop rule 1 did not fire; no reading correction for G0-e or G0-f).
- The fixture margin holds: the largest null p is 0.756 against 0.95, so moderation did not bring any null near 0.95.

## Acceptance criteria

| # | State | Evidence |
|---|---|---|
| A1 | met | G0 45590869970 (table above) |
| A2 | met: every listed row fails by assertion at G0's value (R1; E7 (b) on R2's base); all pass on the fix | R1 45592793437, R2 45601895895 |
| A3 | **not met: stop rules 3 and 4.** Every file `same` except the new file (`fixed`) and `test-workflow-saved-da.R` (`message-changed`, a stored reference that encodes the defect); check identical | R2 45601895895 |
| A4 | not yet (the merge train's or checkpoint's full suite) | |
| A5 | NEWS and ledger row written (82470db1); `NAMESPACE`, `DESCRIPTION`, `man/` unchanged; change-report section below | 82470db1 |

## Decision needed (OD-2, stop rule 4)

The KIT-112 W1 (a) and W2 frozen captures in `tests/testthat/test-workflow-saved-da.R` (`savedDaW1Frozen()`) move on
head in exactly the leaves of AUDIT-KIT-249-impl-20261125 I01 (every move traced to the six fallback comparisons;
none where qvalue() succeeded; limpa does not move): they froze the empty q-value column. Put to the owner as parked row 36;
this branch does not edit that file. Re-recording them is the owner's word by file
(VAL-130 (4)); the ticket's section 3 and section 7 prediction ("none moves") missed this file's W1 captures. Whether
to accept the move (and re-record, or leave the laptop-recorded values failing as they already do on runners for most
fields) is OD-2. One correction pass remains in the budget; it is not used, because the result is valid, not an
instrument defect.

## R2 shard projection (R1-01)

`r2-shards.py` (stdlib Python) projects each shard from the per-file seconds of full suite 45573760967 (and train 3's
45573711299 for its three new files), 12 minutes set-up per job, files heaviest first, round-robin: 6 shards per side,
largest projected 45.6 minutes (`r2-shards.txt`). `r2-files.txt`: the 77 files recounted at the fix's tree (70 by section 9's grep, including the new file and KIT-218's
and KIT-220's, + 7 named); 6 shards per side by the ticket's rule (r1.2), largest projected 45.7 minutes;
`r2-shards.txt`'s "chosen k = 5" is the script's own rule (largest shard at most 50 minutes), not the ticket's.

## For the change report

- **What was wrong.** In the proteomics differential analysis (`runTestsContrasts`), when Storey's q-value could not be
  estimated for a comparison (qvalue 2.44.0: no p-value at or above 0.95), the Benjamini-Hochberg fallback was computed
  inside the error handler and assigned to a local copy, so the q-value column (`fdr_qvalue`) stayed empty for every
  protein of that comparison: nothing significant, an empty static volcano, no proteins for enrichment. The DA tab's
  warning and the analysis record said BH values were used. On `main` (and the tags v0.3.6.3, v0.5.0-rc1) the DA tab
  showed a modal naming the failed comparisons and stating "The analysis used Benjamini-Hochberg FDR correction
  (p.adjust) instead" (`main:R/mod_prot_da.R:1439-1441`), which was false in the same way; there was no analysis record.
  No notice to `main` users is drafted or promised.
- **Who was affected.** Proteomics DIA (iq), label-free and TMT on the default route (in process, managed, replay), the
  scripted wrappers, and the lipidomics/metabolomics exported S4 `differentialAbundanceAnalysis`, in comparisons with
  few unchanged proteins or where nearly every protein changes strongly (at most 0.95^m0 per comparison under
  independent uniform null p-values). Not affected: DIA with limpa, the declared RUV designs, lipid/metab in the
  application, transcriptomics.
- **Since.** a71b8618 (2025-11-12); on `main` at `R/func_prot_da.R:2642`, `:2763`; in v0.9.0 at
  `R/func_prot_da_model_stats.R:640`, `:762`; at the base 16112917 `R/func_prot_da_model_stats.R:686`, `:809`.
- **Before a71b8618** the same comparison stopped the analysis with qvalue's error; after it, the analysis completed
  with no q-values for that comparison. On `main` the DIA-with-limpa route has no handler and still stops with
  qvalue's error on the same trigger (`main:R/func_prot_limpa.R:1284`).
- **How an affected result is recognised** (every version since a71b8618): a comparison whose `fdr_qvalue` column is
  empty in every row while `raw_pvalue` is present; on `main` and the tags v0.3.6.3 and v0.5.0-rc1 also the modal
  naming it and the console line "q_values_all has 0 non-NA values"; on staging builds from v0.9.0, `fdr_method` "BH"
  for that comparison in `analysis_record.json`, `qvalue_warnings` in the DA products and, since KIT-101
  (2026-11-07), `Interactive_Volcano_Plots/<comparison>.unavailable.txt` in place of the interactive volcano.
- **Seen earlier and misread.** KIT-101 (F-08, 2026-11-07) found every committed proteomics fixture of 5 or 6
  proteins with an all-missing FDR and read it as q-values that "cannot be fitted on so few p-values"; its fix
  (DAS-066) let the publication complete with a 'not available' note. The KIT-112 frozen captures (2026-11-10) and
  KIT-220's runner values for them were recorded with the defect; KIT-218's W2 count (61 'Not significant') and
  KIT-204's B-A fixture showed it too.
- **After the fix, statistically.** Benjamini-Hochberg is Storey's q-value with pi0 = 1. It is used exactly when pi0
  is likely small, so there it is conservative by about 1/pi0, and the switch from Storey's q-value is discontinuous.
  Comparisons of one analysis can carry different methods; the record states the method per comparison.
- **The fix.** 82470db1, `R/func_prot_da_model_stats.R:674-691` and `:799-816`: the value of the inner `tryCatch` is
  assigned (as `adjustProtDaResultTable` does); `qvalue()` stays unqualified.
- **What changed in results.** Only comparisons where qvalue() fails: `fdr_qvalue` now equals `p.adjust(raw_pvalue, "BH")`
  over the valid p-values. On G0's constructed comparison (30 proteins, 6 planted +2 log2 effects): before, 0 of 6
  planted significant on every route; after, 6 of 6 and no null (R2's new file). Every other comparison is unchanged
  (75 files `same`). Older records recompute on replay with BH values (OD-1 (a), VAL-237).
- **Evidence.** G0 45590869970 (before, every route), R1 45592793437 (tests at the before values), R2 45601895895
  (after; the KIT-112 frozen captures move, OD-2). Limits: synthetic data only; frequency statement inferred, not
  measured; no real data read.

## Owed elsewhere (not done here)

- I10: the F112-03 comment `tests/testthat/test-workflow-saved-da.R:1155-1156` ("At q 0.05 without TREAT this fixture
  publishes no Glimma page") becomes false with the fix; it is NOT edited now (the file's captures are the owner's,
  parked row 36) and is owed at the re-record. `diaInteractiveVolcanoUnavailableText`'s "(the FDR could not be
  estimated, for example on too few proteins)" goes to the owner's wording (`wording.md` W4).
- I12 (conventions divergence): the fix's two new comment lines are 116 characters and keep the function's 2-space
  idiom; the new test file has lines over 100 characters and `function(` lambdas (R.md: `\(x)`), new code that
  diverges from `conventions/R.md`; recorded, not reformatted.
- I13: A4 (the next checkpoint's full suite) must read `test-integration-model.R`, `test-integration-mofa-model.R` and
  `test-prot-missing-value-counts.R`, which replay records and were outside R2; the four e2e proteomics lanes ran 0
  tests in R2 (skipped) and are no evidence.

## Not checked

- A4 (the merge train's or checkpoint's full suite).
- `daAnalysisWrapperFunction` in G0 (KIT-218 merged after G0; E7 (b) covers it in R2).
- The managed path and the record were not part of G0 (G0-d in process only); E5 and E6 cover them in R1.
- G0-d's published outputs (Glimma pages, figures) were produced but not read.
- Real data: not read.
