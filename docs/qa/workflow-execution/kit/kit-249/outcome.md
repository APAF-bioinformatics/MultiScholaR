# KIT-249 outcome (part 1: G0, tests-only commit, R1)

Ticket: `tickets/kit/KIT-249.md` r1.2. Branch `kit-249` from `v0.9-staging` at **16112917329a63ab062655f04db5bf1dd7df9729**
(merge train 3). G0 branch `kit-249-g0` (not merged). **The fix commit is not written: it waits for the owner's OD-1
(older records on replay), which sets its file list.** R2, A3 to A5 and the change-report evidence follow the fix.

## What was run

| Run | Commit | What | Job minutes |
|---|---|---|---|
| G0 45590869970 | dcbce79f (16112917 + `tools/kit249/g0.R`) | scripts job, installed package, no package change | 10.2 |
| R1 45592793437 | 5f072f85 (tests-only: `tests/testthat/test-kit249-prot-da-qvalue-fallback.R`) | the new file alone | 12.1 |

Runner: ubuntu-22.04, R 4.6.1, qvalue 2.44.0, limma 3.68.4, reference BLAS (`g0-45590869970/provenance.txt`).

## Observed

**G0** (`g0-45590869970/g0-routes.tsv`, per protein `g0-proteins.tsv`, the matrix `g0-matrix.tsv`; c shift 0, so
the fixture rule of section 8 was not needed):

| Route | B_vs_A max p | qvalue() / bootstrap | non-NA fdr_qvalue | planted sig. (of 6) | flagged | C_vs_A |
|---|---|---|---|---|---|---|
| G0-a runTestsContrasts topTable | 0.756 | both error | 0 of 30 | 0 | yes | Storey, 30 of 30, 6 of 6, 0 nulls |
| G0-b runTestsContrasts topTreat (lfc 0) | 0.754 | both error | 0 of 30 | 0 | yes | the same |
| G0-c S4 method, LFQ object | 0.754 | both error | 0 of 30 | 0 | yes; count table "Not significant: 30" | the same; "Significant and Up: 6" |
| G0-d application kernel in process, DIA owner (publish on) | 0.754 | both error | 0 of 30 (da_proteins_long 0 of 30) | 0 | manifest names B_vs_A; record `fdr_method` "BH"; notification "Some contrasts used the existing BH fallback after q-value estimation failed." | the same; record "qvalue" |
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

## Inferred

- The defect is as section 2.1 reads it on every claimed route, including the scripted wrapper and the lipid/metab
  exported S4 method (stop rule 1 did not fire; no reading correction for G0-e or G0-f).
- The fixture margin holds: the largest null p is 0.756 against 0.95, so moderation did not bring any null near 0.95.

## Acceptance criteria

| # | State | Evidence |
|---|---|---|
| A1 | met | G0 45590869970 (table above) |
| A2 | first half met (every listed row fails by assertion at G0's value; guards pass); the second half (all pass on the fix) waits for the fix | R1 45592793437 |
| A3 | not yet (R2 follows the fix; 6 shards per side, section 9 r1.2, `r2-shards.py`, `r2-shards.txt`, `r2-files.txt`) | |
| A4 | not yet (the merge train's or checkpoint's full suite) | |
| A5 | not yet (NEWS, ledger row and the change-report section come with the fix) | |

## R2 shard projection (R1-01)

`r2-shards.py` (stdlib Python) projects each shard from the per-file seconds of full suite 45573760967 (and train 3's
45573711299 for its three new files), 12 minutes set-up per job, files heaviest first, round-robin: 6 shards per side,
largest projected 45.6 minutes (`r2-shards.txt`). `r2-files.txt`: the 74 files at 16112917 (67 by section 9's grep + 7 named).

## Not checked

- The fix and everything after it (R2, A3 to A5, the "For the change report" section): waiting for OD-1.
- `daAnalysisWrapperFunction` (KIT-218 not merged).
- The managed path and the record were not part of G0 (G0-d in process only); E5 and E6 cover them in R1.
- G0-d's published outputs (Glimma pages, figures) were produced but not read.
- Real data: not read.
