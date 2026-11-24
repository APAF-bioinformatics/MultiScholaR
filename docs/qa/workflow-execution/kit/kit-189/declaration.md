# KIT-189 declaration (before the counted run S)

Fixed before any counted row exists. Ticket `tickets/kit/KIT-189.md` revision 1.3 (sections 5, 7, 10), the owner's answers in VAL-86 (D1, D2 with D2a to D2d, D3, D4) and VAL-101 (the budget after G2). L1 (runs 45275836889, 45281701515) was blind: its rows and scoring tables were never uploaded (`lifecycle.md`).

## Bindings

- Package code: `v0.9-staging` at 3ff60190448150fa13e06d624378cad1481d3820; this branch adds only `docs/qa/workflow-execution/kit/kit-189/`. Every S, D and SC request binds `package_commit` to the commit that adds this file.
- Scripts (SHA-256): `sim-k189.R` 2fc0b3397c2250529a8a664cf4b457d71a33d7e2077d24f62140482684e3008a; `score-k189.R` 0764541fd8fed0f41400f77268c9731863586671ccf80d8c5b53074e4f805321; `e7-k189.R` f4d8ed869712a08c8cd64bcbe919efe1f51bfe0e18681286a0dd26797b15dbca.
- B2's entries: `b2-entries.tsv` SHA-256 664fecce8a6c169c5f1bc630e08c5b6a944314b0d6fa70b795e8756b51820d82 (2,000 accessions; `sim-k189.R` stops if its draw differs).
- Library: the run-tests snapshot of the bound commit (R 4.6.1; clusterProfiler 4.20.0, enrichit 0.2.0, qvalue 2.44.0, GO.db 3.23.1), reference BLAS, ubuntu-22.04 runners.

## Design

- **Cells** (index c): null, lists of 25, 126, 400 from the common universe U: B1 (E. coli full, 4,531 entries) 1 to 3, B2 (2,000-entry E. coli subset) 4 to 6, B3 (yeast, 6,733 entries) 7 to 9, R = 2,000 each. Planted, lists of 126, strata small (10 to 30 in U), mid (31 to 100), large (101 to 500): B1 10 to 12, B2 15 to 17, B3 20 to 22, R = 400; mid at 25 and 400: B1 13, 14; B2 18, 19; B3 23, 24, R = 200. Each list is run under `"none"` and `"is_a_part_of_regulates"` (paired).
- **Generator**: seed `189000000 + 100000 c + r` (Mersenne-Twister, Inversion, Rejection), r = 1 to R; null: L proteins uniform from U; planted: T uniform among the stratum's BP terms of the propagated sets (sizes in U), s = min(floor(L/2), floor(|T|/2)) members of T and L − s from U \ T.
- **Truth**: ER_V = E_V / (L |V| / |U|) with E_V = s |V ∩ T| / |T| + (L − s) |V \ T| / (|U| − |T|), sizes in U, compared in exact integer arithmetic; null term ER ≤ 1, non-null ER > 1, material ER ≥ 1.5; every term null under the null scenario.
- **Route**: `processEnrichments(da, "999999", q_cutoff = 0.05, pathway_dir = <new tempfile() directory>, go_annotations, exclude_iea = FALSE, protein_id_column = "Protein.IDs", contrast_names = "effect", go_propagation)`, the up result's `@result` and `@geneSets`; plot saving replaced by no-ops. A call is `crashed` (an error, or the route's own "Error processing up-regulated genes" warning), `empty` (no result: zero calls, valid) or `ok`.
- **Arms** (strict `< 0.05`; an NA q-value is not a call): A-Q `qvalue`; A-BH `p.adjust`; A-BHf BH over `c(p, rep(1, m − m'))`; A-Qf `enrichit:::calculate_qvalue` on the same vector; A-BY, A-Holm on the kept p; A-raw. Deciding: A-Q, A-BH, A-BHf.

## Estimands, intervals, line, classification

As ticket 5.5 and 5.6, computed by `score-k189.R`: (a) any-call rate (null; deciding), (b) false terms per list, (c) share with 10 or more calls, (d) mean false discovery proportion (planted; deciding at L = 126), (d′) at ER 1.5, (e) power (planted term called among lists where tested; material terms called; any non-null call), (f) agreement, (g) the mechanism record. Rates: Clopper-Pearson; means: t on R − 1 degrees of freedom. **Line 0.075** (D1). Per (cell, propagation, arm): **controls** when the one-sided 95 % upper limit is below 0.075; **fails** when the one-sided lower limit at 0.05 / 108 is above 0.075; otherwise **not determined**; more than 5 % crashed calls: not determined. 108 classifications: 18 null cells for (a) and 18 planted cells (L = 126) for (d), × 3 deciding arms; the `"none"` cells decide (D2c).

**Material power gain**: A-Q minus A-BH in planted-term recovery, paired, one-sided lower limit at 0.05 / 18 at or above 0.05 in at least one deciding planted cell. **Large-event condition** (D2d): in every one of the 18 null cells, the one-sided 95 % Clopper-Pearson upper limit of the share of lists (q-values not all NA) on which A-Q calls 10 or more terms and A-BH fewer than 10 is below 0.01.

## Outcomes and readings (D2 as answered, VAL-86; first matching row)

| Outcome | Condition | Default and options |
|---|---|---|
| Stop | E0, E1 or E2 failed (did not: L1) | — |
| O3 | BH fails in a deciding cell | BHf controls everywhere: BHf default, the q-value and BH offered; otherwise no default (owner, further ticket), all three stated |
| O1 | Q fails, BH fails in none | BH by rule F default; the q-value offered with its measured failure (not withdrawn) |
| O2 | Q, BH, BHf control everywhere | no material gain: BH by rule F default (BHf), the q-value and BH offered; material gain and the large-event condition: the q-value default, BH by rule F offered; material gain without the condition: BH by rule F default, the q-value offered with its large-event rate |
| O5 | BHf fails; Q and BH fail in none | BHf not recommended; Q and BH by O2's branches if both control everywhere, else O4's, with BH over the kept terms |
| O4 | otherwise | the q-value default only if Q controls everywhere, the gain is material and the large-event condition holds; otherwise BH by rule F; the other offered; every not-determined cell stated with its upper limit |

**Rule F** (D2a): BHf if it controls in every deciding cell, else BH over the kept terms. **Text T** for every option: generated by `score-k189.R` (`options.tsv`) in the form of ticket 5.7.

## Runs and attempts

24 S jobs as ticket 10 revision 1.3 (G2, VAL-101), D (replicates 1 to 10 of cells 2 and 11, compared with S's rows without `part` and `seconds`: E5), SC once on the committed S rows. A timed-out job's unwritten replicates run once as a continuation from a committed list (S′). Instrument defects: VAL-54. A run with a result is never re-run. A defect found in `score-k189.R` after scoring is fixed in a dated section and re-scored once.
