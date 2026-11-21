# RUV options at association 0.6: results (2026-11-21)

Declaration: `confirm06-declaration.md` (written and pushed before any job). Package under test: `kit-129-integration`
44b5f2195e355bb6e6efe1dae21905e30b81c26c, installed in every job. Scored once with `aggregate.R` sha256 `089d852c...`
(the confirmation round's), named in `score.json` before any result was read. Tables in `confirm06-gha/` (layout of
`confirm-gha/`); the raw per-run tables are not committed (sha256 in `confirm06-gha/not-committed.sha256`).

## Runs

| Step | Run | Outcome | Core-hours (job time x 4 cores) |
|---|---|---|---|
| Inert branch push (placeholders) | 44790983069, 44790983023 | both stopped at their first step, as intended | 0 |
| Lifecycle `ruv06-lc` | 44791180548 | passed (declaration section 10) | about 8.7 (one job waited 30 min in a runner library step) |
| Lifecycle scoring | 44799035290 | success | under 0.1 |
| Counted `ruv06` | 44799202191 | 5 simulator jobs, all success; gate status 0; 1 h 40 min wall | 22.1 (18.9 of simulator compute) |
| Scoring `confirm06-gha` | 44814537401 | success | under 0.1 |
| **Total** | | no correction needed | **about 31 of the 50 declared** |

## Checks

- **Gate:** S4's S3 agrees with the committed laptop rows in every column. S1, S2 and S3 have no reference for this cell
  (declared).
- **Arms 1 and 2 (reported, not gating):** in every simulator the cell rows of arms 1 and 2 (both subsets where earlier
  rows exist) are identical to the earlier scored rows of the same seeds (`final-20261118/s1`, `s2`; `confirm-gha/s3`,
  `s4`) in runs, k, features, null rates, mean FDP, share above 10 %, power, bias and rmse.
- **Refusals and crashes:** none. 14,500 rows, every one `ok`; the refusal share is 0 % for every figure below (in
  particular the lipid/metab standards never reached k = 5, so no assay ran without RUV).

## The 13 figures (mean false discovery proportion at a 5 % target, subset `full`, 500 runs, extra-hard cell S3)

"Rules failing at 0.6" lists the pass rules this arm misses in the extra-hard cell (they do not fail a design; rule 1
null rate, 2 mean FDP above 0.075, 3 share of data sets above 10 % FDP against the true-factor analysis, 4 reduction
against no correction). Neighbours: the omic's recommended design and no correction at 0.6, from the current table.

| Omic | Option (arm) | At 0.6 (SE) | Data sets > 10 % FDP | Power | Rules failing at 0.6 | Refused | Same option at 0.3 (current) | Recommended at 0.6 | No correction at 0.6 |
|---|---|---|---|---|---|---|---|---|---|
| Classic | recommended (34) | **6.3 %** (0.12) | 8.0 % | 84 % | none | 0 % | 5.8 % | (this row) | 59.3 % (new) |
| Classic | supplied list (35) | **6.0 %** (0.13) | 7.6 % | 84 % | none | 0 % | 5.4 % | 6.3 % | 59.3 % |
| Classic | first-pass nulls, separation rule (40) | **28.0 %** (0.35) | 98.8 % | 88 % | 1, 2, 3 | 0 % | 7.9 % | 6.3 % | 59.3 % |
| Classic | no correction, M30 test (37) | **59.3 %** (0.25), shown in all 3 classic check sentences | 100 % | 76 % | 1, 2, 3 | 0 % | | 6.3 % | (this row) |
| DIA limpa | recommended (34) | **6.6 %** (0.13) | 9.4 % | 79 % | 3 | 0 % | 5.3 % | (this row) | 53.3 % |
| DIA limpa | supplied list (35) | **7.3 %** (0.21) | 17.8 % | 78 % | 3 | 0 % | 5.5 % | 6.6 % | 53.3 % |
| DIA limpa | first-pass nulls, separation rule (40) | **30.2 %** (0.63) | 89.8 % | 86 % | 1, 2, 3 | 0 % | 15.1 % | 6.6 % | 53.3 % |
| Transcriptomics | supplied list (35) | **4.8 %** (0.06) | 0.4 % | 54 % | none | 0 % | 4.8 % | 5.2 % | 45.0 % |
| Transcriptomics | first-pass nulls, separation rule (40) | **9.6 %** (0.78) | 41.6 % | 27 % | 1, 2, 3 | 0 % | 6.4 % | 5.2 % | 45.0 % |
| Lipid/metab | supplied list, standards (35) | **13.6 %** (0.91) | 41.6 % | 12 % | 1, 2, 3 | 0 % | 7.9 % | 14.8 % | 64.2 % |
| Lipid/metab | first-pass nulls, separation rule (40) | **56.2 %** (0.42) | 100 % | 48 % | 1, 2, 3 | 0 % | 26.1 % | 14.8 % | 64.2 % |

SE in percentage points. For reference, the analysis given the true factor (arm 2) at 0.6: classic 5.0 %, DIA 4.7 %,
transcriptomics 4.8 %, lipid/metab 7.2 %. Every figure is within 1 point of the earlier-method figure the table
recorded as "was", except lipid/metab supplied list (12.0 -> 13.6 %) and first-pass nulls (53.7 -> 56.2 %), which
move as the confirmation round moved the other lipid/metab figures. Successor table:
`docs/qa/workflow-execution/kit/kit-129/final-figures-06.tsv` (these 13 rows filled; the other 143 rows unchanged).

## What contradicts a text shown in the app (for the user; nothing changed)

1. **DIA with limpa, recommended approach.** The preset note says "In every scenario tested it did as well as an
   analysis that was given the true unwanted factor (5.3% false discoveries with unwanted variation tied to the groups,
   correlation 0.3)" (`R/utils_omic_gui_config.R`, DIA `groups_all_features` note). At 0.6 it does not: 6.6 % against
   4.7 % with the true factor, and 9.4 % of data sets above 10 % false discoveries against 0.6 %, so it misses rule 3
   in this cell. The note's "every scenario tested" was true when 0.6 had not been tested for this arm; it is no longer.
   Decision for the user: reword the note, or keep it with the scenario named.

No other figure contradicts a current recommendation or a note:
- Classic recommended holds at 0.6 (no rule fails).
- First-pass nulls get worse in every omic (28 to 56 %), which agrees with the cautions. These still quote the 0.3
  figure, and classic's caution says "more" where 28 % at 0.6 is far more (a wording point, not a contradiction).
- The lipid/metab standards (13.6 %) are not clearly better than the recommended design (14.8 %; the difference is
  about one standard error) and have about half its power (12 % against 21 %). The recommendation is not contradicted.

## Not done here

Nothing in `R/` changed. Updating `ruvGroupAssociationConsequences()` (the `NA_real_` at 0.6 and classic
`uncorrected_06`) and its tests, the ledger row, and any wording change are for the main session and the user.
