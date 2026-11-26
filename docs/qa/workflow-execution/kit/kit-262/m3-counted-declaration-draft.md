# KIT-262 M3: counted declaration, DRAFT for the auditor's read (nothing is sent from it; stage 1 and the counted jobs wait for the coordinator's word)

**Status.** L3r2 (run 45906333805) met every declared check (`l3r2/`; `runs.tsv`). This draft is the declaration of stage 1 of the re-estimation and of the counted jobs. It is to be committed alone, with the final script hashes, only after the auditor has read it and the open items of section 9 are closed. Owner decisions in force: VAL-302 (O-1 (ii): P and C1 only; O-2 (a); O-3: stress cells and skip rule; O-4 (a)), VAL-381 (stress cell unchanged, KIT-172's floor order), VAL-390 (ceiling 250.32), VAL-396 (M3 runs without the sample-correlation filter), VAL-400 (the zero-spread rule). No package code.

## 1. Binding and scripts

`package_commit` of every request is the commit that adds the final declaration to `kit-262`: `origin/v0.9-staging` a3c69901ea9b50b8de875d79ca89fb872912c092 plus `docs/qa/workflow-execution/kit/kit-262/` only (package diff empty, `m3-binding-a7.md`, including "A7 for E9's reference"). Scripts as of L3r2 (SHA-256 recomputed at the declaring commit): `m3-k262.R` 08eea38e..., `rules-k262.R` 9547d350..., `score-m3-k262.R` 0998cc59..., `reestimate-m3-k262.R` b6d18e54..., `e5hand-m3-expected.tsv` 3e78da07.... Environment in every request: `OMP_NUM_THREADS=1`, `K262_M3_WORKERS=4`.

## 2. Design (ticket section 6, as reduced by the owner's decisions)

Arms on the same data set: **P**, the shipped Automatic RUV step in KIT-172's journey; **C1**, R5 judged at min(k, cap), ranking unchanged, run as the Manual journey at C1's percentage and k where its choice differs, otherwise P's row. Cells and replicates: N-r 100, N-s 100, N-w 250, E-r 200, E-s 200, E-p 300, E-g 200, E-w 200 (VAL-235's seeds `190000000 + 10000 x idx + r`, KIT-190's generator, 1,550 data sets) and the two 24-sample stress cells S-N (null) and S-E (10 % of features at +-1 log2), 2 x 12 samples, 300 features, MCAR 10 %, KIT-172's floor order `after`, no `IS_` rows (ITSD off), seeds `262200000 + 10000 x idx + r` (idx 9, 10); R from stage 1 (minimum 250 and 200). Checks as in L3r2: E1, E7 (reduced), E9 on the eight cells (derived from `m3-binding-a7.md`; a difference not explained by it is F3: no counted job), E10, E11 (5 % of re-used units, r %% 20 = 7), E8 (replicates 1 to 5), F2 (stress cells: OBSERVE, reported as `F2_STRESS`).

## 3. Stage 1 (declared and sent first; blind)

One request, `m3-k262.R` with `K262_M3_STAGE1=1`, `K262_M3_REPS=1:50`, `K262_M3_E8=0` (cells S-N, S-E; rows private; `out` holds `reestimate.tsv` and `provenance.json`). Then `reestimate-m3-k262.R plan` gives R per stress cell (floors 250 and 200, multiples of 50, at most 1,000, target 0.80^(1/2), one-sided 0.05). **Route fixed now: if any R exceeds 400 the plan flags `above_400`, nothing is sent, the projected cost and the ceiling line go to the coordinator and, per VAL-390, to the owner if it passes 250.32; R is never reduced.** The final R is committed before any counted job; stage 2 re-runs seeds 1 to 50 and checks `reestimate.tsv` byte for byte.

## 4. Plan from L3r2's measured seconds (G2) and the running line

Measured (4 workers side by side, `l3r2/s3/timing.tsv`, `s5`): P journey 28 to 31 s on the eight cells and 21 s on the stress cells; a unit with all arms 43 s (eight cells) and 44 s (stress, with E10); L3's forced arms are an upper bound on a counted unit. Plan = units x seconds x 1.93 / 4 workers; jobs cut at 150 minutes of plan; `timeout_minutes` = 2 x plan, rounded up to 5.

| Job | Cells, replicates | Units (+ E8 units) | Plan (min at 1.93) | Timeout | Expected (min at 1.0) |
|---|---|---|---|---|---|
| Stage 1 | S-N, S-E, 1:50 (rows private) | 100 | 35 | 75 | 18 |
| J1 | E-p 1:300 | 300 | 104 | 210 | 54 |
| J2 | N-w 1:250 | 250 | 86 | 175 | 45 |
| J3 | E-r, E-s 1:200 | 400 (+10) | 142 | 285 | 73 |
| J4 | E-g, E-w 1:200 | 400 (+10) | 142 | 285 | 73 |
| J5 | N-r, N-s 1:100 | 200 (+10) | 73 | 145 | 38 |
| J6 | S-N 1:R (placeholder 250) | 250 (+5) | 90 | 180 | 47 |
| J7 | S-E 1:R (placeholder 200) | 200 (+5) | 73 | 145 | 38 |

Jobs J6 and J7 are recomputed from the final R. Each job: `K262_M3_CELLS`, `K262_M3_REPS` as listed (one replicate range per job; cells in a job share R), `K262_M3_E8=1`, `K262_M3_E8_MAX=5`.

Core-hours (job-minutes x 4 / 60, 10 minutes of set-up per job): stage 1 expected 1.89, 3.03 at spread 1.93; counted jobs expected 29.16 (437 job-minutes with 7 jobs), at spread 1.93 51.94, and with the stress cells at R doubled (two further jobs) 63.89; scoring 1.87 left of its line. **Running line against 250.32 (VAL-390):** used 87.24 + about 4.1 (L3, L3r, L3r2) = about 91.3; after stage 1 94.4; after the counted jobs at R doubled 158.3; after scoring 160.1; with both reserves on paper (correction 8.77, instrument 7.43 less the L3 overrun of about 0.1) 176.3: **about 74 core-hours of headroom**. The stress R at which the line would pass 250.32 is about 1,800 per stress cell (about 3,100 units beyond the doubled R); a stage-1 R above 400 is routed as in section 3 whatever the headroom. If the plan recomputed after stage 1 passes 250.32, it returns to the owner (VAL-390). Stage 1 is charged to M3's line.

## 5. Decision rules (ticket 6.2 and 7.1 reduced by O-1 (ii); one-sided alpha 0.05, no split)

Null cells (N-r, N-s, N-w, S-N): paired any-call, Tango's score test at margin 0.10, Berger-Boos beside it. Effect cells (E-r, E-s, E-p, E-g, E-w, S-E): paired mean FDP difference, sign-flip test (100,000 vectors, seed 262000002) at margin 0.05, Johnson's t when |skewness| > 1. **C1 is "no worse" (NW) when every one of the ten cells rejects.** Only "worse" is tested. A not-NW cell is "shown worse" when its one-sided lower bound exceeds 0, otherwise "inconclusive"; a cell with a missing unit after the continuation is "not determined (did not complete)", never tested. Outcome: Stop (E1, E7, E9 or F2 failing in a counted data set; E11 has its own stop), K1 (all ten NW and M2's gain `a`, read from `results/m2-verdict.tsv`; `a` also holds with the AD cells from v2), otherwise K4 naming the inconclusive and not-determined cells. O-2 (a): "no worse" is judged on the cell average; the changed-data-set figures, calibration by abundance third, power, k and residual degrees of freedom (6.3, 6.4) are reported beside it. M3's own P arm gives the baselines: P's calls and FDP at the binding differ from the pre-KIT-224 figures of the ticket and the decision table shows M3's.

**VAL-400's rule (declared blind, before any counted data).** Effect cells: if the standard deviation of a cell's paired FDP differences is below 10^-14, the differences are treated as equal and the cell is decided by the rule already declared for a standard deviation of exactly zero (KIT-190 `decide.R`, S4 (c)): no test statistic; "no worse" if and only if their mean is below the margin 0.05. Reason: FDPs are ratios of counts with denominators of at most 515, so two different paired differences differ by at least 1.4 x 10^-11, which gives a standard deviation of at least about 4.5 x 10^-13 even at R = 1,000 with a single differing data set; the stored rows carry 15 significant digits, so representation error is about 10^-15 at most. Under such noise Johnson's correction can return "no worse" for differences far above the margin. The change cannot alter a cell whose differences genuinely vary. **`p_signflip`**, the sign-flip p, is a column of `m3-verdicts.tsv` beside Johnson's verdict, **for reporting only**: it decides nothing.

## 6. Limits stated with the result

- **VAL-396 (the owner):** "measured without the sample-correlation filter; data sets where the filter would remove a sample are not covered."
- **The guard's limit (the auditor's words):** the partial-match guard covers the scripts' own code, and package or dependency code that is not under a warning-muffling wrapper. A partial match inside `suppressWarnings`, or inside a handler that muffles warnings (6 uses of `suppressWarnings` in `utils_workflow_omics_declarations.R`; the M2-form helpers `m2Prepare` and `m2Sweep`), is honoured silently. That is exactly what the shipped product does, because users run without `warnPartialMatchDollar`. No warning-catching handler exists in the route files (a grep for `warning = function` in 12 route files found none; dependencies not checked), so the option cannot divert the product's path there. It does not matter for what M3 claims: every field M3 compares is read with `[[` or from trace columns outside any muffling wrapper, and a partial match under a muffling wrapper is the shipped behaviour and acts on both arms alike.
- Ticket section 8: simulated lipid data (metabolomics by the identity check E8), one generator with one unwanted factor, no ITSD acting, two groups of equal size, error rates only at n = 12 and 24 (caps 2 and 3), legacy route only, one filter setting; no statement about real data.

## 7. Continuation, errors, instrument rows

`m3-continuation.md` applies: missing units run once by tag (`score-m3-k262.R missing`), no partial re-run, each data set scored once, timeouts at 2 x plan, the continuation charged to the correction reserve. An error row (a stopped journey of the product) is never a denominator and never re-run.

## 8. Attempts and stop rules

Stage 1 once; each counted job once; one continuation per job; scoring once; one correction pass at most; an unfavourable or inconclusive valid result is never re-run (VAL-54). A failed counted check with a mechanics failure and no result judged is an instrument defect (at most three runs for the item, one re-send used so far on L3 and two on its successors: re-sends 1 and 2 of 3 were L3r and L3r2, the third remains).

## 9. Open before this can be committed (script edits not made; they need the auditor's read)

1. **Error counts by arm** (coordinator, item 11): the scorer counts error rows without the arm (the row's error text begins "P journey" or "C1 Manual journey"). To add: error counts by arm per cell in `m3-units.tsv`, and any C1-only error row (P ok, C1 Manual journey stopped) named beside the verdict and handed to the coordinator, because dropping the data set from both arms makes C1 look no worse where it fails. **Decision needed:** declare or remove `owner_asked` (excluded > 5, inherited from KIT-190's `decide.R`, not declared in the ticket); proposal: remove it and use the C1-only naming.
2. **Instrument-row re-run path (VAL-54):** `loadUnits` stops on any instrument row, and a re-run unit in another source would be "in two sources". To add: a declared way to supersede an instrument row of a unit that is re-run under VAL-54's four conditions (the re-run's complete row replaces the instrument row; both stay in the records), with the continuation tag list including instrument units on request. **Decision needed** on the mechanism before the declaration.
3. The cost of E11's remedy (if E11 fails in the counted run: every re-used eight-cell unit also runs a Manual journey: at most about 1,200 units x 17.8 s x 1.93 / 240 = about 2.5 core-hours of CPU plus set-up, inside the headroom).
4. The final script hashes and the stage-1 R.
