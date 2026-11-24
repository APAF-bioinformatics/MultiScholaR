# Using the validation harness `msrval` (KIT-159)

For the authors of the validation tickets (KIT-160 to KIT-181). The contract is the ticket
`tickets/kit/KIT-159.md`; this page says how to use what is built, and what is not built yet.
Nothing here runs on a laptop: every computation is a `run-tests` scripts request (or, for parity
and Part B, the lanes of sections 7.9 and 7.10).

## State of the harness (keep this section current)

| Part | Built | Verified on a runner |
|---|---|---|
| Scoring, metrics, intervals, classifier, baseline, coverage machinery (group a) | yes | yes (T1a correction; T1w1 to T1w3) |
| Records, licence guard, declarations, provenance, determinism, units and chunks (group b) | yes | yes (T1w1 to T1w3) |
| Driver parts, helper extraction, chain reader, kernel identity, hook, result tables (group c) | yes | yes (T1w1 to T1w3) |
| Driver core: action scripts, `valRunWorkflow`, arm builders, standard tables, replay, custom arms, parity comparisons, Part B runner (group d) | yes | yes (T1w3, 168 of 169; T-TAB in the VAL-60 run) |
| `val-run.R` phases `record`, `pre`, `check` | yes | no |
| Coverage regimes T-COV (T2) | yes | yes (offered regimes below) |
| `val-run.R` phases that run cells (`lifecycle`, `counted`, `score`, `determinism`, `correction`, `amended`), `val-parity.R`, the `distributable.yml` stage `parity` | yes (2b18b1e3) | no; the cell-to-verdict path is repaired after the 2026-11-24 audit (I01 to I20, at 8f37e3e1) and its end-to-end tests (`test-e2e.R`) are not yet run |
| Ticket coverage regimes (`val-tests.R --ticket <dir> --declaration <json>`, pinned in the declaration as `ticket_coverage: {path, sha256}`) | yes (8f37e3e1) | no |
| Replay of declared cells (`replay`) inside a cell phase | no (recorded `not_built` in provenance) | no |
| Known pending change: KIT-202 (lipid and metabolite intensity filter counts the design's samples only) | the five metabolomics corpora are MS-DIAL tables whose integer `Alignment ID` is counted as a sample by today's filter, so every stored metabolomics value after the filter step moves when KIT-202 merges: S0's committed `s0/s0r1/` metabolomics rows (`routes.tsv`, `steps.tsv`, `alternatives.tsv`, ledgers of `metab_lc`, `metab_gc`, `metab_combined`, `metab_lc-export_load`) (the ledgers hold the settings each step saw, so they move only where a default is read from the data; the row counts in `routes.tsv` and `steps.tsv` move); the harness stores no `bindings.json` (the `bindings.json` files under `docs/qa/workflow-execution/` belong to other tickets); lipidomics (LipidSearch) is unaffected unless its tables carry a numeric id. No gate compares a new run against them (S0 gates nothing; G1(a) compares against the product's GUI-entry oracles, not these; parity and the lifecycle run fresh); no test reads them. After KIT-202 merges: re-run S0's three metabolomics routes and the `metab_lc` export-then-load alternative (S0 run 1 took 2,425 s for nine routes; about 1 core-hour) and the T-TAB metabolomics cases in the next whole-suite run | no |

## A ticket's evidence folder

`docs/qa/workflow-execution/kit/<ticket>/` holds the two declaration stages (section 5.3), each
committed alone: start from `harness/templates/` (`declaration.json`, `declaration.md`,
`draft-inputs.tsv`, `inputs.tsv`, `units.tsv`, `prepare.R`, `scoring.R`). Record the inputs first
(`val-run.R record`), then write the declaration from the recorded `inputs.tsv`.

## Arms

| Kind | How to declare it | What the harness does |
|---|---|---|
| `workflow` | a route and a builder: `valArmShipped()`, `valArmRuvLegacy()`, `valArmRuvCalibrated()`, `valArmRuvOff()`, `valArmRuvDesign(design, k_rule)`, `valArmRuvSupplied(role)` | `valRunWorkflow` with the route's accepted script (`valActionScript`, `action-scripts.md`); the builder's inputs per route from `valArmInputs` |
| `kernel` | a function in the ticket's pinned code | `valKernelIdentity` on seeded identity units and every refusal step seen |
| `baseline` | `normalise`, optional `batch` | `valBaselineArm` (section 5.5) |
| `published` | a function and an `anchor` check, `anchor_on_fail` | `valRunCustomArm`; labelled "published processing re-run" |
| `custom` | a function | `valRunCustomArm`; labelled "not the harness baseline" |

What a default user does on each route, and what the driver may set, is section 7.4; the driver
refuses (`val_spec_invalid`) any input the script, the declaration's `forced` settings or the arm's
overrides do not account for, and records every value with its source in `ledger-<route>.tsv`.

## Standard tables of a workflow cell

`valReadArmTables` writes, in the cell's work folder: `da.tsv` (`feature_id, assay, contrast,
log2fc, p, s_own, call_statistic, called_own, mean_expression`), `matrix_pre_ruv.tsv`,
`matrix.tsv` (with `assay`), `features.tsv`, `ruv.tsv` when RUV ran, `itsd.tsv` on lipidomics
and metabolomics, and `matrix_display.tsv` ("display; not what is tested") where RUV factors enter
the model. Which state each table is read from is the rule of section 7.3 (`valMatrixStages`,
`valCellRuv`). The column maps of each route's own DA table are provisional until G3 records them.
The post-cell hook receives these tables' paths and, by named formals only, `unit`, `arm`,
`design_path`, `out_dir`.

## Offered interval regimes (T-COV, T2; provisional ledger row H2)

From `t2/coverage-table.tsv` (T2 requests A 45149480750, B 45149529939, C 45149578690 at 266d9dc8;
2,000 simulations per cell, targets from 20,000). The rule is section 10.2's, applied as written to
each row's recorded status, with no row re-classified: a method failing its band in a regime is not
offered there. A regime is method x statistic class x unit form (and design and transform for C3) x
unit count x level; `valClassify` refuses an interval with no offered regime (`val_spec_invalid`)
unless the consumer registers a passing ticket coverage regime (section 4.4).

Cell legend: `Y` offered (pass); `Yc` offered, reported `conservative` (T2's rows: a run-level
method above min(level + 0.04, 1), the band as run at 266d9dc8; the rank rows of T2R: r1.13's rule
on the miss rate); `Y/Yc` offered, `pass` under T2's band as run and `conservative` under r1.13's
miss-rate rule (O5 of VAL-109: both labels printed; 117 cells in T2, of which the 19 C4 cells were
re-pinned by the I10 diagnostic, leaving 98: C3 accuracy 45, C3 threshold 53; `valHarnessRegimes()` carries the second label as `status_r113`, and
`results.md` prints both beside a verdict on such a row). A `conservative` row with no observed miss
(C4e g = 14 rank, coverage 1.000) licenses only the statistic it simulated (`pauc_diff`): declare
it as `interval.statistic`, or the interval is refused (O4). `Yb` `clopper_pearson`: offered (conservative by construction; its expected
status, not a flag); `Dg` `degenerate`, not offered (more than 1 % zero-width intervals; r1.13); `N` not offered (fails
its band); `Fl` `flagged`, not offered unless the declaration sets `accept_flagged` naming the
owner's ledger row; `-` no interval at that count and level (`exact_median` below its smallest n).
The numbers are the simulated coverage, the lowest and highest over the row's data conditions and
statistics (one number where they agree). `coverage-table.tsv` carries the same ends
(`coverage_min`, `coverage_max`) and, for a row not offered, the condition that failed and on which
side (`failed_condition`, `failed_side`: `low`, `high`) (5.7, r1.12; computed from the recorded
detail columns, no rerun). By VAL-69 (r1.13) `exact_median` at n = 10 and 98.33 % and above keeps `fail` (an
exact binomial tail is a candidate rule only). **The four `rank` rows (C4 g = 4, 6; C4e g = 5, 14)
come from the rank re-run, request R (run 45227824000 at 1827d420; `t2/t2r-*`)**, with the
generator re-pinned (group effects scaled by 0.6, 0.5, 0.35 and 0.25 so that the plain arm's mean
normalised pAUC at FPR 0.01 is in [0.4, 0.8]; second arm the unmoderated ranking), the `degenerate`
rule (no cell degenerate: zero-width share 0, targets 0.27, 0.18, 0.024, 0.049) and run-level
`conservative` read on the miss rate (r1.13); every other row is T2's. Levels in per cent: 95, 97.5, 98.33 (m = 3), 98.75 (m =
4), 99.17 (m = 6), 99.375 (m = 8), 99.58 (m = 12), 99.67 (m = 15). C3 designs: `pairs5`, `pairs8`,
`pairs21` (pairs from 5, 8, 21 batches), `r8_2x2` (2 + 2 of 3 and 5), `r21_2x2` (2 + 2 of 7 and
14), `en3x3` (3 + 3 of 7 and 14), `labs8` (the laboratory grouping); n is the batch count.

| Regime | Method | Class | Form | Design | Transform | n | 95 | 97.5 | 98.33 | 98.75 | 99.17 | 99.375 | 99.58 | 99.67 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C1 | exact_median |  |  |  |  | 6 | Y 0.9665-0.9725 | - | - | - | - | - | - | - |
| C1 | exact_median |  |  |  |  | 7 | Y 0.9860-0.9875 | Y 0.9860-0.9875 | Y 0.9860-0.9875 | - | - | - | - | - |
| C1 | exact_median |  |  |  |  | 8 | Y 0.9910-0.9925 | Y 0.9910-0.9925 | Y 0.9910-0.9925 | Y 0.9910-0.9925 | Y 0.9910-0.9925 | - | - | - |
| C1 | exact_median |  |  |  |  | 9 | Y 0.9590-0.9640 | Y 0.9950-0.9975 | Y 0.9950-0.9975 | Y 0.9950-0.9975 | Y 0.9950-0.9975 | Y 0.9950-0.9975 | Y 0.9950-0.9975 | - |
| C1 | exact_median |  |  |  |  | 10 | Y 0.9700-0.9770 | Y 0.9700-0.9770 | N 0.9950-0.9970 | N 0.9950-0.9970 | N 0.9950-0.9970 | N 0.9950-0.9970 | N 0.9950-0.9970 | N 0.9950-0.9970 |
| C1 | exact_median |  |  |  |  | 13 | Y 0.9760-0.9785 | Y 0.9760-0.9785 | Y 0.9960-0.9965 | Y 0.9960-0.9965 | Y 0.9960-0.9965 | Y 0.9960-0.9965 | Y 0.9960-0.9965 | Y 0.9995-1.0000 |
| C1 | exact_median |  |  |  |  | 21 | Y 0.9730-0.9750 | Y 0.9905-0.9920 | Y 0.9905-0.9920 | Y 0.9905-0.9920 | Y 0.9905-0.9920 | Y 0.9975-0.9985 | Y 0.9975-0.9985 | Y 0.9975-0.9985 |
| C2 | paired_t |  |  |  |  | 4 | Fl 0.8120-0.9510 | Fl 0.8585-0.9765 | Fl 0.8815-0.9870 | Fl 0.9010-0.9905 | Fl 0.9220-0.9915 | Fl 0.9305-0.9930 | Fl 0.9445-0.9960 | Fl 0.9510-0.9960 |
| C2 | paired_t |  |  |  |  | 5 | Fl 0.8350-0.9475 | Fl 0.8790-0.9740 | Fl 0.8945-0.9840 | Fl 0.9105-0.9880 | Fl 0.9310-0.9915 | Fl 0.9410-0.9940 | Fl 0.9540-0.9960 | Fl 0.9590-0.9960 |
| C2 | paired_t |  |  |  |  | 6 | Fl 0.8115-0.9545 | Fl 0.8485-0.9715 | Fl 0.8720-0.9765 | Fl 0.8885-0.9825 | Fl 0.9065-0.9880 | Fl 0.9145-0.9900 | Fl 0.9270-0.9950 | Fl 0.9355-0.9965 |
| C2 | paired_t |  |  |  |  | 8 | Fl 0.8270-0.9515 | Fl 0.8610-0.9755 | Fl 0.8825-0.9840 | Fl 0.8930-0.9865 | Fl 0.9080-0.9925 | Fl 0.9175-0.9945 | Fl 0.9285-0.9960 | Fl 0.9345-0.9970 |
| C2 | paired_t |  |  |  |  | 13 | Fl 0.8380-0.9495 | Fl 0.8775-0.9775 | Fl 0.8905-0.9840 | Fl 0.8985-0.9885 | Fl 0.9070-0.9935 | Fl 0.9180-0.9950 | Fl 0.9300-0.9965 | Fl 0.9340-0.9970 |
| C2 | paired_t |  |  |  |  | 21 | Fl 0.8660-0.9560 | Fl 0.8950-0.9810 | Fl 0.9065-0.9870 | Fl 0.9140-0.9905 | Fl 0.9245-0.9945 | Fl 0.9305-0.9960 | Fl 0.9405-0.9970 | Fl 0.9430-0.9990 |
| C2 | paired_t |  |  |  |  | 100 | Fl 0.9105-0.9435 | Fl 0.9470-0.9675 | Fl 0.9580-0.9815 | Fl 0.9605-0.9855 | Fl 0.9685-0.9900 | Fl 0.9705-0.9915 | Fl 0.9745-0.9935 | Y 0.9790-0.9945 |
| C2 | paired_t |  |  |  |  | 800 | Y 0.9425-0.9550 | Y 0.9735-0.9750 | Y 0.9800-0.9815 | Y 0.9830-0.9855 | Y 0.9895-0.9915 | Y 0.9915-0.9935 | Y 0.9930-0.9945 | Y 0.9945-0.9960 |
| C3 | jackknife | accuracy | unit | pairs5 | identity | 5 | Y 0.9325-0.9755 | Y 0.9600-0.9905 | Y 0.9710-0.9945 | Y 0.9755-0.9955 | Y 0.9825-0.9970 | Y 0.9895-0.9975 | Y 0.9905-0.9980 | Y 0.9930-0.9985 |
| C3 | jackknife | accuracy | unit | pairs5 | log | 5 | Y 0.9595-0.9780 | Y 0.9790-0.9905 | Y 0.9860-0.9945 | Y/Yc 0.9900-0.9970 | Y 0.9925-0.9975 | Y/Yc 0.9940-0.9985 | Y 0.9960-0.9985 | Y 0.9970-0.9985 |
| C3 | jackknife | accuracy | unit | labs8 | identity | 8 | Yc 0.9695-0.9935 | Y/Yc 0.9830-0.9975 | Y/Yc 0.9895-0.9985 | Y/Yc 0.9915-0.9990 | Y/Yc 0.9935-0.9990 | Y/Yc 0.9945-0.9995 | Y/Yc 0.9945-0.9995 | Y/Yc 0.9955-0.9995 |
| C3 | jackknife | accuracy | unit | labs8 | log | 8 | Yc 0.9680-0.9920 | Y/Yc 0.9825-0.9965 | Y/Yc 0.9865-0.9980 | Y/Yc 0.9895-0.9980 | Y/Yc 0.9920-0.9985 | Y/Yc 0.9930-0.9995 | Y/Yc 0.9940-1.0000 | Y/Yc 0.9955-1.0000 |
| C3 | jackknife | accuracy | unit | pairs8 | identity | 8 | N 0.9190-0.9715 | N 0.9465-0.9870 | N 0.9575-0.9915 | N 0.9625-0.9950 | N 0.9690-0.9965 | Y 0.9740-0.9970 | Y 0.9820-0.9980 | Y 0.9860-0.9980 |
| C3 | jackknife | accuracy | unit | pairs8 | log | 8 | Y 0.9340-0.9770 | Y 0.9630-0.9895 | Y 0.9735-0.9925 | Y 0.9795-0.9950 | Y 0.9870-0.9970 | Y/Yc 0.9890-0.9985 | Y/Yc 0.9910-0.9995 | Y/Yc 0.9930-0.9995 |
| C3 | jackknife | accuracy | unit | r8_2x2 | identity | 8 | Yc 0.9530-0.9920 | Y/Yc 0.9780-0.9950 | Y/Yc 0.9825-0.9970 | Y/Yc 0.9855-0.9975 | Y/Yc 0.9885-0.9980 | Y/Yc 0.9930-0.9990 | Y/Yc 0.9950-0.9995 | Y/Yc 0.9950-0.9995 |
| C3 | jackknife | accuracy | unit | r8_2x2 | log | 8 | Yc 0.9575-0.9925 | Y/Yc 0.9730-0.9970 | Y/Yc 0.9820-0.9985 | Y/Yc 0.9855-0.9985 | Y/Yc 0.9900-0.9995 | Y/Yc 0.9915-1.0000 | Y/Yc 0.9955-1.0000 | Y/Yc 0.9960-1.0000 |
| C3 | jackknife | accuracy | unit | en3x3 | identity | 21 | Y 0.9450-0.9815 | Y 0.9650-0.9920 | Y 0.9720-0.9940 | Y 0.9755-0.9955 | Y 0.9820-0.9965 | Y/Yc 0.9860-0.9985 | Y/Yc 0.9895-0.9990 | Y 0.9920-0.9990 |
| C3 | jackknife | accuracy | unit | en3x3 | log | 21 | Y 0.9520-0.9855 | Y/Yc 0.9760-0.9955 | Y/Yc 0.9845-0.9980 | Y/Yc 0.9860-0.9985 | Y/Yc 0.9895-0.9995 | Y/Yc 0.9910-0.9995 | Y/Yc 0.9935-1.0000 | Y/Yc 0.9940-1.0000 |
| C3 | jackknife | accuracy | unit | pairs21 | identity | 21 | N 0.9250-0.9610 | Y 0.9585-0.9835 | Y 0.9650-0.9905 | Y 0.9705-0.9935 | Y 0.9735-0.9960 | Y 0.9775-0.9975 | Y 0.9815-0.9975 | Y 0.9830-0.9975 |
| C3 | jackknife | accuracy | unit | pairs21 | log | 21 | Y 0.9330-0.9620 | Y 0.9575-0.9815 | Y 0.9690-0.9895 | Y 0.9760-0.9930 | Y 0.9820-0.9970 | Y/Yc 0.9840-0.9990 | Y/Yc 0.9880-0.9995 | Y/Yc 0.9885-0.9995 |
| C3 | jackknife | accuracy | unit | r21_2x2 | identity | 21 | Y 0.9335-0.9720 | Y 0.9560-0.9865 | Y 0.9655-0.9930 | Y 0.9720-0.9955 | Y 0.9745-0.9960 | Y 0.9785-0.9965 | Y 0.9835-0.9985 | Y 0.9860-0.9985 |
| C3 | jackknife | accuracy | unit | r21_2x2 | log | 21 | Y 0.9435-0.9705 | Y 0.9655-0.9880 | Y 0.9775-0.9920 | Y 0.9810-0.9930 | Y 0.9850-0.9955 | Y 0.9875-0.9970 | Y 0.9915-0.9970 | Y 0.9920-0.9980 |
| C3 | jackknife | threshold_count | unit | pairs5 | identity | 5 | N 0.9195-0.9730 | N 0.9530-0.9855 | Y 0.9670-0.9890 | Y 0.9735-0.9925 | Y 0.9795-0.9945 | Y 0.9820-0.9960 | Y 0.9860-0.9980 | Y 0.9885-0.9990 |
| C3 | jackknife | threshold_count | unit | pairs5 | log | 5 | Y 0.9495-0.9820 | Y 0.9755-0.9920 | Y 0.9810-0.9950 | Y 0.9840-0.9960 | Y 0.9860-0.9975 | Y/Yc 0.9890-0.9990 | Y/Yc 0.9920-0.9995 | Y/Yc 0.9930-1.0000 |
| C3 | jackknife | threshold_count | unit | labs8 | identity | 8 | Yc 0.9670-0.9940 | Y/Yc 0.9830-0.9975 | Y/Yc 0.9870-0.9985 | Y/Yc 0.9880-0.9985 | Y/Yc 0.9915-0.9990 | Y/Yc 0.9930-0.9995 | Y/Yc 0.9955-1.0000 | Y/Yc 0.9965-1.0000 |
| C3 | jackknife | threshold_count | unit | labs8 | log | 8 | Yc 0.9745-0.9940 | Y/Yc 0.9870-0.9985 | Y/Yc 0.9915-0.9995 | Y/Yc 0.9930-1.0000 | Y/Yc 0.9945-1.0000 | Y/Yc 0.9955-1.0000 | Y/Yc 0.9975-1.0000 | Y/Yc 0.9975-1.0000 |
| C3 | jackknife | threshold_count | unit | pairs8 | identity | 8 | N 0.9060-0.9740 | N 0.9390-0.9885 | N 0.9495-0.9910 | N 0.9590-0.9935 | N 0.9640-0.9965 | N 0.9710-0.9975 | Y 0.9795-0.9985 | Y/Yc 0.9815-0.9995 |
| C3 | jackknife | threshold_count | unit | pairs8 | log | 8 | Y 0.9415-0.9805 | Y 0.9680-0.9895 | Y 0.9815-0.9930 | Y 0.9870-0.9950 | Y 0.9885-0.9965 | Y 0.9895-0.9965 | Y 0.9940-0.9975 | Y 0.9950-0.9975 |
| C3 | jackknife | threshold_count | unit | r8_2x2 | identity | 8 | Yc 0.9375-0.9925 | Y/Yc 0.9640-0.9965 | Y/Yc 0.9740-0.9985 | Y/Yc 0.9795-0.9985 | Y/Yc 0.9855-0.9990 | Y/Yc 0.9875-0.9990 | Y/Yc 0.9890-0.9990 | Y/Yc 0.9910-0.9995 |
| C3 | jackknife | threshold_count | unit | r8_2x2 | log | 8 | Yc 0.9645-0.9920 | Y/Yc 0.9810-0.9965 | Y/Yc 0.9865-0.9975 | Y/Yc 0.9890-0.9985 | Y/Yc 0.9920-0.9990 | Y/Yc 0.9935-0.9995 | Y/Yc 0.9945-1.0000 | Y/Yc 0.9955-1.0000 |
| C3 | jackknife | threshold_count | unit | en3x3 | identity | 21 | Y/Yc 0.9415-0.9880 | Y/Yc 0.9605-0.9945 | Y/Yc 0.9690-0.9960 | Y/Yc 0.9745-0.9985 | Y/Yc 0.9815-0.9995 | Y/Yc 0.9855-0.9995 | Y/Yc 0.9875-1.0000 | Y/Yc 0.9885-1.0000 |
| C3 | jackknife | threshold_count | unit | en3x3 | log | 21 | Y/Yc 0.9585-0.9895 | Y/Yc 0.9780-0.9950 | Y/Yc 0.9820-0.9990 | Y/Yc 0.9855-0.9995 | Y/Yc 0.9890-1.0000 | Y/Yc 0.9915-1.0000 | Y/Yc 0.9950-1.0000 | Y/Yc 0.9950-1.0000 |
| C3 | jackknife | threshold_count | unit | pairs21 | identity | 21 | N 0.9210-0.9555 | N 0.9480-0.9750 | N 0.9570-0.9840 | N 0.9655-0.9870 | N 0.9710-0.9895 | Y 0.9745-0.9940 | Y 0.9790-0.9955 | Y 0.9805-0.9965 |
| C3 | jackknife | threshold_count | unit | pairs21 | log | 21 | Y 0.9320-0.9590 | Y 0.9595-0.9785 | Y 0.9700-0.9865 | Y 0.9775-0.9910 | Y 0.9825-0.9940 | Y 0.9865-0.9955 | Y 0.9910-0.9965 | Y 0.9920-0.9965 |
| C3 | jackknife | threshold_count | unit | r21_2x2 | identity | 21 | Y 0.9360-0.9775 | Y 0.9585-0.9915 | Y 0.9685-0.9955 | Y/Yc 0.9730-0.9970 | Y 0.9780-0.9975 | Y 0.9795-0.9980 | Y 0.9820-0.9985 | Y 0.9825-0.9990 |
| C3 | jackknife | threshold_count | unit | r21_2x2 | log | 21 | Y 0.9525-0.9750 | Y 0.9730-0.9900 | Y 0.9810-0.9925 | Y 0.9845-0.9950 | Y/Yc 0.9855-0.9990 | Y/Yc 0.9880-0.9990 | Y/Yc 0.9930-0.9995 | Y/Yc 0.9930-1.0000 |
| C5 | wilson |  |  |  |  | 21 | Y 0.9430-0.9790 | Y 0.9725-0.9790 | Y 0.9790-0.9840 | Y 0.9790-0.9840 | Y 0.9790-0.9955 | Y 0.9790-0.9955 | Y 0.9790-0.9955 | Y 0.9790-0.9955 |
| C5 | wilson |  |  |  |  | 48 | Y 0.9440-0.9675 | Y 0.9660-0.9675 | Y 0.9870-0.9890 | Y 0.9870-0.9890 | Y 0.9870-0.9890 | Y 0.9890-0.9935 | Y 0.9890-0.9945 | Y 0.9890-0.9945 |
| C6 | newcombe_paired |  |  |  |  | 21 | N 0.9765-0.9925 | Y 0.9910-0.9985 | Y 0.9985-0.9990 | Y 1.0000 | Y 1.0000 | Y 1.0000 | Y 1.0000 | Y 1.0000 |
| C6 | newcombe_paired |  |  |  |  | 48 | Y 0.9615-0.9720 | Y 0.9860-0.9925 | Y 0.9925-0.9965 | Y 0.9945-0.9970 | Y 0.9970-0.9985 | Y 0.9990-1.0000 | Y 0.9990-1.0000 | Y 0.9995-1.0000 |
| C7 | bootstrap_within_run | accuracy |  |  |  | 200 | Y 0.9735 | Y 0.9920 | Y 0.9950 | Y 0.9970 | Y 0.9990 | Y 0.9995 | Y 0.9995 | Y 0.9995 |
| C8 | clopper_pearson |  |  |  |  | 10 | Yb 0.9835-0.9910 | Yb 0.9835-0.9910 | Yb 0.9990-1.0000 | Yb 0.9990-1.0000 | Yb 0.9990-1.0000 | Yb 0.9990-1.0000 | Yb 0.9990-1.0000 | Yb 0.9990-1.0000 |
| C8 | clopper_pearson |  |  |  |  | 21 | Yb 0.9770-0.9865 | Yb 0.9915-0.9980 | Yb 0.9915-0.9980 | Yb 0.9915-0.9980 | Yb 0.9915-0.9980 | Yb 0.9985-1.0000 | Yb 0.9985-1.0000 | Yb 0.9985-1.0000 |
| C8 | clopper_pearson |  |  |  |  | 250 | Yb 0.9440-0.9770 | Yb 0.9745-0.9885 | Yb 0.9805-0.9885 | Yb 0.9860-0.9915 | Yb 0.9900-0.9955 | Yb 0.9930-0.9975 | Yb 0.9955-0.9975 | Yb 0.9965-0.9975 |
| C4e | jackknife | rank | block |  |  | 14 | Yc 1.0000 | Yc 1.0000 | Yc 1.0000 | Yc 1.0000 | Yc 1.0000 | Yc 1.0000 | Yc 1.0000 | Yc 1.0000 |
| C7 | bootstrap_within_run | threshold_count |  |  |  | 200 | Y 0.9485 | Y 0.9735 | Y 0.9830 | Y 0.9875 | Y 0.9915 | Y 0.9925 | Y 0.9940 | Y 0.9950 |
| C4 | jackknife | accuracy | replicate |  |  | 4 | Y 0.9810-0.9825 | Y 0.9895-0.9900 | Y 0.9915-0.9930 | Y 0.9945-0.9950 | Y 0.9965-0.9970 | Y 0.9970-0.9980 | Yc 0.9980-0.9995 | Yc 0.9980-0.9995 |
| C4 | jackknife | accuracy | replicate |  |  | 6 | Yc 0.9950-0.9960 | Yc 0.9980-0.9985 | Yc 0.9990-0.9995 | Yc 0.9995 | Yc 0.9995-1.0000 | Yc 1.0000 | Yc 1.0000 | Yc 1.0000 |
| C4 | jackknife | rank | replicate |  |  | 4 | Y 0.9685 | Y 0.9885 | Y 0.9925 | Y 0.9945 | Y 0.9960 | Y 0.9965 | Y 0.9980 | Y 0.9985 |
| C4 | jackknife | rank | replicate |  |  | 6 | Y 0.9855 | Yc 0.9950 | Yc 0.9975 | Yc 0.9985 | Yc 0.9985 | Yc 0.9985 | Yc 0.9995 | Yc 0.9995 |
| C4 | jackknife | threshold_count | replicate |  |  | 4 | Yc 0.9685-0.9975 | Yc 0.9815-0.9975 | Yc 0.9875-0.9980 | Yc 0.9910-0.9980 | Yc 0.9920-0.9990 | Yc 0.9935-0.9990 | Yc 0.9940-0.9990 | Y 0.9940-0.9990 |
| C4 | jackknife | threshold_count | replicate |  |  | 6 | Dg 0.7265-0.9990 | Dg 0.7265-0.9990 | Dg 0.7265-0.9995 | Dg 0.7265-0.9995 | Dg 0.7265-0.9995 | Dg 0.7265-1.0000 | Dg 0.7265-1.0000 | Dg 0.7265-1.0000 |
| C4e | jackknife | accuracy | block |  |  | 5 | Yc 0.9605-0.9890 | Yc 0.9785-0.9955 | Yc 0.9860-0.9980 | Yc 0.9885-0.9990 | Yc 0.9910-0.9995 | Yc 0.9925-0.9995 | Yc 0.9940-1.0000 | Yc 0.9965-1.0000 |
| C4e | jackknife | accuracy | block |  |  | 14 | Yc 0.9825-1.0000 | Yc 0.9905-1.0000 | Yc 0.9960-1.0000 | Yc 0.9965-1.0000 | Yc 0.9970-1.0000 | Yc 0.9985-1.0000 | Yc 0.9990-1.0000 | Yc 0.9990-1.0000 |
| C4e | jackknife | rank | block |  |  | 5 | Yc 0.9900 | Y 0.9935 | Yc 0.9960 | Yc 0.9970 | Yc 0.9980 | Yc 0.9985 | Yc 0.9995 | Yc 0.9995 |
| C4e | jackknife | threshold_count | block |  |  | 5 | Y 0.9835 | Y 0.9930 | Yc 0.9960 | Yc 0.9970 | Yc 0.9980 | Yc 0.9985 | Yc 0.9990 | Y 0.9990 |
| C4e | jackknife | threshold_count | block |  |  | 14 | Yc 0.9995 | Yc 1.0000 | Yc 1.0000 | Yc 1.0000 | Yc 1.0000 | Yc 1.0000 | Yc 1.0000 | Yc 1.0000 |

What the table means for the consumers (read off it; nothing decided here):

- `paired_t` (C2) is `flagged` at every count from 4 to 100 and every level, except n = 100 at
  99.67 % (offered), and offered at n = 800 at every level: at the small counts it is usable only
  with `accept_flagged`.
- `jackknife` form (iii), delete-one-block (C4e): `accuracy` and `threshold_count` are not offered at
  g = 5 (KIT-162) or g = 14 (KIT-166) at any level; `rank` is offered at both (re-run R): at g = 5
  `conservative` at 95 % (coverage 0.990) and from 98.33 %, `pass` at 97.5 %; at g = 14
  `conservative` at every level (coverage 1.000: no miss in 2,000).
- `jackknife` form (ii), delete-one-replicate-per-group (C4): `accuracy` and `threshold_count` are not
  offered at 95 % or 97.5 % at g = 4 (KIT-160, KIT-161), nor below 99.17 % (accuracy) or 98.75 %
  (threshold) at g = 6 (KIT-170); `rank` is offered at both counts.
- `jackknife` form (i) (C3): offered on every design at every level except the identity transform on
  `pairs8` (both classes, below 99.375 % or 99.58 %), on `pairs21` (accuracy at 95 %; threshold below
  99.375 %) and on `pairs5` (threshold at 95 % and 97.5 %); the log transform is offered on every
  design at every level.
- `newcombe_paired` (C6) is not offered at 21 pairs at 95 % (within-run, over its upper band at
  (0.9, 0.9)); every other C5, C6, C7 and C8 row is offered.
- A consumer whose count or level has no offered regime is told before its declaration (4.4); its
  verdict there is descriptive unless a ticket coverage regime passes.

## G1(a) classes and the DIA-with-limpa link

`limpa_da_last_digits` is not a class and cannot be admitted (D19: VAL-109 (5); the KIT-199 plan
audit's ruling C1 (a)). KIT-199's re-freeze of the DIA-with-limpa reference is read as
`declared_change (KIT-199)`: `valClassifyLiveDiffs(diffs, lane_steps, declared = <lane, step, key,
ticket>)` classes a declared key first and the step's other differences among themselves, so one
step can carry a declared change and `integration_artifacts`; `valReadChain` link (i) accepts the
link's own declared change only where its record lists the step and key.

**Counting unit (KIT-199 re-check 1, Q01).** Class counts are per (step, class), never per
difference row: `g1a-reference.tsv` has one row per lane and differing step with its class (G1(a)),
so a step with seven changed integration keys is one `integration_artifacts` unit. A control record
(`gui-entry-control.tsv` of KIT-172, 177, 178, 180, 198, 199) may list one row per key (`lane`,
`step`, `field`, `key`, `class`); the harness reduces them to (step, class) before comparing. A
re-freezing ticket's declared keys are listed one row per key (`lane`, `step`, `key`, `expected`;
KIT-199: `kit-199/declared-keys.tsv`); only listed keys are accepted as its `declared_change`.

Owner's decisions VAL-133 (to build, tests first; not built yet): G13 compares the bundle's
numbers under the SWD-083 tolerance rule (not 10 significant digits); G1(a) at the release commit
keeps identical masked captures as its first check and, when that fails, compares values (keys,
order, missing values and calls exactly; numbers under SWD-083) and reports `certified` with the
difference printed or `not_certified`; the Linux bundle's GUI-entry lane's six DIA-with-limpa DA
steps are excluded by name from installer acceptance (no harness stage reads that lane).

Not built yet (I15 remaining; KIT-199 implementation audit): the loader of `gui-entry-control.tsv`
and `declared-keys.tsv`, link (i) checking declared digest values, reference-row adoption by key
sets, and the G1(a) part (1) and (2) runners, the loader of committed
`gui-entry-control.tsv` records (to be built against KIT-172's, 177's and 180's records as
fixtures), and the release-commit DA-value comparison of the DIA-with-limpa lane (VAL-120 (2);
KIT-199 section 12.3).

## Transcriptomics RUV grid summary (KIT-192)

After KIT-192 merges, the transcriptomics RUV grid summary file is written only when a legacy grid
ran (VAL-102 D1 (a)); it does not appear on declared-design arms (calibrated, off,
`valArmRuvDesign`). The harness reads the trx DA tables from `de_results_*.csv` and expects no
summary file on any arm.

**I10 re-pin (2026-11-24, VAL-109 (2)).** The C4 and C4e `accuracy` and `threshold_count` rows
(C4 g = 4, 6; C4e g = 5, 14) come from the I10 diagnostic (run in `t2/di10d-*`, request D at
169aee9f): T2's generator drew the feature baselines anew in every simulation with one-directional
changes, so median normalisation shifted every estimate by a feature-level amount the unit
jackknife cannot see; the fixed per-cell feature draw with changes split +1/-1 (`fixed_symmetric`,
ticket 10.2 r1.16) removes it (C4e targets: bias 0.0008 and -0.0091, null rejection 0.051 and
0.051, against -0.21/-0.19 and 0.37/0.63 in T2). Every one of the 11 rows that failed at 95 % now
covers at or above 0.9605; no row fails. C4 g = 6 `threshold_count` is `degenerate` (its
sensitivity statistic gives zero-width intervals in more than 1 % of simulations): not offered.
`half_width_ratio` is printed for these rows (about 1.4 to 5.3: the jackknife intervals are
conservative) and for the two consumer-used rank rows in `t2/di10h-*` (200 simulations: C4 g = 6
rank 2.0 to 2.7, C4e g = 5 rank 2.1 to 3.2; not in the pinned table, which keeps T2R's rank rows).

## Within-run intervals beside a run-level verdict (4.4; O8 of VAL-109)

- Jackknife forms (ii) and (iii) (`full` plus `deletion` units): the run-level interval is widened,
  bound by bound, to the `full` unit's within-run interval of the same view; mandatory below 6
  units, the declarable default at 6 or more (`widen_to_within_run = false` with its reason).
- Independent (or `draw`) units: nothing is widened at any unit count. Each unit's within-run
  interval is reported beside the verdict (`within_run_reported=<unit>:[low,high]` in `flags`), with
  `within_run_widening=report_only(O8:below_6)` below 6 units (the owner's ruling) or
  `report_only(O8:6_or_more_open)` at 6 or more: the audit register leaves open which unit's
  interval would widen there, so the harness reports only; `widen_to_within_run` has no effect on
  independent units, and `results.md` says so beside the verdict. A widening rule for independent
  units is the owner's to make.

## Not yet usable

- The cell phases are built but have not run end to end on a runner (the repair's runs are listed
  in `outcome.md`, "Resume here"). A request's `out` for a cell phase is `val-out` and its
  `out_pattern` must admit the run folder's files: `\\.(tsv|txt|md|json|sha256|gz)$` (the files of
  5.7 and `allowed-paths.txt`; ledgers are `ledger-<route>-<arm>-<unit>.tsv`).
- The ticket's `prepareInputs()` returns `features` (`assay`, `feature_id`: every feature of the
  declared input role) beside `units`; `featureKey(feature_id, arm)` returns `key`, `status`
  (`ok`, `unmapped`, `ambiguous`) and `reason`. The post-cell hook gets the proteomics DA results
  object as `results` and `NULL` on lipidomics, metabolomics and transcriptomics, and writes only
  under `out_dir` (`runs/<label>/ticket/`).
- After a refused unit a verdict needs a regime for n - k units (O7 (b), VAL-109), else
  `not_evaluable`.
- The I10 diagnostic of the C4 and C4e 500-feature cells is not yet run; until it is, the C4 and
  C4e `accuracy` and `threshold_count` rows marked fail are not evidence about the method.
- Parity (G13) needs `val-parity.R` and the `distributable.yml` stage on top of KIT-176's merged
  workflow.
