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
| `val-run.R` phases that run cells (`lifecycle`, `counted`, `score`, `determinism`, `correction`, `amended`), `val-parity.R`, the `distributable.yml` stage `parity` | no | no |

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

Cell legend: `Y` offered (pass); `Yc` offered, reported `conservative` (a run-level method above
min(level + 0.04, 1), or `clopper_pearson`, conservative by construction); `N` not offered (fails
its band); `Fl` `flagged`, not offered unless the declaration sets `accept_flagged` naming the
owner's ledger row; `-` no interval at that count and level (`exact_median` below its smallest n).
The numbers are the simulated coverage, the lowest and highest over the row's data conditions and
statistics (one number where they agree). Levels in per cent: 95, 97.5, 98.33 (m = 3), 98.75 (m =
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
| C3 | jackknife | accuracy | unit | pairs5 | log | 5 | Y 0.9595-0.9780 | Y 0.9790-0.9905 | Y 0.9860-0.9945 | Y 0.9900-0.9970 | Y 0.9925-0.9975 | Y 0.9940-0.9985 | Y 0.9960-0.9985 | Y 0.9970-0.9985 |
| C3 | jackknife | accuracy | unit | labs8 | identity | 8 | Yc 0.9695-0.9935 | Y 0.9830-0.9975 | Y 0.9895-0.9985 | Y 0.9915-0.9990 | Y 0.9935-0.9990 | Y 0.9945-0.9995 | Y 0.9945-0.9995 | Y 0.9955-0.9995 |
| C3 | jackknife | accuracy | unit | labs8 | log | 8 | Yc 0.9680-0.9920 | Y 0.9825-0.9965 | Y 0.9865-0.9980 | Y 0.9895-0.9980 | Y 0.9920-0.9985 | Y 0.9930-0.9995 | Y 0.9940-1.0000 | Y 0.9955-1.0000 |
| C3 | jackknife | accuracy | unit | pairs8 | identity | 8 | N 0.9190-0.9715 | N 0.9465-0.9870 | N 0.9575-0.9915 | N 0.9625-0.9950 | N 0.9690-0.9965 | Y 0.9740-0.9970 | Y 0.9820-0.9980 | Y 0.9860-0.9980 |
| C3 | jackknife | accuracy | unit | pairs8 | log | 8 | Y 0.9340-0.9770 | Y 0.9630-0.9895 | Y 0.9735-0.9925 | Y 0.9795-0.9950 | Y 0.9870-0.9970 | Y 0.9890-0.9985 | Y 0.9910-0.9995 | Y 0.9930-0.9995 |
| C3 | jackknife | accuracy | unit | r8_2x2 | identity | 8 | Yc 0.9530-0.9920 | Y 0.9780-0.9950 | Y 0.9825-0.9970 | Y 0.9855-0.9975 | Y 0.9885-0.9980 | Y 0.9930-0.9990 | Y 0.9950-0.9995 | Y 0.9950-0.9995 |
| C3 | jackknife | accuracy | unit | r8_2x2 | log | 8 | Yc 0.9575-0.9925 | Y 0.9730-0.9970 | Y 0.9820-0.9985 | Y 0.9855-0.9985 | Y 0.9900-0.9995 | Y 0.9915-1.0000 | Y 0.9955-1.0000 | Y 0.9960-1.0000 |
| C3 | jackknife | accuracy | unit | en3x3 | identity | 21 | Y 0.9450-0.9815 | Y 0.9650-0.9920 | Y 0.9720-0.9940 | Y 0.9755-0.9955 | Y 0.9820-0.9965 | Y 0.9860-0.9985 | Y 0.9895-0.9990 | Y 0.9920-0.9990 |
| C3 | jackknife | accuracy | unit | en3x3 | log | 21 | Y 0.9520-0.9855 | Y 0.9760-0.9955 | Y 0.9845-0.9980 | Y 0.9860-0.9985 | Y 0.9895-0.9995 | Y 0.9910-0.9995 | Y 0.9935-1.0000 | Y 0.9940-1.0000 |
| C3 | jackknife | accuracy | unit | pairs21 | identity | 21 | N 0.9250-0.9610 | Y 0.9585-0.9835 | Y 0.9650-0.9905 | Y 0.9705-0.9935 | Y 0.9735-0.9960 | Y 0.9775-0.9975 | Y 0.9815-0.9975 | Y 0.9830-0.9975 |
| C3 | jackknife | accuracy | unit | pairs21 | log | 21 | Y 0.9330-0.9620 | Y 0.9575-0.9815 | Y 0.9690-0.9895 | Y 0.9760-0.9930 | Y 0.9820-0.9970 | Y 0.9840-0.9990 | Y 0.9880-0.9995 | Y 0.9885-0.9995 |
| C3 | jackknife | accuracy | unit | r21_2x2 | identity | 21 | Y 0.9335-0.9720 | Y 0.9560-0.9865 | Y 0.9655-0.9930 | Y 0.9720-0.9955 | Y 0.9745-0.9960 | Y 0.9785-0.9965 | Y 0.9835-0.9985 | Y 0.9860-0.9985 |
| C3 | jackknife | accuracy | unit | r21_2x2 | log | 21 | Y 0.9435-0.9705 | Y 0.9655-0.9880 | Y 0.9775-0.9920 | Y 0.9810-0.9930 | Y 0.9850-0.9955 | Y 0.9875-0.9970 | Y 0.9915-0.9970 | Y 0.9920-0.9980 |
| C3 | jackknife | threshold_count | unit | pairs5 | identity | 5 | N 0.9195-0.9730 | N 0.9530-0.9855 | Y 0.9670-0.9890 | Y 0.9735-0.9925 | Y 0.9795-0.9945 | Y 0.9820-0.9960 | Y 0.9860-0.9980 | Y 0.9885-0.9990 |
| C3 | jackknife | threshold_count | unit | pairs5 | log | 5 | Y 0.9495-0.9820 | Y 0.9755-0.9920 | Y 0.9810-0.9950 | Y 0.9840-0.9960 | Y 0.9860-0.9975 | Y 0.9890-0.9990 | Y 0.9920-0.9995 | Y 0.9930-1.0000 |
| C3 | jackknife | threshold_count | unit | labs8 | identity | 8 | Yc 0.9670-0.9940 | Y 0.9830-0.9975 | Y 0.9870-0.9985 | Y 0.9880-0.9985 | Y 0.9915-0.9990 | Y 0.9930-0.9995 | Y 0.9955-1.0000 | Y 0.9965-1.0000 |
| C3 | jackknife | threshold_count | unit | labs8 | log | 8 | Yc 0.9745-0.9940 | Y 0.9870-0.9985 | Y 0.9915-0.9995 | Y 0.9930-1.0000 | Y 0.9945-1.0000 | Y 0.9955-1.0000 | Y 0.9975-1.0000 | Y 0.9975-1.0000 |
| C3 | jackknife | threshold_count | unit | pairs8 | identity | 8 | N 0.9060-0.9740 | N 0.9390-0.9885 | N 0.9495-0.9910 | N 0.9590-0.9935 | N 0.9640-0.9965 | N 0.9710-0.9975 | Y 0.9795-0.9985 | Y 0.9815-0.9995 |
| C3 | jackknife | threshold_count | unit | pairs8 | log | 8 | Y 0.9415-0.9805 | Y 0.9680-0.9895 | Y 0.9815-0.9930 | Y 0.9870-0.9950 | Y 0.9885-0.9965 | Y 0.9895-0.9965 | Y 0.9940-0.9975 | Y 0.9950-0.9975 |
| C3 | jackknife | threshold_count | unit | r8_2x2 | identity | 8 | Yc 0.9375-0.9925 | Y 0.9640-0.9965 | Y 0.9740-0.9985 | Y 0.9795-0.9985 | Y 0.9855-0.9990 | Y 0.9875-0.9990 | Y 0.9890-0.9990 | Y 0.9910-0.9995 |
| C3 | jackknife | threshold_count | unit | r8_2x2 | log | 8 | Yc 0.9645-0.9920 | Y 0.9810-0.9965 | Y 0.9865-0.9975 | Y 0.9890-0.9985 | Y 0.9920-0.9990 | Y 0.9935-0.9995 | Y 0.9945-1.0000 | Y 0.9955-1.0000 |
| C3 | jackknife | threshold_count | unit | en3x3 | identity | 21 | Y 0.9415-0.9880 | Y 0.9605-0.9945 | Y 0.9690-0.9960 | Y 0.9745-0.9985 | Y 0.9815-0.9995 | Y 0.9855-0.9995 | Y 0.9875-1.0000 | Y 0.9885-1.0000 |
| C3 | jackknife | threshold_count | unit | en3x3 | log | 21 | Y 0.9585-0.9895 | Y 0.9780-0.9950 | Y 0.9820-0.9990 | Y 0.9855-0.9995 | Y 0.9890-1.0000 | Y 0.9915-1.0000 | Y 0.9950-1.0000 | Y 0.9950-1.0000 |
| C3 | jackknife | threshold_count | unit | pairs21 | identity | 21 | N 0.9210-0.9555 | N 0.9480-0.9750 | N 0.9570-0.9840 | N 0.9655-0.9870 | N 0.9710-0.9895 | Y 0.9745-0.9940 | Y 0.9790-0.9955 | Y 0.9805-0.9965 |
| C3 | jackknife | threshold_count | unit | pairs21 | log | 21 | Y 0.9320-0.9590 | Y 0.9595-0.9785 | Y 0.9700-0.9865 | Y 0.9775-0.9910 | Y 0.9825-0.9940 | Y 0.9865-0.9955 | Y 0.9910-0.9965 | Y 0.9920-0.9965 |
| C3 | jackknife | threshold_count | unit | r21_2x2 | identity | 21 | Y 0.9360-0.9775 | Y 0.9585-0.9915 | Y 0.9685-0.9955 | Y 0.9730-0.9970 | Y 0.9780-0.9975 | Y 0.9795-0.9980 | Y 0.9820-0.9985 | Y 0.9825-0.9990 |
| C3 | jackknife | threshold_count | unit | r21_2x2 | log | 21 | Y 0.9525-0.9750 | Y 0.9730-0.9900 | Y 0.9810-0.9925 | Y 0.9845-0.9950 | Y 0.9855-0.9990 | Y 0.9880-0.9990 | Y 0.9930-0.9995 | Y 0.9930-1.0000 |
| C5 | wilson |  |  |  |  | 21 | Y 0.9430-0.9790 | Y 0.9725-0.9790 | Y 0.9790-0.9840 | Y 0.9790-0.9840 | Y 0.9790-0.9955 | Y 0.9790-0.9955 | Y 0.9790-0.9955 | Y 0.9790-0.9955 |
| C5 | wilson |  |  |  |  | 48 | Y 0.9440-0.9675 | Y 0.9660-0.9675 | Y 0.9870-0.9890 | Y 0.9870-0.9890 | Y 0.9870-0.9890 | Y 0.9890-0.9935 | Y 0.9890-0.9945 | Y 0.9890-0.9945 |
| C6 | newcombe_paired |  |  |  |  | 21 | N 0.9765-0.9925 | Y 0.9910-0.9985 | Y 0.9985-0.9990 | Y 1.0000 | Y 1.0000 | Y 1.0000 | Y 1.0000 | Y 1.0000 |
| C6 | newcombe_paired |  |  |  |  | 48 | Y 0.9615-0.9720 | Y 0.9860-0.9925 | Y 0.9925-0.9965 | Y 0.9945-0.9970 | Y 0.9970-0.9985 | Y 0.9990-1.0000 | Y 0.9990-1.0000 | Y 0.9995-1.0000 |
| C7 | bootstrap_within_run | accuracy |  |  |  | 200 | Y 0.9735 | Y 0.9920 | Y 0.9950 | Y 0.9970 | Y 0.9990 | Y 0.9995 | Y 0.9995 | Y 0.9995 |
| C8 | clopper_pearson |  |  |  |  | 10 | Yc 0.9835-0.9910 | Yc 0.9835-0.9910 | Yc 0.9990-1.0000 | Yc 0.9990-1.0000 | Yc 0.9990-1.0000 | Yc 0.9990-1.0000 | Yc 0.9990-1.0000 | Yc 0.9990-1.0000 |
| C8 | clopper_pearson |  |  |  |  | 21 | Yc 0.9770-0.9865 | Yc 0.9915-0.9980 | Yc 0.9915-0.9980 | Yc 0.9915-0.9980 | Yc 0.9915-0.9980 | Yc 0.9985-1.0000 | Yc 0.9985-1.0000 | Yc 0.9985-1.0000 |
| C8 | clopper_pearson |  |  |  |  | 250 | Yc 0.9440-0.9770 | Yc 0.9745-0.9885 | Yc 0.9805-0.9885 | Yc 0.9860-0.9915 | Yc 0.9900-0.9955 | Yc 0.9930-0.9975 | Yc 0.9955-0.9975 | Yc 0.9965-0.9975 |
| C4e | jackknife | rank | block |  |  | 14 | Yc 1.0000 | Y 1.0000 | Y 1.0000 | Y 1.0000 | Y 1.0000 | Y 1.0000 | Y 1.0000 | Y 1.0000 |
| C7 | bootstrap_within_run | threshold_count |  |  |  | 200 | Y 0.9485 | Y 0.9735 | Y 0.9830 | Y 0.9875 | Y 0.9915 | Y 0.9925 | Y 0.9940 | Y 0.9950 |
| C4 | jackknife | accuracy | replicate |  |  | 4 | N 0.9180-0.9775 | N 0.9535-0.9875 | Y 0.9695-0.9940 | Y 0.9765-0.9955 | Y 0.9855-0.9980 | Y 0.9880-0.9985 | Y 0.9920-0.9990 | Y 0.9935-1.0000 |
| C4 | jackknife | accuracy | replicate |  |  | 6 | N 0.8980-0.9905 | N 0.9395-0.9960 | N 0.9560-0.9975 | N 0.9655-0.9985 | Y 0.9765-0.9990 | Y 0.9830-0.9990 | Y 0.9865-0.9995 | Y 0.9890-1.0000 |
| C4 | jackknife | rank | replicate |  |  | 4 | Y 0.9840 | Y 0.9920 | Y 0.9955 | Y 0.9970 | Y 0.9985 | Y 0.9985 | Y 0.9990 | Y 0.9995 |
| C4 | jackknife | rank | replicate |  |  | 6 | Yc 0.9975 | Y 0.9980 | Y 0.9990 | Y 0.9990 | Y 0.9995 | Y 0.9995 | Y 0.9995 | Y 0.9995 |
| C4 | jackknife | threshold_count | replicate |  |  | 4 | N 0.9080-0.9855 | N 0.9490-0.9950 | Y 0.9640-0.9975 | Y 0.9735-0.9980 | Y 0.9815-0.9980 | Y 0.9865-0.9985 | Y 0.9885-0.9995 | Y 0.9925-0.9995 |
| C4 | jackknife | threshold_count | replicate |  |  | 6 | N 0.8995-0.9990 | N 0.9405-1.0000 | N 0.9575-1.0000 | Y 0.9675-1.0000 | Y 0.9780-1.0000 | Y 0.9845-1.0000 | Y 0.9870-1.0000 | Y 0.9900-1.0000 |
| C4e | jackknife | accuracy | block |  |  | 5 | N 0.7450-0.7565 | N 0.8230-0.8370 | N 0.8685-0.8815 | N 0.8925-0.9015 | N 0.9170-0.9310 | N 0.9350-0.9445 | N 0.9515-0.9610 | N 0.9575-0.9630 |
| C4e | jackknife | accuracy | block |  |  | 14 | N 0.7065-0.7070 | N 0.7810-0.7820 | N 0.8145 | N 0.8330 | N 0.8555 | N 0.8720 | N 0.8880-0.8885 | N 0.8960-0.8965 |
| C4e | jackknife | rank | block |  |  | 5 | Yc 0.9940 | Y 0.9985 | Y 0.9990 | Y 0.9990 | Y 0.9990 | Y 0.9995 | Y 0.9995 | Y 0.9995 |
| C4e | jackknife | threshold_count | block |  |  | 5 | N 0.7065 | N 0.8005 | N 0.8385 | N 0.8660 | N 0.8950 | N 0.9130 | N 0.9325 | N 0.9445 |
| C4e | jackknife | threshold_count | block |  |  | 14 | N 0.7260 | N 0.7895 | N 0.8235 | N 0.8460 | N 0.8750 | N 0.8865 | N 0.9035 | N 0.9160 |

What the table means for the consumers (read off it; nothing decided here):

- `paired_t` (C2) is `flagged` at every count from 4 to 100 and every level, except n = 100 at
  99.67 % (offered), and offered at n = 800 at every level: at the small counts it is usable only
  with `accept_flagged`.
- `jackknife` form (iii), delete-one-block (C4e): `accuracy` and `threshold_count` are not offered at
  g = 5 (KIT-162) or g = 14 (KIT-166) at any level; `rank` is offered at both.
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

## Not yet usable

- Cells cannot yet be run through `val-run.R` (its cell phases are not written); S0 and the
  lifecycle need them.
- Parity (G13) needs `val-parity.R` and the `distributable.yml` stage on top of KIT-176's merged
  workflow.
