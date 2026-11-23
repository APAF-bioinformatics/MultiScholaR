# Using the validation harness `msrval` (KIT-159)

For the authors of the validation tickets (KIT-160 to KIT-181). The contract is the ticket
`tickets/kit/KIT-159.md`; this page says how to use what is built, and what is not built yet.
Nothing here runs on a laptop: every computation is a `run-tests` scripts request (or, for parity
and Part B, the lanes of sections 7.9 and 7.10).

## State of the harness (keep this section current)

| Part | Built | Verified on a runner |
|---|---|---|
| Scoring, metrics, intervals, classifier, baseline, coverage machinery (group a) | yes | yes (T1a correction) |
| Records, licence guard, declarations, provenance, determinism, units and chunks (group b) | yes | all but T-LOG (fix unrun) |
| Driver parts, helper extraction, chain reader, kernel identity, hook, result tables (group c) | yes | no (unrun) |
| Driver core: action scripts, `valRunWorkflow`, arm builders, standard tables, replay, custom arms, parity comparisons, Part B runner (group d) | yes | no (unrun; T0d showed its tests failing first) |
| `val-run.R` phases `record`, `pre`, `check` | yes | no |
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

## Not yet usable

- Cells cannot yet be run through `val-run.R` (its cell phases are not written); S0 and the
  lifecycle need them.
- Parity (G13) needs `val-parity.R` and the `distributable.yml` stage on top of KIT-176's merged
  workflow.
- Coverage regimes (T-COV, `coverage-table.tsv`) are not run yet (T2), so no interval method is
  offered for scoring until it is.
