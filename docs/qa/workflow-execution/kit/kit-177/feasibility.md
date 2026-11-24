# KIT-177 feasibility (ticket section 8, F1 to F4)

Ticket revision 1.6. Read at `v0.9-staging` fc48b6ec (KIT-180 merged). F1 and F2 are readings; F3 and F4 ran in R0
(`runs.tsv`).

## F1 consumers of post-log values (reading, no run)

(a) Searched `R/func_prot_*`, `R/func_general_plotting_*`, `R/func_general_nipals.R`, `R/func_omics_ruv_*`,
`R/utils_workflow_prot_declarations.R` for comparisons with 0, `abs`, `sqrt`, a log, a constant fill of missing
values, `2^x - 1`, and sums or products of log values. Classified:

| Site | What it does on log values | Class |
|---|---|---|
| `R/func_general_plotting_pca_rle_helpers.R:86-103, 356-375` | PCA input filter divides the SD by the mean log value; a mean within 1e-6 of 0 is "no variation" | changes a displayed figure (known, threat a; listed for the owner, not changed: "part of the recorded science") |
| `R/func_prot_s4_accession_methods.R:192-197` (and the `islogged` aggregation `log2(sum(2^x))` at `:267`) | duplicate accessions aggregated on log values (sum or mean) | changes a saved value only where duplicate accessions exist (known; F177 has none) |
| `R/func_prot_s4_norm_methods.R:80` (`scale` = limma `normalizeMedianAbsValues`) | divides by the median absolute log value | changes a saved value when log values are below 0 (known, D6: warning added) |
| `R/func_prot_da_heatmap.R:181` | fills missing with 0 after row scaling | display only (known) |
| `R/func_general_nipals.R:36, 108` | missing set to 0 after centring | display only (centred values; unchanged by a shift) |
| `R/func_prot_qc_diagnostic_helpers.R:348` | missing correlations set to 0 | display only |
| `R/func_prot_s4_missingness.R:160, 600` | `2^x - 1` inverses | not on the label-free or TMT routes (2.3) |
| `R/func_prot_norm_optimization_helpers.R:327` | `values_fill = 0` | design matrix indicator, not a quantity |
| volcano and summary colour rules (`abs(log2FC) >= 1`) | fold changes, not log intensities | unchanged by a shift |

No site beyond the known four changes a saved or displayed number because the log values may now be below 0 or
because zeros become missing. Stop rule not reached.

(b) to (e): answered in the ticket (2.1, 2.2, 2.3).

(f) The record lists every `fixed` provenance entry as a parameter with provenance `fixed`
(`R/utils_workflow_analysis_record.R:421-432`); the validator checks each parameter's name, source and provenance
vocabulary only (`:607-636`). The two added fields `nonpositive_set_to_na` and `nonfinite_set_to_na` are therefore
accepted and rendered with the label "fixed" (ticket 6, "Records"). Shown on the runner by T3b and T4.

(g) Both sites can read the restored state's provenance:
- the design tab's server start, `mod_prot_design_server` (`R/mod_prot_design.R:141`), holds `workflow_data$state_manager`;
- the load-session handler, `da_server_load_session_handler` (`R/mod_prot_da_handlers_load_session.R`), holds the
  restored owner as `manager` after the session's states are restored.
On an on-disk owner (`DIAStateAdapter`) the design step's provenance is `manager$.context$snapshot$step_provenance$protein_s4_initial`
(`operation = "design_builder"`, `fixed$transform`). The notifier `notifyProtDesignEarlierTransform` reads exactly this.
The texts table: `omicGuiProtLog2Texts()` in `R/utils_omic_gui_config.R` (as KIT-172's `omicGuiLog2Texts()`).

## F2 Proteome Discoverer's "Reporter Abundance Based On" (public Thermo documentation)

Read 2026-11-24: Proteome Discoverer 3.1 User Guide, "Using signal-to-noise values as quantification channel values",
https://docs.thermofisher.com/r/Proteome-Discoverer-3.1-User-Guide/en-US1467223947v1. Found: the Reporter Ions
Quantifier node's parameter "Reporter Abundance Based On" chooses between the reporter peaks' signal-to-noise values
and their intensities; "with an Orbitrap analyzer, the intensity of a signal is not proportional to the number of
ions", whereas "the number of ions is proportional to the S/N value of a peak". The guide refers to the node's Help
topic for the parameter's options; its default ("Automatic") and what it selects on Orbitrap data were **not found** in
the pages read. So NEWS says "some exports" (signal-to-noise exports have small values), not "the default export".

## F3 replay of a label-free and a TMT record (R0, run 45222712083, commit 025c9542, R/ = base fc48b6ec)

Script `replay-check.R`; result `r0/replay.tsv`, `r0/replay-slots.tsv`. Known before the run: KIT-180's run R0d
(45149789668) replayed GUI-made records in order with no replayed state `identical()` (AUDIT-KIT-180-impl I06).

Observed, on `prot_lfq_maxquant` (LFQ) and `prot_tmt_pd` (TMT) alike, records made in process by the package's own
import and design step (base: `log2(x + 1)`) and the route's declared operations:
- record valid; steps import, design_builder, the four protein QC steps, normalization-ruv, pearson-selection; the
  replay ran every later step and the DA (`DA-contrast-output`); replayed history equal to the record's;
- 5 of 6 replayed states `identical()` to the first run's; `correlation_filtered` not identical, in one slot:
  `@args$daAnalysisUI` (not path-bearing);
- the DA long table not `identical()`, and its column names differ between the two reads;
- the kit fingerprint (`kitFingerprintStore`, `kitReproduceMismatches`): no mismatch on either lane.

**Outcome: (3) "replays in order, not identical"** on both lanes. Under the ticket's stop rule (revision 1.6) this
stops the ticket before R1 and goes to the owner. Nothing was re-run.

Inferred, not tested: `daAnalysisUI` (with its timestamp) is committed onto the head state by the DA's request half
(`commitDIADaArguments`, `R/utils_workflow_dia_da.R`); the script read the first run's `correlation_filtered` after the
DA had run, but captured the replayed one in `on_step` at the correlation step, before the replayed DA. If so, the
state difference comes from the moment of capture, not from the replay. The DA table's differing column names are
not explained by this reading. Whether this is an instrument defect or a result is the owner's decision (VAL-54: in
doubt it is a result).

## F4 the control's mechanism from a script (R0, same run; S2 `gui-entry-control.R` on `prot_lfq`)

Result `r0/f4-verdict.tsv`, `r0/f4-tables.tsv`, `r0/f4-steps.tsv`. Head = 025c9542, base = fc48b6ec (R/ and tests/
equal), so every comparison is expected empty. Observed: both child processes under C collation; the base side's
`test_path()` inside its worktree; the head side's oracle directory redirected into `out`; base and head records and
the frozen file written into `out`; the checkout's oracles unchanged; step lists equal; frozen against base 0 rows;
base against head 0 rows; head record against its frozen file 0 rows; 8 state tables exported on each side, all cells
identical. The one table classed `unexplained (design rule)` is `protein_s4_initial`, because with no change the
design-rule check (base minus head = log2(1 + 1/v)) cannot hold; on the fix commit it is the check intended.
**Answer: yes**, the state tables can be exported without editing `helper-gui-entry.R` and S2's mechanism works from a
script (R4f not needed).

Scope completed by the implementation audit (AUDIT-KIT-177-impl I11): `R/mod_prot_*`, `R/utils_workflow_dia_*`, the
other `R/func_general_*`, `R/utils_workflow_saved_exports.R`, the integration writers, the two report templates and the
proteomics workbooks scanned for the same patterns; nothing further on the label-free and TMT routes.
