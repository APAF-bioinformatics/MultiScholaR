# KIT-095 remediation (2026-11-06)

**Mode and authority.** Mode `audit_and_repair`. The user authorised the repairs on 2026-11-06 ("fix all recommendations"), together with five decisions:

1. RUVSeq and EDASeq are installed in the kit library by KIT-040.
2. The shared DA saving is extended to accept an edgeR fit. KIT-095 owns this.
3. Existing code is parameterised, not copied.
4. RUVIII is an export after DA, as a product only.
5. The bindings hash is KIT-098's, under protocol v1.8.

**What changed.** Only `tickets/kit/KIT-095.md`, from revision 1.0 to 1.1. The audit report is not rewritten. The ticket validator passes with 0 errors and 0 warnings, and the status is consistent.

## Index and graph changes (for the coordinator)

- **Title:** the file title is unchanged. The index title still differs from it ("…partitioned where exact…" instead of "…partitioned variants where exact…", F15). Sync the index to the file title.
- **Dependencies:** unchanged (`["KIT-094"]`). KIT-093 is reached through KIT-094, so the renderer dependency (F02) needs no new edge.
- **`time_estimate`:** 8h to 10h, because KIT-095 now extends the shared DA saving.
- **`files_to_modify`:** these were added:
  - `R/utils_workflow_saved_entities.R`
  - `tests/testthat/test-workflow-operation-registry.R`
  - `tests/testthat/test-workflow-omics-declarations.R`
  - `tests/testthat/test-workflow-omics-worker.R`
  - `tests/testthat/test-workflow-saved-consumers.R`

## Finding register

### Findings from the KIT-095 audit

| ID | Severity | Disposition | Where in KIT-095 1.1 | Closure evidence (after implementation) |
|---|---|---|---|---|
| KIT095-F01 | blocking | **Closed in ticket** (decision 2). `newSavedWorkflowDaFit` accepts a `DGEGLM` (`engine = "edgeR_QL"`); the `MArrayLM` envelope is unchanged. `fixed`, `da_apply` and the manifest `model` become binding fields. The edgeR model record is set out. The realised formula is recorded, and replay refuses a mismatch. Multiple contrasts are declared as an extension of the notebook. | "DA: the edgeR contract…"; acceptance; backward compatibility | Envelope test (fails first). Manifest `model` reads edgeR. Replay exact. Lipid/metab dump identical. |
| KIT095-F02 | blocking | **Routed to KIT-093**, which owns `R/func_trx_publication.R`. KIT-095 registers the operations and adaptors that call it, mapped per output. | Ownership note; publication table | KIT-093 renderer tests; KIT-095 worker-harness runs of each publication operation. |
| KIT095-F03 | should-fix | **Closed in ticket.** The global statistics are restated: effective library size with the unfiltered TMM factors, the median and CPM cutoff, `MinSampleSize`, the defaults, and the zero-variance step; for logCPM, the filtered library sizes × factors. The kernel facts are KIT-093's, and the declaration follows them. | "Global statistics" | Step-map rows with edgeR 4.10.1 line evidence. Partitioned equals whole on a corpus whose TMM factors are not 1. |
| KIT095-F04 | should-fix | **Closed in ticket** (state/product/dependency table). The core changes it needs (replay switch, dependency parameters, filtering names) are **routed to KIT-094**, and KIT-095 uses its fields. | "States, products and dependencies" | Journey replays every step, including product-only ones, exactly. |
| KIT095-F05 | should-fix | **Closed in ticket.** A TRX `on_template` with no edgeR call on the zero-row template. | "Worker steps and partitioned variants" | Partitioned test, plus an assertion that the template path never calls edgeR. |
| KIT095-F06 | should-fix | **Closed in ticket.** A binding `integration_export` adaptor on the shared handler: lipid/metab keep the S4 RDS; TRX writes the `newY` TSV, samples by genes, with a `Run` column. | Publication table and notes | The TSV hash equals the KIT-040 reference; lipid/metab dump unchanged. |
| KIT095-F07 | should-fix | **Routed to KIT-094** (per-binding publication list). The TRX list includes `trx_ruv_publication`; lipid/metab keep their nine. | Publication section | Registry test changes only by the `trx_*` names. |
| KIT095-F08 | should-fix | **Closed in ticket.** Every notebook output is mapped to a product and operation, or dropped with a reason. PNG/PDF/XLSX are compared by their inputs. | Publication table | Worker-harness outputs; reference hashes for CSV/TSV. |
| KIT095-F09 | should-fix | **Closed in ticket.** The test files are added. Defaults are `fixed`/`no_shipped_default`, so `config.ini` is not edited and the legacy sections are not reused. | Front matter; Registry | Targeted tests pass, and the lipid/metab assertions are unchanged. |
| KIT095-F10 | should-fix | **Closed in ticket.** The step map comes first and any flag stops the ticket for the user. The bindings hash is **routed to KIT-098** (v1.8, decision 5), and the acceptance criterion is reworded to match. | Planning gate; Simulation Setup; acceptance | Map committed before the declarations, with `open_flags: 0` or the answers recorded. |
| KIT095-F11 | should-fix | **Closed in ticket.** The RUV mode, the grid settings, the manual override, `ruviii_m_source`, the record outcomes (including `controls = "inferred"`) and the report-source names are declared here for KIT-097. | "RUV selection"; Registry | Registry and record tests; manual- and skip-mode acceptance. |
| KIT095-F12 | should-fix | **Closed in ticket.** New criteria: the worker harness for every step and publication operation, products compared alongside states, and the registry-shape and refusal tests. | Acceptance criteria | The tests. |
| KIT095-F13 | minor | **Closed in ticket.** A determinism line in the step map's `kernel_observations`. | Step map | The line is present. |
| KIT095-F14 | minor | **Closed in ticket.** "Prepared once" now means the transpose and the ranking; no imputation; no `W1` precomputation. | "Worker steps and partitioned variants" | The ticket text. |
| KIT095-F15 | minor | **Closed in ticket.** Units use `whole_matrix`. The default variant is `whole` until KIT-099; the variant field itself is **routed to KIT-094**. The section-map content is stated. The index title is for the coordinator to sync. | Tables; Registry | Ticket text; index sync. |
| KIT095-F16 | minor | **Closed in ticket.** `gene_qc_dir` and `da_genes` are named, with the notebook's file names. The RUVSeq/EDASeq precondition is per decision 1 (installed by KIT-040), and tests skip when they are absent. | Publication section; Experiment budget | `packageVersion` recorded in `kit-095/`. |

### Findings routed in from other audits

| ID | Severity | Disposition | Where in KIT-095 1.1 | Closure evidence (after implementation) |
|---|---|---|---|---|
| KIT096-F04 | should-fix | **Closed in ticket** (handoff). State `trx_raw_data_s4`; provenance through `recordWorkflowOmicsDesignProvenance(<trx binding>)`; type `TRANSCRIPTOMICS` (KIT-094 field); roles `trx_counts_<n>` and `design_file`; annotation as the name column, with a separate file only under the `trx_annotation` role, never as an assay. | "Inputs, annotation and handoff to the GUI" | KIT-096's test that the GUI design commit equals the journey set-up. |
| KIT096-F09 | should-fix | **Closed in ticket.** Count `filtering_metrics`: genes before, after `filterByExpr` and after zero-variance; library size and detected genes per sample. It is the filter's display product. | Publication table | Worker-harness filter run returns the metrics. |
| KIT096-F08 | minor (routed to KIT-093/095) | **Covered by F03.** | "Global statistics" | As F03. |
| KIT097-F04 | should-fix | **Closed in ticket.** Automatic, manual or skip mode; grid settings; chosen values and scores as parameters and record outcomes; a manual `k` above the cap is refused (`trx_ruv_k_above_cap`); skip means no W and an unavailable RUVIII. | "RUV selection" | Manual- and skip-mode acceptance. |
| KIT097-F06 | should-fix | **Closed in ticket.** The notebook-output to product and operation map, paired with KIT-093's renderers. KIT-097 maps these products to tabs, the summary and the report. | Publication table | KIT-097's e2e test asserts each mapped output. |
| KIT097-F08 | should-fix | **Closed in ticket.** `trx.ruviii` `newY` TSV (the notebook format, two paths, `M` recorded) through the parameterised shared integration export. It is the declared input KIT-061 will consume; the integration workbook's file name matches. | Publication table and notes | The TSV is identical to the reference `newY`. |

## Not closable here

- **Repairs that depend on other tickets.** These are closed only when the owning tickets land:
  - KIT-093: the renderers (F02) and the kernel pre-pass arguments (F03);
  - KIT-094: the dependency, replay, publication-list, filtering-names, default-variant and workflow-type binding fields (F04, F07, F15, KIT096-F04);
  - KIT-098: the bindings hash (F10).

  KIT-095 1.1 consumes those fields by role. Their exact names are fixed in KIT-094's plan.
- **Decisions taken in this repair that the user may reverse:**
  - A manual `k` above the cap is **refused**. The audit offered "refused or warned".
  - The skip mode makes the RUVIII export unavailable.
