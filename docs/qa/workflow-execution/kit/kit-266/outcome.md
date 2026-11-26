# KIT-266 outcome (branch kit-266, from origin/v0.9-staging 7b4ae7da)

Work item (AGENTS.md, stated before editing): the DIA-with-limpa route builds Figure 1 (VAL-308; answers VAL-316; ceiling
VAL-333). Ticket KIT-266 r1.5. Files: the ticket's section 12 list. Checks: R1 (tests-only, failing first), R2 (base against
head, 55 files and R CMD check), R3 (render), M1 (controller memory). Exit: sections 10 and 13 of the ticket. One correction
pass; instrument defects per VAL-54 (at most three re-sends across R1 and M1).

## Step 0: the ticket's lines read against the merged tree (7b4ae7da)

Read: AGENTS.md, conventions (4), IMPLEMENTER-BRIEF, KIT-266 r1.5, AUDIT-KIT-266-plan-20261126 (Re-check 2 in full), ledger rows
VAL-308, 316, 333, 334, 342, 343, 346, 348, 349.

G0's reuse condition (ticket section 11): `git diff 894b9cc8 7b4ae7da` over the files G0 observed, by hunk. No observed hunk
changed: the tab entry (`R/mod_prot_norm_workflow_helpers.R` 374-458), the gate and the aesthetics path
(`R/mod_prot_norm_qc_support_helpers.R` 116-155), the pre generator and source rule (`R/mod_prot_norm_qc_generation_helpers.R`
1-160), the pass-through (`R/utils_workflow_dia_norm.R`, now 745-769), the Reset (`R/mod_prot_norm_session_helpers.R`), Copy's
composite entries and the harness (`tests/testthat/helper-gui-entry.R`) are unchanged in the observed hunks. G0's premise
stands; G0 is not repeated.

What moved (lines now, at 7b4ae7da; the ticket's are at 176f9d3f / origin/kit-259b):

| Item | Ticket | Now |
|---|---|---|
| KIT-259 merged (train 23): `protNormFigureTitle(filename, adjusted)` | `origin/kit-259b:338-445` | `R/utils_workflow_dia_norm.R:339`; `writeProtNormFigurePlaceholder` :360; `drawProtNormQcFigure` :396; `writeProtNormCompositePlaceholders(ruvMode, qcDir, omicType)` :424 (no `resolveInputsFn` yet); `withProtNormFigureNotice` :464 |
| resolver / builder / pass-through | `:617`, `:673`, `:702-725` | :617, :673, :745-769 (+43 above them); `diaLimpaPassThroughSource` :777; `removeProtNormOwnFiles` :850 |
| KIT-271 merged (train 26): `resolveProtNormCorrelationTableNames` | `:742-768` | :822; Copy's limpa branch `R/func_general_filemgmt_results_summary_helpers.R:386-400` |
| composite call and list call in the publication | `:1107`, `:1110` (kit-259b) | `R/utils_workflow_dia_products.R:1107-1110` |
| `computeWorkflowDiaQcFiguresAction`, `diaQcFiguresParameters`, twin `saveWorkflowDiaCompositeFigure`, `writeDIACompositeFigureFiles` | `:2090`, `:2068`, `:2177-2221`, `:2232` | :2095, :2073, :2182-2226, :2237 |
| module helpers: covariates `finish` | `:100-101` | :100-101 (the observer's call is now wrapped in `withProtNormFigureNotice`, :107-113); summary text :180; `dispatchProtNormPassThrough` :268-283; `dispatchProtNormRuvCovariates` :293-325 (`publish` wrapped :300-310); `runProtNormNormalizationWorkflow` limpa branch :504-533; Apply :685-704; `dispatchProtNormCorrelation` :733-765; Skip :785; `dispatchProtNormNormalization` `preQc` wrapper :656-658 |
| limpa template | `:586-631`, `:1429-1436`, `:1476`, `:1486` | Result Section :577-640, Main results call :1429-1436, `fig.cap` :1476, `QC_metrics_placeholders` chunk :1486-1517 (`describe()` maps `ruv_corrected_cancor` only :1502-1503) |

Nothing contradicts the ticket's design. Ticket predictions restated on top of KIT-271:
- `prot_dia_limpa.json` moves as the ticket's section 8 says, plus KIT-271's 13 lines (the table is named
  `normalised_results_cln_with_replicates.tsv` / `Publication_tables/normalised_results.tsv`); not edited here.
- `test-gui-entry-proteomics.R` is 92/101 on staging, with the four file-name key lines already failing at the ten steps from
  `correlation` on (VAL-349 (c)). R2's prediction: head failures = base failures plus the limpa steps from `norm` on that pass on
  base (they gain the composite and `pre_norm_*` keys), compared with the four file-name key lines masked, as merge train 26 did.
- `test-prot-norm-core.R` stays 1547/307 on base (stored values already failing); its declared edits add base-only failures
  (section 8 item 6) and no head-only failure.
