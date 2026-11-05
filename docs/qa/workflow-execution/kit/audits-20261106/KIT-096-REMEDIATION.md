# KIT-096 remediation (audit_and_repair, 2026-11-06)

- **Authority:** the user, 2026-11-06, "fix all recommendations". The F01 decision, verbatim: "reusing existing code and extending it to accommodate a new omic type would be ideal" (option (b)).
- **Audit:** [KIT-096-AUDIT.md](KIT-096-AUDIT.md), unchanged.
- **Edited:** `tickets/kit/KIT-096.md` only, revision 1.0 to 1.1. No code, other ticket, index or plan was touched. Checkboxes are left unchecked.

## Index sync (for the coordinator)

| Field | Change |
|---|---|
| title | unchanged |
| dependencies | unchanged (`KIT-095`) |
| time_estimate | **8h → 12h** (parameterising the lipid import, design and QC modules under unchanged lipid CI, plus the module-CI infrastructure and the residency and probe extensions) |
| description (ticket frontmatter) | now says the tabs are built by extending the lipid modules with an omic profile; the index description needs no change |

## Finding closure

| ID | Severity | Disposition | Repair in KIT-096 rev 1.1 | Closure evidence (at implementation) |
|---|---|---|---|---|
| F01 | blocking | Closed by the user's decision | Option (b). The lipid import, design and QC modules take a `profile` whose default is today's lipid behaviour. Transcriptomics uses them, with no copies. Metab is untouched. "Unchanged" is defined as the lipid and metab module CI, the declaration tests and the frozen dump passing unmodified. The false "shared design builder" citation is removed, and the genuinely shared helpers are listed. `files_to_modify` drops `R/mod_trx_import_*.R` and `R/mod_trx_design_*.R` and adds `R/mod_lipid_import_*.R`, `R/mod_lipid_design_*.R`, `R/mod_lipid_qc.R` and `R/mod_lipid_qc_s4.R`. | The lipid and metab module-CI lanes, the declaration tests and the frozen dump pass with their test files unmodified. |
| F02 | should-fix | Closed in ticket | Added to `files_to_modify`: `helper-module-ci-manifest.R`, `helper-module-ci-fixtures.R`, `test-module-ci-manifest.R`, `tests/testdata/mci/fixture_packs.json`, `tests/testdata/mci/fx/transcriptomics/` and `generate-fixtures.R`. Named the MCI ids (from MCI-026) and the lane `transcriptomics-unit`. Fixtures come from the KIT-040 generator. | `test-module-ci-manifest.R`, `test-module-ci-fixture-integrity.R` and `test-module-ci-trx-*.R` pass. |
| F03 | should-fix | Closed in ticket | The acceptance is reworded: the KIT-075 residency walk with `TranscriptomicsAssayData` finds nothing held; the KIT-082 probe passes on the QC tab; import and design reading in process for one action is the stated exception. Positive control in Simulation Setup. The residency and probe test files are added to scope. | The walk and probe fail on the holding variant and pass on the tabs. |
| F04 | should-fix | **Routed to KIT-095** (the hand-off); consumer side closed in KIT-096 | The design commit uses the state name, workflow type `TRANSCRIPTOMICS` and input roles that KIT-095 hands over, with provenance through `recordWorkflowOmicsDesignProvenance(<trx binding>)`. There is no separate annotation file, so there is no stray assay role. The acceptance is reworded: the design commit is identical to the setup of KIT-095's journey, the QC tab commits KIT-095's filter and finalise states, and the record is valid. | KIT-095 names the hand-off; a KIT-096 test compares state metadata and provenance. |
| F05 | should-fix | Closed in ticket | One format, the notebook count table, with the id and name columns chosen by the user. The annotation is the name column. featureCounts and salmon/RSEM are out of scope. The KIT-093 count policy runs at import and its refusal is shown by class, with no rounding option. | Import module-CI accepts the screen corpus and refuses the three invalid fixtures. |
| F06 | should-fix | Closed in ticket | Under the transcriptomics profile, sanitisation is off and run order follows the counts file (not `mixedsort`). An acceptance case requires names and order identical to the corpus design. | That module-CI case; later KIT-098's identity check. |
| F07 | should-fix | **Routed to KIT-093** (kernel rule, including more than two groups); tab side closed | The design tab records the group level order and passes it to the kernels, and it applies KIT-093's multi-group rule. | Design-case test that the level order is recorded and passed. |
| F08 | minor | **Routed to KIT-093** (TMM-scaled pre-pass); tab side closed | The QC tab passes the group factor and level order and no library-size argument. The zero-variance step is shown. The acceptance requires the mask to be identical to the KIT-040 reference. | QC mask identity test. |
| F09 | should-fix | **Routed to KIT-095** (the `filtering_metrics` product for counts); consumer side closed | The QC tab shows the counts metrics (genes before, after `filterByExpr` and after zero-variance; library size and detected genes per sample), the filtering figure, and the raw-logCPM PCA and RLE. TMM and RUVs figures are KIT-097's. Under a managed session these come from the worker. | Managed module-CI case reads them from the worker result; the probe sees no assay read. |
| F10 | minor | Closed in ticket | `tech_rep_group` is committed as `technical_replicate_id`; a group with fewer than two samples is refused at save. | Design module-CI cases. |
| F11 | minor | Closed in ticket | Stable ids fixed: the module, tabs and sub-module ids, the builder's shared input ids and completion text, and the `trx-tab-*` test ids. | Test asserting the ids. |
| F12 | minor | Closed in ticket | The divergence (shinydashboard shell, managed runner in place of `ExtendedTask`, `fileInput` fallback) is recorded under the user's instruction and the scope rule in `R.md`. Functions stay at about 50 lines or fewer, and the plan declares a line budget for each new file. | Outcome records the divergence; plan carries the budget. |
| F13 | minor | Closed in ticket; **e2e part routed to KIT-097** | Wiring narrowed to defining `mod_transcriptomics_ui`/`_server` and removing the stale comment in `app_ui.R`. `R/app_server.R` is dropped from `files_to_modify`. The `.E2E_VALID_OMICS` and "Unsupported omic" (`test-e2e-browser-harness.R:120`) change moves to KIT-097. | Test that selecting transcriptomics renders the module; KIT-097 carries the e2e files. |
| F14 | minor | Closed in ticket; **session Open routed to KIT-097** | No legacy adoption, since transcriptomics has no older build. Session Open and export are KIT-097's. Import and design help text state the whole-number policy and matching by name only. | Ticket text; UI help strings. |
| F15 | minor | Closed in ticket | Acceptance cases: a mismatched design is refused with no positional mapping; unassigned, duplicate, unknown-level, empty and singleton groups are refused by class. | Design module-CI cases. |

**Counts:** 15 findings.
- 10 are closed in the ticket: F01 (by the user's decision), F02, F03, F05, F06, F10, F11, F12, F13 and F15. The e2e part of F13 goes to KIT-097.
- 5 are closed on the consumer side and routed for the producer side: F04 and F09 to KIT-095, F07 and F08 to KIT-093, and session Open (F14) to KIT-097.
- None is unclosable.

**Validation (2026-11-06, after the edit):** `OPENAI_TICKETS_DIR=…/tickets/kit next-openai-ticket.sh KIT-096` gives `validation` valid (no errors, no warnings), `openai_validation` valid, and `status_consistent` true.
