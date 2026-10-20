# KIT-042 audit remediation (2026-10-20, user: "fix findings including medium low")

Register: [AUDIT-20261020.md](AUDIT-20261020.md) (unchanged). Repairs at `c1189ef9`; evidence reruns on that commit; the install replaces `kit_042`. The three decisions the audit left open were taken as proposed: F-01 marks the entries, F-02 refuses by name, F-03 rewords the criterion.

| ID | Repair | Closure evidence | State |
|---|---|---|---|
| F-01 | `diaStepProvenanceEntry` carries `adopted = TRUE`; the ingest writes it on every legacy entry; the record builder treats an `adopted` entry as a backfill (software and defaults not attributed). | Unit test: both adopted steps' record `software$sha` and `defaults_sha256` are NA; the KIT-018 DIA test asserts `adopted` on every adopted state; adoption check `legacy-at-ingest` requires the mark on all four cases. | closed |
| F-02 | The ingest refuses a state of another class by name (`workflow_session_unsupported`, naming the state and class, pointing the user to export before enrichment) before the self-contained check. | Unit test: a session with an `EnrichmentResults` state is refused with the class and the state name; adoption check `enriched`: refused, owner unchanged. | closed |
| F-03 | Criterion 3 reworded (ticket revision 1.3): the archive's `manifest.json` names the type; the latest pointer keeps the unchanged native tar schema. | Ticket; the unit test asserts the pointer's `archive`. | closed |
| F-04 | A `<TYPE>-native` owner (adopted from a session without a recorded format) accepts any import of its type and refuses another type; the dead `tolower()` removed. | Unit test (TMT session without `data_format`): `TMT-native`, `pd_tmt` re-import returns the manager, `maxquant` refused. | closed |
| F-05 | `ruv_normalised_for_da_analysis_obj` released with `data_tbl`/`data_cln` when an owner is replaced; the Norm tab's `norm_data` references handed to KIT-048 (ticket section "Carried from KIT-042"). | Code; `tickets/kit/KIT-048.md`. | closed (handed off) |
| F-06 | Outcome reworded (the three deliberate DIA deltas named); the KIT-018 test asserts the DIA provenance delta. | Outcome; test. | closed |
| F-07 | Two distinct hydrated objects; the trivial post-refusal identity replaced by "no adapter attached and no new staging directory"; `latest$archive` asserted; a TMT unit case added. | `test-workflow-native-session.R` 606 pass. | closed |
| F-08 | Handed to KIT-046 (ticket section "Carried from KIT-042"). | `tickets/kit/KIT-046.md`. | closed (handed off) |
| F-09 | Handed to KIT-045 (naming). | `tickets/kit/KIT-045.md`. | closed (handed off) |

## Reruns at `c1189ef9`

- Tests: `test-workflow-native-session.R` 606 (two KIT-042 tests), `test-workflow-analysis-record.R` 177, `test-workflow-checkpoint-dia-state.R` 236, `test-prot-07b-da-handlers-compat.R` 55, `test-workflow-saved-exports.R` 276; 0 failures, 0 errors. (One commit was made while an over-strict staging-directory assertion still failed; the assertion was corrected and folded into the same commit before it was pushed.)
- Adoption check: `head_commit = c1189ef9`, clean tree; the four adoption cases as before (component, identical, `legacy` marked `adopted`, DA 3/3 equal to the control, typed export and reopen); `limpa` and `enriched` refused with the named class, the owner unchanged.
- `R CMD check --no-manual --no-tests` at `c1189ef9` (clean worktree): 3 WARNINGs, 3 NOTEs, the baseline ([log](check-c1189ef9/r-cmd-check.log)).
- Install `kit_042` = `core/browser-smoke.YfqkNV` (source commit `c1189ef9` verified), replacing `browser-smoke.8xsfwt`; bindings regenerated.

## States after remediation

`implementation_verified`: all items. `execution_complete`: adoption check (both phases, six sessions), tests, check, install. `decision_supported`: LFQ/TMT sessions saved at correlation filtering adopt, hydrate identically, reproduce the control build's DA and export natively; sessions saved after enrichment are refused by name (adoption of enrichment results is outside this ticket). `promotion_authorized`: none sought. `claim_scope`: as the outcome states.
