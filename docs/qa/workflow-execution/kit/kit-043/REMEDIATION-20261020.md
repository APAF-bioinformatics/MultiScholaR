# KIT-043 audit remediation (2026-10-20, user: "fix findings")

Register: [AUDIT-20261020.md](AUDIT-20261020.md) (unchanged). Repairs at `c78bd4af`; evidence reruns on that commit; the install replaces `kit_043`. The two decisions the audit left open were taken as proposed: F-01 repaired now; F-02 declared (the record is the authority, KIT-033's rule), not scoped away.

| ID | Repair | Closure evidence | State |
|---|---|---|---|
| F-01 | The duplicate `max_groups_percentage_cutoff` entry deleted; the test asserts the exact note text and `anyDuplicated(names(workflowRecordReportSources)) == 0`. | `test-workflow-analysis-record.R` 201 pass. | closed |
| F-02 | Declared in the source map's comment, the outcome and the ticket (revision 1.2): with a record the DIA-NN DA text takes the record's q-value and fold-change cutoff; a DIA-NN-shaped parameter list in the test shows the values taken and no new note. | Test; outcome; ticket. | closed (declared) |
| F-03 | The outcome states the rendered "2st" text; the hard-coded ordinal suffix in the four templates handed to KIT-050 (fixing it here would touch the DIA-NN templates this ticket must leave untouched). | Outcome; `tickets/kit/KIT-050.md` "Carried from KIT-043". | closed (handed off) |
| F-04 | A writer test on a protein-only owner (`LFQ` store, `design_builder` provenance, a `protein_intensity_filter` step with a user-set percentile): `study_parameters.txt` carries the record section, no truncation, every recorded parameter, the user-set value. | `test-workflow-analysis-record.R` 201 pass. | closed |
| F-05 | `min_reps_per_group` and `min_groups` entries dropped (the templates never print them); the test asserts the variable is absent. | Test; outcome reworded. | closed |
| F-06 | Exact-note assertion; the DIA-NN-shaped case. | Test. | closed |
| F-07 | Handed to KIT-050 with F-03. | `tickets/kit/KIT-050.md`. | closed (handed off) |
| F-08, F-09 | Claim limits, already stated in the outcome. | Outcome. | closed |

## Reruns at `c78bd4af`

- Tests: `test-workflow-analysis-record.R` 201, `test-general-filemgmt-report-helpers-shared.R` 84; dependents `test-module-ci-prot-summary.R` 110, `test-prot-12` 446, `test-prot-12a` 62, `test-module-ci-prot-da.R` 245, `test-e2e-test-mode-integration.R` 63, `test-metab-05` 160; 0 failures, 0 errors.
- Render check on the installed library: `head_commit = c78bd4af`, clean tree, library commit `c78bd4af`; both fixtures as before (record values printed, software from the record, no-record renders equal to the baseline and the control, DIA-NN templates untouched, `study_parameters.txt` complete).
- `R CMD check --no-manual --no-tests` at `c78bd4af` (clean worktree): 3 WARNINGs, 3 NOTEs, the baseline ([log](check-c78bd4af/r-cmd-check.log)).
- Install `kit_043` = `core/browser-smoke.AuaIVL` (source commit `c78bd4af` verified), replacing `browser-smoke.fFZHS8`; bindings regenerated.

## States after remediation

`implementation_verified`: the source map holds each variable once; the two templates read the record as the DIA-NN template does. `execution_complete`: renders (installed library), tests, check, install. `decision_supported`: LFQ/TMT reports print the record's values and software and fall back unchanged without it; the DIA-NN DA text follows the record (declared). `promotion_authorized`: none sought. `claim_scope`: as the outcome states (five-protein fixtures, RUV skipped).
