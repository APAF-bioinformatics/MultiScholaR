# KIT-014 post-hoc audit remediation (2026-10-13, user: "fix findings")

Register: `KIT-014-POSTHOC-AUDIT-20261013.md` (unchanged). Stage `kit-014` evidence and decision untouched (the decision was `accepted` with no refusal, so the PH-01 rule change does not alter it).

| ID | Closure | Evidence |
|---|---|---|
| PH-01 | Repaired. The verifier's journey-oracle branch names an inexact state in its log (`oracle mismatch: <fields>`); `kitSingleDecision` labels a member `rejected` only when its verifier refusal names an oracle mismatch (failing rule `oracle_mismatch`), and `invalid` for a verifier refusal of any other kind. | `test-workflow-kit-harness.R`: the single-member test covers an oracle mismatch (`rejected`) and a binding refusal (`invalid`); bindings regenerated (harness hashes). |
| PH-02 | Recorded as a claim limit in the ticket's acceptance evidence and handed to KIT-017's closure checklist (declaration without admission). No code change (bounding the streamed path is a mechanism of its own). | `tickets/kit/KIT-014.md` acceptance line; `tickets/kit/KIT-017.md` hand-off. |
| PH-03 | Hand-off to KIT-017 (protocol v1.4 line is a user decision). | `tickets/kit/KIT-017.md` hand-off. |
| PH-04 | Hand-off to KIT-017 (the closure journey pair observes the five modules under a real browser). | `tickets/kit/KIT-017.md` hand-off. |
| PH-05 | Note only; no change. | — |
| PH-06 | No change (revision 1.1 is the record). | — |
