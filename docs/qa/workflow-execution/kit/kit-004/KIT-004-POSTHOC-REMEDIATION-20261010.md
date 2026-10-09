# KIT-004 audit remediation (2026-10-10)

Register: `KIT-004-POSTHOC-AUDIT-20261010.md` (unchanged). Repairs at `3bf3ff91`; user: "fix findings".

| ID | Repair | Closure evidence |
|---|---|---|
| PC-01 | The callr child loads the source tree with `pkgload` only when `DESCRIPTION` exists at the root, otherwise the installed package, and attaches the namespace so the unexported helpers resolve either way. | `test-workflow-operations.R` 66 passed from the source tree. The installed branch could not be exercised here: the only installed build (`kit_001`) predates the runner; it runs with the next installed suite. |
| PC-02 | Post pass rerun on the committed code. | `fingerprints-post.json` now records `head_commit 3bf3ff91`; `compare-pre-a-vs-post.txt` regenerated: same two exclusions, nothing else. The uncommitted-tree pass kept as `fingerprints-post-uncommitted-tree.json`. |
| PC-03 | The reference branch of `runProtNormApplyCorrelationObserver` and `runProtNormSkipCorrelationObserver` runs the declaration and the completion inside `withDIAStateScope(getDIAActiveManager(), …)`, as the legacy branch does. | The committed post pass drives the apply observer through that path: `correlation_filtered` state, products, holders identical to baseline. |
| PC-04 | New test in `test-workflow-saved-qc.R`: the module edge writes `config_list$normaliseBetweenSamples$method` and `config_list$ruvParameters` (k, percentage, control count) on the computed pass and again on the identity hit; the completion notification fires on both. | `test-workflow-saved-qc.R` 120 passed (was 108). |
| PC-05 | Hand-off line in the outcome: the adapter delegates to `resolveProtRuvParameters` when KIT-011/012 retire the legacy path. | outcome.md |
| PC-06, PC-07, PC-08 | No change required; recorded. | — |

Closeout after remediation: implementation_verified true; execution_complete true; decision_supported true (same claim scope); promotion_authorized false. Open: PC-01's installed branch is verified by the next installed suite run.
