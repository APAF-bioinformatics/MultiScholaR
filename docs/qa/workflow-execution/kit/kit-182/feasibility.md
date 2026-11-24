# KIT-182 feasibility gates F1 to F3 (ticket section 8)

Answered by the new test file on both commits of R1 and its instrument reruns (runs 45294678011, 45302634386,
45306797628). Values are the `KIT-182 F1` and `KIT-182 F2` log lines of the test logs, identical on base and head.

## F1: the lineage predicate

| Case | What | Expected | Observed (both commits) | Clauses that fired | Shipped predicate (head) |
|---|---|---|---|---|---|
| (a) | Protein Intensity Filter on PRF1 | FALSE | FALSE | none | equal |
| (b) | the filter on head `normalized` | TRUE | TRUE | i, ii | equal |
| (c) | the replicate filter (`rows`) on head `normalized` | TRUE | TRUE | i, ii | equal |
| (d) | manual, correlation applied, filter on `correlation_filtered` | TRUE | TRUE | ii | equal |
| (e) | as (d), correlation skipped | TRUE | TRUE | ii | equal |
| (f) | as (d) after `commitDIADaArguments` re-commits `correlation_filtered` | TRUE | TRUE | ii | equal |
| (g) b | (b) with the QC step on the worker path | TRUE | TRUE | i, ii | equal |
| (g) d | (d) with the QC step on the worker path | TRUE | TRUE | ii | equal |
| (h) | (d) after a session export, a reopen into a new owner and the tab's restore (R1b repeat) | TRUE | TRUE (head side of the R1b repeat; the file is absent on its staging base) | ii | equal |
| (i) | QC from PRF1 after a normalisation and a Revert | FALSE | FALSE | none | equal |
| (j) | adopted head saved while `ruv_corrected` is in the history | TRUE | TRUE | iii only | equal |
| (k) | design re-saved after RUV, and a QC change on it | FALSE, FALSE | FALSE, FALSE | none | equal |

Clause (i) alone gives FALSE in (d) to (h) and (j), as the ticket states (asserted in (d) to (f)). Every expectation was
met; F1's stop rule did not apply.

Limitations (observed): (g) the worker's dispatch, after the controller's acceptance committed the state, stopped in the
helper's display read ("selective tables require readSavedWorkflowTable", logged); the predicate was evaluated on the
committed head.

(h), first run (R1, 45294678011): `runProtNormExportSessionWorkflow` stopped in `describeDIANativeFile` →
`workflowTableFileIdentity` ("saved table changed during read", which that function raises for a missing, linked or
directory path; `R/utils_workflow_saved_table_codec.R:486-488`): a declared dependency was absent at the owner's root. The
test passed a temporary `sourceDir` other than the owner's `paths$source_dir` (inferred cause). The two instrument repairs
that followed first reordered and then removed the export; that narrowed the declared case (implementation audit I01)
and the table row above for (h) is from that narrowed case (restore only, same owner). Restored as declared in the
correction pass (R1b, 45326674826): with `fixture$paths$source_dir` the export and `openDIANativeSession` completed;
`adoptDIANativeSession` then stopped in `validateProtDaFilteredSession` with "Filtered session data does not contain a
valid contrasts table." (the fixture's workflow carries no `contrasts_tbl`; the application's design step sets one).
Stopped for the coordinator; the owner allowed one test line (VAL-138): the fixture's workflow carries a
`contrasts_tbl` before the export, as `test-workflow-native-session.R` sets it. In the R1b repeat (45334768790) F1 (h)
ran to its end on the merge 1644014f: export with the project's source folder, `openDIANativeSession`,
`adoptDIANativeSession` into a new owner, `restoreProtNormTabState`, lineage TRUE (clause ii, shipped predicate equal),
and "Run Normalization & RUV" on the reopened project refused with class `prot_norm_input_after_normalisation`.
F1 (h) is answered as declared.

## F2: the QC change removes proteins

Protein Intensity Filter (`strict_mode = FALSE`, `min_reps_per_group = 2`, `min_groups = 2`, percentile 10) on PRF1:
64 to 58 proteins, k = 6 (inferred 4 to 7). On the RUV-corrected tables of the S1 cases (63 proteins: the post-RUV clean-up removes one, as KIT-123 records) 63 to 58.
The fallback (drop the first protein row) was never used.

## F3: the DIA fixture as LFQ and TMT

T1 ran in full on all three routes (in process and managed) after `manager$setWorkflowType()`; no failure on either
commit, so the request-only fallback was not used.

## T8 managed: the drive

The KIT-122 drive produced the image on this fixture: the full T8 managed assertions apply (no limited claim). Mechanics
repaired once: the managed chain's component-read counter (`helper-workflow-omics-worker.R`) passes `fields = NULL` on,
which reads no assay table; the scripted QC change in the application now reads through the package's own function.
