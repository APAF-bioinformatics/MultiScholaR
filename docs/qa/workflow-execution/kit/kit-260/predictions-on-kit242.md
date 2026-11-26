# KIT-260 on KIT-242: predictions, written before the run

Run: base = origin/kit-242 8162ea52 (KIT-242 Part B on staging aa119308), head = kit-260 after the single merge of
origin/kit-242 (666ace2e) and these records; the same 130 files; R CMD check both sides; env
MULTISCHOLAR_TRX_MANIFEST_NOTE=1 (the list reaches a transcriptomics file). The comparison shows KIT-260's change alone.

No test or expectation of KIT-260 changed in the merge; failing-first by value stands from R1 45618785756 (base
tests-only 63e29c04) and R1b 45635735448 (base 7af2d51e). CC01's chunk edit (c545f24e) runs here for the first time.

| File | Predicted verdict | Why |
|---|---|---|
| `test-kit260-itsd-default-pattern.R` | `only-on-head`, every test passes | KIT-242's ITSD changes leave these objects' results as before: both assays use standards (E3, E4, E6, E8, I02, I10); G4's assay A2 matches no standard and is skipped with KIT-242's notice, the test reads A1 only; I05's older record runs under the earlier token for both tickets (its product holds only `final`); the mirror replays exactly. The report tests need CC01's chunk (`ReportKnown()`) |
| `test-omic-norm-module.R` | fails on both sides (kit-242 583/56); `message-changed` | the ITSD captures' `log`, `notes`, `outputs` and `provenance` keys hold neither branch's combined values: base fails them with KIT-242's lines alone, head with both tickets' lines (KIT-242's count line, KIT-260's per-assay rows lines and notification). The keys that move are listed by file and line after the run, not edited (VAL-303 and VAL-352 cover the ITSD marker's digests) |
| `test-omics-norm-core.R` | `count-changed` (one more passing expectation on head) | KIT-260's edited expectation (its notification first, the error second) holds on the merged tree: KIT-242 raises no notice on that run |
| `test-integration-gui.R` | `same` | VAL-275's line and the fact found by its text (513ae741) |
| `test-metab-00-import-detection-characterization.R`, `test-metab-import-detection-direct-shared.R`, `test-workflow-omics-import.R` | `same` | base and head each pass their own expectations (VAL-241 (i), (ii)) |
| `test-workflow-lipid-declarations.R`, `test-workflow-metab-declarations.R` | `same` | the ITSD declaration's `fixed` is identical on both sides (one shared marker) |
| `test-workflow-omics-declarations.R`, `test-workflow-omics-binding-steps.R` (KIT-073, KIT-094 laptop fixtures) | fail on both; `same` or `message-changed` | KIT-260's record outcomes (`itsd_features_per_assay`) may appear in the characterisation; deferred to the re-record (VAL-241 (iii), VAL-181, VAL-189) |
| `test-gui-entry-lipidomics.R`, `test-gui-entry-metabolomics.R` | `same` (19/16 and 57/48 on both) | unchanged by KIT-260 in every earlier comparison |
| `test-prot-norm-core.R` | `same` (1547/307 on both) | proteomics; untouched |
| every other file | `same` | — |

R CMD check: findings identical on both sides; nothing only on head.

Budget: before this run 1,753.5 job-minutes (116.9 core-hours) of 3,020 (201.33). This run's ceiling 1,100 job-minutes;
expected about 490, so the total stays below the ticket's ceiling even at this run's ceiling (2,853.5). One correction
pass remains possible only if its run fits the remainder.
