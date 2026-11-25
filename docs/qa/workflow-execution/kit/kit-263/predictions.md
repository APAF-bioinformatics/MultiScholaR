# KIT-263: predictions and the branch-point re-check (written before any run)

Branch `kit-263` from `origin/v0.9-staging` 127011429f5160e6a1ed57fd0552dc16a90d5aa9. Ticket r1.1.

## r1.2 note: re-read of every cited location at 12701142

Between 0a774e0b (the plan audit) and 12701142 only `docs/methods-decision-ledger.md` changed (`git diff --stat`).
No cited line moved: the writer's signature `:12-17`, the config_list resolution `:83-91`, the merge `:151-167`,
the Summary core's call `:264-266`, the recorder `test-omic-summary-module.R:451-458` (mocked at `:487`), the Rd usage
`:6-17`. The filter goes after `:91` (the `config_list` log line), before the merge at `:151`; nothing between reads
`s4_params` or `config_list` except `:643` (`buildWorkflowRuvOptimizationLines`, RUV sections, none foreign).
The ticket's grep recount: 44 files (41 test files + 3 helpers), the same 41 test files as the ticket's list; with
`test-gui-entry-proteomics.R` and the new file, 43 files.

## Predictions

R1 (tests-only commit, base = 12701142):
- `test-kit263-prot-study-parameters-sections.R`: T1 fails (6 foreign lines against 0), T5 fails (`NULL` against
  "transcriptomics"), T4 errors (`studyParameterForeignSections` not found), T2 and T3 pass.
  Stop rule 1: T1 at 0 on base.

R2 (fix against tests-only commit; 43 files + R CMD check; env MULTISCHOLAR_TRX_MANIFEST_NOTE=1):
- every file `same`, except the new file (`only-on-head`; all five tests pass on head);
- `test-omic-summary-module.R`: `same`, its failures identical on both sides (merge train 11: 1675 passed, 28 failed);
- `test-gui-entry-proteomics.R`: `same` (failures, if any, identical on both sides);
- R CMD check: identical on both sides (the Rd matches the roxygen).
